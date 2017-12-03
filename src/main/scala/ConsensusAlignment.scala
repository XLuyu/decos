import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.Set
import scala.collection.mutable.Map
import scala.io.Source
import java.io.RandomAccessFile

import scala.collection.mutable
import org.apache.spark._

/**
  * Created by workshop on 04-Jul-17.
  */
object MappingRead {
  def apply(fileno: Int, offset: Long, pos: Int, complement: Boolean, file: RandomAccessFile): MappingRead = {
    file.seek(offset)
    val head = file.readLine().length + 1
    val seq = file.readLine()
    file.readLine()
    val qual = file.readLine()
    new MappingRead(fileno, offset + head, pos, seq, qual, complement)
  }
}

class MappingRead(fileNumber: Int, offset: Long, pos: Int, rawSeq: String, quality: String, complement: Boolean = false) {
  val fileno = fileNumber
  val id = offset
  var complemented = complement
  var kmeridx = pos
  var seq = if (complement) Util.reverseComplement(rawSeq) else rawSeq
  val qual = if (complement) quality.reverse else quality
  var column = Array.fill(seq.length)(-1)
}

class ConsensusSequence(init: String) {
  var columnID = init.indices.toBuffer
  var vote = ArrayBuffer.fill(init.length)(mutable.Map('A' -> 0, 'C' -> 0, 'G' -> 0, 'T' -> 0, 'N'-> 0 ,'-' -> 0))
  for (i <- 0 until init.length) vote(i)(init(i)) += 1

  def charAt(idx: Int): Char = vote(idx).maxBy(x => x._2 * 1000 - x._1)._1

  def addCount(idx: Int, c: Char): Unit = {
    vote(idx)(c) += 1
  }

  def sequence(keepGap: Boolean = false): String = {
    val consensusWithGap = for (i <- vote.indices) yield charAt(i)
    if (keepGap) consensusWithGap.mkString else consensusWithGap.filter(_ != '-').mkString
  }

  def insert(idx: Int, newChars: Array[Char]): Unit = {
    columnID.insertAll(idx, columnID.size until (columnID.size + newChars.length))
    val newColumn = for (c <- newChars) yield mutable.Map('A' -> 0, 'C' -> 0, 'G' -> 0, 'T' -> 0, 'N'-> 0 ,'-' -> 0) + (c -> 1)
    vote.insertAll(idx, newColumn)
  }

  def checkGap(idx: Int, newChars: Array[Char]): Unit = {
    var gapCount = 0
    var i = idx
    while (i < vote.length && gapCount < newChars.length && vote(i).maxBy(_._2)._2 == vote(i)('-')) {
      vote(i)(newChars(gapCount)) += 1
      i += 1
      gapCount += 1
    }
    if (gapCount < newChars.length) insert(i, newChars.slice(gapCount, newChars.length))
    for (j <- 0 until newChars.length - gapCount)
      vote(i + j) += ('-' -> (vote(idx - 1).values.sum - 1))
  }
}

object ConsensusAlignment {

  val K = Settings.K
  val NWaligner = new NeedlemanWunschAligner(Settings.MAX_READ_LEN * 2, Settings.MAX_INDEL)

  def alignOneEnd(groupSeq: String, readSeq: String, readST: SuffixTree): (Int, Int, Array[Int]) = {
    // return mapping array like: -2 -2 0 1 2 -1 3 4 6 -3 -3 -3
    def sharp(a: (Int, Int, Int), b: (Int, Int, Int)): (Int, Int, Int) = { // remove overlap part from a
      val (a1, a2, a3, a4) = (a._1, a._1 + a._3, a._2, a._2 + a._3)
      val (b1, b2, b3, b4) = (b._1, b._1 + b._3, b._2, b._2 + b._3)
      if (Math.abs(a1 - b1 - a3 + b3) > Settings.MAX_INDEL) return null
      if (a1 < b1 && a2 <= b2 && a3 < b3 && a4 <= b4)
        return (a1, a3, a._3 - Math.max(Math.max(a2 - b1, a4 - b3), 0))
      if (b1 <= a1 && b2 < a2 && b3 <= a3 && b4 < a4) {
        val ovlp = Math.max(Math.max(b2 - a1, b4 - a3), 0)
        return (a1 + ovlp, a3 + ovlp, a._3 - ovlp)
      }
      null
    }

    def compactLength(a: (Int, Int, Int)): Int = {
      var i = a._2
      var j = a._2 + a._3 - 1
      while (i < a._2 + a._3 - 1 && groupSeq(i) == groupSeq(i + 1)) i += 1
      while (j > a._2 && groupSeq(j) == groupSeq(j - 1)) j -= 1
      Math.max(j - i + 1, 1)
    }

    val mapping = Array.fill[Int](groupSeq.length)(-1)
    val lengthThreshold = Math.log(readSeq.length) / Math.log(2)
    //    println(readST.pairwiseLCS(groupSeq))
    val LCS = readST.pairwiseLCS(groupSeq).filter(_._3 > K).filter(compactLength(_) > 2).sortBy(-_._3)
    //
    //        println(groupSeq)
    //        println(LCS)
    var matchSeg = List[(Int, Int, Int)]()
    for (m <- LCS) {
      var sharped = m
      for (nm <- matchSeg)
        if (sharped != null) sharped = sharp(sharped, nm)
      if (sharped != null) matchSeg = sharped :: matchSeg
    }
    matchSeg = matchSeg.sortBy(_._2)
    if (matchSeg.isEmpty) return (0, 0, null)
    var matchCount = 0
    var spanCount = 0
    val unmatchedReadHead = Math.min(matchSeg.head._1, matchSeg.head._2 + 2)
    val unmatchedGroupHead = Math.min(matchSeg.head._2, matchSeg.head._1 + 2)
    var lastR = matchSeg.head._1
    var lastG = matchSeg.head._2
    val HeadMapping = NWaligner.align(groupSeq.substring(lastG - unmatchedGroupHead, lastG).reverse,
      readSeq.substring(lastR - unmatchedReadHead, lastR).reverse, tailFlex = true)
    for (i <- HeadMapping.indices)
      mapping(lastG - 1 - i) = if (HeadMapping(i) == -1) -1 else lastR - 1 - HeadMapping(i)
    for (i <- lastG - unmatchedGroupHead until lastG if mapping(i) != -1 && groupSeq(i) == readSeq(mapping(i))) matchCount += 1
    spanCount += Math.min(unmatchedReadHead, unmatchedGroupHead)
    for ((rs, gs, len) <- matchSeg) {
      if (gs - lastG == rs - lastR) {
        for (i <- 0 until (gs - lastG)) mapping(lastG + i) = lastR + i
      } else {
        val regionMapping = NWaligner.align(groupSeq.substring(lastG, gs), readSeq.substring(lastR, rs))
        for (i <- regionMapping.indices)
          mapping(lastG + i) = if (regionMapping(i) == -1) -1 else lastR + regionMapping(i)
      }
      for (i <- lastG until gs if mapping(i) != -1 && groupSeq(i) == readSeq(mapping(i))) matchCount += 1
      matchCount += len
      spanCount += Math.max(gs - lastG, rs - lastR) + len
      for (i <- 0 until len) mapping(gs + i) = rs + i
      lastG = gs + len
      lastR = rs + len
    }
    val unmatchedReadTail = Math.min(readSeq.length - lastR, groupSeq.length - lastG + 2)
    val unmatchedGroupTail = Math.min(groupSeq.length - lastG, readSeq.length - lastR + 2)
    val tailMapping = NWaligner.align(groupSeq.substring(lastG, lastG + unmatchedGroupTail),
      readSeq.substring(lastR, lastR + unmatchedReadTail), tailFlex = true)
    for (i <- tailMapping.indices)
      mapping(lastG + i) = if (tailMapping(i) == -1) -1 else lastR + tailMapping(i)
    for (i <- lastG until lastG + unmatchedGroupTail if mapping(i) != -1 && groupSeq(i) == readSeq(mapping(i))) matchCount += 1
    spanCount += Math.min(unmatchedReadTail, unmatchedGroupTail)
    var i = 0
    while (mapping(i) == -1) {
      mapping(i) = -2
      i += 1
    }
    i = mapping.length - 1
    while (mapping(i) == -1) {
      mapping(i) = -3
      i -= 1
    }
    //    for (i <- mapping) printf("%d,",i)
    //    println()
    (matchCount, spanCount, mapping)
  }

  def main(args: Array[String]): Unit = {
    val readArray = ArrayBuffer[MappingRead]()
    var n = 0
    var ca: ConsensusAlignment = null
    var r: (Int, Array[Int]) = null
    for (line <- Source.fromFile("Demo.txt").getLines) {
      n += 1
      readArray += new MappingRead(0, n, 1, line, line, true)
      if (n == 1) {
        ca = new ConsensusAlignment(readArray(0))
      } else {
        r = ca.align(readArray.last, new SuffixTree(readArray.last.seq))
        ca.joinAndUpdate(readArray.last, r._2)
        println(ca.printPileup())
      }
    }
    val result = ArrayBuffer[List[Long]]()
    ca.reportAllEdgesTuple("CCTCTGCTGGCGCTGGTTGCC", result)
    println(ca.printPileup())
    println(result)
  }
}

class ConsensusAlignment(read: MappingRead) extends ArrayBuffer[MappingRead]() {
  var consensus = new ConsensusSequence(read.seq)
  this += read
  read.column = consensus.columnID.toArray
  val ufs = new UnionFindSet[Int]()

  def updateConsensus(consensus: ConsensusSequence, mapping: Array[Int], seq: String, readColumn: Array[Int]) {
    var inserted = 0
    if (mapping(0) > 0) {
      consensus.insert(0, seq.substring(0, mapping(0)).toCharArray)
      for (i <- 0 until mapping(0)) readColumn(i) = consensus.columnID(i)
      inserted += mapping(0)
    }
    var last = -1
    for (i <- mapping.indices) {
      if (mapping(i) > -1) {
        if (last != -1 && last + 1 != mapping(i)) {
          consensus.checkGap(inserted + i, seq.substring(last + 1, mapping(i)).toCharArray)
          for (j <- last + 1 until mapping(i)) readColumn(j) = consensus.columnID(inserted + i + j - last - 1)
          inserted += mapping(i) - last - 1
        }
        while (consensus.charAt(inserted + i) == '-') {
          consensus.addCount(inserted + i, '-')
          inserted += 1
        }
        consensus.addCount(inserted + i, seq(mapping(i)))
        readColumn(mapping(i)) = consensus.columnID(inserted + i)
        last = mapping(i)
      } else if (mapping(i) == -1) {
        while (consensus.charAt(inserted + i) == '-') {
          consensus.addCount(inserted + i, '-')
          inserted += 1
        }
        consensus.addCount(inserted + i, '-')
      } else {
        while (consensus.charAt(inserted + i) == '-') inserted += 1
      }
    }
    if (last + 1 != seq.length) {
      consensus.insert(inserted + mapping.length, seq.substring(last + 1).toCharArray)
      for (i <- last + 1 until seq.length) readColumn(i) = consensus.columnID(inserted + mapping.length + i - last - 1)
      inserted += mapping(0)
    }
    if (readColumn.contains(-1)) {
      println("Exception at Line 365!")
      for (i <- this) println(i.seq)
      println(consensus.sequence())
      for (k <- consensus.columnID) print(k + ",")
      println()
      for (k <- mapping) print(k + ",")
      println()
      for (k <- readColumn) print(k + ",")
      println()
      println(seq)
    }
  }

  def joinAndUpdate(read: MappingRead, mapping: Array[Int]) {
    this += read
    updateConsensus(consensus, mapping, read.seq, read.column)
  }

  def align(read: MappingRead, readST: SuffixTree): (Int, Array[Int]) = {
    val consensusSeq = consensus.sequence()
    val (count, span, mapping) = ConsensusAlignment.alignOneEnd(consensusSeq, read.seq, readST)
    if (mapping == null || count < span * Settings.MATCH_RATE) return null
    //shifting indel alignment to right-most place(or left-most for complemented read)
//    var last = -1
//    if (!read.complemented){
//      for ( i <- mapping.indices) {
//        if (mapping(i) == -1) {
//          var j = i
//          while (j+1 < consensusSeq.length && consensusSeq(j + 1) == consensusSeq(i) && mapping(j + 1) == -1) j += 1
//          if (j+1 < consensusSeq.length && consensusSeq(j + 1) == consensusSeq(i) && mapping(j + 1) != -1) {
//            mapping(i) = mapping(j+1)
//            mapping(j+1) = -1
//          }
//        }
//        if (mapping(i) >= 0) {
//          if (last != -1 && mapping(i) > last + 1)
//            while (mapping(i)-1 > last && read.seq(mapping(i)-1) == read.seq(mapping(i))) mapping(i) -= 1
//          last = mapping(i)
//        }
//      }
//    } else { // complemented
//      for ( i <- mapping.length-1 to 0 by -1) {
//        if (mapping(i) == -1) {
//          var j = i
//          while (j>0 && consensusSeq(j-1)==consensusSeq(i) && mapping(j-1) == -1) j -= 1
//          if (j>0 && consensusSeq(j-1)==consensusSeq(i) && mapping(j-1) != -1) {
//            mapping(i) = mapping(j-1)
//            mapping(j-1) = -1
//          }
//        }
//        if (mapping(i) >= 0) {
//          if (last != -1 && mapping(i) < last - 1)
//            while (mapping(i)+1 < last && read.seq(mapping(i)+1) == read.seq(mapping(i))) mapping(i) += 1
//          last = mapping(i)
//        }
//      }
//    }
    (2 * count - span, mapping)
  }

  def reportAllEdgesTuple(kmer: String, report: ArrayBuffer[List[Long]]): Unit = {
    def baseEncode(fileno: Int, position: Long, gapCount: Int, baseCode: Int, quality: Int): Long = {
      // tried bitwise operation, but it gathers hashCode of records and skew the partition
      var id = (gapCount * 64 + quality) * 1000000000000L + position
      (id * 6 + baseCode) * 2 + fileno
    }

//    val t1 = System.currentTimeMillis()
    val KmerSet = new Array[mutable.Set[Int]](this.size)
    val isAnchoredByMinKmer = Array.ofDim[Boolean](this.size, this.size)
    val thishash = kmer.hashCode
    for (i <- this.indices) {
      val seq = this (i).seq
      KmerSet(i) = mutable.Set[Int]()
      for (j <- 0 to seq.length - ConsensusAlignment.K) {
        val kmerhash = Util.minKmerByRC(seq.substring(j, j + ConsensusAlignment.K)).hashCode
        if (kmerhash < thishash) KmerSet(i) += kmerhash
      }
      for (j <- 0 until i)
        isAnchoredByMinKmer(i)(j) = (KmerSet(i) intersect KmerSet(j)).isEmpty
    }
//    val t2 = System.currentTimeMillis()
    val table = mutable.Map[Char, Int]('A' -> 0, 'C' -> 1, 'G' -> 2, 'T' -> 3, 'N' -> 4, '-' -> 5)
    val columns = Array.fill[ArrayBuffer[(Int, Long)]](consensus.columnID.size)(ArrayBuffer[(Int, Long)]())
    val columnCount = Array.fill[mutable.Map[Char, Long]](consensus.columnID.size)(mutable.Map[Char, Long]('A' -> 0, 'C' -> 0, 'G' -> 0, 'T' -> 0, 'N' -> 0, '-' -> 0))
    for (idx <- this.indices) {
      val read = this (idx)
      var refidx = 0
      var first = true
      for (i <- read.column.indices) {
        var gapCount = 0
        while (consensus.columnID(refidx + gapCount) != read.column(i)) {
          gapCount += 1
        }
        if (!first) {
          for (g <- 1 to gapCount) {
            val baseCode = if (read.complemented)
              baseEncode(read.fileno, read.id + read.seq.length - i - 1, gapCount + 1 - g, table('-'), (read.qual(i - 1) + read.qual(i)) / 2 - 33)
            else
              baseEncode(read.fileno, read.id + i - 1, g, table('-'), (read.qual(i - 1) + read.qual(i)) / 2 - 33)
            columns(refidx + g - 1) += ((idx, baseCode))
            columnCount(refidx + g - 1)('-') += (read.qual(i - 1) + read.qual(i)) / 2 - 33
          }
        }
        refidx += gapCount
        val baseCode = if (read.complemented)
          baseEncode(read.fileno, read.id + read.seq.length - 1 - i, 0, table(Util.complement(read.seq(i))), read.qual(i) - 33)
        else
          baseEncode(read.fileno, read.id + i, 0, table(read.seq(i)), read.qual(i) - 33)
        columns(refidx) += ((idx, baseCode))
        columnCount(refidx)(read.seq(i)) += read.qual(i) - 33
        refidx += 1
        first = false
      }
    }
//    val t3 = System.currentTimeMillis()
    val rtable = Array('A', 'C', 'G', 'T', 'N', '-')
    for (i <- columns.indices if columnCount(i).count(_._2>Settings.CUTTHRESHOLD)>1) {
      val prev = mutable.Map[Char,Int]()
      for ((idx,baseCode) <- columns(i)) {
        var base = rtable((baseCode/2%6).toInt)
        if (this (idx).complemented) base = Util.complement(base)
        if (columnCount(i)(base)>Settings.CUTTHRESHOLD){
          if (!prev.contains(base)) prev(base) = idx
          else ufs.union(idx,prev(base))
        }
      }
    }
    val defaultValue = if (ufs.isEmpty()) 0 else ufs.getAnyKey()
    val father = new Array[Int](this.size)
    for (idx <- this.indices){
      father(idx) = if (!ufs.contains(idx)) defaultValue else ufs.find(idx)
    }
//    val t4 = System.currentTimeMillis()
    for (column <- columns) {
      val node = mutable.Map[Int,List[Long]]()
      var prev = mutable.Map[Int,List[Int]]()
      for ((idx,baseCode) <- column) {
        val s = father(idx)
        if (!node.contains(s)) { node(s) = List[Long](); prev(s) = List[Int]()}
        val k = prev(s).iterator
        if (k.isEmpty || isAnchoredByMinKmer(idx)(k.next))
          node(s) ::= baseCode * (if (this (idx).complemented) -1 else 1)
        prev(s) ::= idx
      }
      for ( (_,nodelist) <- node if nodelist.size>1) report += nodelist
    }
//    val t5 = System.currentTimeMillis()
//    println("In report:",t2-t1,t3-t2,t4-t3,t5-t4)
  }

  def printPileup(): String = {
    var pileup = f"${"Consensus"}%20s " + consensus.sequence(true) + "\n"
    for (idx <- this.indices) {
      val read = this(idx)
      pileup += f"${ufs.find(idx)}%5s${read.id}%12s(${read.fileno}${if (read.complemented) "-" else "+"})"
      var i = 0
      for (column <- consensus.columnID if i < read.column.length) {
        pileup += (if (column == read.column(i)) read.seq(i) else " ")
        if (column == read.column(i)) i += 1
      }
      pileup += "\n"
    }
    pileup + "-----\n"
  }
}