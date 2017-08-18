import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.Set
import scala.collection.mutable.Map
import java.io.RandomAccessFile

import scala.collection.mutable
import org.apache.spark._
import org.apache.spark.graphx._

/**
  * Created by workshop on 04-Jul-17.
  */
object MappingRead{
  def apply(idx: Long, pos: Int, file1: RandomAccessFile, file2: RandomAccessFile): MappingRead = {
    file1.seek(idx)
    file2.seek(idx)
    val seq1 = file1.readLine()
    val seq2 = file2.readLine()
    new MappingRead(idx, pos, seq1, seq2)
  }
}
class MappingRead(idx: Long, pos: Int, file1: String, file2: String) {
  val id: Long = idx
  var reversed = false
  var kmeridx: Int = pos
  var seq1 = file1
  var seq2 = file2
  var column1 = Array.fill(seq1.length)(-1)
  var column2 = Array.fill(seq2.length)(-1)
  if (kmeridx < 0) {
    reversed = !reversed
    kmeridx = -kmeridx
    val tmp = seq1
    seq1 = seq2
    seq2 = tmp
  }
}

class ConsensusSequence(init: String) {
  var columnID = init.indices.toBuffer
  var vote = ArrayBuffer.fill(init.length)(mutable.Map('A' -> 0, 'C' -> 0, 'G' -> 0, 'T' -> 0, '-' -> 0))
  for (i <- 0 until init.length) vote(i)(init(i)) += 1

  def charAt(idx: Int): Char = vote(idx).maxBy(x => x._2 * 1000 - x._1)._1

  def addCount(idx: Int, c: Char): Unit = {
    if (c != 'N') vote(idx)(c) += 1
  }

  def sequence(keepGap: Boolean = false): String = {
    val consensusWithGap = for (i <- vote.indices) yield charAt(i)
    if (keepGap) consensusWithGap.mkString else consensusWithGap.filter(_ != '-').mkString
  }

  def insert(idx: Int, newChars: Array[Char]): Unit = {
    columnID.insertAll(idx, columnID.size until (columnID.size + newChars.length))
    val newColumn = for (c <- newChars) yield mutable.Map('A' -> 0, 'C' -> 0, 'G' -> 0, 'T' -> 0, '-' -> 0) + (c -> 1)
    vote.insertAll(idx, newColumn)
  }

  def checkGap(idx: Int, newChars: Array[Char]): Unit = {
    var gapCount = 0
    var i = idx
    while (i < vote.length && vote(i).maxBy(_._2)._2 == vote(i)('-')) {
      gapCount += 1
      i += 1
    }
    if (gapCount < newChars.length) insert(i, newChars.slice(gapCount, newChars.length))
    for (j <- 0 until newChars.length - gapCount)
      vote(i + j) += ('-' -> (vote(idx - 1).values.sum - 1))
  }
}

object ConsensusAlignment {
  val K = Settings.K
  val NWaligner = new NeedlemanWunschAligner(Settings.MAX_READ_LEN*2)

  def alignOneEnd(groupSeq: String, readSeq: String, readST: SuffixTree): (Int, Int, Array[Int]) = {
    // return mapping array like: -2 -2 0 1 2 -1 3 4 6 -3 -3 -3

    def sharp(a: (Int, Int, Int), b: (Int, Int, Int)): (Int, Int, Int) = { // remove overlap part from a
      val (a1, a2, a3, a4) = (a._1, a._1 + a._3, a._2, a._2 + a._3)
      val (b1, b2, b3, b4) = (b._1, b._1 + b._3, b._2, b._2 + b._3)
      if (Math.abs(a1-b1-a3+b3)>Settings.MAX_INDEL) return null
      if (a1 < b1 && a2 <= b2 && a3 < b3 && a4 <= b4)
        return (a1, a3, a._3 - Math.max(Math.max(a2 - b1, a4 - b3), 0))
      if (b1 <= a1 && b2 < a2 && b3 <= a3 && b4 < a4) {
        val ovlp = Math.max(Math.max(b2 - a1, b4 - a3), 0)
        return (a1 + ovlp, a3 + ovlp, a._3 - ovlp)
      }
      null
    }
    val mapping = Array.fill[Int](groupSeq.length)(-1)
    val LCS = readST.pairwiseLCS(groupSeq, K).filter(_._3 > K).sortBy(-_._3)
    //    println(groupSeq)
    //    println(matches)
    var matchSeg = List[(Int, Int, Int)]()
    for (m <- LCS) {
      var sharped = m
      for (nm <- matchSeg)
        if (sharped != null) sharped = sharp(sharped, nm)
      if (sharped != null) matchSeg = sharped :: matchSeg
    }
    matchSeg = matchSeg.sortBy(_._2)
    if (matchSeg.isEmpty) return (0, 0, null)
    val unmatchedHeadLen = Math.min(matchSeg.head._1, matchSeg.head._2)
    var lastR = matchSeg.head._1 - unmatchedHeadLen
    var lastG = matchSeg.head._2 - unmatchedHeadLen
    var matchCount = 0
    var spanCount = 0
    for ((rs, gs, len) <- matchSeg) {
      if (gs - lastG == rs - lastR) {
        for (i <- 0 until (gs - lastG)) mapping(lastG + i) = lastR + i
      } else {
        val regionMapping = ConsensusAlignment.NWaligner.align(groupSeq.substring(lastG, gs), readSeq.substring(lastR, rs))
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
    val unmatchedTailLen = Math.min(groupSeq.length - lastG, readSeq.length - lastR)
    for (i <- 0 until unmatchedTailLen) {
      mapping(lastG + i) = lastR + i
      if (groupSeq(lastG + i) == readSeq(mapping(lastG + i))) matchCount += 1
    }
    spanCount += unmatchedTailLen
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
    (matchCount, spanCount, mapping)
  }

  def main(args: Array[String]) {
    val read2 = new MappingRead(1, 1, "CCCACAAAGTCCAGCGTACCATAAACGCAAGCCTCAACGCAGCGACGAGCACGAGAGCGGTCAGTAGCAATCCAAACTTAGTTACTCGTCAGAAAATCGG",
      "CCCACAAAGTCCAGCGTACCATAAACGCAAGCCTCAACGCAGCGACGAGCACGAGAGCGGTCAGTAGCAATCCAAACTTAGTTACTCGTCAGAAAATCGG")
    val read1 = new MappingRead(2, 1, "ATCCGGTACCCACAAAGTCCAGCGTACCATAAACGCAAGCCTCAACGCAGCGACGAGCACGAGAGCGGTCAGTAGCAATCCAAACTTAGTTACTCGTCAG",
      "AGCGTACCATAAACGCAAGCCTCAACGCAGCGACGAGCACGAGGAGCGGTCAGTAGCAATCCAAACTTAGTTACTCGTCAGAAAATCGGATCCTACCAAA")
    val read3 = new MappingRead(3, 1, "TCCAGCGTACCATAAACGCAAGCCTCAACGCAGCGACGAGCACGAGAGCGGTCAGTAGCAATCCAAACTTAGTTACTCGTCAGAAATCGGATGGCTGAAT",
      "ATGGGTCCAAGATCCCACAAAGTCCAGCGTACCATAAACGCAAGCCTCAACGCAGCGACGAGCACGAGAGCGGTCAGTAGCAATCCAAACTTAGTTACTC")
    val read4 = new MappingRead(4, 1, "AAAGTCCAGCGTACCATAAACGCAAGCCTCAACGCAGCGACGAGCACGAGAGCGGTCAGTAGCAATCCAAACTTAGTTACTCGTCAGAAAATCGGATGGC",
      "GTCCAAGATCCCACAAAGTCCAGCGTACCATAAACGCAAGCCTCAACGCAGCGACGAGCACGAGAGCGGTCAGTAGCAATCCAAACTTAGTTACTCGTCA")
    val read5 = new MappingRead(5, 1, "CAACGCAGCGACGAGCACGAGAGCGGTCAGTAGCAATCCAAACTTAGTTACTCGTCAGAAAATCGGATGGCTGAATGAGTACCCAATGATCAGTACCAGT",
      "ACGTTTGACTAGTACGGTACGTAGCATGGGTCCAAGATCCCACAAAGTCCAGCGTACCATAAACGCAAGCCTCAACGCAGCGACGAGCACGAGAGCGGTC")
    val ca = new ConsensusAlignment(read1)
    var r: (Int, Array[Int], Array[Int]) = null
    r = ca.align(read2, new SuffixTree(read2.seq1), new SuffixTree(read2.seq2))
    ca.joinAndUpdate(read2, r._2, r._3)
    r = ca.align(read3, new SuffixTree(read3.seq1), new SuffixTree(read3.seq2))
    ca.joinAndUpdate(read3, r._2, r._3)
    r = ca.align(read4, new SuffixTree(read4.seq1), new SuffixTree(read4.seq2))
    ca.joinAndUpdate(read4, r._2, r._3)
    r = ca.align(read5, new SuffixTree(read5.seq1), new SuffixTree(read5.seq2))
    ca.joinAndUpdate(read5, r._2, r._3)
    println(ca.printPileup())
    val result = ArrayBuffer[List[Long]]()
    ca.reportAllEdgesTuple("GAGCGGTCAGTAGC",result)
    println(result)
  }
}

class ConsensusAlignment(read: MappingRead) extends ArrayBuffer[MappingRead]() {
  var consensus1 = new ConsensusSequence(read.seq1)
  var consensus2 = new ConsensusSequence(read.seq2)
  this += read
  read.column1 = consensus1.columnID.toArray
  read.column2 = consensus2.columnID.toArray

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
      }
    }
    if (last + 1 != seq.length) {
      consensus.insert(inserted + mapping.length, seq.substring(last + 1).toCharArray)
      for (i <- last + 1 until seq.length) readColumn(i) = consensus.columnID(inserted + mapping.length + i - last - 1)
      inserted += mapping(0)
    }
  }

  def joinAndUpdate(read: MappingRead, mapping1: Array[Int], mapping2: Array[Int]) {
    this += read
    updateConsensus(consensus1, mapping1, read.seq1, read.column1)
    updateConsensus(consensus2, mapping2, read.seq2, read.column2)
  }

  def align(read: MappingRead, readST1: SuffixTree, readST2: SuffixTree): (Int, Array[Int], Array[Int]) = {
    val (count1, span1, mapping1) = ConsensusAlignment.alignOneEnd(consensus1.sequence(), read.seq1, readST1)
    if (mapping1 == null || count1 < span1 * Settings.MATCH_RATE) return null
    val (count2, span2, mapping2) = ConsensusAlignment.alignOneEnd(consensus2.sequence(), read.seq2, readST2)
    if (mapping2 == null || count2 < span2 * Settings.MATCH_RATE) return null
    (2 * count1 - span1 + 2 * count2 - span2, mapping1, mapping2)
  }

  def reportAllEdgesTuple(kmer: String, report: ArrayBuffer[List[Long]]): Unit = {
    val seq1KmerSet = new Array[mutable.Set[Int]](this.size)
    val minCommonKmer = Array.ofDim[Boolean](this.size, this.size)
    val thishash = kmer.hashCode
    for (i <- this.indices) {
      val seq1 = this (i).seq1
      seq1KmerSet(i) = mutable.Set[Int]()
      for (j <- 0 to seq1.length - ConsensusAlignment.K) {
        val kmerhash = seq1.substring(j, j + ConsensusAlignment.K).hashCode
        if (kmerhash <= thishash) seq1KmerSet(i) += kmerhash
      }
      for (j <- 0 until i)
        minCommonKmer(i)(j) = (seq1KmerSet(i) intersect seq1KmerSet(j)).min == thishash
    }
    val table = mutable.Map[Char, Int]('A' -> 0, 'C' -> 1, 'G' -> 2, 'T' -> 3, 'N' -> 4, '-' -> 5)
    val columns1 = Array.fill[ArrayBuffer[(Int, Long)]](consensus1.columnID.size)(ArrayBuffer[(Int, Long)]())
    val columns2 = Array.fill[ArrayBuffer[(Int, Long)]](consensus2.columnID.size)(ArrayBuffer[(Int, Long)]())
    for (idx <- this.indices) {
      val read = this (idx)
      var refidx = 0
      var first = true
      for (i <- read.column1.indices) {
        var gapCount = 0
        while (consensus1.columnID(refidx) != read.column1(i)) {
          if (!first) {
            // 48-bit identifier(base offset in input file) | 8 bit alignment gap offset | 3 bit base code
            val baseCode = ((((read.id + i) << 8) + gapCount) << 3) + 5
            columns1(refidx) += ((idx, baseCode * (if (read.reversed) -1 else 1)))
          }
          refidx += 1
          gapCount += 1
        }
        val baseCode = ((read.id + i) << 11) + table(read.seq1(i))
        columns1(refidx) += ((idx, baseCode * (if (read.reversed) -1 else 1)))
        refidx += 1
        first = false
      }
      refidx = 0
      first = true
      for (i <- read.column2.indices) {
        var gapCount = 0
        while (consensus2.columnID(refidx) != read.column2(i)) {
          if (!first) {
            val baseCode = ((((read.id + i) << 8) + gapCount) << 3) + 5
            columns2(refidx) += ((idx, baseCode * (if (read.reversed) 1 else -1)))
          }
          refidx += 1
          gapCount += 1
        }
        val baseCode = ((read.id + i) << 11) + table(read.seq2(i))
        columns2(refidx) += ((idx, baseCode * (if (read.reversed) 1 else -1)))
        refidx += 1
        first = false
      }
    }
    for (column <- columns1) {
      var edge = List[Long]()
      var prev = List[Int]()
      for (j <- column) {
        var allCommonMin = true
        for (k <- prev) allCommonMin &= minCommonKmer(j._1)(k)
        if (allCommonMin) edge ::= j._2
        prev ::= j._1
      }
      if (edge.size>1) report += edge
    }
    for (column <- columns2) {
      var edge = List[Long]()
      var prev = List[Int]()
      for (j <- column) {
        var allCommonMin = true
        for (k <- prev) allCommonMin &= minCommonKmer(j._1)(k)
        if (allCommonMin) edge ::= j._2
        prev ::= j._1
      }
      if (edge.size>1) report += edge
    }
  }

  def printPileup(): String = {
    var pileup = f"${"Consensus"}%20s " + consensus1.sequence(true) + "|" + consensus2.sequence(true) + "\n"
    for (read <- this) {
      pileup += f"${read.id}%20s "
      var refidx = 0
      for (i <- read.column1.indices) {
        while (consensus1.columnID(refidx) != read.column1(i)) {
          pileup += " "
          refidx += 1
        }
        pileup += read.seq1(i)
        refidx += 1
      }
      pileup += " " * (consensus1.columnID.length - refidx) + "|"
      refidx = 0
      for (i <- read.column2.indices) {
        while (consensus2.columnID(refidx) != read.column2(i)) {
          pileup += " "
          refidx += 1
        }
        pileup += read.seq2(i)
        refidx += 1
      }
      pileup += "\n"
    }
    pileup + "-----\n"
  }
}