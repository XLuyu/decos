import scala.collection.mutable.ArrayBuffer
import scala.io.Source

import scala.collection.mutable
import org.apache.spark._

/**
  * Created by workshop on 04-Jul-17.
  */

class MappingRead(fileNumber: Int, offset: Long, pos: Int, rawSeq: String, quality: IndexedSeq[Int], complement: Boolean = false) {
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

  def addCount(idx: Int, c: Char) = vote(idx)(c) += 1

  def sequence(keepGap: Boolean = false): String = {
    val consensusWithGap = for (i <- vote.indices) yield charAt(i)
    if (keepGap) consensusWithGap.mkString else consensusWithGap.filter(_ != '-').mkString
  }

  def insert(idx: Int, newChars: String): Unit = {
    columnID.insertAll(idx, columnID.size until (columnID.size + newChars.length))
    val newColumn = for (c <- newChars) yield mutable.Map('A' -> 0, 'C' -> 0, 'G' -> 0, 'T' -> 0, 'N'-> 0 ,'-' -> 0) + (c -> 1)
    vote.insertAll(idx, newColumn)
  }

  def checkGap(idx: Int, newChars: String): Unit = {
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
  val NWaligner = new NeedlemanWunschAligner(Settings.MAX_READ_LEN * 2, Settings.MAX_INDEL)
  def main(args: Array[String]): Unit = {
    val readArray = ArrayBuffer[MappingRead]()
    var n = 0
    var ca: ConsensusAlignment = null
    var r: (Int, Array[Int]) = null
    for (line <- Source.fromFile("Demo.txt").getLines) {
      n += 1
      readArray += new MappingRead(0, n, 1, line, line.map(_/1), true)
      if (n == 1) {
        ca = new ConsensusAlignment(readArray(0))
      } else {
        r = ca.align(readArray.last, new SuffixTree(readArray.last.seq))
        ca.joinAndUpdate(readArray.last, r._2)
        println(ca.printPileup())
      }
    }
    val result = ArrayBuffer[List[Long]]()
    ca.reportAllEdgesTuple("CCTCTGCTGGCGCTGGTTGCC",Int.MaxValue, result)
    println(ca.printPileup())
    println(result)
  }
}

class ConsensusAlignment(read: MappingRead) extends ArrayBuffer[MappingRead]() {
  var consensus = new ConsensusSequence(read.seq)
  this += read
  read.column = consensus.columnID.toArray
  val ufs = new UnionFindSet[Int]()
  var father:Array[Int] = _

  def alignOneEnd(groupSeq: String, read: MappingRead, readST: SuffixTree): (Int, Int, Array[Int]) = {
    // return mapping array like: -2 -2 0 1 2 -1 3 4 6 -3 -3 -3
    def sharp(a: (Int, Int, Int), b: (Int, Int, Int)): (Int, Int, Int) = { // remove overlap part from a
      val (a1, a2, a3, a4) = (a._1, a._1 + a._3, a._2, a._2 + a._3)
      val (b1, b2, b3, b4) = (b._1, b._1 + b._3, b._2, b._2 + b._3)
      if (Math.abs(a1 - b1 - a3 + b3) > Settings.MAX_INDEL) return null
      if (a1 < b1 && a3 < b3) {
        val ovlp = Math.max(Math.max(a2 - b1, a4 - b3), 0)
        return (a1, a3, a._3 - ovlp)
      }
      if (b2 < a2 && b4 < a4) {
        val ovlp = Math.max(Math.max(b2 - a1, b4 - a3), 0)
        return (a1 + ovlp, a3 + ovlp, a._3 - ovlp)
      }
      null
    }
    val readSeq = read.seq
    def compactLength(a: (Int, Int, Int)): Int = {
      var i = a._2
      var j = a._2 + a._3 - 1
      while (i < a._2 + a._3 - 1 && groupSeq(i) == groupSeq(i + 1)) i += 1
      while (j > a._2 && groupSeq(j) == groupSeq(j - 1)) j -= 1
      Math.max(j - i + 1, 1)
    }

    val mapping = Array.fill[Int](groupSeq.length)(-1)
    //  val lengthThreshold = Math.log(readSeq.length) / Math.log(2)
    //    println(readST.pairwiseLCS(groupSeq))
    var LCS = readST.pairwiseLCS(groupSeq).filter(_._3 > Settings.K).filter(compactLength(_) > 2).sortBy(-_._3)
    // LCS = List[(read, group, length)]
    //        println(groupSeq)
    //        println(LCS)

    val offset = mutable.Map[Int,Int]().withDefaultValue(0)
    LCS.foreach( x => offset(x._1-x._2) += x._3 )
    if (offset.isEmpty) return (0, 0, null)
    val maxOffset = offset.keys.max
    val minOffset = offset.keys.min
    var core = (maxOffset,offset(maxOffset))
    for ( i <- minOffset until maxOffset) {
      var sum = 0
      for ( j <- 0 to Settings.MAX_INDEL) sum += offset(i-j)
      if (sum>core._2) core = (i,sum)
    }
    LCS = LCS.filter(x => 0 <= x._1-x._2-core._1 && x._1-x._2-core._1 <= Settings.MAX_INDEL).sortBy( x => -offset(x._1-x._2)*1000+x._2-x._1)
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
    val HeadMapping = ConsensusAlignment.NWaligner.align(groupSeq.substring(lastG - unmatchedGroupHead, lastG).reverse,
      readSeq.substring(lastR - unmatchedReadHead, lastR).reverse, read.qual.slice(lastR - unmatchedReadHead, lastR).reverse, tailFlex = true)
    for (i <- HeadMapping.indices)
      mapping(lastG - 1 - i) = if (HeadMapping(i) == -1) -1 else lastR - 1 - HeadMapping(i)
    for (i <- lastG - unmatchedGroupHead until lastG if mapping(i) != -1 && groupSeq(i) == readSeq(mapping(i))) matchCount += 1
    spanCount += Math.min(unmatchedReadHead, unmatchedGroupHead)
    for ((rs, gs, len) <- matchSeg) {
      if (gs - lastG == rs - lastR) {
        for (i <- 0 until (gs - lastG)) mapping(lastG + i) = lastR + i
      } else {
        val regionMapping = ConsensusAlignment.NWaligner.align(groupSeq.substring(lastG, gs), readSeq.substring(lastR, rs), read.qual.slice(lastR, rs))
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
    val tailMapping = ConsensusAlignment.NWaligner.align(groupSeq.substring(lastG, lastG + unmatchedGroupTail),
      readSeq.substring(lastR, lastR + unmatchedReadTail), read.qual.slice(lastR, lastR + unmatchedReadTail), tailFlex = true)
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
  def updateConsensus(consensus: ConsensusSequence, mapping: Array[Int], seq: String, readColumn: Array[Int]) {
    var inserted = 0
    if (mapping(0) > 0) {
      consensus.insert(0, seq.substring(0, mapping(0)))
      for (i <- 0 until mapping(0)) readColumn(i) = consensus.columnID(i)
      inserted += mapping(0)
    }
    var last = -1
    for (i <- mapping.indices) {
      if (mapping(i) > -1) {
        if (last != -1 && last + 1 != mapping(i)) {
          consensus.checkGap(inserted + i, seq.substring(last + 1, mapping(i)))
          for (j <- last + 1 until mapping(i)) readColumn(j) = consensus.columnID(inserted + i + j - last - 1)
          inserted += mapping(i) - last - 1
        }
        while (consensus.charAt(inserted + i) == '-') {
          if (last != -1) consensus.addCount(inserted + i, '-')
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
      consensus.insert(inserted + mapping.length, seq.substring(last + 1))
      for (i <- last + 1 until seq.length) readColumn(i) = consensus.columnID(inserted + mapping.length + i - last - 1)
    }
  }

  def joinAndUpdate(read: MappingRead, mapping: Array[Int]) {
    this += read
    updateConsensus(consensus, mapping, read.seq, read.column)
  }

  def align(read: MappingRead, readST: SuffixTree): (Int, Array[Int]) = {
    val consensusSeq = consensus.sequence()
    val (_, _, mapping) = alignOneEnd(consensusSeq, read, readST)
    var last = -1
    var tot = 0
    var pm = 0
    if (mapping != null){
      for (i <- mapping.indices) {
        if (mapping(i) > -1) {
          if (last != -1 && last + 1 != mapping(i)) {
            tot += read.qual.slice(last+1,mapping(i)).sum
          }
          if (consensusSeq(i)==read.seq(mapping(i))) pm += read.qual(mapping(i))
          tot += read.qual(mapping(i))
          last = mapping(i)
        } else if (mapping(i) == -1) {
          tot += (read.qual(last)+read.qual(last+1))/2
        }
      }
    }
    if (mapping == null || pm < tot * 0.90) return null
    //    if (mapping == null || count < span * Settings.MATCH_RATE) return null
    //shifting indel alignment to right-most place(or left-most for complemented read)
    var i = 1
    var leftBound = Math.max(mapping.indexWhere(_ >= 0),1)
    var rightBound = Math.min(mapping.lastIndexWhere(_ >= 0),mapping.length-2)
    while ( i < mapping.length-1) {
      if (mapping(i) == -1) {
        val left = (i-1 to leftBound by -1).find(consensusSeq(_)!=consensusSeq(i)).getOrElse(leftBound-1)+1
        val right = (i+1 to rightBound).find(consensusSeq(_)!=consensusSeq(i)).getOrElse(rightBound+1)-1
        val (start,end,step) = if (consensusSeq(i)=='A' || consensusSeq(i)=='C') (left,right,1) else (right,left,-1)
        var nextEmptyPos = start
        for ( j <- start to end by step if mapping(j) != -1){
          mapping(nextEmptyPos) = mapping(j)
          nextEmptyPos += step
        }
        for ( j <- nextEmptyPos to end by step) mapping(j) = -1
        i = right
      }
      i += 1
    }
    val readMapping = Array.fill[Int](read.seq.length)(-1)
    for ( i <- mapping.indices if mapping(i)>=0) readMapping(mapping(i)) = i
    leftBound = 0
    while ( readMapping(leftBound) == -1) {readMapping(leftBound) = -2; leftBound += 1}
    if (leftBound==0) leftBound = 1
    rightBound = readMapping.length-1
    while ( readMapping(rightBound) == -1) {readMapping(rightBound) = -3; rightBound -= 1}
    if (rightBound==readMapping.length-1) rightBound = readMapping.length-2
    i = 1
    while ( i < readMapping.length-1){
      if (readMapping(i) == -1) {
        val left = (i-1 to leftBound by -1).find(read.seq(_)!=read.seq(i)).getOrElse(leftBound-1)+1
        val right = (i+1 to rightBound).find(read.seq(_)!=read.seq(i)).getOrElse(rightBound+1)-1
        val (start,end,step) = if (read.seq(i)=='A' || read.seq(i)=='C') (left,right,1) else (right,left,-1)
        var nextEmptyPos = start
        for ( j <- start to end by step if readMapping(j) != -1){
          mapping(readMapping(j)) = nextEmptyPos
          readMapping(nextEmptyPos) = readMapping(j)
          nextEmptyPos += step
        }
        i = right
      }
      i += 1
    }
    (2 * pm - tot, mapping)
  }

  def reportAllEdgesTuple(kmer: String, kmerCov:Int, report: ArrayBuffer[List[Long]]): Unit = {
    def baseEncode(read:MappingRead, position: Long, gapCount: Int, baseCode: Int, quality: Int): Long = {
      // tried bitwise operation, but it gathers hashCode of records and skew the partition
      //      if (read.complemented)
      val id = (gapCount * 64 + quality) * 1000000000000L + position
      (id * 6 + baseCode) * 2 + read.fileno
    }
    class FingerPrint(s:StringBuilder, q:ArrayBuffer[Int]){
      var seq = s.clone()
      var qual = q.clone()
      def dissimilarity(other:FingerPrint) = (for (i <- s.indices if seq(i)!=' ' && other.seq(i)!=' ' && seq(i)!=other.seq(i)) yield Math.min(qual(i),other.qual(i))).sum
      def similarity(other:FingerPrint) = (for (i <- s.indices if seq(i)!=' ' && seq(i)==other.seq(i)) yield Math.min(qual(i),other.qual(i))).sum
      def simAndDissim(other:FingerPrint) = (similarity(other),dissimilarity(other))
      def += (other:FingerPrint): Unit ={
        for ( i <- s.indices ) {
          if (seq(i)==' ') seq(i) = other.seq(i)
          qual(i) += other.qual(i)
        }
      }
    }
    //    val t1 = System.currentTimeMillis()
    //compute pair-wise k-mer intersection
    val KmerSet = new Array[mutable.Set[Int]](this.size)
//    val isAnchoredByMinKmer = Array.ofDim[Boolean](this.size, this.size)
    val isAnchoredByMinKmer = new Array[mutable.Map[Int,Boolean]](this.size)
    val thishash = kmer.hashCode
    for (i <- this.indices) {
      val seq = this (i).seq
      KmerSet(i) = mutable.Set[Int]()
      for (j <- 0 to seq.length - Settings.K) {
        val kmerhash = Util.minKmerByRC(seq.substring(j, j + Settings.K)).hashCode
        if (kmerhash < thishash) KmerSet(i) += kmerhash
      }
      isAnchoredByMinKmer(i) = mutable.Map[Int,Boolean]()
//      for (j <- 0 until i)
//        isAnchoredByMinKmer(i)(j) = (KmerSet(i) intersect KmerSet(j)).isEmpty //TODO: compute by use
    }
    //    val t2 = System.currentTimeMillis()
    // compute column linked list
    val table = mutable.Map[Char, Int]('A' -> 0, 'C' -> 1, 'G' -> 2, 'T' -> 3, 'N' -> 4, '-' -> 5)
    val columns = Array.fill[ArrayBuffer[(Int, Long)]](consensus.columnID.size)(ArrayBuffer[(Int, Long)]())
    val columnCount = Array.fill[mutable.Map[Char, Long]](consensus.columnID.size)(mutable.Map[Char, Long]('A' -> 0, 'C' -> 0, 'G' -> 0, 'T' -> 0, 'N' -> 0, '-' -> 0))
    for (idx <- this.indices) {
      val read = this (idx)
      var refidx = 0
      for (i <- read.column.indices) {
        var gapCount = consensus.columnID.indexWhere(_ == read.column(i),refidx) - refidx
        if (i>0) {
          for (g <- 1 to gapCount) {
            val baseCode = if (read.complemented)
              baseEncode(read, read.id + read.seq.length - i - 1, gapCount + 1 - g, table('-'), (read.qual(i - 1) + read.qual(i)) / 2)
            else
              baseEncode(read, read.id + i - 1, g, table('-'), (read.qual(i - 1) + read.qual(i)) / 2)
            columns(refidx + g - 1) += ((idx, baseCode))
            columnCount(refidx + g - 1)('-') += (read.qual(i - 1) + read.qual(i)) / 2
          }
        }
        refidx += gapCount
        val baseCode = if (read.complemented)
          baseEncode(read, read.id + read.seq.length - 1 - i, 0, table(Util.complement(read.seq(i))), read.qual(i))
        else
          baseEncode(read, read.id + i, 0, table(read.seq(i)), read.qual(i))
        columns(refidx) += ((idx, baseCode))
        columnCount(refidx)(read.seq(i)) += read.qual(i)
        refidx += 1
      }
    }
    //    val t3 = System.currentTimeMillis()
    val rtable = Array('A', 'C', 'G', 'T', 'N', '-')
    val fingerprint = Array.fill[StringBuilder](this.size)(new StringBuilder())
    val fpQuality = Array.fill[ArrayBuffer[Int]](this.size)(new ArrayBuffer[Int]())
    for (i <- columns.indices) {
      val baseList = columns(i).map(x=>{
        var base = rtable((x._2/2%6).toInt)
        if (this (x._1).complemented) base = Util.complement(base)
        (base,(x._2/12000000000000L%64).toInt)
      })
      val errorList = Util.errorTest(baseList,kmerCov)
      if (errorList.count(_._2==false)>1){
        var last = 0
        for ((idx,baseCode) <- columns(i)) {
          var base = rtable((baseCode/2%6).toInt)
          if (this (idx).complemented) base = Util.complement(base)
          fingerprint(idx).+=(if (!errorList(base)) base else ' ')
          fpQuality(idx).+=((baseCode/12000000000000L%64).toInt)
          while (last!=idx) {fingerprint(last) += ' '; fpQuality(last) += 40 ; last += 1}
          last += 1
        }
        while (last<this.size) {fingerprint(last) += ' '; fpQuality(last) += 40 ; last += 1}
      }
    }
    val fplen = fingerprint.head.length
    father = new Array[Int](this.size)
    if (fplen==0){
      for ( i <- father.indices) father(i) = 0
    } else {
      val fp2idx = mutable.Map[StringBuilder, Int]()
      val fpGroup = ArrayBuffer[FingerPrint]()
      for (i <- this.indices){
        val fp = new FingerPrint(fingerprint(i),fpQuality(i))
        if (fp2idx.contains(fp.seq))
          fpGroup(fp2idx(fp.seq)) += fp
        else{
          fp2idx(fp.seq.clone()) = fpGroup.size
          fpGroup.append(fp)
        }
      }
      val fpCluster = new Array[Int](fpGroup.size)
      val solid = fpGroup.map(_.qual.min>Settings.CUTTHRESHOLD)
      for ( i <- fpGroup.indices if solid(i)){
        var best = (-1,-1)
        for ( j <- 0 until i if solid(j) ) {
          val (sim,dissim) = fpGroup(j).simAndDissim(fpGroup(i))
          if (dissim==0 && (best._1 == -1 || best._2 < sim)) best = (j,sim)
        }
        if (best._1 == -1) fpCluster(i) = i else {
          fpGroup(best._1) += fpGroup(i)
          fpCluster(i) = best._1
        }
      }
      val solidlist = solid.indices.filter(solid)
      for ( i <- fpGroup.indices if !solid(i)){
        fpCluster(i) = if (solidlist.isEmpty) i else {
          val simlist = solidlist.filter(x=>fpGroup(i).dissimilarity(fpGroup(x))<=Settings.CUTTHRESHOLD)
          if (simlist.isEmpty) i else simlist.maxBy(x=>fpGroup(i).similarity(fpGroup(x)))
        }
      }
      for (idx <- this.indices){
        father(idx) = fpCluster(fp2idx(fingerprint(idx)))
      }
    }
    //    val t4 = System.currentTimeMillis()
    for (column <- columns) {
      val node = mutable.Map[Int,List[Long]]()
      val prev = mutable.Map[Int,List[Int]]()
      for ((idx,baseCode) <- column) {
        val s = father(idx)
        if (!node.contains(s)) { node(s) = List[Long](); prev(s) = List[Int]()}
        if (prev(s).isEmpty)
          node(s) ::= baseCode * (if (this (idx).complemented) -1 else 1)
        else {
          val k = prev(s).head
          if (!isAnchoredByMinKmer(idx).contains(k))
            isAnchoredByMinKmer(idx)(k) = (KmerSet(idx) intersect KmerSet(k)).isEmpty
          if (isAnchoredByMinKmer(idx)(k))
            node(s) ::= baseCode * (if (this (idx).complemented) -1 else 1)
        }
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
      pileup += f"${father(idx)}%5s${read.id}%12s(${read.fileno}${if (read.complemented) "-" else "+"})"
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