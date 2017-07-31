import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.Set
import java.io.RandomAccessFile

import scala.collection.mutable
import org.apache.spark._
import org.apache.spark.graphx._

/**
  * Created by workshop on 04-Jul-17.
  */
object anchoredRead {
  final val OO = Settings.OO
  val MATCH_RATE = Settings.MATCH_RATE
  val K = Settings.K

  def matchOneEnd(seq1: String, seq2: String, offset: Int, heuristic: Boolean): Int = {
    // return match-mismatch. return -OO if similarity is low
    // offset==3
    // seq1    012345
    // seq2 012345
    // seq1[i]==seq2[i+offset]
    if (offset < 0) return matchOneEnd(seq2, seq1, -offset, heuristic)
    val overlap = Math.min(seq1.length, seq2.length - offset)
    var matchCount = 0
    for (i <- 0 until overlap) {
      if (seq1.charAt(i) == seq2.charAt(i + offset)) matchCount += 1
      if (heuristic && i == 4 && matchCount < 4) return -OO
    }
    if (matchCount.toDouble / overlap < MATCH_RATE || matchCount < K) return -OO
    matchCount - (overlap - matchCount)
  }
}

class anchoredRead(idx: Long, pos: Int, file1:RandomAccessFile,file2:RandomAccessFile) {
  final val OO = Settings.OO
  val id:Long = idx
  var reversed = false
  var kmeridx:Int = pos
  file1.seek(id)
  file2.seek(id)
  var seq1 = file1.readLine()
  var seq2 = file2.readLine()
  var offset1,offset2 = 0
  if (kmeridx < 0) {
    reversed = !reversed
    kmeridx = -kmeridx
    val tmp = seq1
    seq1 = seq2
    seq2 = tmp
  }
  def matchRead(read:anchoredRead,offset:Int):(Double,Int,Int) = {
    // return best pair-alignment with score,x,y
    val x = this.kmeridx - read.kmeridx
    val anchoredScore = anchoredRead.matchOneEnd(read.seq1,this.seq1,x,false)
    if (anchoredScore == -OO) return null
    var bestMatch:(Double,Int,Int) = null
    val iterateRange = if (offset== -OO) 1-read.seq2.length until this.seq2.length else List(offset)
    for ( y <- iterateRange){
      var result = anchoredRead.matchOneEnd(read.seq2,this.seq2,y,true)
      if (result!= -OO){
        result += anchoredScore
        if (bestMatch==null || bestMatch._1 < result)
          bestMatch = (result,x+this.offset1,y+this.offset2)
      }
    }
    bestMatch
  }
}

class anchorAlignment(read:anchoredRead) extends ArrayBuffer[anchoredRead](){
  val K = Settings.K
  def anchor(read:anchoredRead, leftOffset:Int, rightOffset:Int){
    this += read
    read.offset1 = leftOffset
    read.offset2 = rightOffset
  }
  this.anchor(read, 0, 0)
  def matchRead(read:anchoredRead,offset:Int= -Settings.OO) = this.last.matchRead(read,offset)
  def reportDoubt(kmer:String,colID:Int, report:ArrayBuffer[(String,(Char,Long))]):Int = {
    var colIDv = colID
    var seq1Left,seq2Left,seq1Right,seq2Right = 0
    for (i <- this){
      seq1Left = Math.min(seq1Left,i.offset1)
      seq1Right = Math.max(seq1Right, i.offset1 + i.seq1.length)
      seq2Left = Math.min(seq2Left, i.offset2)
      seq2Right = Math.max(seq2Right, i.offset2 + i.seq2.length)
    }
    for (colIdx <- seq1Left until seq1Right){
      val bases = for ( i <- this if 0 <= colIdx - i.offset1 && colIdx - i.offset1 < i.seq1.length)
        yield i.seq1.charAt(colIdx - i.offset1)
      if (bases.max!=bases.min){
        for ( i <- this if 0 <= colIdx - i.offset1 && colIdx - i.offset1 < i.seq1.length){
          report += ((kmer+colIDv,(i.seq1.charAt(colIdx-i.offset1),(i.id+colIdx-i.offset1)*(if (i.reversed) -1 else 1))))
        }
        colIDv += 1
      }
    }
    for (colIdx <- seq2Left until seq2Right){
      val bases = for ( i <- this if 0 <= colIdx - i.offset2 && colIdx - i.offset2 < i.seq2.length)
        yield i.seq2.charAt(colIdx - i.offset2)
      if (bases.max!=bases.min){
        for ( i <- this if 0 <= colIdx - i.offset2 && colIdx - i.offset2 < i.seq2.length){
          report += ((kmer+colIDv,(i.seq2.charAt(colIdx-i.offset2),(i.id+colIdx-i.offset2)*(if (i.reversed) 1 else -1))))
        }
        colIDv += 1
      }
    }
    colIDv
  }
  def reportAll(kmer:String,colID:Int, report:ArrayBuffer[(String,(Char,Long))]):Int = {
    var colIDv = colID
    var seq1Left,seq2Left,seq1Right,seq2Right = 0
    for (i <- this){
      seq1Left = Math.min(seq1Left,i.offset1)
      seq1Right = Math.max(seq1Right, i.offset1 + i.seq1.length)
      seq2Left = Math.min(seq2Left, i.offset2)
      seq2Right = Math.max(seq2Right, i.offset2 + i.seq2.length)
    }
    val seq1KmerSet = new Array[mutable.Set[String]](this.size)
    val seq2KmerSet = new Array[mutable.Set[String]](this.size)
    for ( i <- this.indices){
      val seq1 = this(i).seq1
      seq1KmerSet(i) = mutable.Set[String]()
      for ( j <- 0 to seq1.length-K){ seq1KmerSet(i) += seq1.substring(j,j+K)}
      val seq2 = this(i).seq2
      seq2KmerSet(i) = mutable.Set[String]()
      for ( j <- 0 to seq2.length-K){ seq2KmerSet(i) += seq2.substring(j,j+K)}
    }
    for (colIdx <- seq1Left until seq1Right){
      val kmerSet = mutable.Set[String]()
      var firstRead = -1
      var reportedSth = false
      for ( i <- this.indices if 0 <= colIdx - this(i).offset1 && colIdx - this(i).offset1 < this(i).seq1.length){
        if (kmerSet.isEmpty){
          if (firstRead== -1) firstRead = i
        } else if (kmerSet.intersect(seq1KmerSet(i)).min==kmer) {
          val read = this (i)
          reportedSth = true
          report += ((kmer + colIDv, (read.seq1.charAt(colIdx - read.offset1), (read.id + colIdx - read.offset1) * (if (read.reversed) -1 else 1))))
        }
        kmerSet ++= seq1KmerSet(i)
      }
      if (reportedSth) {
        report += ((kmer + colIDv, (this (firstRead).seq1.charAt(colIdx - this (firstRead).offset1),
                                  (this (firstRead).id + colIdx - this (firstRead).offset1) * (if (this (firstRead).reversed) -1 else 1))))
      }
      colIDv += 1
    }
    for (colIdx <- seq2Left until seq2Right){
      val kmerSet = mutable.Set[String]()
      var firstRead = -1
      var reportedSth = false
      for ( i <- this.indices if 0 <= colIdx - this(i).offset2 && colIdx - this(i).offset2 < this(i).seq2.length){
        if (kmerSet.isEmpty){
          if (firstRead== -1) firstRead = i
        } else {
          val intersect = kmerSet.intersect(seq2KmerSet(i))
          if (intersect.isEmpty || intersect.min==kmer) {
            val read = this (i)
            reportedSth = true
            report += ((kmer + colIDv, (read.seq2.charAt(colIdx - read.offset2), (read.id + colIdx - read.offset2) * (if (read.reversed) -1 else 1))))
          }
        }
        kmerSet ++= seq2KmerSet(i)
      }
      if (reportedSth) {
        report += ((kmer + colIDv, (this (firstRead).seq2.charAt(colIdx - this (firstRead).offset2),
                                  (this (firstRead).id + colIdx - this (firstRead).offset2) * (if (this (firstRead).reversed) -1 else 1))))
      }
      colIDv += 1
    }
    colIDv
  }
  def reportEdgesTuple(kmer:String, colID:Int, report:ArrayBuffer[List[Long]]):Int = {
    var colIDv = colID
    var seq1Left,seq2Left,seq1Right,seq2Right = 0
    for (i <- this){
      seq1Left = Math.min(seq1Left,i.offset1)
      seq1Right = Math.max(seq1Right, i.offset1 + i.seq1.length)
      seq2Left = Math.min(seq2Left, i.offset2)
      seq2Right = Math.max(seq2Right, i.offset2 + i.seq2.length)
    }
    for (colIdx <- seq1Left until seq1Right){
      val bases = for ( i <- this if 0 <= colIdx - i.offset1 && colIdx - i.offset1 < i.seq1.length)
        yield i.seq1.charAt(colIdx - i.offset1)
      if (bases.max!=bases.min){
        var first_last:(Long,Long) = null
        for ( i <- this if 0 <= colIdx - i.offset1 && colIdx - i.offset1 < i.seq1.length){
          var vid = (i.id+colIdx-i.offset1)*(if (i.reversed) -1 else 1)*5
          if (i.seq1.charAt(colIdx-i.offset1)=='A') vid += 0
          if (i.seq1.charAt(colIdx-i.offset1)=='C') vid += 1
          if (i.seq1.charAt(colIdx-i.offset1)=='G') vid += 2
          if (i.seq1.charAt(colIdx-i.offset1)=='T') vid += 3
          if (i.seq1.charAt(colIdx-i.offset1)=='N') vid += 4
          if (first_last == null) first_last = (vid,vid) else {
            report += List(vid,first_last._2)
            first_last = (first_last._1,vid)
          }
        }
        if (first_last!=null && first_last._1!=first_last._2)
          report += List(first_last._2,first_last._1)
        colIDv += 1
      }
    }
    for (colIdx <- seq2Left until seq2Right){
      val bases = for ( i <- this if 0 <= colIdx - i.offset2 && colIdx - i.offset2 < i.seq2.length)
        yield i.seq2.charAt(colIdx - i.offset2)
      if (bases.max!=bases.min){
        var first_last:(Long,Long) = null
        for ( i <- this if 0 <= colIdx - i.offset2 && colIdx - i.offset2 < i.seq2.length){
          var vid = (i.id+colIdx-i.offset2)*(if (i.reversed) 1 else -1)*5
          if (i.seq2.charAt(colIdx-i.offset2)=='A') vid += 0
          if (i.seq2.charAt(colIdx-i.offset2)=='C') vid += 1
          if (i.seq2.charAt(colIdx-i.offset2)=='G') vid += 2
          if (i.seq2.charAt(colIdx-i.offset2)=='T') vid += 3
          if (i.seq2.charAt(colIdx-i.offset2)=='N') vid += 4
          if (first_last == null) first_last = (vid,vid) else {
            report += List(vid,first_last._2)
            first_last = (first_last._1,vid)
          }
        }
        if (first_last!=null && first_last._1!=first_last._2)
          report += List(first_last._2,first_last._1)
        colIDv += 1
      }
    }
    colIDv
  }
  def reportAllEdgesTuple(kmer:String, colID:Int, report:ArrayBuffer[List[Long]]):Int = {
    var seq1Left,seq2Left,seq1Right,seq2Right = 0
    for (i <- this){
      seq1Left = Math.min(seq1Left,i.offset1)
      seq1Right = Math.max(seq1Right, i.offset1 + i.seq1.length)
      seq2Left = Math.min(seq2Left, i.offset2)
      seq2Right = Math.max(seq2Right, i.offset2 + i.seq2.length)
    }
    val seq1KmerSet = new Array[mutable.Set[Int]](this.size)
//    val seq2KmerSet = new Array[mutable.Set[Int]](this.size)
    val minCommonKmer = Array.ofDim[Boolean](this.size,this.size)
//    val hasCommonKmer = Array.ofDim[Boolean](this.size,this.size)
    val thishash = kmer.hashCode
    for ( i <- this.indices){
      val seq1 = this(i).seq1
      seq1KmerSet(i) = mutable.Set[Int]()
      for ( j <- 0 to seq1.length-K){
        val kmerhash = seq1.substring(j,j+K).hashCode
        if (kmerhash<=thishash) seq1KmerSet(i) += kmerhash
      }
      for ( j <- 0 until i )
        minCommonKmer(i)(j) = (seq1KmerSet(i) intersect seq1KmerSet(j)).min == thishash
//      val seq2 = this(i).seq2
//      seq2KmerSet(i) = mutable.Set[Int]()
//      for ( j <- 0 to seq2.length-K){
//        seq2KmerSet(i) += seq2.substring(j,j+K).hashCode
//      }
//      for ( j <- 0 until i )
//        hasCommonKmer(i)(j) = !(seq2KmerSet(i) intersect seq2KmerSet(j)).isEmpty
    }
    val table = Map[Char,Int]('A'->0,'C'->1,'G'->2,'T'->3, 'N'->4)
    for (colIdx <- seq1Left until seq1Right){
      var clique = List[Long]()
      var readInColumn = List[Int]()
      for ( (i,idx) <- this.zipWithIndex if 0 <= colIdx - i.offset1 && colIdx - i.offset1 < i.seq1.length){
        var vid = (i.id+colIdx-i.offset1)*(if (i.reversed) -1 else 1)*5+table(i.seq1.charAt(colIdx-i.offset1))
        if (clique.isEmpty || !(for ( prev <- readInColumn ) yield minCommonKmer(idx)(prev)).contains(false))
          clique ::= vid
        readInColumn ::= idx
      }
      if (clique.size>1) report += clique
    }
    for (colIdx <- seq2Left until seq2Right){
      var clique = List[Long]()
      var readInColumn = List[Int]()
      for ( (i,idx) <- this.zipWithIndex if 0 <= colIdx - i.offset2 && colIdx - i.offset2 < i.seq2.length){
        var vid = (i.id+colIdx-i.offset2)*(if (i.reversed) 1 else -1)*5+table(i.seq2.charAt(colIdx-i.offset2))
        if (clique.isEmpty || !(for ( prev <- readInColumn ) yield minCommonKmer(idx)(prev)).contains(false))
          clique ::= vid
        readInColumn ::= idx
      }
      if (clique.size>1) report += clique
    }
    colID
  }
  def printPileup():String = {
    var pileup = ""
    var seq1Left,seq2Left,seq1Right,seq2Right = 0
    for (i <- this){
      seq1Left = Math.min(seq1Left,i.offset1)
      seq1Right = Math.max(seq1Right, i.offset1 + i.seq1.length)
      seq2Left = Math.min(seq2Left, i.offset2)
      seq2Right = Math.max(seq2Right, i.offset2 + i.seq2.length)
    }
    for ( i <- this){
      pileup += f"${i.id}%10d" + "." * (i.offset1 - seq1Left)
      pileup += i.seq1 + "." * (seq1Right - i.offset1 - i.seq1.length)
      pileup += "..." + "." * (i.offset2 - seq2Left) + i.seq2
      pileup += "." * (seq2Right - i.offset2 - i.seq2.length) +"\n"
    }
    pileup+ "-----\n"
  }
}