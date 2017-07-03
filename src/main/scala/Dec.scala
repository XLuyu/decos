import org.apache.spark.SparkContext
import org.apache.spark.SparkConf
import org.apache.spark.rdd._
import org.apache.spark.storage.StorageLevel._
import org.apache.spark.HashPartitioner
import java.io.RandomAccessFile
import java.util.ArrayList
import java.io.FileWriter
import scala.collection.mutable.ArrayBuffer

object Dec {
  final val OO = 2000000000
  final val CUTTHRESHOLD = 2
  val K = 14
  val MATCH_RATE = 0.9
//  val refFilePath = ("/temp/refread/Phix1.fq", "/temp/refread/Phix2.fq")
  val refFilePath = ("/temp/refread/Ecoli1.fq", "/temp/refread/Ecoli2.fq")
  
  def scanKmer(read: String): TraversableOnce[(String, String)] = {
    val id_read1_read2 = read.split("\t")
    val id = id_read1_read2(0)
    val read1 = id_read1_read2(1)
    val read2 = id_read1_read2(2)
		// s is starting index of kmer (1-based), t is ending index (0-based)
		// 1-base enable negative position to represent second end of read
    (
    for ( t <- K to read1.length() if !read1.substring(t-K,t).contains('N'))
      yield (read1.substring(t-K,t),id + "+" + (t-K+1))
    )++(
    for ( t <- K to read1.length() if !read2.substring(t-K,t).contains('N'))
      yield (read2.substring(t-K,t),id + "-" + (t-K+1))
    )
  }
  
  def alignByAnchor(keyValues: (String, Iterable[String])): TraversableOnce[(String, (Char,Int))] = {
    val readFile1 = new RandomAccessFile(refFilePath._1, "r")
    val readFile2 = new RandomAccessFile(refFilePath._2, "r")
    
    def matchOneEnd(seq1: String, seq2: String, offset: Int, heuristic: Boolean): Int = {
      // return match-mismatch. return -OO if similarity is low
      // offset==3
      // seq1    012345
      // seq2 012345
      // seq1[i]==seq2[i+offset]
      if (offset < 0) return matchOneEnd(seq2, seq1, -offset, heuristic)
      val overlap = Math.min(seq1.length(), seq2.length() - offset)
      var matchCount = 0
      for (i <- 0 until overlap) {
        if (seq1.charAt(i) == seq2.charAt(i + offset)) matchCount += 1
        if (heuristic && i == 4 && matchCount < 4) return -OO
      }
      if (matchCount.toDouble / overlap < MATCH_RATE || matchCount < K) return -OO
      matchCount - (overlap - matchCount)
    }
    
    class anchoredRead(idx: Int, pos: Int) {
      val id:Int = idx
      var reversed = false
      var kmeridx:Int = pos
      readFile1.seek(id)
      readFile2.seek(id)
      var seq1 = readFile1.readLine()
      var seq2 = readFile2.readLine()
      var offset1,offset2 = 0
      if (kmeridx < 0) {
        reversed = !reversed
        kmeridx = -kmeridx
        val tmp = seq1
        seq1 = seq2
        seq2 = tmp
      }
      def matchRead(read:anchoredRead):(Double,Int,Int) = {
        // return best pair-alignment with score,x,y
        val x = this.kmeridx - read.kmeridx
        val anchoredScore = matchOneEnd(read.seq1,this.seq1,x,false)
        if (anchoredScore == -OO) return null
        var bestMatch:(Double,Int,Int) = null
        for ( y <- 1-read.seq2.length until this.seq2.length){
          var result = matchOneEnd(read.seq2,this.seq2,y,true)
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
      def anchor(read:anchoredRead, leftOffset:Int, rightOffset:Int){
        this += read
        read.offset1 = leftOffset
        read.offset2 = rightOffset
      }
      this.anchor(read, 0, 0)
      def matchRead(read:anchoredRead) = this.last.matchRead(read)
      def reportDoubt(kmer:String,colID:Int, report:ArrayBuffer[(String,(Char,Int))]):Int = {
        var colIDv = colID
        var seq1Left,seq2Left,seq1Right,seq2Right = 0
        for (i <- this){
				  seq1Left = Math.min(seq1Left,i.offset1)
				  seq1Right = Math.max(seq1Right, i.offset1 + i.seq1.length())
				  seq2Left = Math.min(seq2Left, i.offset2)
				  seq2Right = Math.max(seq2Right, i.offset2 + i.seq2.length())
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
      def printPileup():String = {
        var pileup = ""
        var seq1Left,seq2Left,seq1Right,seq2Right = 0
        for (i <- this){
				  seq1Left = Math.min(seq1Left,i.offset1)
				  seq1Right = Math.max(seq1Right, i.offset1 + i.seq1.length())
				  seq2Left = Math.min(seq2Left, i.offset2)
				  seq2Right = Math.max(seq2Right, i.offset2 + i.seq2.length())
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
    val Ast = new SuffixTree()
    val Ust = new SuffixTree()
    val alignmentGroup = ArrayBuffer[anchorAlignment]()
    var Values:Array[String] = keyValues._2.toArray.sortBy(-_.split("[+-]")(1).toInt)
    for (value <- Values){
      val tokens = value.split(raw"(?=[-+])")
      val read = new anchoredRead(tokens(0).toInt,tokens(1).toInt)
      var best = -1.0 // max score
      var bestidx = (-1,0,0) // best group to match, left offset, right offset
//      var Avote = Ast.tidVote(read.seq1)
//      var Uvote = Ust.tidVote(read.seq2)
      for (i <- alignmentGroup.indices){
//        val vote = Avote(i)._2+Uvote(i)._2
//        if (vote>best){
//          best = vote
//          bestidx = (i,Avote(i)._1,Uvote(i)._1)
//        }
        val result = alignmentGroup(i).matchRead(read)
        if (result!=null && (best < 0 || best < result._1)){
          best = result._1
          bestidx = (i,result._2,result._3)
        }
      }
//      if (best>=0 && (matchOneEnd(read.seq1,alignmentGroup(bestidx._1)(0).seq1,bestidx._2,false)== -OO ||
//                     matchOneEnd(read.seq2,alignmentGroup(bestidx._1)(0).seq2,bestidx._3,false)== -OO)){
//        best = -1
//      }
      if (best<0){
//        Ast.append(read.seq1)
//        Ust.append(read.seq2)
        alignmentGroup += new anchorAlignment(read)
      } else
        alignmentGroup(bestidx._1).anchor(read, bestidx._2, bestidx._3)
    }
    var col = 0
    val report = ArrayBuffer[(String,(Char,Int))]()
    for (alignment <- alignmentGroup ){
      col = alignment.reportDoubt(keyValues._1, col, report)
    }
//    var logFile = new FileWriter("/home/x/xieluyu/output/log" + keyValues._1);
//    for (alignment <- alignmentGroup ){
//      logFile.write(alignment.printPileup())
//    }
//    logFile.close()
    readFile1.close()
    readFile2.close()
    report
  }
  
  def mergeWithinPartition(KVS:Iterator[(String, (Char,Int))]):Iterator[(String, (Char,Int))] = {
    val father = collection.mutable.Map[String,String]()
    def find(k:String):String = if (father(k)==k) k else {father(k) = find(father(k)); father(k)}
    val basedict = collection.mutable.Map[(Char,Int),String]()
    for ((v,k) <- KVS){
      if (!father.contains(v)) father(v) = v
      if (!basedict.contains(k))
        basedict(k) = v
      else
        father(find(v)) = find(basedict(k))
    }
    for ((k,v) <- basedict.iterator) yield (find(v),k)
  }
  def judgeColBases(KVs:(String, Iterable[(Char,Int)])): TraversableOnce[String] = {
    val values = KVs._2
    val ACGT =  collection.mutable.Map[Char,Int]('A'->0,'C'->0,'G'->0,'T'->0)
    for (value <- values if value._1!='N') ACGT(value._1) += 1
    val (ref,refcnt) = ACGT.maxBy(_._2)
    if (refcnt <= CUTTHRESHOLD) return List()
    for (value <- values if ACGT(value._1)<=CUTTHRESHOLD && value._1!='N') //TODO: value._1=='N' || ACGT(value._1)<=CUTTHRESHOLD
      yield value._2 + "\t"+ ref
  }
  def run(sc: SparkContext, ifile: String = "/home/x/xieluyu/reads/lineno_seq1_seq2.txt", odir: String = "/home/x/xieluyu/output") {
    val javaRuntime = Runtime.getRuntime
    javaRuntime.exec("rm -r "+odir)
    val readsFile = sc.textFile(ifile,64)
    val P1 = readsFile.flatMap(Dec.scanKmer)
    val P2 = P1.groupByKey(1024).filter(_._2.size >1).flatMap(Dec.alignByAnchor)
    var partitionsNum = P2.getNumPartitions
    var P3 = P2.mapPartitions(Dec.mergeWithinPartition, true)
    while(partitionsNum!=1){
      partitionsNum = (partitionsNum - 1) / 4 + 1
      P3 = P3.partitionBy(new HashPartitioner(partitionsNum)).mapPartitions(Dec.mergeWithinPartition, true)
    }
    val P4 = P3.groupByKey().flatMap(Dec.judgeColBases)
    P4.saveAsTextFile("file://" + odir)
  }
  def main(args: Array[String]) {
    val conf = new SparkConf().setAppName("decos")
    val sc = new SparkContext(conf)
    run(sc)
  }
}