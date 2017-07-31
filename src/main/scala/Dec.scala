import org.apache.spark.SparkContext
import org.apache.spark.SparkConf
import org.apache.spark.HashPartitioner
import java.io.{FileWriter, RandomAccessFile}

import scala.collection.mutable.ArrayBuffer
import org.apache.spark._
import org.apache.spark.graphx._
import com.kwartile.lib.cc.ConnectedComponent

object Dec {
  final val OO = Settings.OO
  val K = Settings.K
  val MATCH_RATE = Settings.MATCH_RATE

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
  
  def alignByAnchor(keyValues: (String, Iterable[String])): TraversableOnce[(String, (Char,Long))] = {
    val readFile1 = new RandomAccessFile(Settings.refFilePath._1, "r")
    val readFile2 = new RandomAccessFile(Settings.refFilePath._2, "r")
//    val Ast = new SuffixTree()
//    val Ust = new SuffixTree()
    val alignmentGroup = ArrayBuffer[anchorAlignment]()
    val Values:Array[String] = keyValues._2.toArray.sortBy(-_.split("[+-]")(1).toInt)
    for (value <- Values){
      val tokens = value.split(raw"(?=[-+])")
      val read = new anchoredRead(tokens(0).toLong,tokens(1).toInt,readFile1,readFile2)
//      val readST = new SuffixTree
//      readST.append(read.seq2)
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
//        val (offset,matchLen) = readST.tidMaxVote(alignmentGroup(i).last.seq2)
//        val result = if (matchLen>=K) alignmentGroup(i).matchRead(read,-offset) else alignmentGroup(i).matchRead(read)
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
    val report = ArrayBuffer[(String,(Char,Long))]()
    for (alignment <- alignmentGroup ){
      col = alignment.reportDoubt(keyValues._1, col, report)
    }
    if (Settings.printAlignment){
      val logFile = new FileWriter("/home/x/xieluyu/output/log" + keyValues._1)
      for (alignment <- alignmentGroup ){
        logFile.write(alignment.printPileup())
      }
      logFile.close()
    }
    readFile1.close()
    readFile2.close()
    report
  }
  def alignByAnchorEdgeTuple(keyValues: (String, Iterable[String])): TraversableOnce[List[Long]] = {
    val readFile1 = new RandomAccessFile(Settings.refFilePath._1, "r")
    val readFile2 = new RandomAccessFile(Settings.refFilePath._2, "r")
    //    val Ast = new SuffixTree()
    //    val Ust = new SuffixTree()
    val alignmentGroup = ArrayBuffer[anchorAlignment]()
    val Values:Array[String] = keyValues._2.toArray.sortBy(-_.split("[+-]")(1).toInt)
    for (value <- Values){
      val tokens = value.split(raw"(?=[-+])")
      val read = new anchoredRead(tokens(0).toLong,tokens(1).toInt,readFile1,readFile2)
      //      val readST = new SuffixTree
      //      readST.append(read.seq2)
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
        //        val (offset,matchLen) = readST.tidMaxVote(alignmentGroup(i).last.seq2)
        //        val result = if (matchLen>=K) alignmentGroup(i).matchRead(read,-offset) else alignmentGroup(i).matchRead(read)
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
    val report = ArrayBuffer[List[Long]]()
    for (alignment <- alignmentGroup ){
      col = alignment.reportAllEdgesTuple(keyValues._1, col, report)
    }
    if (Settings.printAlignment){
      val logFile = new FileWriter("/home/x/xieluyu/output/log" + keyValues._1)
      for (alignment <- alignmentGroup ){
        logFile.write(alignment.printPileup())
      }
      logFile.close()
    }
    readFile1.close()
    readFile2.close()
    report
  }
  def mergeWithinPartition(KVS:Iterator[(String, (Char,Long))]):Iterator[(String, (Char,Long))] = {
    val father = collection.mutable.Map[String,String]()
    def find(k:String):String = if (father(k)==k) k else {father(k) = find(father(k)); father(k)}
    val basedict = collection.mutable.Map[(Char,Long),String]()
    for ((v,k) <- KVS){
      if (!father.contains(v)) father(v) = v
      if (!basedict.contains(k))
        basedict(k) = v
      else
        father(find(v)) = find(basedict(k))
    }
    for ((k,v) <- basedict.iterator) yield (find(v),k)
  }
  def judgeColBases(KVs:(String, Iterable[(Char,Long)])): TraversableOnce[String] = {
    val values = KVs._2
    val ACGT =  collection.mutable.Map[Char,Long]('A'->0,'C'->0,'G'->0,'T'->0)
    for (value <- values if value._1!='N') ACGT(value._1) += 1
    val (ref,refcnt) = ACGT.maxBy(_._2)
    if (refcnt <= Settings.CUTTHRESHOLD) return List()
    for (value <- values if value._1=='N' || ACGT(value._1)<=Settings.CUTTHRESHOLD)
      yield value._2 + "\t"+ ref
  }
  def judgeColBasesK(KVs:(Long, Iterable[Long])): TraversableOnce[String] = {
    val table = Array('A','C','G','T','N')
    val values = KVs._2.map(x => (table((x%5+5).toInt%5),(x-(if (x>=0) 0 else 4))/5))
    val ACGT =  collection.mutable.Map[Char,Long]('A'->0,'C'->0,'G'->0,'T'->0)
    for (value <- values if value._1!='N') ACGT(value._1) += 1
    val (ref,refcnt) = ACGT.maxBy(_._2)
    if (refcnt <= Settings.CUTTHRESHOLD) return List()
    for (value <- values if value._1=='N' || ACGT(value._1)<=Settings.CUTTHRESHOLD)
      yield value._2 + "\t"+ ref
  }
  def judgeColBasesM(KVs:(Long, Iterable[Long])): TraversableOnce[String] = {
    val table = Array('A','C','G','T','N')
    val values = (KVs._1 :: KVs._2.toList).map(x => (table((x%5+5).toInt%5),(x-(if (x>=0) 0 else 4))/5))
    val ACGT =  collection.mutable.Map[Char,Long]('A'->0,'C'->0,'G'->0,'T'->0)
    for (value <- values if value._1!='N') ACGT(value._1) += 1
    val (ref,refcnt) = ACGT.maxBy(_._2)
    if (refcnt <= Settings.CUTTHRESHOLD) return List()
    for (value <- values if value._1=='N' || ACGT(value._1)<=Settings.CUTTHRESHOLD)
      yield value._2 + "\t"+ ref
  }
  def runK(sc: SparkContext, ifile: String = "/home/x/xieluyu/reads/lineno_seq1_seq2.txt", odir: String = "/home/x/xieluyu/output") {
    val javaRuntime = Runtime.getRuntime
    javaRuntime.exec("rm -r "+odir)
    val readsFile = sc.textFile(ifile,64)
    val P1 = readsFile.flatMap(Dec.scanKmer)
    val P2 = P1.groupByKey(1024).filter(_._2.size >1).flatMap(Dec.alignByAnchorEdgeTuple)
    val (id_clique,converged,iterCount) = ConnectedComponent.run(sc,P2,Int.MaxValue)
    val P3 = id_clique.map(kv=>(kv._2,kv._1)).groupByKey()
    val P4 = P3.flatMap(judgeColBasesK)
    P4.saveAsTextFile("file://" + odir)
  }
  def runM(sc: SparkContext, ifile: String = "/home/x/xieluyu/reads/lineno_seq1_seq2.txt", odir: String = "/home/x/xieluyu/output") {
    val javaRuntime = Runtime.getRuntime
    javaRuntime.exec("rm -r "+odir)
    val readsFile = sc.textFile(ifile,64)
    val P1 = readsFile.flatMap(Dec.scanKmer)
    val P2 = P1.groupByKey(1024).filter(_._2.size >1).flatMap(Dec.alignByAnchorEdgeTuple)
    val (id_clique,converged,iterCount) = ConnectedComponent.runG(sc,P2,Int.MaxValue)
    val P3 = id_clique.map(kv=>(kv._2,kv._1)).groupByKey()
    val P4 = P3.flatMap(judgeColBasesM)
    P4.saveAsTextFile("file://" + odir)
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
    runM(sc)
  }
}