import java.io.{FileWriter, RandomAccessFile}
import scala.collection.mutable.ArrayBuffer

import org.apache.spark.SparkContext
import org.apache.spark.SparkConf
import org.apache.spark.HashPartitioner
import org.apache.spark._
import org.apache.spark.graphx._

import ConnectedComponent._

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

  def alignByAnchor(keyValues: (String, Iterable[String])): TraversableOnce[List[Long]] = {
    val readFile1 = new RandomAccessFile(Settings.refFilePath._1, "r")
    val readFile2 = new RandomAccessFile(Settings.refFilePath._2, "r")
    val alignmentGroup = ArrayBuffer[ConsensusAlignment]()
    val Values:Array[String] = keyValues._2.toArray.sortBy(-_.split("[+-]")(1).toInt)
    var cnt = 0
    for (value <- Values){
      println("------- read",cnt)
      cnt += 1
      val tokens = value.split(raw"(?=[-+])")
      val read = MappingRead(tokens(0).toLong,tokens(1).toInt,readFile1,readFile2)
      val readST1 = new SuffixTree(read.seq1)
      val readST2 = new SuffixTree(read.seq2)
      var best = -1
      var bestAlign:(Int, Array[Int], Array[Int]) = null
      println("--- read seq1",read.seq1)
      for (i <- alignmentGroup.indices){
        println(alignmentGroup(i).consensus1.sequence())
        val align = alignmentGroup(i).align(read,readST1,readST2)
//        if (align!=null) {
//          println(readST1.pairwiseLCS(alignmentGroup(i).consensus1.sequence(), K))
//          for ( i <- align._2) {
//            printf("%d,", i)
//          }
//          println("----")
//          for ( i <- align._3) {
//            printf("%d,", i)
//          }
//          println("----")
//        }
        if (align!=null && (bestAlign==null || align._1 > bestAlign._1)){
          best = i
          bestAlign = align
        }
      }
      if (bestAlign==null){
        alignmentGroup += new ConsensusAlignment(read)
//        println("new AG")
      } else {
//        println("best",best)
//        println(readST1.pairwiseLCS(alignmentGroup(best).consensus1.sequence(), K))
        alignmentGroup(best).joinAndUpdate(read, bestAlign._2, bestAlign._3)
      }
    }
    val report = ArrayBuffer[List[Long]]()
    for (alignment <- alignmentGroup ){
      alignment.reportAllEdgesTuple(keyValues._1, report)
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
  def judgeColBasesM(KVs:(Long, Iterable[Long])): TraversableOnce[String] = {
    val table = Array('A','C','G','T','N','-')
    val values = (KVs._1 :: KVs._2.toList).map(x => (table((x%5+5).toInt%5),(x-(if (x>=0) 0 else 4))/5))
    val ACGT = collection.mutable.Map[Char,Long]('A'->0,'C'->0,'G'->0,'T'->0)
    for (value <- values if value._1!='N') ACGT(value._1) += 1
    val (ref,refcnt) = ACGT.maxBy(_._2)
    if (refcnt <= Settings.CUTTHRESHOLD) return List()
    for (value <- values if value._1=='N' || ACGT(value._1)<=Settings.CUTTHRESHOLD)
      yield value._2 + "\t"+ ref
  }
  def judgeColBases(KVs:(Long, Iterable[Long])): TraversableOnce[String] = {
    val table = Array('A','C','G','T','N','-')
    val values = (KVs._1 :: KVs._2.toList).filter(x=>(Math.abs(x)>>3&0xff)==0).filter(x=>(Math.abs(x)&0x7)<5).map(
      x => (table(Math.abs(x).toInt&7),(Math.abs(x)>>11)*(if (x>=0) 1 else -1)))
    val ACGT = collection.mutable.Map[Char,Long]('A'->0,'C'->0,'G'->0,'T'->0)
    for (value <- values if value._1!='N') ACGT(value._1) += 1
    val (ref,refcnt) = ACGT.maxBy(_._2)
    if (refcnt <= Settings.CUTTHRESHOLD) return List()
    for (value <- values if value._1=='N' || ACGT(value._1)<=Settings.CUTTHRESHOLD)
      yield value._2 + "\t"+ ref
  }
  def runST(sc: SparkContext, ifile: String = "/home/x/xieluyu/reads/lineno_seq1_seq2.txt", odir: String = "/home/x/xieluyu/output") {
    val javaRuntime = Runtime.getRuntime
    javaRuntime.exec("rm -r "+odir)
    val readsFile = sc.textFile(ifile,64)
    val P1 = readsFile.flatMap(Dec.scanKmer)
    val P2 = P1.groupByKey(1024).filter(_._2.size >1).flatMap(Dec.alignByAnchor)
    val id_clique = ConnectedComponent.runByNodeList(sc,P2)
    val P3 = id_clique.map(kv=>(kv._2,kv._1)).groupByKey()
    val P4 = P3.flatMap(judgeColBases)
    P4.saveAsTextFile("file://" + odir)
  }
  def run(sc: SparkContext, ifile: String = "/home/x/xieluyu/reads/lineno_seq1_seq2.txt", odir: String = "/home/x/xieluyu/output") {
    val javaRuntime = Runtime.getRuntime
    javaRuntime.exec("rm -r "+odir)
    val readsFile = sc.textFile(ifile,64)
    val P1 = readsFile.flatMap(Dec.scanKmer)
    val P2 = P1.groupByKey(1024).filter(_._2.size >1).flatMap(Dec.alignByAnchor)
    val id_clique = ConnectedComponent.runByNodeList(sc,P2)
    val P3 = id_clique.map(kv=>(kv._2,kv._1)).groupByKey()
    val P4 = P3.flatMap(judgeColBasesM)
    P4.saveAsTextFile("file://" + odir)
  }
  def main(args: Array[String]) {
    val conf = new SparkConf().setAppName("decos")
    val sc = new SparkContext(conf)
    runST(sc)
  }
}