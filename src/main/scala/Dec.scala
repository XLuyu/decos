import java.io._

import scala.collection.mutable.ArrayBuffer
import org.apache.hadoop.io.{LongWritable, Text}
import org.apache.spark.SparkContext
import org.apache.spark.SparkConf
import org.apache.spark.HashPartitioner
import org.apache.spark.rdd._
import org.apache.spark._

import scala.io.Source
import ConnectedComponent._
import org.apache.spark.storage.StorageLevel

object Dec {
  final val OO = Settings.OO
  val K = Settings.K

  def decomposeKmer(read: (Int, Long, String)): TraversableOnce[(String, (Int, Long, Int, Boolean, String))] = {
    // return: (kmer,(fileno,offset,kmerPos,Complement))
    val (fileno, offset, record) = read
    val lines = record.split("\n")
    val head = lines(0).length + 1
    val sequence = lines(1)
    val quality = lines(3)
    val compressedSeq = (for (i <- sequence.indices) yield ((quality(i) - 33) * 5 + sequence(i) % 5).toChar).mkString
    for (t <- K to sequence.length() if !sequence.substring(t - K, t).contains('N')) yield
      if ("AC" contains sequence(t - (K + 1) / 2))
        (sequence.substring(t - K, t), (fileno, offset + head, t - K, false, compressedSeq))
      else
        (Util.reverseComplement(sequence.substring(t - K, t)), (fileno, offset + head, sequence.length() - t, true, compressedSeq))
  }

  def alignByAnchorST(keyValues: (String, Iterable[(Int, Long, Int, Boolean, String)])): TraversableOnce[List[Long]] = {
    //    val fileHandles = for ( filename <- Settings.inputFilePath) yield new RandomAccessFile(filename, "r")
    val compressTable = Array('A', 'G', 'C', 'N', 'T')
    val alignmentGroup = ArrayBuffer[ConsensusAlignment]()
    val Values = keyValues._2.toArray.sortBy(x => (-x._3, x._2))
    //    val t1 = System.currentTimeMillis()
    for ((fileno, offset, pos, compl, cseq) <- Values) {
      val seq = cseq.map(x => compressTable(x % 5))
      val qual = cseq.map(_ / 5)
      val read = new MappingRead(fileno, offset, pos, seq, qual, compl)
      val readST = new SuffixTree(read.seq)
      var best = -1
      var bestAlign: (Int, Array[Int]) = null
      for (i <- alignmentGroup.indices) {
        val align = alignmentGroup(i).align(read, readST)
        if (align != null && (bestAlign == null || align._1 > bestAlign._1)) {
          best = i
          bestAlign = align
        }
      }
      if (bestAlign == null) {
        alignmentGroup += new ConsensusAlignment(read)
      } else {
        alignmentGroup(best).joinAndUpdate(read, bestAlign._2)
      }
    }
    //    val t2 = System.currentTimeMillis()
    val report = ArrayBuffer[List[Long]]()
    alignmentGroup.foreach(_.reportAllEdgesTuple(keyValues._1, report))
    //    if (report.nonEmpty && Settings.debugPrint) {
    //      val logFile = new FileWriter("/home/x/xieluyu/log/" + keyValues._1)
    //      for (alignment <- alignmentGroup) {
    //        logFile.write(alignment.printPileup())
    //      }
    //      logFile.close()
    //    }
    //    val t3 = System.currentTimeMillis()
    //    println("alignment:",t2-t1)
    //    println("report:",t3-t2)
    //    fileHandles.foreach(_.close())
    report
  }

  def judgeColBases(KVs: (Long, Iterable[Long])): TraversableOnce[(Int, Long, Int, Char, Int)] = {
    val table = Array('A', 'C', 'G', 'T', 'N', '-')
    val values = (KVs._1 :: KVs._2.toList).map(x => {
      val complemented = if (x > 0) false else true
      val abs = Math.abs(x)
      val gap = abs / (12000000000000L * 64)
      val quality = abs / 12000000000000L % 64
      val pos = abs / 12 % 1000000000000L
      val base = if (complemented) Util.complement(table((abs / 2 % 6).toInt)) else table((abs / 2 % 6).toInt)
      val fileno = abs % 2
      (fileno, pos, gap, base, quality, complemented)
    })
    // .filter(x => x._3==0 && x._4!='-')
    val ACGT = collection.mutable.Map[Char, Long]('A' -> 0, 'C' -> 0, 'G' -> 0, 'T' -> 0, '-' -> 0)
    for (value <- values if "ACGT-" contains value._4) ACGT(value._4) += value._5
    //    for (value <- values if value._1!='N') println(KVs._1.toString+value.toString)
    val (ref, refcnt) = ACGT.maxBy(_._2)
    if (refcnt <= Settings.CUTTHRESHOLD) return List()
    for (value <- values if value._4 == 'N' || ACGT(value._4) <= Settings.CUTTHRESHOLD) yield
      // value._2*(-2*value._1+1) + "\t" + (if (value._6) Util.complement(ref) else ref)
      (value._1.toInt, value._2, value._3.toInt, if (value._6) Util.complement(ref) else ref, value._5.toInt)
  }

  def printGenomeColumns(data: Iterator[(Long, Iterable[Long])]): Unit = {
    val table = Array('A', 'C', 'G', 'T', 'N', '-')
    val logFile = new FileWriter("/home/x/xieluyu/log/" + Math.abs(data.hashCode))
    for ((k, v) <- data) {
      val values = (k :: v.toList).map(x => {
        val complemented = if (x > 0) false else true
        val abs = Math.abs(x)
        val gap = abs / (12000000000000L * 64)
        val quality = abs / 12000000000000L % 64
        val pos = abs / 12 % 1000000000000L
        val base = if (complemented) Util.complement(table((abs / 2 % 6).toInt)) else table((abs / 2 % 6).toInt)
        val fileno = abs % 2
        (fileno, pos, gap, base, quality, complemented)
      }) //.filter(x => x._3==0 && "ACGTN".contains(x._4))
      logFile.write(values.mkString(" ") + "\n")
    }
    logFile.close()
  }

  def readFastqFiles(sc: SparkContext) = {
    val filelist = Settings.inputFilePath
    val s = for (i <- filelist.indices) yield
      sc.newAPIHadoopFile(filelist(i), classOf[FastqInputFormat], classOf[LongWritable], classOf[Text]).map(x => (i, x._1.get, x._2.toString))
    s.reduce((a, b) => a ++ b)
  }

  def writeFastqFile(KVs: (Int, Iterable[(Int, Long, Int, Char, Int)])): Unit = {

    val bis = new BufferedInputStream(new FileInputStream(new File(Settings.inputFilePath(KVs._1))))
    val in = new BufferedReader(new InputStreamReader(bis), 10 * 1024 * 1024)
    while (true) {
      val name = in.readLine()
      var seq = in.readLine()
      val name2 = in.readLine()
      var qual = in.readLine()
    }
    val inputFile = Source.fromFile(Settings.inputFilePath(KVs._1))
    val outputFile = new FileWriter(Settings.inputFilePath(KVs._1) + ".corr")
    for (line <- inputFile.getLines())
      for ((_, offset, gap, base, quality) <- KVs._2.toArray.sortBy(_._2)) {

      }
  }

  def runST(sc: SparkContext, odir: String = "/home/x/xieluyu/output") {
    val javaRuntime = Runtime.getRuntime
    javaRuntime.exec("rm -r " + odir)
    val readsFile = Dec.readFastqFiles(sc)
    val P1 = readsFile.flatMap(Dec.decomposeKmer)
//    P1.persist(StorageLevel.DISK_ONLY)
//        P1.cache()
        val KmerCount = P1.mapValues(x=>1).reduceByKey(_+_).map(x => (x._2,1)).reduceByKey(_ + _).collect().sorted
        val peak = for ( i <- 1 until KmerCount.size-1 if KmerCount(i-1)._2<=KmerCount(i)._2 && KmerCount(i)._2>=KmerCount(i+1)._2) yield KmerCount(i)
        val cov = peak.maxBy(_._2)._1*10
        println(cov)
//    val cov = 1000000000
    val P2 = P1.groupByKey(10007).filter(x => (x._2.size > 1) && (x._2.size < cov)).flatMap(Dec.alignByAnchorST)
    val id_clique = ConnectedComponent.runByNodeList(sc, P2, P1)
    val P3 = id_clique.map(kv => if (kv._2 < 0) (-kv._2, -kv._1) else (kv._2, kv._1)).groupByKey()
    //    if (Settings.debugPrint) P3.foreachPartition(printGenomeColumns)
    val P4 = P3.flatMap(judgeColBases)
    P4.saveAsTextFile("file://" + odir)
    //    P4.groupBy(_._1).foreach(writeFastqFile)
  }

  def main(args: Array[String]) {
    val conf = new SparkConf().setAppName("decos")
      .registerKryoClasses(Array[Class[_]](
        Class.forName("org.apache.hadoop.io.LongWritable"),
        Class.forName("org.apache.hadoop.io.Text")))
    val sc = new SparkContext(conf)
    sc.setLogLevel("WARN")
    sc.hadoopConfiguration.setInt("mapred.max.split.size", 4 * 1024 * 1024) //mapreduce.input.fileinputformat.split.maxsize
    runST(sc)
    println("[ DECOS ] finished successfully!")
  }
}