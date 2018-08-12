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
import org.apache.spark.util.{DoubleAccumulator, LongAccumulator}

object Dec {
  final val OO = Settings.OO
  val K = Settings.K

  def decomposeKmerWithoutRead(totalBase:LongAccumulator,totalKmer:LongAccumulator,totalQual:DoubleAccumulator)
                   (read: (Int, Long, String)): TraversableOnce[(String,Int)] = {
    // return: (kmer,(fileno,offset,kmerPos,Complement))
    val (_, _, record) = read
    val lines = record.split("\n")
    val sequence = lines(1)
    val quality = lines(3)
    totalBase.add(sequence.length)
    totalKmer.add(sequence.length-K+1)
    totalQual.add(quality.map(x=>math.pow(10,(33-x)/10.0)).sum)
    for (t <- K to sequence.length() if !sequence.substring(t - K, t).contains('N')) yield
      if ("AC" contains sequence(t - (K + 1) / 2)) (sequence.substring(t - K, t),1)
      else (Util.reverseComplement(sequence.substring(t - K, t)),1)
  }
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

  def alignByAnchorST(kmerCov:Int)(keyValues: (String, Iterable[(Int, Long, Int, Boolean, String)])): TraversableOnce[List[Long]] = {
    val compressTable = Array('A', 'G', 'C', 'N', 'T')
    val alignmentGroup = ArrayBuffer[ConsensusAlignment]()
    val Values = keyValues._2.toArray.sortBy(x => (-x._3, x._2))
    val report = ArrayBuffer[List[Long]]()
    //    val t1 = System.currentTimeMillis()
    val reads = ArrayBuffer[MappingRead]()
    for ((fileno, offset, pos, compl, cseq) <- Values) {
      val seq = cseq.map(x => compressTable(x % 5))
      val qual = cseq.map(_ / 5)
      val read = new MappingRead(keyValues._1.hashCode, fileno, offset, pos, seq, qual, compl)
      reads.append(read)
    }
    if ((1 until reads.length).map(i=>(reads(i).kmerset intersect reads(i-1).kmerset).isEmpty).indexOf(true)== -1) return report
    for (read <- reads){
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
    alignmentGroup.foreach(_.reportAllEdgesTuple(keyValues._1, kmerCov, report))
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

  def judgeColBases(kmerCov:Int)(KVs: (Long, Iterable[Long])): TraversableOnce[(Int, Long, Int, Char, Int)] = {
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
    val ACGT = collection.mutable.Map[Char, Long]('A' -> 0, 'C' -> 0, 'G' -> 0, 'T' -> 0, '-' -> 0, 'N' -> 0)
    for (value <- values if "ACGT-" contains value._4) ACGT(value._4) += value._5
    val (ref, refcnt) = ACGT.maxBy(_._2)
    val Error = Util.errorTest(values.map(x=>(x._4,x._5.toInt)),kmerCov*2,1)
    for (value <- values if Error(value._4)) yield
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
      })
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
    val totalBase = sc.longAccumulator("Base")
    val totalKmer = sc.longAccumulator("Kmer")
    val totalQual = sc.doubleAccumulator("Qual")
    val readsFile = Dec.readFastqFiles(sc)
    val P1 = readsFile.flatMap(decomposeKmer)
    val parallel = Settings.prime.find(_>P1.getNumPartitions*5).getOrElse(Settings.prime.last)
//    P1.persist(StorageLevel.DISK_ONLY)
//    P1.cache()
    val KmerCount = readsFile.flatMap(decomposeKmerWithoutRead(totalBase,totalKmer,totalQual))
                    .reduceByKey(_+_,parallel).map(x => (x._2,1)).reduceByKey(_ + _).collect().sorted
    val trough = (1 until KmerCount.length-1).find(i=>KmerCount(i-1)._2>=KmerCount(i)._2 && KmerCount(i)._2<=KmerCount(i+1)._2).get
    val peak = for (i <- trough until KmerCount.length-1 if KmerCount(i-1)._2<=KmerCount(i)._2 && KmerCount(i)._2>=KmerCount(i+1)._2) yield KmerCount(i)
    val kmerCov = peak.maxBy(_._2)._1
    val avgErrRate = totalQual.value/totalBase.value
    val cov = kmerCov/math.pow(1-avgErrRate,K)*totalBase.value/totalKmer.value
    println(s"[K-mer Spectrum] Trough         = $trough")
    println(s"[K-mer Spectrum] K-mer Coverage = $kmerCov")
    println(s"[Dataset]        Avg Err Rate   = $avgErrRate")
    println(s"[Dataset]        Base  Coverage = $cov")
    println(s"[Dataset]        Genome Size    = ${totalBase.value/cov}")
    val P2 = P1.groupByKey(parallel).filter(x => (x._2.size > 1) && (x._2.size < kmerCov*20)).flatMap(Dec.alignByAnchorST(kmerCov))
    val id_clique = ConnectedComponent.runByNodeList(sc, P2, P1)
    val P3 = id_clique.map(kv => if (kv._2 < 0) (-kv._2, -kv._1) else (kv._2, kv._1)).groupByKey()
//    if (Settings.debugPrint) P3.foreachPartition(printGenomeColumns)
    val P4 = P3.flatMap(judgeColBases(kmerCov))
    P4.saveAsTextFile("file://" + odir)
    javaRuntime.exec(f"cat $odir/part* > $odir/result")
  }

  def main(args: Array[String]) {
    new Arg(args)
    val conf = new SparkConf().setAppName("decos").registerKryoClasses(Array[Class[_]](
        Class.forName("org.apache.hadoop.io.LongWritable"),Class.forName("org.apache.hadoop.io.Text")))
    val sc = new SparkContext(conf)
    sc.setLogLevel("WARN")
    sc.hadoopConfiguration.setInt("mapred.max.split.size", 4 * 1024 * 1024) //mapreduce.input.fileinputformat.split.maxsize
    runST(sc)
    println("\n[DECOS] finished successfully!\n")
  }
}