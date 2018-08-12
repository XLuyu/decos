import org.rogach.scallop._

class Arg(arguments: Seq[String]) extends ScallopConf(arguments) {
  val inputFilePath = trailArg[List[String]](required = false, default = Some(List("/home/x/xieluyu/reads/Ecoli1.fq", "/home/x/xieluyu/reads/Ecoli2.fq")))
  val K = opt[Int](default=Some(21))
  val MAX_READ_LEN = opt[Int](default=Some(150))
  val MAX_INDEL = opt[Int](default=Some(2))
  val CUTTHRESHOLD = opt[Int](default=Some(60))
  val MATCH_RATE = opt[Double](default=Some(0.9))
  val debugPrint = opt[Boolean]()
  verify()
  Settings.inputFilePath = inputFilePath().toArray
  Settings.K = K()
  Settings.MAX_READ_LEN = MAX_READ_LEN()
  Settings.MAX_INDEL = MAX_INDEL()
  Settings.CUTTHRESHOLD = CUTTHRESHOLD()
  Settings.MATCH_RATE = MATCH_RATE()
  Settings.debugPrint = debugPrint()
}

object Settings {
  val prime = Array(2003,4001,8009,16001,30011,60013)
  val OO = 2000000000
  var debugPrint = true
  var MAX_READ_LEN = 150
  var MAX_INDEL = 2
  var CUTTHRESHOLD = 60
  var K = 21
  var MATCH_RATE = 0.9
  var inputFilePath = Array("/home/x/xieluyu/reads/Ecoli1.fq", "/home/x/xieluyu/reads/Ecoli2.fq")
//  val inputFilePath = Array("/research/wongls-group/xieluyu/benchmark/simulated/fv/D03.fq")
//  val inputFilePath = Array("/research/wongls-group/xieluyu/benchmark/simulated/Ecoli/Ecoli.fq")
//    val inputFilePath = Array("/research/wongls-group/xieluyu/benchmark/simulated/Ecoli100/Ecoli100.fq")
//  val inputFilePath = Array("/research/wongls-group/xieluyu/benchmark/illumina2014/Hiseq/L.pneumophila_SRR801797.fastq")
//  val inputFilePath = Array("/research/wongls-group/xieluyu/benchmark/illumina2014/Hiseq/E.coli_SRR490124.fastq")
//    val inputFilePath = Array("/research/wongls-group/xieluyu/benchmark/illumina2014/Hiseq/S.cerevisiae_ERR422544.fastq")
//  val inputFilePath = Array("/research/wongls-group/xieluyu/benchmark/illumina2014/Hiseq/Human_ERX069715.fastq")
}