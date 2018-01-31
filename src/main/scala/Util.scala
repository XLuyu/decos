import scala.collection.mutable

/**
  * Created by workshop on 19-Sep-17.
  */
object Util {
  val complement = Map('A' -> 'T', 'T' -> 'A', 'C' -> 'G', 'G' -> 'C', 'N' -> 'N', '-' -> '-')
  val Q2E = (0 to 50).map(x=>math.pow(10,-x/10.0)).toArray

  def reverseComplement(seq: String): String = {
    for (i <- seq.reverse) yield complement(i)
  }

  def minKmerByRC(kmer: String): String = {
    // for kmer with odd length, check the middle base
    // if it's 'G'/'T', return its reverse complement, otherwise return itself
    if ("GT".contains(kmer(kmer.length / 2))) reverseComplement(kmer) else kmer
  }

  //  def minHashInWindow(seq:String, window:Int = 10, K:Int = Settings.K) = {
  //    val kmerlist = for ( i <- K to seq.length) yield minKmerByRC(seq.substring(i-K,i)).hashCode
  //    (for ( i <- window to kmerlist.size) yield kmerlist.slice(i-window,i).min).toSet
  //    kmerlist.sorted.slice(0,10)
  //  }
  def binomialCDF(n: Int, p: Double): Double = {
    var cdf, cur = math.pow(p, n)
    for (i <- 1 to n) {
      cur = cur * (n + 1 - i) / i * (1 - p) / p
      cdf += cur
      if (cdf > 0.9999) return i - (cdf - 0.9999) / cur
    }
    0
  }

  def errorTestWithMargin(bases: Iterable[(Char, Int)], coverageUpperbound: Int = Int.MaxValue): Map[Char,Boolean] = {
    val count = bases.groupBy(_._1).mapValues(_.size)
    val qualSum = bases.groupBy(_._1).mapValues(x => if (x.head._1=='N') 0 else x.map(_._2).sum)
    val ref = qualSum.maxBy(_._2)._1
    val e = bases.map(x=>Q2E(x._2)).sum/bases.size
    val p = 1-e
    val q = -10*math.log10(e)
    count.map(x => (x._1, x._1!=ref && qualSum(x._1)<=q+Util.binomialCDF(x._2+math.min(count(ref),coverageUpperbound),p)*q))
  }
  def errorTest(bases: Iterable[(Char, Int)], coverageUpperbound: Int = Int.MaxValue): Map[Char,Boolean] = {
    val count = bases.groupBy(_._1).mapValues(_.size)
    val qualSum = bases.groupBy(_._1).mapValues(x => if (x.head._1=='N') 0 else x.map(_._2).sum)
    val ref = qualSum.maxBy(_._2)._1
    val e = bases.map(x=>Q2E(x._2)).sum/bases.size
    val p = 1-e
    val q = -10*math.log10(e)
    count.map(x => (x._1, x._1!=ref && qualSum(x._1)<=Util.binomialCDF(x._2+math.min(count(ref),coverageUpperbound),p)*q))
  }
}
