import scala.collection.mutable
/**
  * Created by workshop on 19-Sep-17.
  */
object Util {
  val complement = Map('A'->'T','T'->'A','C'->'G','G'->'C','N'->'N','-'->'-')
  def reverseComplement(seq:String):String={
    for ( i <- seq.reverse ) yield complement(i)
  }
  def minKmerByRC(kmer:String):String={
    // for kmer with odd length, check the middle base
    // if it's 'G'/'T', return its reverse complement, otherwise return itself
    if ("GT".contains(kmer(kmer.length/2))) reverseComplement(kmer) else kmer
  }
//  def minHashInWindow(seq:String, window:Int = 10, K:Int = Settings.K) = {
//    val kmerlist = for ( i <- K to seq.length) yield minKmerByRC(seq.substring(i-K,i)).hashCode
//    (for ( i <- window to kmerlist.size) yield kmerlist.slice(i-window,i).min).toSet
//    kmerlist.sorted.slice(0,10)
//  }
}
