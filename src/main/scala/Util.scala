/**
  * Created by workshop on 19-Sep-17.
  */
object Util {
  val complement = Map('A'->'T','T'->'A','C'->'G','G'->'C')
  def reverseComplement(seq:String):String={
    for ( i <- seq.reverse ) yield complement(i)
  }
  def minKmerByRC(kmer:String):String={
    // for kmer with odd length, check the middle base
    // if it's 'G'/'T', return its reverse complement, otherwise return itself
    if ("GT".contains(kmer(kmer.length/2))) reverseComplement(kmer) else kmer
  }
}
