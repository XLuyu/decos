class NeedlemanWunschAligner(size: Int, mismatchPenalty: Int = 1, gapPenalty: Int = 100) {
  val f = Array.ofDim[Int](size + 1, size + 1)
  // f(i,j) = min score after finish alignment of t[0..i-1] and s[0..j-1]

  /**
    * Needleman-Wunsch Algorithm to align 2 string
    *
    * @param t String on top as reference
    * @param s String on bottom
    * @return a array of t.length to indicate what position in s each char of t align to
    */
  def align(t: String, s: String): Array[Int] = {
    val mapping = Array.fill[Int](t.length)(-1)
    // init f
    for (i <- 0 to t.length) f(i)(0) = gapPenalty * i
    for (i <- 0 to s.length) f(0)(i) = gapPenalty * i
    // recursion
    for (i <- 1 to t.length; j <- 1 to s.length)
      f(i)(j) = Math.min(Math.min(f(i-1)(j)+gapPenalty,
                                  f(i)(j-1)+gapPenalty),
                                  f(i-1)(j-1)+(if (t(i-1)==s(j-1)) 0 else mismatchPenalty))
    var (i,j) = (t.length,s.length)
    while (i!=0 && j!=0){
      if (f(i)(j)==f(i-1)(j)+gapPenalty) i -= 1
      else if (f(i)(j)==f(i)(j-1)+gapPenalty) j -= 1
      else {
        i -= 1
        j -= 1
        mapping(i) = j
      }
    }
    mapping
  }
}