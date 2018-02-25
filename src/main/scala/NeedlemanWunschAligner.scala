class NeedlemanWunschAligner(size: Int, offsetLimit: Int = Settings.OO, matchBonus: Int = 1, mismatchPenalty: Int = 1, gapPenalty: Int = 50) {
  val f = Array.ofDim[Int](size + 1, size + 1)
  val r = Array.ofDim[(Int,Int)](size + 1, size + 1)
  var tmp = 0
  // f(i,j) = min score after finish alignment of t[0..i-1] and s[0..j-1]

  /**
    * Needleman-Wunsch Algorithm to align 2 string
    *
    * @param t String on top as reference
    * @param s String on bottom
    * @return a array of t.length to indicate what position in s each char of t align to
    */
  def align(t: String, s: String, q:IndexedSeq[Int], matchBonus: Int = matchBonus, mismatchPenalty: Int = mismatchPenalty, gapPenalty: Int = gapPenalty,
            offsetLimit: Int = offsetLimit, tailFlex: Boolean = false): Array[Int] = {
    val mapping = Array.fill[Int](t.length)(-1)
    if (s.length==0) return mapping
    this.synchronized {
      // init f
      for (i <- 0 to t.length) f(i)(0) = if (i > offsetLimit) Settings.OO else (q(0)+50) * i
      for (i <- 1 to s.length) f(0)(i) = if (i > offsetLimit) Settings.OO else f(0)(i-1)+q(i-1)+50
      // recursion
      for (i <- 1 to t.length; j <- 1 to s.length) {
        if (Math.abs(i - j) > offsetLimit) f(i)(j) = Settings.OO else {
          f(i)(j) = f(i - 1)(j) + (if (j == s.length && tailFlex) 0 else (if (j==s.length) q(j-1) else (q(j-1)+q(j))/2)+50)
          r(i)(j) = (i - 1, j)
          tmp = f(i)(j - 1) + (if (i == t.length && tailFlex) 0 else q(j-1)+50)
          if (tmp < f(i)(j)) { f(i)(j) = tmp; r(i)(j) = (i, j - 1) }
          if (t(i - 1) != s(j - 1)) {
            tmp = f(i - 1)(j - 1) + q(j-1)
            if (tmp < f(i)(j)) { f(i)(j) = tmp; r(i)(j) = (i - 1, j - 1) }
          } else {
            var k = 1
            while (k<2 && k <= i && k <= j && t(i - k) == s(j - k)) {
              tmp = f(i - k)(j - k)
              if (tmp < f(i)(j)) { f(i)(j) = tmp; r(i)(j) = (i - k, j - k) }
              k += 1
            }
          }
        }
      }
      var (i, j) = (t.length, s.length)
      while (i != 0 && j != 0) {
        val next = r(i)(j)
        if (next._1==i || next._2==j) { i = next._1; j = next._2}
        else while (i != next._1 && j != next._2){
          mapping(i - 1) = j - 1
          i -= 1
          j -= 1
        }
      }
    }
    mapping
  }
}