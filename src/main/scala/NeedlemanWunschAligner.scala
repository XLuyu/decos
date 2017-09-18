class NeedlemanWunschAligner(size: Int, offsetLimit: Int = Settings.OO, matchBonus: Int = 1, mismatchPenalty: Int = 1, gapPenalty: Int = 50) {
  val f = Array.ofDim[Int](size + 1, size + 1)
  // f(i,j) = min score after finish alignment of t[0..i-1] and s[0..j-1]

  /**
    * Needleman-Wunsch Algorithm to align 2 string
    *
    * @param t String on top as reference
    * @param s String on bottom
    * @return a array of t.length to indicate what position in s each char of t align to
    */
  def align(t: String, s: String, matchBonus: Int = matchBonus, mismatchPenalty: Int = mismatchPenalty, gapPenalty: Int = gapPenalty,
            offsetLimit: Int = offsetLimit, headFlex: Boolean = false, tailFlex: Boolean = false): Array[Int] = {
    val mapping = Array.fill[Int](t.length)(-1)
    this.synchronized {
      // init f
      for (i <- 0 to t.length) f(i)(0) = if (i > offsetLimit) Settings.OO else if (headFlex) 0 else gapPenalty * i
      for (i <- 0 to s.length) f(0)(i) = if (i > offsetLimit) Settings.OO else if (headFlex) 0 else gapPenalty * i
      // recursion
      for (i <- 1 to t.length; j <- 1 to s.length) {
        if (Math.abs(i - j) > offsetLimit) f(i)(j) = Settings.OO else {
          f(i)(j) = Math.min(
            f(i - 1)(j) + (if (j == s.length && tailFlex) 0 else gapPenalty),
            f(i)(j - 1) + (if (i == t.length && tailFlex) 0 else gapPenalty))
          if (t(i - 1) != s(j - 1)) f(i)(j) = Math.min(f(i)(j), f(i - 1)(j - 1) + mismatchPenalty) else {
            var k = 1
            while (k <= i && k <= j && t(i - k) == s(j - k)) {
              f(i)(j) = Math.min(f(i)(j), f(i - k)(j - k) - k * k * matchBonus)
              k += 1
            }
          }
        }
      }
      var (i, j) = (t.length, s.length)
      while (i != 0 && j != 0) {
        if (f(i)(j) == f(i - 1)(j) + (if (j == s.length && tailFlex) 0 else gapPenalty)) i -= 1
        else if (f(i)(j) == f(i)(j - 1) + (if (i == t.length && tailFlex) 0 else gapPenalty)) j -= 1
        else if (t(i - 1) != s(j - 1)) {
          i -= 1
          j -= 1
          mapping(i) = j
        } else {
          var k = 1
          while (f(i)(j) != f(i - k)(j - k) - k * k * matchBonus) {
            mapping(i - k) = j - k
            k += 1
          }
          mapping(i - k) = j - k
          i -= k
          j -= k
        }
      }
    }
    mapping
  }
}