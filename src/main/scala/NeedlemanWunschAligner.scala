import scala.collection.mutable
class NeedlemanWunschAligner(size: Int, offsetLimit: Int = Settings.OO, matchBonus: Int = 1, mismatchPenalty: Int = 1, gapPenalty: Int = 50) {
  val f = Array.ofDim[Int](size + 1, size + 1)
  val route = Array.ofDim[(Int, Int)](size + 1, size + 1)
  val combo = Array.ofDim[Int](size + 1, size + 1)
  // f(i,j) = min score after finish alignment of t[0..i-1] and s[0..j-1]

  /**
    * Needleman-Wunsch Algorithm to align 2 string
    *
    * @param t String on top as reference
    * @param s String on bottom
    * @return a array of t.length to indicate what position in s each char of t align to
    */
  def align(t: String, tq: IndexedSeq[mutable.Map[Char,Int]], s: String, sq: IndexedSeq[Int],
            matchBonus: Int = matchBonus, mismatchPenalty: Int = mismatchPenalty, gapPenalty: Int = gapPenalty,
            offsetLimit: Int = offsetLimit, tailFlex: Boolean = false): Array[Int] = {
    val mapping = Array.fill[Int](t.length)(-1)
    if (s.length == 0 || t.length == 0) return mapping
    this.synchronized {
      // init f
      f(0)(0) = 0; route(0)(0) = (-1, -1); combo(0)(0) = 0
      for (i <- 1 to t.length) {
        f(i)(0) = if (i > offsetLimit) Settings.OO else f(i - 1)(0) + Math.min(tq(i-1)(t(i-1))-tq(i - 1)('-'), sq(0)) * gapPenalty
        route(i)(0) = (i - 1, 0)
        combo(i)(0) = 0
      }
      for (i <- 1 to s.length) {
        f(0)(i) = if (i > offsetLimit) Settings.OO else f(0)(i - 1) + Math.min(tq(0)(t(0)), sq(i - 1)) * gapPenalty
        route(0)(i) = (0, i - 1)
        combo(0)(i) = 0
      }
      // recursion
      for (i <- 1 to t.length; j <- 1 to s.length) {
        if (Math.abs(i - j) > offsetLimit) f(i)(j) = Settings.OO else {
          combo(i)(j) = if (t(i - 1) != s(j - 1)) 0 else combo(i-1)(j-1)+Math.min(tq(i - 1)(t(i-1)), sq(j - 1))
          f(i)(j) = f(i - 1)(j) + (if (j == s.length && tailFlex) 0 else Math.min(tq(i - 1)(t(i-1))-tq(i-1)('-'), sq(j - 1)) * gapPenalty)
          route(i)(j) = (i - 1, j)
          val tmp = f(i)(j - 1) + (if (i == t.length && tailFlex) 0 else Math.min(tq(i - 1)(t(i-1)), sq(j - 1)) * gapPenalty)
          if (tmp < f(i)(j)) { f(i)(j) = tmp; route(i)(j) = (i, j - 1) }
          if (t(i - 1) != s(j - 1)) {
            val tmp = f(i - 1)(j - 1) + Math.min(tq(i - 1)(t(i-1))-tq(i-1)(s(j-1)), sq(j - 1)) * mismatchPenalty
            if (tmp < f(i)(j)) { f(i)(j) = tmp; route(i)(j) = (i - 1, j - 1) }
          } else {
//            val tmp = f(i - 1)(j - 1) - combo(i)(j) * matchBonus
//            if (tmp < f(i)(j)) { f(i)(j) = tmp; route(i)(j) = (i - 1, j - 1) }
            var k = 1
            var score = 0
            while (k <= i && k <= j && t(i - k) == s(j - k)) {
              score += Math.min(tq(i - k)(t(i-k)), sq(j - k))
              val tmp = f(i - k)(j - k) - k * score * matchBonus
              if (tmp < f(i)(j)) { f(i)(j) = tmp; route(i)(j) = (i - k, j - k) }
              k += 1
            }
          }
        }
      }
      var (i, j) = (t.length, s.length)
      while (i != 0 && j != 0) {
        val next = route(i)(j)
        if (next._1 == i || next._2 == j)
          { i = next._1; j = next._2 }
        else
        while (i != next._1 && j != next._2) {
          mapping(i - 1) = j - 1
          i -= 1
          j -= 1
        }
      }
    }
    mapping
  }
}