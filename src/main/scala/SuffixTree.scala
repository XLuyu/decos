import scala.collection.mutable
import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.Map
class Edge(s: Int, e: Int, d: Node, t: (Int, Int)) {
  var start = s
  var end = e
  var dest = d
  var tag = t
}

class Node() {
  var edge = Array[Edge](null, null, null, null, null)
  var suffix: Node = _
}

class SuffixTree() {
  var tid, start_idx = 0
  var T = ""
  var root = new Node()
  var falsum = new Node()
  var toRoot = new Edge(-1, 0, root, null)
  falsum.edge = Array[Edge](toRoot, toRoot, toRoot, toRoot, toRoot)
  root.suffix = falsum

  private def atoi(c: Char): Int = c match {
    case 'A' => 0
    case 'C' => 1
    case 'G' => 2
    case 'T' => 3
    case _ => 4
  }

  private def test_and_split(s: Node, k: Int, p: Int, t: Char): (Boolean, Node) = {
    if (k < p) {
      val current = s.edge(atoi(T(k)))
      if (t == T(current.start + p - k)) return (true, s)
      val r = new Node()
      s.edge(atoi(T(k))) = new Edge(current.start, current.start + p - k, r, null)
      r.edge(atoi(T(current.start + p - k))) = current
      current.start += p - k
      (false, r)
    } else (s.edge(atoi(t)) != null, s)
  }
  private def canonize(s: Node, k: Int, p: Int): (Node, Int) = {
    var flag = true
    var (s_, k_) = (s, k)
    while (k_ < p && flag) {
      val current = s_.edge(atoi(T(k_)))
      if (current.end - current.start > p - k_) { flag = false } else {
        k_ += current.end - current.start
        s_ = current.dest
      }
    }
    (s_, k_)
  }
  private def update(s: Node, k: Int, i: Int): (Node, Int) = {
    var oldr = root
    var end_r = test_and_split(s, k, i - 1, T(i - 1))
    var sk = (s, k)
    while (!end_r._1) {
      end_r._2.edge(atoi(T(i - 1))) = new Edge(i - 1, T.length, null, (tid, start_idx))
      if (oldr != root) oldr.suffix = end_r._2
      oldr = end_r._2
      sk = canonize(sk._1.suffix, sk._2, i - 1)
      start_idx += 1
      end_r = test_and_split(sk._1, sk._2, i - 1, T(i - 1))
    }
    if (oldr != root) oldr.suffix = sk._1
    sk
  }
  def append(template: String) {
    start_idx = 0
    val oldlen = T.length
    T += template + '$'
    var sk = (root, oldlen)
    for (i <- oldlen until T.length) {
      sk = update(sk._1, sk._2, i + 1)
      sk = canonize(sk._1, sk._2, i + 1)
    }
    tid += 1
  }
  def matchPatternSuffix(pattern: String): List[(Int, Int, Int, Int)] = {
    // return list of (append tid, template start, pattern start, pattern end), pattern end not inclusive
    var ret = List[(Int, Int, Int, Int)]()
    var s = root
    var k, p, pstart_idx = 0
    var current: Edge = null
    for (p <- 0 until pattern.length) {
      var reported = false
      var find = false
      while ((k < p || s != root) && !find) {
        if (k == p) {
          if (s.edge(atoi(pattern(p))) != null) find = true
        } else {
          current = s.edge(atoi(pattern(k)))
          if (current.start + p - k < T.length && T(current.start + p - k) == pattern(p)) find = true
        }
        if (!find) {
          if (s == root) k += 1
          else {
            if (!reported) {
              if (k < p && current.tag != null)
                ret = (current.tag._1, current.tag._2, pstart_idx, p) :: ret
              reported = true
            }
            s = s.suffix
            pstart_idx += 1
            if (k < p) current = s.edge(atoi(pattern(k)))
            while (k < p && current.end - current.start <= p - k) {
              s = current.dest
              k += current.end - current.start
              if (k < p) current = s.edge(atoi(pattern(k)))
            }
          }
        }
      }
      if (k == p && s == root && s.edge(atoi(pattern(p))) == null)
        k += 1
      else {
        if (k == p) current = s.edge(atoi(pattern(p)))
        if (k <= p && current.end - current.start <= p - k + 1) {
          s = current.dest
          k += current.end - current.start
        }
      }
    }
    if (k < pattern.length && s.edge(atoi(pattern(k))).tag != null)
      ret = (s.edge(atoi(pattern(k))).tag._1, s.edge(atoi(pattern(k))).tag._2, pstart_idx, pattern.length) :: ret
    ret
  }
  def tidVote(pattern: String): Array[(Int, Int)] = {
    val voteCount = mutable.Map[(Int, Int), Int]()
    val matches = matchPatternSuffix(pattern)
    for (record <- matches) {
      val v = voteCount.getOrElse((record._1, record._2 - record._3), 0) + record._4 - record._3
      voteCount((record._1, record._2 - record._3)) = v
    }
    val vote = Array.fill[(Int, Int)](tid)((0, 0))
    for ((k, v) <- voteCount)
      if (vote(k._1)._2 < v) vote(k._1) = (k._2, v)
    vote
  }
  def tidMaxVote(pattern: String): (Int, Int) = {
    val voteCount = mutable.Map[(Int, Int), Int]()
    val matches = matchPatternSuffix(pattern)
    if (matches.isEmpty) return (0,0)
    val maxVote = matches.maxBy(x=>x._4-x._3)
    (maxVote._2-maxVote._3,maxVote._4 - maxVote._3)
  }
}
object SuffixTree {
  def main(args: Array[String]) {
    var st = new SuffixTree()
//    st.append("TCTGATGAGTCGAAAAATTATCTTGATAAAGCAGGAATTACTACTGCTTGTTTACGAATTAAATCGAAGTGGACTGCTGGCGGAAAATGAGAAAATTCGA")
//    st.append("ATTAAATCGAAGTGGACTGCTGGCGGAAAATGAGAAAATTCGACCTATCCTTGCGCAGCTCGAGAAGCTCTTACTTTGCGACCTTTCGCCATCAACTAAC")
//    st.append("CCATAAACGCAAGCCTCAACGCAGCGACGAGCACGAGAGCGGTCAGTAGCAATCCAAACTTTGTTACTCGTCAGAAAATCGAAATCATCTTCGGTTAAAT")
//    st.append("CCCACAAAGTCCAGCGTACCATAAACGCAAGCCTCAACGCAGCGACGAGCACGAGAGCGGTCAGTAGCAATCCAAACTTAGTTACTCGTCAGAAAATCGG")
//    st.append("AAACTCAACAGGAGCAGGAAAGCGAGGGTATCCCACAAAGTCCAGCGTACCATAAACGCAAGCCTCAACGCAGCGACGAGCACGAGAGCGGTCAGTAGCA")
    st.append("GTCGAAAAATTATCTTGATAAAGCAGGAATTACTACTGCTTGTTTACGAATTAAATCGAAGTGGACTGCTGGCGGAAAATGAGAAAATTCGACCTATCCT")
    println(st.matchPatternSuffix("AACGCAGCGACGAGCACGAGAGCGGTCAGTAGCTATCCAAACTTTGTTACTCGTCAGAAAATCGAAATCATCTTCGGTTAAATCCAAAACGGCAGAAGCC"))
    for ((k,v) <- st.tidVote("AACGCAGCGACGAGCACGAGAGCGGTCAGTAGCTATCCAAACTTTGTTACTCGTCAGAAAATCGAAATCATCTTCGGTTAAATCCAAAACGGCAGAAGCC"))
      println(k,v)
  }
}