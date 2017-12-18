import org.apache.spark.SparkContext
import org.apache.spark.rdd.RDD

object ConnectedComponent extends Serializable {
  /**
    * MapReduce-based connected Component algorithm
    * @param sc sparkContext
    * @param nodePair is a RDD representing edges in graph
    * @return a RDD[(Long,Long)] where first element is node and the second is representative of connected component.
    *         Note that representative itself is missing from RDD (i.e. No (Long,Long)'s first element is a representative)
    *         This ensures there is no loop in new graph.(no multiple edges as well)
    */
  def run(sc: SparkContext, nodePair:RDD[(Long,Long)]): RDD[(Long, Long)] = {
    val accum = sc.longAccumulator("Converged")
    def largeStarReduce(uT:(Long,Iterable[Long])):TraversableOnce[(Long,Long)] = {
      val u = uT._1
      val T = uT._2.toSet
      val m = (T+u).minBy(Math.abs)
      val newEdge = for ( node <- T if Math.abs(node) > u) yield if (node<0) (-node,-m) else (node,m)
      if (m!=u && newEdge.nonEmpty) accum.add(1)
      newEdge
    }
    def smallStarReduce(uN:(Long,Iterable[Long])):TraversableOnce[(Long,Long)] = {
      val u = uN._1
      val N = uN._2.toSet
      val m = N.minBy(Math.abs)
      if (N.size > 1) accum.add(1)
      (for ( node <- N if node != m) yield (node,m)) + Tuple2(u,m)
    }
    def duplicateEdge(e:(Long,Long)) = List(
      if (e._1<0) (-e._1,-e._2) else e,
      if (e._2<0) (-e._2,-e._1) else (e._2,e._1)
    )
    var edges = nodePair
    var oldEdges:RDD[(Long,Long)] = null
    do {
      accum.reset()
      edges = edges.flatMap(duplicateEdge).groupByKey().flatMap(largeStarReduce)
      edges = edges.groupByKey().flatMap(smallStarReduce)
//      edges.cache()
      edges.count()
      if (oldEdges!=null) oldEdges.unpersist()
//      oldEdges = edges
    } while (accum.value!=0)
    edges
  }
  def runByNodeList(sc: SparkContext, nodeList:RDD[List[Long]]): RDD[(Long, Long)] = {
    def listToPair(nodes:List[Long]) : List[(Long, Long)] = {
      val m = nodes.minBy(Math.abs)
      for ( i <- nodes if i!=m) yield (i,m)
    }
    val nodePair = nodeList.flatMap(listToPair)
    run(sc,nodePair)
  }
}


