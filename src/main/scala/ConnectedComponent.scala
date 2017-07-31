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
      val m = Math.min(u,T.min)
      val newEdge = for ( node <- T if node > u) yield (node,m)
      if (m!=u && !newEdge.isEmpty) accum.add(1)
      newEdge
    }
    def smallStarReduce(uN:(Long,Iterable[Long])):TraversableOnce[(Long,Long)] = {
      val u = uN._1
      val N = uN._2.toSet
      val m = N.min
      if (N.size > 1) accum.add(1)
      (for ( node <- N if node != m) yield (node,m)) + Tuple2(u,m)
    }

    var edges = nodePair
    do {
      accum.reset()
      edges = edges.flatMap(x => List(x,(x._2,x._1))).groupByKey().flatMap(largeStarReduce)
      edges = edges.groupByKey().flatMap(smallStarReduce)
      edges.count()
    } while (accum.value!=0)
    edges
  }
  def runByNodeList(sc: SparkContext, nodeList:RDD[List[Long]]): RDD[(Long, Long)] = {
    def listToPair(nodes:List[Long]) : List[(Long, Long)] = {
      val m = nodes.min
      for ( i <- nodes if i!=m) yield (i,m)
    }
    val nodePair = nodeList.flatMap(listToPair)
    run(sc,nodePair)
  }
}


