import scala.collection.mutable.Map
/**
  * Created by workshop on 02-Nov-17.
  */
class UnionFindSet[DataType] {
  private val father = Map[DataType,DataType]()
  def find(t: DataType):DataType = {
    if (!father.contains(t)) {father(t)=t; t} else iterFind(t)
  }
  def iterFind(t: DataType):DataType = {
    if (father(t)!=t) father(t) = iterFind(father(t))
    father(t)
  }

  def union(a: DataType, b: DataType) = {
    val A = find(a)
    val B = find(b)
    if (A!=B) father(A) = B
  }

  def isEmpty() = father.isEmpty
  def contains(elem:DataType) = father.contains(elem)
  def getAnyKey() = father.last._1

  def printDebugInfo() = {
    for (entry <- father) {
      println(entry)
    }
  }
}
