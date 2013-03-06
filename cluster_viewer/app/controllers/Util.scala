package controllers
import scala.collection.Map

object Formatter {
  private def consensusString(values: Array[Array[Float]]) = {
    val Alphabet = Array('A', 'C', 'G', 'T')
    val result = new StringBuilder
    for (row <- 0 until values.length) {
      var maxcol = 0
      var maxval = 0.0
      for (col <- 0 until values(row).length) {
        if (values(row)(col) > maxval) {
          maxcol = col
          maxval = values(row)(col)
        }
      }
      result.append(Alphabet(maxcol))
    }
    result.toString
  }

  def formatMotifInfos(motifInfos: Map[String, Map[Int, Seq[MotifInfo]]],
                       cluster: Int): String = {
    val comps = motifInfos.keys.map { seqType =>
      val clusterMotifInfos = motifInfos(seqType)
      val out = if (clusterMotifInfos.contains(cluster)) {
        clusterMotifInfos(cluster).map { motifInfo => 
          "%s (%.2e)".format(consensusString(motifInfo.pssm), motifInfo.evalue)
        }.mkString(", ")
      }
      "(%s -> [%s])".format(seqType, out)
    }
    comps.mkString("<br>")
  }
}
