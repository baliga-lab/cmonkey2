package controllers

object Formatter {
  def formatEvalues(clusterMotifInfos: Map[Int, Map[String, Array[MotifInfo]]], cluster: Int): String = {
    if (clusterMotifInfos.contains(cluster)) {
      val comps = clusterMotifInfos(cluster).keys.map { seqType =>
        val evals = clusterMotifInfos(cluster)(seqType).map { motifInfo => 
          "%.8f".format(motifInfo.evalue)
        }.mkString(",")
        "(%s -> [%s])".format(seqType, evals)
      }
      comps.mkString("<br>")
    } else "-"
  }
}
