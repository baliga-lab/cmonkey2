package controllers

import play.api.libs.json._
import java.io._
import java.util.regex._
import scala.collection.mutable.HashMap

case class ClusterStats(numRows: Int, numColumns: Int, residual: Double)
case class IterationStats(clusters: Map[Int, ClusterStats], medianResidual: Double,
                          motifPValue: Double, networkScores: Map[String, Double])

object StatsReader {
  val FilePattern = Pattern.compile("(\\d+)-stats.json")
}

class StatsReader {
  import StatsReader._

  implicit object StatsFormat extends Format[IterationStats] {
    def reads(json: JsValue): IterationStats = {
      val residual = (json \ "median_residual").asInstanceOf[JsNumber].value.doubleValue
      val clusters = (json \ "cluster").as[JsObject]
      val motifPValue = (json \ "motif-pvalue").asInstanceOf[JsNumber].value.doubleValue
      val networkScoresJson = (json \ "network-scores").as[JsObject]

      val clusterStats = new HashMap[Int, ClusterStats]
      for (field <- clusters.fields) {
        val cstats = field._2
        clusterStats(field._1.toInt) = ClusterStats((cstats \ "num_rows").asInstanceOf[JsNumber].value.intValue,
                                                    (cstats \ "num_columns").asInstanceOf[JsNumber].value.intValue,
                                                    (cstats \ "residual").asInstanceOf[JsNumber].value.doubleValue)
      }
      val networkScores = new HashMap[String, Double]
      for (field <- networkScoresJson.fields) {
        networkScores(field._1) =
          field._2.asInstanceOf[JsNumber].value.doubleValue
      }
      IterationStats(clusterStats.toMap, residual, motifPValue, networkScores.toMap)
    }
    def writes(stats: IterationStats): JsValue = JsUndefined("TODO")
  }

  def readStats(OutDirectory: File) = {
    val files = OutDirectory.listFiles(new FilenameFilter {
      def accept(dir: File, name: String) = FilePattern.matcher(name).matches
    })
    val stats = new HashMap[Int, IterationStats]

    for (i <- 0 until files.length) {
      val file = files(i)
      val in = new BufferedReader(new FileReader(file))
      val buffer = new StringBuilder
      var line = in.readLine
      while (line != null) {
        buffer.append(line)
        line = in.readLine
      }
      in.close

      val statsOption = Some(play.api.libs.json.Json.parse(buffer.toString).as[IterationStats])      
      if (statsOption != None) {
        val matcher = FilePattern.matcher(file.getName)
        matcher.matches
        val iteration = matcher.group(1).toInt
        stats(iteration) = statsOption.get
      }
    }
    stats
  }
}
