package controllers

import play.api.libs.json._
import java.io._
import java.util.regex._
import scala.collection.mutable.HashMap

case class ClusterStats(numRows: Int, numColumns: Int, residual: Double)
case class IterationStats(clusters: Map[Int, ClusterStats], medianResidual: Double)

object StatsReader {
  val StatsFilePattern = Pattern.compile("(\\d+)-stats.json")
}

class StatsReader(OutFile: File) {
  import StatsReader._

  implicit object StatsFormat extends Format[IterationStats] {
    def reads(json: JsValue): IterationStats = {
      val residual = (json \ "median_residual").asInstanceOf[JsNumber].value.doubleValue
      val clusters = (json \ "cluster").as[JsObject]
      val clusterStats = new HashMap[Int, ClusterStats]
      for (field <- clusters.fields) {
        val cstats = field._2.as[JsObject]
        clusterStats(field._1.toInt) = ClusterStats(cstats.value("num_rows").asInstanceOf[JsNumber].value.intValue,
                                                    cstats.value("num_columns").asInstanceOf[JsNumber].value.intValue,
                                                    cstats.value("residual").asInstanceOf[JsNumber].value.doubleValue)
      }
      IterationStats(clusterStats.toMap, residual)
    }
    def writes(snapshot: IterationStats): JsValue = JsUndefined("TODO")
  }

  def readStats(OutDirectory: File) = {
    val files = OutDirectory.listFiles(new FilenameFilter {
      def accept(dir: File, name: String) = StatsFilePattern.matcher(name).matches
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
        val matcher = StatsFilePattern.matcher(file.getName)
        matcher.matches
        val iteration = matcher.group(1).toInt
        stats(iteration) = statsOption.get
      }
    }
    stats
  }
}
