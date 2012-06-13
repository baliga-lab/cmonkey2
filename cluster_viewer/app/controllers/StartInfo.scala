package controllers

import play.api.libs.json._
import java.io._

case class StartInfo(startTime: String, numIterations: Int, organismCode: String,
                     species: String, numRows: Int, numConditions: Int)


class StartInfoReader {

  implicit object StartInfoFormat extends Format[StartInfo] {
    def reads(json: JsValue): StartInfo = {
      val jsonObj = json.as[JsObject]
      val startTime = (jsonObj \ "start_time").as[String]
      val numIterations = (jsonObj \ "num_iterations").as[Int]
      val organismCode = (jsonObj \ "organism-code").as[String]
      val species = (jsonObj \ "species").as[String]
      StartInfo(startTime, numIterations, organismCode, species,
                (jsonObj \ "num_rows").as[Int],
                (jsonObj \ "num_columns").as[Int])
    }
    def writes(log: StartInfo): JsValue = JsUndefined("TODO")
  }

  def readStartInfo(OutDirectory: File): Option[StartInfo] = {
    val file = new File(OutDirectory, "start.json")
    if (file.exists) {
      val in = new BufferedReader(new FileReader(file))
      val buffer = new StringBuilder
      var line = in.readLine
      while (line != null) {
        buffer.append(line)
        line = in.readLine
      }
      in.close
      Some(play.api.libs.json.Json.parse(buffer.toString).as[StartInfo])
    } else None
  }
}
