package controllers

import play.api.libs.json._
import java.io._

case class FinishInfo(finishTime: String)

class FinishInfoReader {

  implicit object StopInfoFormat extends Format[FinishInfo] {
    def reads(json: JsValue): FinishInfo = {
      val jsonObj = json.as[JsObject]
      val finishTime = (jsonObj \ "finish_time").as[String]
      FinishInfo(finishTime)
    }
    def writes(log: FinishInfo): JsValue = JsUndefined("TODO")
  }

  def readFinishInfo(OutDirectory: File): Option[FinishInfo] = {
    val file = new File(OutDirectory, "finish.json")
    if (file.exists) {
      val in = new BufferedReader(new FileReader(file))
      val buffer = new StringBuilder
      var line = in.readLine
      while (line != null) {
        buffer.append(line)
        line = in.readLine
      }
      in.close
      Some(play.api.libs.json.Json.parse(buffer.toString).as[FinishInfo])
    } else None
  }
}
