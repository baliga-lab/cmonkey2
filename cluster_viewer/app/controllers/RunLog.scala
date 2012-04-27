package controllers

import play.api.libs.json._
import java.io._
import java.util.regex._
import scala.collection.mutable.HashMap

case class RunLog(functionName: String, active: Array[Boolean], scaling: Array[Float]) {
  override def toString = {
    "RunLog('%s', %d active, %d scaling)".format(functionName, active.length, scaling.length)
  }
}

object RunLogReader {
  val FilePattern = Pattern.compile("(\\d+)-runlog.json")
}

class RunLogReader {
  import RunLogReader._

  implicit object RunLogFormat extends Format[Array[RunLog]] {
    def reads(json: JsValue): Array[RunLog] = {
      val result = Array.ofDim[RunLog](json.asInstanceOf[JsArray].value.length)
      var i = 0
      json.asInstanceOf[JsArray].value.foreach { log =>
        val logObj = log.as[JsObject]
        result(i) = RunLog((logObj \ "name").asInstanceOf[JsString].value,
                           (logObj \ "active").asInstanceOf[JsArray].value.map(_.asInstanceOf[JsBoolean].value).toArray,
                           (logObj \ "scaling").asInstanceOf[JsArray].value.map(_.asInstanceOf[JsNumber].value.toFloat).toArray)
        i += 1
      }
      result
    }
    def writes(log: Array[RunLog]): JsValue = JsUndefined("TODO")
  }

  def readLogs(OutDirectory: File): Option[Array[RunLog]] = {
    val files = OutDirectory.listFiles(new FilenameFilter {
      def accept(dir: File, name: String) = FilePattern.matcher(name).matches
    })
    var max = 0
    var maxFile: File = null
    for (file <- files) {
      val matcher = FilePattern.matcher(file.getName)
      if (matcher.matches) {
        val num = matcher.group(1).toInt
        if (num > max) {
          max = num
          maxFile = file
        }
      }
    }
    if (maxFile == null) None
    else {
      val in = new BufferedReader(new FileReader(maxFile))
      val buffer = new StringBuilder
      var line = in.readLine
      while (line != null) {
        buffer.append(line)
        line = in.readLine
      }
      in.close
      Some(play.api.libs.json.Json.parse(buffer.toString).as[Array[RunLog]])
    }
  }
}
