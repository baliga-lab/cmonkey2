package controllers

import java.io._
import java.util.regex._
import scala.collection.mutable.HashMap
import scala.collection.mutable.ArrayBuffer

case class RunLog(functionName: String, scaling: Array[Float],
                  active: Array[Boolean]) {
  override def toString = {
    "RunLog('%s', %d scaling)".format(functionName, scaling.length)
  }
}

object RunLogReader {
  val FilePattern = Pattern.compile("(.+)\\.runlog")
}

class RunLogReader {
  import RunLogReader._

  def readLogs(OutDirectory: File): Array[RunLog] = {
    val files = OutDirectory.listFiles(new FilenameFilter {
      def accept(dir: File, name: String) = FilePattern.matcher(name).matches
    })
    files.map { file =>
      val inreader = new BufferedReader(new FileReader(file))
      var line = inreader.readLine
      val activeVals = new ArrayBuffer[Boolean]
      val scalingVals = new ArrayBuffer[Float]
      while (line != null) {
        val comps = line.split(":")
        activeVals += (if (comps(1) == "1") true else false)
        scalingVals += comps(2).toFloat
        line = inreader.readLine
      }
      RunLog(file.getName.replace(".runlog", ""), scalingVals.toArray,
             activeVals.toArray)
    }.toArray
  }
}
