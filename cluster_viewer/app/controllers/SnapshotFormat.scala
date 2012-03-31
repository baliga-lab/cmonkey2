package controllers

import play.api.libs.json._
import java.io._
import java.util.regex._
import scala.collection.mutable.HashMap

case class MotifInfo(motifNum: Int, evalue: Double, pssm: Array[Array[Float]])
case class Snapshot(rows: Map[Int, List[String]], columns: Map[Int, List[String]],
                    motifs: Map[Int, Map[String, Array[MotifInfo]]]) {
  def clusters: Seq[Int] = {
    rows.keys.toSeq.sorted
  }
}

object SnapshotReader {
  val BaseResultsFileName = "%d-results.json"
  val JsonFilePattern = Pattern.compile("(\\d+)-results.json")
}
class SnapshotReader(OutDirectory: File, Synonyms: SynonymsMap) {
  import SnapshotReader._

  implicit object SnapshotFormat extends Format[Snapshot] {
    def reads(json: JsValue): Snapshot = {
      val rows = (json \ "rows").as[JsObject]
      val cols = (json \ "columns").as[JsObject]
      val motifsVal = (json \ "motifs")
      val clusterRows = new HashMap[Int, List[String]]
      for (field <- rows.fields) {
        clusterRows(field._1.toInt) = field._2.as[List[String]].map(str => Synonyms(str))
        //println(clusterRows(field._1.toInt))
      }

      val clusterCols = new HashMap[Int, List[String]]
      for (field <- cols.fields) {
        clusterCols(field._1.toInt) = field._2.as[List[String]]
      }

      val clusterMotifs = new HashMap[Int, Map[String, Array[MotifInfo]]]
      try {
        val motifs = motifsVal.as[JsObject]        
        for (field <- motifs.fields) {
          val cluster = field._1
          val seqTypeObj = field._2.as[JsObject]

          val seqTypeMotifs = new HashMap[String, Array[MotifInfo]]

          // iterate over the sequence types, which are keys
          for (stfield <- seqTypeObj.fields) {
            val seqType = stfield._1
            // an array of motif objects (motif_num, evalue, annotations, pssm)
            // annotations are triples of (gene, position, pvalue)
            val stMotifs = stfield._2.asInstanceOf[JsArray].value
            val motifInfos = new java.util.ArrayList[MotifInfo]
            for (motif <- stMotifs) {
              val motifObj = motif.asInstanceOf[JsObject]
              val pssm = motifObj.value("pssm").as[Array[Array[Float]]]
              val evalue = if (motifObj.keys.contains("evalue")) {
                motifObj.value("evalue").asInstanceOf[JsNumber].value.doubleValue
              } else 0.0
              val motifNum = motifObj.value("motif_num").asInstanceOf[JsNumber].value.intValue
              motifInfos.add(MotifInfo(motifNum, evalue, pssm))
            }
            seqTypeMotifs(seqType) = motifInfos.toArray(new Array[MotifInfo](0))
          }
          clusterMotifs(field._1.toInt) = seqTypeMotifs.toMap
        }
      } catch {
        case _ => println("\nNo motifs found !!!")
      }
      Snapshot(clusterRows.toMap, clusterCols.toMap, clusterMotifs.toMap)
    }

    def writes(snapshot: Snapshot): JsValue = JsUndefined("TODO")
  }

  def readSnapshot(iteration: Int) : Option[Snapshot] = {
    val pathname = (OutDirectory + "/" + BaseResultsFileName).format(iteration)
    printf("Reading snapshot: %s\n", pathname)
    val infile = new File(pathname)
    if (infile.exists) {
      val in = new BufferedReader(new FileReader(infile))
      val buffer = new StringBuilder
      var line = in.readLine
      while (line != null) {
        buffer.append(line)
        line = in.readLine
      }
      in.close
      Some(play.api.libs.json.Json.parse(buffer.toString).as[Snapshot])
    } else {
      printf("File '%s' does not exist !\n", infile.getName)
      None
    }
  }
}
