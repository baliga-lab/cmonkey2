package controllers

import play.api._
import play.api.mvc._
import play.api.libs.json._
import java.io._
import java.util.regex._
import scala.collection.mutable.HashMap
import scala.collection.mutable.ArrayBuffer
import scala.collection.JavaConversions._


case class MotifInfo(motifNum: Int, evalue: Double, pssm: Array[Array[Float]])
case class Snapshot(rows: Map[Int, List[String]], columns: Map[Int, List[String]],
                    motifs: Map[Int, Map[String, Array[MotifInfo]]]) {
  def clusters: Seq[Int] = {
    rows.keys.toSeq.sorted
  }
}

case class ClusterStats(numRows: Int, numColumns: Int, residual: Double)
case class IterationStats(clusters: Map[Int, ClusterStats], medianResidual: Double)

object Application extends Controller {

  val JsonFilePattern = Pattern.compile("(\\d+)-results.json")
  val StatsFilePattern = Pattern.compile("(\\d+)-stats.json")
  val AppConfig = Play.current.configuration
  val ProjectConfigFile = new File("project.conf")
  val ProjectConfig = new java.util.Properties
  ProjectConfig.load(new FileReader(ProjectConfigFile))

  val OutDirectory = new File(ProjectConfig.getProperty("cmonkey.out.directory"))
  val BaseResultsFileName = "%d-results.json"
  val Synonyms = SynonymsFactory.getSynonyms(
    ProjectConfig.getProperty("cmonkey.synonyms.format"),
    ProjectConfig.getProperty("cmonkey.synonyms.file"))

  val ratiosFile = if (ProjectConfig.getProperty("cmonkey.ratios.file") != null) {
    new File(ProjectConfig.getProperty("cmonkey.ratios.file"))
  } else {
    new File(OutDirectory, "ratios.tsv")
  }
  val RatiosFactory = new RatioMatrixFactory(ratiosFile, Synonyms)

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

  private def readSnapshot(iteration: Int): Option[Snapshot] = {
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

  private def availableIterations = {
    val fileNames = OutDirectory.list(new FilenameFilter {
      def accept(dir: File, name: String) = JsonFilePattern.matcher(name).matches
    })
    val result = new ArrayBuffer[Int]
    for (name <- fileNames) {
      val matcher = JsonFilePattern.matcher(name)
      matcher.matches
      result += matcher.group(1).toInt
    }
    result.sortWith((v1: Int, v2: Int) => v1 < v2)
  }

  private def readStats = {
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

  // **********************************************************************
  // ****** VIEWS
  // **********************************************************************

  def index = index2(1)

  def index2(iteration: Int) = Action {
    val iterations = availableIterations
    println("# Iterations found: " + iterations.length)

    // create sorted results
    val stats = readStats
    val statsIterations = stats.keySet.toArray
    java.util.Arrays.sort(statsIterations)
    val meanResiduals = new Array[Double](stats.size)
    var i = 0
    for (key <- statsIterations) {
      meanResiduals(i) = stats(key).medianResidual
      i += 1
    }

    Ok(views.html.index(readSnapshot(iteration), iterations, iteration,
                        statsIterations, meanResiduals, stats.toMap))
  }

  def cluster(iteration: Int, cluster: Int) = Action {
    val snapshot = readSnapshot(iteration)
    val ratios = RatiosFactory.readRatios(snapshot.get.rows(cluster).toArray)
    val rows    = snapshot.get.rows(cluster)
    val columns = snapshot.get.columns(cluster)

    if (snapshot.get.motifs.contains(cluster)) {
      val motifInfos = new java.util.ArrayList[MotifInfo]
      val motifMap = snapshot.get.motifs(cluster)
      for (seqType <- motifMap.keys) {
        val motifMapInfos = motifMap(seqType)
        for (info <- motifMapInfos) motifInfos.add(info)
      }

      // NOTE: there is no differentiation between sequence types yet !!!!
      val pssm = if (snapshot.get.motifs.size > 0) {
        toJsonPssm(snapshot.get.motifs(cluster))
      } else {
        new Array[String](0)
      }
      Ok(views.html.cluster(iteration, cluster, rows, columns, ratios,
                            motifInfos.toArray(new Array[MotifInfo](0)), pssm))
    } else {
      Ok(views.html.cluster(iteration, cluster, rows, columns, ratios,
                            new Array[MotifInfo](0), new Array[String](0)))
    }
  }

  private def toJsonPssm(motifMap: Map[String, Array[MotifInfo]]): Array[String] = {
    val result = new java.util.ArrayList[String]
    for (seqType <- motifMap.keys) {
      val motifInfos = motifMap(seqType)
      printf("SEQ TYPE: %s, # PSSMs: %d\n", seqType, motifInfos.length)
      for (i <- 0 until motifInfos.length) {
        result.add(Json.stringify(JsObject(List("alphabet" -> Json.toJson(Array("A", "C", "G", "T")),
                                                "values" -> Json.toJson(motifInfos(i).pssm)))))
      }
    }
    println("# PSSMS: " + result.length)
    result.toArray(new Array[String](0))
  }
}
