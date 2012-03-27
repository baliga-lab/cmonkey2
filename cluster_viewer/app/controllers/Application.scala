package controllers

import play.api._
import play.api.mvc._
import java.io._
import scala.collection.mutable.HashMap
import scala.collection.mutable.ArrayBuffer
import scala.collection.JavaConversions._
import play.api.libs.json._

object Application extends Controller {

  val AppConfig = Play.current.configuration
  val ProjectConfigFile = new File("project.conf")
  val ProjectConfig = new java.util.Properties
  ProjectConfig.load(new FileReader(ProjectConfigFile))

  val OutDirectory = new File(ProjectConfig.getProperty("cmonkey.out.directory"))
  val Synonyms = SynonymsFactory.getSynonyms(
    ProjectConfig.getProperty("cmonkey.synonyms.format"),
    ProjectConfig.getProperty("cmonkey.synonyms.file"))
  val snapshotReader = new SnapshotReader(OutDirectory, Synonyms)

  val ratiosFile = if (ProjectConfig.getProperty("cmonkey.ratios.file") != null) {
    new File(ProjectConfig.getProperty("cmonkey.ratios.file"))
  } else {
    new File(OutDirectory, "ratios.tsv")
  }
  val RatiosFactory = new RatioMatrixFactory(ratiosFile, Synonyms)


  private def availableIterations = {
    val fileNames = OutDirectory.list(new FilenameFilter {
      def accept(dir: File, name: String) = SnapshotReader.JsonFilePattern.matcher(name).matches
    })
    val result = new ArrayBuffer[Int]
    for (name <- fileNames) {
      val matcher = SnapshotReader.JsonFilePattern.matcher(name)
      matcher.matches
      result += matcher.group(1).toInt
    }
    result.sortWith((v1: Int, v2: Int) => v1 < v2)
  }

  // **********************************************************************
  // ****** VIEWS
  // **********************************************************************

  def index = index2(1)

  def index2(iteration: Int) = Action {
    val iterations = availableIterations
    println("# Iterations found: " + iterations.length)

    // create sorted results
    val stats = StatsReader.readStats(OutDirectory)
    val statsIterations = stats.keySet.toArray
    java.util.Arrays.sort(statsIterations)
    val meanResiduals = new Array[Double](stats.size)
    var i = 0
    for (key <- statsIterations) {
      meanResiduals(i) = stats(key).medianResidual
      i += 1
    }

    Ok(views.html.index(snapshotReader.readSnapshot(iteration), iterations, iteration,
                        statsIterations, meanResiduals, stats.toMap))
  }

  def cluster(iteration: Int, cluster: Int) = Action {
    val snapshot = snapshotReader.readSnapshot(iteration)
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
