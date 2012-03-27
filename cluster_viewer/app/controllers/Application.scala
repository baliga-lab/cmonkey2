package controllers

import play.api._
import play.api.mvc._
import java.io._
import scala.collection.mutable.HashMap
import scala.collection.mutable.ArrayBuffer
import scala.collection.JavaConversions._
import play.api.libs.json._

case class IntHistogram(xvalues: Array[Int], yvalues: Array[Int])
case class ResidualHistogram(xvalues: Array[String], yvalues: Array[Int])

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
  val statsReader = new StatsReader(OutDirectory)

  val ratiosFile = if (ProjectConfig.getProperty("cmonkey.ratios.file") != null) {
    new File(ProjectConfig.getProperty("cmonkey.ratios.file"))
  } else {
    new File(OutDirectory, "ratios.tsv")
  }
  val RatiosFactory = new RatioMatrixFactory(ratiosFile, Synonyms)

  private def snapshotIterations = {
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
    // sort keys ascending by iteration number
    val stats = statsReader.readStats(OutDirectory).toMap
    val statsIterations = stats.keySet.toArray
    java.util.Arrays.sort(statsIterations)

    makeRowStats(stats)
    Ok(views.html.index(snapshotReader.readSnapshot(iteration),
                        snapshotIterations,
                        iteration,
                        statsIterations,
                        makeMeanResiduals(statsIterations, stats),
                        stats,
                        makeRowStats(stats),
                        makeColumnStats(stats),
                        makeResidualHistogram(stats)))
  }

  private def makeRowStats(stats: Map[Int, IterationStats]) = {
    makeIntHistogram(stats, (cs: ClusterStats) => cs.numRows)
  }
  private def makeColumnStats(stats: Map[Int, IterationStats]) = {
    makeIntHistogram(stats, (cs: ClusterStats) => cs.numColumns)
  }

  private def makeIntHistogram(stats: Map[Int, IterationStats], getKey: ClusterStats => Int) = {
    val histogram = new HashMap[Int, Int]   // # rows -> # clusters
    val maxIteration = stats.keys.max
    for (cluster <- stats(maxIteration).clusters.keys) {
      val key = getKey(stats(maxIteration).clusters(cluster))
      if (!histogram.contains(key)) histogram(key) = 0
      histogram(key) += 1
    }

    val sortedKeys = histogram.keySet.toArray
    java.util.Arrays.sort(sortedKeys)
    val values = new Array[Int](sortedKeys.length)
    var i = 0
    for (key <- sortedKeys) {
      values(i) = histogram(key)
      i += 1
    }
    IntHistogram(sortedKeys, values)
  }

  private def makeResidualHistogram(stats: Map[Int, IterationStats]): ResidualHistogram = {
    val numBuckets = 20
    val maxIteration = stats.keys.max
    var minResidual = 10000000.0f
    var maxResidual = -10000000.0f
    for (cluster <- stats(maxIteration).clusters.keys) {
      val residual = stats(maxIteration).clusters(cluster).residual.asInstanceOf[Float]
      if (residual < minResidual) minResidual = residual
      if (residual > maxResidual) maxResidual = residual
    }

    val xvalues = new Array[String](numBuckets)
    val yvalues = new Array[Int](numBuckets)
    val interval = (maxResidual - minResidual) / numBuckets
    // making labels
    for (i <- 0 until numBuckets) {
      xvalues(i) = "%.2f".format(minResidual + interval * i)
    }
    // dump residuals in the right buckets
    for (cluster <- stats(maxIteration).clusters.keys) {
      val residual = stats(maxIteration).clusters(cluster).residual.asInstanceOf[Float]
      val bucketnum = math.min(numBuckets - 1, ((residual - minResidual) / interval).asInstanceOf[Int])
      yvalues(bucketnum) += 1
    }
    ResidualHistogram(xvalues, yvalues)
  }

  private def makeMeanResiduals(statsIterations: Array[Int], stats: Map[Int, IterationStats]) = {
    val meanResiduals = new Array[Double](stats.size)
    var i = 0
    for (key <- statsIterations) {
      meanResiduals(i) = stats(key).medianResidual
      i += 1
    }
    meanResiduals
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
