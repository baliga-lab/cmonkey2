package controllers

import play.api._
import play.api.mvc._
import java.io._
import scala.collection.mutable.HashMap
import scala.collection.mutable.ArrayBuffer
import scala.collection.JavaConversions._
import scala.util.Sorting
import play.api.libs.json._

case class IntHistogram(xvalues: Array[Int], yvalues: Array[Int])
case class ResidualHistogram(xvalues: Array[String], yvalues: Array[Int])
case class RunInfo(startInfo: StartInfo,
                   finished: Boolean,
                   finishTime: String,
                   currentIteration: Int,
                   clusters: Seq[Int],
                   iterations: Seq[Int],
                   statsIterations: Seq[Int],
                   progress: Double)

object Application extends Controller {

  val AppConfig = Play.current.configuration
  val ProjectConfigFile = new File("project.conf")
  val ProjectConfig = new java.util.Properties
  ProjectConfig.load(new FileReader(ProjectConfigFile))

  val OutDirectory = new File(ProjectConfig.getProperty("cmonkey.out.directory"))
  val Synonyms = SynonymsFactory.getSynonyms(
    ProjectConfig.getProperty("cmonkey.synonyms.format"),
    ProjectConfig.getProperty("cmonkey.synonyms.file"))
  val snapshotReader   = new SnapshotReader(OutDirectory, Synonyms)
  val statsReader      = new StatsReader
  val runlogReader     = new RunLogReader
  val startInfoReader  = new StartInfoReader
  val finishInfoReader = new FinishInfoReader

  val ratiosFile = if (ProjectConfig.getProperty("cmonkey.ratios.file") != null) {
    // explicit naming of ratios file
    new File(ProjectConfig.getProperty("cmonkey.ratios.file"))
  } else {
    // finding the ratios file in the output directory, looks for compressed file
    // first, then falls back to plain tsv
    val gzipfile = new File(OutDirectory, "ratios.tsv.gz")
    val tsvfile = new File(OutDirectory, "ratios.tsv")
    if (gzipfile.exists) gzipfile
    else if (tsvfile.exists) tsvfile
    else throw new FileNotFoundException("could not find ratios file !")
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
    val startInfo = startInfoReader.readStartInfo(OutDirectory)
    val finishInfo = finishInfoReader.readFinishInfo(OutDirectory)
    val stats = statsReader.readStats.toMap
    val runLogs = runlogReader.readLogs

    val statsIterations = stats.keySet.toArray
    java.util.Arrays.sort(statsIterations)

    makeRowStats(stats)
    val snapshotOption = snapshotReader.readSnapshot(iteration)
    val clusters = sortByResidual(snapshotOption)
    val progress = math.min((snapshotIterations(snapshotIterations.length - 1) /
                             startInfo.get.numIterations.toDouble * 100.0), 100.0)
    val runInfo = RunInfo(startInfo.get,
                          finishInfo != None,
                          if (finishInfo != None) finishInfo.get.finishTime else "",
                          iteration, clusters, snapshotIterations, statsIterations,
                          progress)
    Ok(views.html.index(runInfo,
                        snapshotOption,
                        makeMeanResiduals(statsIterations, stats),
                        stats,
                        makeRowStats(stats),
                        makeColumnStats(stats),
                        makeResidualHistogram(stats),
                        runLogs.get))
  }

  private def sortByResidual(snapshotOption: Option[Snapshot]): Seq[Int] = {
    object MyOrdering extends Ordering[(Int, Double)] {
      def compare(a: (Int, Double), b: (Int, Double)) = {
        if (a._2 < b._2) -1
        else if (a._2 > b._2) 1
        else 0
      }
    }
    if (snapshotOption != None) {
      val residuals: Array[(Int, Double)] =
        snapshotOption.get.residuals.toArray
      Sorting.quickSort(residuals)(MyOrdering)
      residuals.map(_._1)
    } else Array[Int]()
  }

  private def makeRowStats(stats: Map[Int, IterationStat]) = {
    makeIntHistogram(stats, (cs: ClusterStat) => cs.numRows)
  }
  private def makeColumnStats(stats: Map[Int, IterationStat]) = {
    makeIntHistogram(stats, (cs: ClusterStat) => cs.numColumns)
  }

  private def makeIntHistogram(stats: Map[Int, IterationStat], getKey: ClusterStat => Int) = {
    val histogram = new HashMap[Int, Int]   // # rows -> # clusters
    if (!stats.isEmpty) {
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
    } else {
      IntHistogram(Array(), Array())
    }
  }

  private def makeResidualHistogram(stats: Map[Int, IterationStat]): ResidualHistogram = {
    if (!stats.isEmpty) {
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
    } else {
      ResidualHistogram(Array(), Array())
    }
  }

  private def makeMeanResiduals(statsIterations: Array[Int], stats: Map[Int, IterationStat]) = {
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
    val motifInfoMap = new HashMap[String, Array[MotifInfo]]
    val pssmMap = new HashMap[String, Array[String]]
    val annotationMap = new HashMap[String, Seq[GeneAnnotations]]

    // motif extraction
    for (seqType <- snapshot.get.motifs.keys) {
      val motifObj = snapshot.get.motifs(seqType) 
      if (motifObj.contains(cluster)) {
        // re-group annotations by gene
        val geneAnnotationMap = new HashMap[String, ArrayBuffer[Annotation]]
        val motifInfos = new java.util.ArrayList[MotifInfo]
        val motifMapInfos = motifObj(cluster)
        for (info <- motifMapInfos) {
          motifInfos.add(info)
          for (annotation <- info.annotations) {
            if (!geneAnnotationMap.containsKey(annotation.gene)) {
              geneAnnotationMap(annotation.gene) = new ArrayBuffer[Annotation]
            }
            geneAnnotationMap(annotation.gene) += annotation
          }
        }
        //println("GENE ANNOTATIONS: " + geneAnnotationMap.keys)
        val geneAnnotationList = new ArrayBuffer[GeneAnnotations]
        for (key <- geneAnnotationMap.keys) {
          // TODO: What is acually plotted is the pvalue of the gene in this cluster
          val geneAnnotations = GeneAnnotations(key, geneAnnotationMap(key).sortWith((a: Annotation, b: Annotation) => a.position < b.position))
          geneAnnotationList += geneAnnotations
        }
        motifInfoMap(seqType) = motifInfos.toArray(Array.ofDim[MotifInfo](0))
        annotationMap(seqType) = geneAnnotationList

        pssmMap(seqType) = if (motifObj(cluster).size > 0) {
          toJsonPssm(motifObj(cluster))
        } else {
          new Array[String](0)
        }
      }
    }
    Ok(views.html.cluster(iteration, cluster, rows, columns, ratios,
                          motifInfoMap.toMap, pssmMap.toMap,
                          annotationMap.toMap))
  }

  private def toJsonPssm(motifInfos: Array[MotifInfo]): Array[String] = {
    val result = new java.util.ArrayList[String]
    for (i <- 0 until motifInfos.length) {
      result.add(Json.stringify(JsObject(List("alphabet" -> Json.toJson(Array("A", "C", "G", "T")),
                                              "values" -> Json.toJson(motifInfos(i).pssm)))))
    }
    println("# PSSMS: " + result.length)
    result.toArray(new Array[String](0))
  }
}
