package controllers

import play.api._
import play.api.data._
import play.api.data.Forms._

import play.api.mvc._
import java.io._
import scala.collection.mutable.HashMap
import scala.collection.mutable.ArrayBuffer
import scala.collection.JavaConversions._
import scala.util.Sorting
import play.api.libs.json._

import org.systemsbiology.DatabaseConfig

case class IntHistogram(xvalues: Array[Int], yvalues: Array[Int])
case class ResidualHistogram(xvalues: Array[String], yvalues: Array[Int])
case class RunInfo(runStatus: RunStatus,
                   currentIteration: Int,
                   numClusters: Int,
                   iterations: Seq[Int],
                   progress: Double) {
  def organismCode = runStatus.organismCode
  def species = runStatus.species
  def finished = runStatus.finished
  def startTime = runStatus.startTime
  def finishTime = runStatus.finishTime
  def numRows = runStatus.numRows
  def numColumns = runStatus.numColumns
  def numIterations = runStatus.numIterations
}
case class RunConfig(organism: String, url: String)

object Application extends Controller {

  val AppConfig = Play.current.configuration
  val OutDirectory = new File((new File(System.getProperty("user.dir"))).getParentFile, "out")
  val Synonyms = SynonymsFactory.getSynonyms(null, null)

  val snapshotReader   = new IterationResultReader(Synonyms)
  val statsReader      = new StatsReader
  val runlogReader     = new RunLogReader
  def ratiosFile = new File(OutDirectory, "ratios.tsv.gz")
  lazy val RatiosFactory = new RatioMatrixFactory(ratiosFile, Synonyms)

  // Define a configuration form
  val configForm = Form(
    mapping(
      "organism" -> text, "url" -> text
    )(RunConfig.apply)(RunConfig.unapply)
  )

  // **********************************************************************
  // ****** VIEWS
  // **********************************************************************

  def index = {
    if (OutDirectory.exists) index2(1)
    else startrun
  }

  def startrun = Action {
    Ok(views.html.startrun(configForm.fill(RunConfig("hal", "(none)"))))
  }
  def submitrun = Action { implicit request =>
    val config = configForm.bindFromRequest.get

    // simply take the uploaded file and copy into a reachable destination
    for (body <- request.body.asMultipartFormData;
         ratios <- body.file("ratios")) {
      val filename = ratios.filename
      if (filename.length > 0) {
        ratios.ref.moveTo(new File(new File("/tmp/play-upload/"), filename))
      }
    }
    Ok(views.html.submitrun())
  }


  def index2(iteration: Int) = Action {

    // sort keys ascending by iteration number
    val start = System.currentTimeMillis
    val runStatus = RunStatusReader.readRunStatus.get
    val stats = statsReader.readStats.toMap
    val runLogs = runlogReader.readLogs(OutDirectory)
    val statsIterations = stats.keySet.toArray
    val lastIteration = runStatus.lastIteration.getOrElse(1)
    java.util.Arrays.sort(statsIterations)

    makeRowStats(stats)
    val progress = math.min((lastIteration / runStatus.numIterations.toDouble * 100.0), 100.0)
    
    val runInfo = RunInfo(runStatus, iteration, runStatus.numClusters, statsIterations,
                          progress)

    val elapsed = System.currentTimeMillis - start
    println("extract index data in " + elapsed + " ms.")
    Ok(views.html.index(runInfo,
                        makeMeanResiduals(statsIterations, stats),
                        stats,
                        makeRowStats(stats),
                        makeColumnStats(stats),
                        makeResidualHistogram(stats),
                        runLogs))
  }

  // TODO: left over from the time data was stored in JSON, we can let the
  // database sort now
  private def sortByResidual(resultOption: Option[IterationResult]): Seq[Int] = {
    object MyOrdering extends Ordering[(Int, Double)] {
      def compare(a: (Int, Double), b: (Int, Double)) = {
        if (a._2 < b._2) -1
        else if (a._2 > b._2) 1
        else 0
      }
    }
    if (resultOption != None) {
      val residuals: Array[(Int, Double)] = resultOption.get.residuals.toArray
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
    val runStatus = RunStatusReader.readRunStatus.get
    val result = snapshotReader.readIterationResult(iteration)
    val ratios = RatiosFactory.readRatios(result.get.rows(cluster).toArray,
                                          result.get.columns(cluster).toArray)
    val rows    = result.get.rows(cluster)
    val columns = result.get.columns(cluster)
    val motifInfoMap = new HashMap[String, Array[MotifInfo]]
    val pssmMap = new HashMap[String, Array[String]]
    val annotationMap = new HashMap[String, Seq[GeneAnnotations]]
    val gagglePSSMs = new HashMap[String, Array[String]]

    // motif extraction
    for (seqType <- result.get.motifs.keys) {
      val motifObj = result.get.motifs(seqType) 
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
        gagglePSSMs(seqType) = if (motifObj(cluster).size > 0) {
          gagglePSSM(motifObj(cluster))
        } else {
          new Array[String](0)
        }
      }
    }
    Ok(views.html.cluster(iteration, cluster, rows, columns, ratios,
                          motifInfoMap.toMap, pssmMap.toMap,
                          annotationMap.toMap,
                          gagglePSSMs.toMap, runStatus.species))
  }

  /**
   * An action that returns the main cluster list of an iteration in JSON format.
   */
  def clusters(iteration: Int) = Action {
    val resultOption =  snapshotReader.readIterationResult(iteration)
    val result = resultOption.get
    val clusters = sortByResidual(resultOption)
    val outlist = Array.ofDim[Array[String]](clusters.length)

    for (i <- 0 until clusters.length) {
      val cluster = clusters(i)
      outlist(i) = Array("%d".format(i + 1),
                         "<a class=\"clusterlink\" id=\"%d\"  href=\"#\">%d</a>".format(cluster, cluster),
                         "%d".format(result.rows(cluster).length),
                         "%d".format(result.columns(cluster).length),
                         "%.2e".format(result.residuals(cluster)),
                         Formatter.formatMotifInfos(result.motifs, cluster))
    }
    Ok(Json.toJson(Map("aaData" -> outlist)))
  }


  private def toJsonPssm(motifInfos: Seq[MotifInfo]): Array[String] = {
    val result = new java.util.ArrayList[String]
    for (i <- 0 until motifInfos.length) {
      result.add(Json.stringify(JsObject(List("alphabet" -> Json.toJson(Array("A", "C", "G", "T")),
                                              "values" -> Json.toJson(motifInfos(i).pssm)))))
    }
    result.toArray(new Array[String](0))
  }

  def gagglePSSM(motifInfos: Seq[MotifInfo]): Array[String] = {
    val result = new java.util.ArrayList[String]
    for (i <- 0 until motifInfos.length) {
      val pssm = motifInfos(i).pssm
      val buffer = new StringBuilder("POSITION\tA\tC\tG\tT\n")
      for (j <- 0 until pssm.length) {
        buffer.append("%d\t".format(j + 1))
        buffer.append(pssm(j).mkString("\t"))
        buffer.append("\n")
      }
      result.add(buffer.toString)
    }
    result.toArray(new Array[String](0))
  }
}
