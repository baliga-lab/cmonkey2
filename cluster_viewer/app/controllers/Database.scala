package controllers

import java.io._
import java.sql.Timestamp
import scala.collection.Map
import java.util.regex._

import language.postfixOps
import scala.collection.mutable.{HashMap, ArrayBuffer}

import play.api.Play.current
import play.api.db._
import scala.slick.session.{Database, Session}
import scala.slick.driver.SQLiteDriver.simple._

import org.systemsbiology.DatabaseConfig

// **********************************************************************
// **** Data access functionality that is visualized is
// **** found here together for consolidation
// **** Except for the run logs, which use a simplistic text
// **** format, all other data is stored in SQLite and accessed through
// **** ScalaQuery/SLICK.
// **** I believe that it is easier to extract data from a relational
// **** database than from JSON in this case.
// **********************************************************************

// **********************************************************************
// **** Table objects
// **********************************************************************
// **********************************************************************
// **** cMonkey results
// **********************************************************************


object RowNames extends Table[(Int, String)]("row_names") {

  def orderNum = column[Int]("order_num", O NotNull)
  def name = column[String]("name", O NotNull)

  def * = orderNum ~ name
}

object ColumnNames extends Table[(Int, String)]("column_names") {

  def orderNum = column[Int]("order_num", O NotNull)
  def name = column[String]("name", O NotNull)

  def * = orderNum ~ name
}

object RowMembers extends Table[(Int, Int, Int)]("row_members") {

  def iteration = column[Int]("iteration", O NotNull)
  def cluster = column[Int]("cluster", O NotNull)
  def orderNum = column[Int]("order_num", O NotNull)
  def * = iteration ~ cluster ~ orderNum
}

object ColumnMembers extends Table[(Int, Int, Int)]("column_members") {

  def iteration = column[Int]("iteration", O NotNull)
  def cluster = column[Int]("cluster", O NotNull)
  def orderNum = column[Int]("order_num", O NotNull)

  def * = iteration ~ cluster ~ orderNum
}

object ClusterResiduals extends Table[(Int, Int, Double)]("cluster_residuals") {

  def iteration = column[Int]("iteration", O NotNull)
  def cluster = column[Int]("cluster", O NotNull)
  def residual = column[Double]("residual", O NotNull)

  def * = iteration ~ cluster ~ residual
}


object MotifInfos extends Table[(Int, Int, Int, String, Int, Double)]("motif_infos") {

  def id = column[Int]("rowid", O NotNull)
  def iteration = column[Int]("iteration", O NotNull)
  def cluster = column[Int]("cluster", O NotNull)
  def seqtype = column[String]("seqtype", O NotNull)
  def motifNum = column[Int]("motif_num", O NotNull)
  def evalue = column[Double]("evalue", O NotNull)

  def * = id ~ iteration ~ cluster ~ seqtype ~ motifNum ~ evalue
}

object MotifPSSMRows
extends Table[(Int, Int, Int, Float, Float, Float, Float)]("motif_pssm_rows") {

  def motifInfoId = column[Int]("motif_info_id", O NotNull)
  def iteration = column[Int]("iteration", O NotNull)
  def row = column[Int]("row", O NotNull)
  def a = column[Float]("a", O NotNull)
  def c = column[Float]("c", O NotNull)
  def g = column[Float]("g", O NotNull)
  def t = column[Float]("t", O NotNull)

  def * = motifInfoId ~ iteration ~ row ~ a ~ c ~ g ~ t
}

object MotifAnnotations
extends Table[(Int, Int, Int, Int, Boolean, Double)]("motif_annotations") {

  def motifInfoId = column[Int]("motif_info_id", O NotNull)
  def iteration = column[Int]("iteration", O NotNull)
  def geneNum = column[Int]("gene_num", O NotNull)
  def position = column[Int]("position", O NotNull)
  def reverse = column[Boolean]("reverse", O NotNull)
  def pvalue = column[Double]("pvalue", O NotNull)

  def * = motifInfoId ~ iteration ~ geneNum ~ position ~ reverse ~ pvalue
}

object ResultQueries {
  lazy val database = Database.forDataSource(DatabaseConfig.datasource)

  def allRowNames = database.withSession { implicit db: Session =>
    Query(RowNames).sortBy(_.name).map(_.name).list.toArray
  }
  def allColumnNames = database.withSession { implicit db: Session =>
    Query(ColumnNames).sortBy(_.name).map(_.name).list.toArray
  }

  def rowMembersForIteration(iteration: Int) = database.withSession { implicit db: Session =>
    Query(RowMembers).filter(_.iteration === iteration).
      sortBy(_.cluster).map(m => (m.cluster, m.orderNum)).list
  }
  def columnMembersForIteration(iteration: Int) = database.withSession { implicit db: Session =>
    Query(ColumnMembers).filter(_.iteration === iteration).
      sortBy(_.cluster).map(m => (m.cluster, m.orderNum)).list
  }

  def clusterResidualsForIteration(iteration: Int) = database.withSession { implicit db: Session =>
    Query(ClusterResiduals).filter(_.iteration === iteration).
      sortBy(_.cluster).map(r => (r.cluster, r.residual)).list
  }

  def motifInfosForIteration(iteration: Int) = database.withSession { implicit db: Session =>
    Query(MotifInfos).filter(_.iteration === iteration).
      sortBy(_.cluster).map(m => (m.id, m.cluster, m.seqtype, m.motifNum, m.evalue)).list
  }

  def pssmRowsForIteration(iteration: Int) = database.withSession { implicit db: Session =>
    Query(MotifPSSMRows).filter(_.iteration === iteration).
      sortBy(_.row).map(r => (r.motifInfoId, r.a, r.c, r.g, r.t)).list
  }

  def motifAnnotationsForIteration(iteration: Int) = database.withSession { implicit db: Session =>
    Query(MotifAnnotations).filter(_.iteration === iteration).
      sortBy(_.geneNum).map(a => (a.motifInfoId, a.geneNum, a.position, a.reverse, a.pvalue)).list
  }
}

object MotifQueries {
  
  lazy val database = Database.forDataSource(DatabaseConfig.datasource)

  def findForIteration(iteration: Int, rowNames: Array[String]) =
    // the result structure, a map from seqtype->cluster->motif info
    database.withSession { implicit db: Session =>

      // motif id -> pssm row
      val pssmMap = new HashMap[Int, ArrayBuffer[Array[Float]]]
      // motif id -> annotation
      val annotMap =(new HashMap[Int, ArrayBuffer[(Int, Int, Boolean, Double)]]).
                          withDefaultValue(new ArrayBuffer[(Int, Int, Boolean, Double)])
      // sequence type -> list((motif id, motif num, evalue))
      val seqtypeMap = new HashMap[String, ArrayBuffer[(Int, Int, Double)]]
      // motif id -> cluster
      val clusterMap = new HashMap[Int, Int]

      val motifInfos = ResultQueries.motifInfosForIteration(iteration)
      val pssmRows = ResultQueries.pssmRowsForIteration(iteration)
      val annotations = ResultQueries.motifAnnotationsForIteration(iteration)

      motifInfos.foreach { info =>
        if (!seqtypeMap.contains(info._3))
          seqtypeMap(info._3) = new ArrayBuffer[(Int, Int, Double)]
        seqtypeMap(info._3) += ((info._1, info._4, info._5))
        clusterMap(info._1) = info._2
      }

      pssmRows.foreach { pssmRow =>
        if (!pssmMap.contains(pssmRow._1)) pssmMap(pssmRow._1) = new ArrayBuffer[Array[Float]]
        pssmMap(pssmRow._1) += Array(pssmRow._2, pssmRow._3, pssmRow._4, pssmRow._5)
      }

      annotations.foreach { annot =>
        if (!annotMap.contains(annot._1)) annotMap(annot._1) =
          new ArrayBuffer[(Int, Int, Boolean, Double)]
        annotMap(annot._1) += ((annot._2, annot._3, annot._4, annot._5))
      }
      val tmpMap = new HashMap[String, HashMap[Int, ArrayBuffer[MotifInfo]]]
      for {
        (seqtype, motifTuples) <- seqtypeMap
         motifTuple <- motifTuples
      } {
        if (!tmpMap.contains(seqtype))
          tmpMap(seqtype) = new HashMap[Int, ArrayBuffer[MotifInfo]]

        val cluster = clusterMap(motifTuple._1)
        if (!tmpMap(seqtype).contains(cluster)) {
          tmpMap(seqtype)(cluster) = new ArrayBuffer[MotifInfo]
        }

        tmpMap(seqtype)(cluster) += MotifInfo(motifTuple._2, motifTuple._3,
                                              pssmMap(motifTuple._1).toArray,
                                              annotMap(motifTuple._1).map { annotTuple =>
                                                Annotation(motifTuple._2,
                                                           annotTuple._3,
                                                           annotTuple._2,
                                                           rowNames(annotTuple._1),
                                                           annotTuple._4)
                                             }.toArray)
      }
      tmpMap
    }
}

// **********************************************************************
// **** cMonkey run information
// **********************************************************************

/**
 * Note that SQLite does not have a timestamp or date type. We only read
 * that value anyways, so String is fine for now
 */
object RunInfos
extends Table[(String, Option[String], Int, Option[Int],
               String, String, Int, Int, Int)]("run_infos") {

  lazy val database = Database.forDataSource(DatabaseConfig.datasource)

  def startTime = column[String]("start_time", O NotNull)
  def finishTime = column[Option[String]]("finish_time")
  def numIterations = column[Int]("num_iterations", O NotNull)
  def lastIteration = column[Option[Int]]("last_iteration")
  def organism = column[String]("organism", O NotNull)
  def species = column[String]("species", O NotNull)
  def numRows = column[Int]("num_rows", O NotNull)
  def numColumns = column[Int]("num_columns", O NotNull)
  def numClusters = column[Int]("num_clusters", O NotNull)

  def * = startTime ~ finishTime ~ numIterations ~ lastIteration ~ organism ~ species ~ numRows ~ numColumns ~ numClusters
  def findAll = database.withSession { implicit db: Session =>
    (for {
      t <- this
     } yield t.startTime ~ t.finishTime ~ t.numIterations ~ t.lastIteration ~ t.organism ~ t.species ~
       t.numRows ~ t.numColumns ~ t.numClusters).list
  }
}

// **********************************************************************
// **** cMonkey run statistics
// **********************************************************************

object IterationStats extends Table[(Int, Double, Double)]("iteration_stats") {

  def iteration = column[Int]("iteration", O NotNull)
  def medianResidual = column[Double]("median_residual", O NotNull)
  def fuzzyCoeff = column[Double]("fuzzy_coeff", O NotNull)

  def * = iteration ~ medianResidual ~ fuzzyCoeff
}

object ClusterStats extends Table[(Int, Int, Int, Int, Double)]("cluster_stats") {

  def iteration = column[Int]("iteration", O NotNull)
  def cluster = column[Int]("cluster", O NotNull)
  def numRows = column[Int]("num_rows", O NotNull)
  def numColumns = column[Int]("num_cols", O NotNull)
  def residual = column[Double]("residual", O NotNull)

  def * = iteration ~ cluster ~ numRows ~ numColumns ~ residual
}

object NetworkStats extends Table[(Int, String, Double)]("network_stats") {

  def iteration = column[Int]("iteration", O NotNull)
  def network = column[String]("network", O NotNull)
  def score = column[Double]("score", O NotNull)

  def * = iteration ~ network ~ score
}

object MotifStats extends Table[(Int, String, Double)]("motif_stats") {

  def iteration = column[Int]("iteration", O NotNull)
  def seqtype = column[String]("seqtype", O NotNull)
  def pval = column[Double]("pval", O NotNull)

  def * = iteration ~ seqtype ~ pval
}

object StatsQueries {
  lazy val database = Database.forDataSource(DatabaseConfig.datasource)
  def iterationStats = database.withSession { implicit db: Session =>
    Query(IterationStats).sortBy(_.iteration).list
  }
  def clusterStats = database.withSession { implicit db: Session =>
    Query(ClusterStats).sortBy(_.cluster).list
  }

  def networkStats = database.withSession { implicit db: Session =>
    Query(NetworkStats).sortBy(_.network).list
  }


  def motifStats = database.withSession { implicit db: Session =>
    Query(MotifStats).sortBy(_.seqtype).list
  }
}


// **********************************************************************
// **********************************************************************
// **********************************************************************
// ***** Access functionality
// **********************************************************************
// **********************************************************************
// **********************************************************************

// **********************************************************************
// **** Result data
// **********************************************************************

case class Annotation(motifNum: Int, reverse: Boolean, position: Int, gene: String, pvalue: Double)
case class GeneAnnotations(gene: String, annotations: Seq[Annotation])
case class MotifInfo(motifNum: Int, evalue: Double, pssm: Array[Array[Float]], annotations: Array[Annotation])
case class IterationResult(rows: Map[Int, List[String]], columns: Map[Int, List[String]],
                           residuals: Map[Int, Double],
                           motifs: Map[String, Map[Int, Seq[MotifInfo]]]) {
}

// **********************************************************************
// **** Result data reader
// **********************************************************************

class IterationResultReader(Synonyms: SynonymsMap) {

  def readIterationResult(iteration: Int) : Option[IterationResult] = {
    val rowNames: Array[String] = ResultQueries.allRowNames
    val columnNames: Array[String] = ResultQueries.allColumnNames
    val rowMembers = ResultQueries.rowMembersForIteration(iteration)
    val colMembers = ResultQueries.columnMembersForIteration(iteration)
    val clusterResiduals = ResultQueries.clusterResidualsForIteration(iteration)

    val rows = new HashMap[Int, List[String]]
    val columns = new HashMap[Int, List[String]]
    val residuals = new HashMap[Int, Double]

    rowMembers.foreach { member =>
      member match {
        case (cluster, orderNum) =>
          if (!rows.contains(cluster)) rows(cluster) = Nil
          rows(cluster) ::= rowNames(orderNum)
      }
    }
    colMembers.foreach { member =>
      member match {
        case (cluster, orderNum) =>
          if (!columns.contains(cluster)) columns(cluster) = Nil
          columns(cluster) ::= columnNames(orderNum)
      }
    }
    clusterResiduals.foreach { clusterResidual =>
      clusterResidual match {
        case (cluster, residual) => residuals(cluster) = residual
      }
    }
    val motifInfos = MotifQueries.findForIteration(iteration, rowNames)

    Some(IterationResult(rows.toMap.withDefaultValue(List()),
                         columns.toMap.withDefaultValue(List()),
                         residuals.toMap.withDefaultValue(0.0),
                         motifInfos.toMap))
  }
}

// **********************************************************************
// **** cMonkey run status information
// **********************************************************************

case class RunStatus(startTime: String, finishTime: Option[String],
                     numIterations: Int, lastIteration: Option[Int],
                     organismCode: String,
                     species: String, numRows: Int, numColumns: Int,
                     numClusters: Int) {
  def finished = finishTime != null
}

object RunStatusReader {

  def readRunStatus: Option[RunStatus] = {
    val runInfos = RunInfos.findAll
    if (runInfos.length > 0) {
      val info = runInfos(0)
      Some(RunStatus(info._1, info._2, info._3, info._4, info._5, info._6, info._7,
                     info._8, info._9))
    } else None
  }
}

// **********************************************************************
// **** cMonkey run log information
// **********************************************************************

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

// **********************************************************************
// **** cMonkey statistics
// **********************************************************************

case class ClusterStat(numRows: Int, numColumns: Int, residual: Double)
case class IterationStat(clusters: Map[Int, ClusterStat], medianResidual: Double,
                         motifPValues: Map[String, Double],
                         networkScores: Map[String, Double],
                         fuzzyCoeff: Double)

class StatsReader {
  def readStats: HashMap[Int, IterationStat] = {

    // for I/O performance reasons, we use the findAll method to build up our stats
    val stats = new HashMap[Int, IterationStat]
    val iterationStats = StatsQueries.iterationStats

    val clusterStats = StatsQueries.clusterStats
    val clusterIterationMap = new HashMap[Int, HashMap[Int, ClusterStat]]
    clusterStats.foreach { cstat =>
      if (!clusterIterationMap.contains(cstat._1)) {
        clusterIterationMap(cstat._1) = new HashMap[Int, ClusterStat]
      }
      clusterIterationMap(cstat._1)(cstat._2) = ClusterStat(cstat._3, cstat._4, cstat._5)
    }

    val motifStats = StatsQueries.motifStats
    val motifIterationMap = (new HashMap[Int, HashMap[String, Double]]).withDefaultValue(new HashMap[String, Double])
    motifStats.foreach { mstat =>
      if (!motifIterationMap.contains(mstat._1)) {
        motifIterationMap(mstat._1) = new HashMap[String, Double]
      }
      motifIterationMap(mstat._1)(mstat._2) = mstat._3
    }

    val networkStats = StatsQueries.networkStats
    val networkIterationMap = (new HashMap[Int, HashMap[String, Double]]).withDefaultValue(new HashMap[String, Double])
    networkStats.foreach { nstat =>
      if (!networkIterationMap.contains(nstat._1)) {
        networkIterationMap(nstat._1) = new HashMap[String, Double]
      }
      networkIterationMap(nstat._1)(nstat._2) = nstat._3
    }


    iterationStats.foreach { stat =>
      val iteration = stat._1
      val networkMap = networkIterationMap(iteration)
      val motifMap = motifIterationMap(iteration)
      val clusterMap = clusterIterationMap(iteration)
      stats(iteration) = IterationStat(clusterMap.toMap, stat._2, motifMap.toMap,
                                       networkMap.toMap, stat._3)
    }
    stats
  }
}
