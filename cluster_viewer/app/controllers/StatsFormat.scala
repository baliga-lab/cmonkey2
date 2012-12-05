package controllers

import play.api.Play.current

import java.io._
import java.util.regex._
import scala.collection.mutable.HashMap

import play.api.db._
import org.scalaquery.ql._
import org.scalaquery.ql.TypeMapper._
import org.scalaquery.ql.extended.{ExtendedTable => Table}
import org.scalaquery.ql.extended.SQLiteDriver.Implicit._
import org.scalaquery.session.{Database, Session}

object IterationStats extends Table[(Int, Double, Double)]("iteration_stats") {

  lazy val database = Database.forDataSource(DB.getDataSource())
  def iteration = column[Int]("iteration", O NotNull)
  def medianResidual = column[Double]("median_residual", O NotNull)
  def fuzzyCoeff = column[Double]("fuzzy_coeff", O NotNull)

  def * = iteration ~ medianResidual ~ fuzzyCoeff
  def findAll = database.withSession { implicit db: Session =>
    (for {
      t <- this
      _ <- Query orderBy(t.iteration)
     } yield t.iteration ~ t.medianResidual ~ fuzzyCoeff).list
  }
}

object ClusterStats extends Table[(Int, Int, Int, Int, Double)]("cluster_stats") {

  lazy val database = Database.forDataSource(DB.getDataSource())
  def iteration = column[Int]("iteration", O NotNull)
  def cluster = column[Int]("cluster", O NotNull)
  def numRows = column[Int]("num_rows", O NotNull)
  def numColumns = column[Int]("num_cols", O NotNull)
  def residual = column[Double]("residual", O NotNull)

  def * = iteration ~ cluster ~ numRows ~ numColumns ~ residual
  def findForIteration(iter: Int) = database.withSession { implicit db: Session =>
    (for {
      t <- this if t.iteration === iter
      _ <- Query orderBy(t.cluster)
     } yield t.cluster ~ t.numRows ~ t.numColumns ~ t.residual).list
  }
}

object NetworkStats extends Table[(Int, String, Double)]("network_stats") {

  lazy val database = Database.forDataSource(DB.getDataSource())
  def iteration = column[Int]("iteration", O NotNull)
  def network = column[String]("network", O NotNull)
  def score = column[Double]("score", O NotNull)

  def * = iteration ~ network ~ score
  def findForIteration(iter: Int) = database.withSession { implicit db: Session =>
    (for {
      t <- this if t.iteration === iter
      _ <- Query orderBy(t.network)
     } yield t.network ~ t.score).list
  }
}

object MotifStats extends Table[(Int, String, Double)]("motif_stats") {

  lazy val database = Database.forDataSource(DB.getDataSource())
  def iteration = column[Int]("iteration", O NotNull)
  def seqtype = column[String]("seqtype", O NotNull)
  def pval = column[Double]("pval", O NotNull)

  def * = iteration ~ seqtype ~ pval
  def findForIteration(iter: Int) = database.withSession { implicit db: Session =>
    (for {
      t <- this if t.iteration === iter
      _ <- Query orderBy(t.seqtype)
     } yield t.seqtype ~ t.pval).list
  }
}

case class ClusterStat(numRows: Int, numColumns: Int, residual: Double)
case class IterationStat(clusters: Map[Int, ClusterStat], medianResidual: Double,
                         motifPValues: Map[String, Double],
                         networkScores: Map[String, Double],
                         fuzzyCoeff: Double)

class StatsReader {
  def readStats: HashMap[Int, IterationStat] = {
    val stats = new HashMap[Int, IterationStat]
    IterationStats.findAll.foreach { stat =>
      val iteration = stat._1
      val clusterMap = new HashMap[Int, ClusterStat]()
      val networkMap = new HashMap[String, Double]()
      val motifMap = new HashMap[String, Double]()
      ClusterStats.findForIteration(iteration).foreach { cstat =>
        clusterMap(cstat._1) = ClusterStat(cstat._2, cstat._3, cstat._4)
      }
      MotifStats.findForIteration(iteration).foreach { mstat =>
        motifMap(mstat._1) = mstat._2
      }
      NetworkStats.findForIteration(iteration).foreach { nstat =>
        networkMap(nstat._1) = nstat._2
      }
      stats(iteration) = IterationStat(clusterMap.toMap, stat._2, motifMap.toMap,
                                       networkMap.toMap, stat._3)
    }
    stats
  }
}
