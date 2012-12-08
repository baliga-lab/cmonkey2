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
  def findAll = database.withSession { implicit db: Session =>
    (for {
      t <- this
      _ <- Query orderBy(t.cluster)
     } yield t.iteration ~ t.cluster ~ t.numRows ~ t.numColumns ~ t.residual).list
  }
}

object NetworkStats extends Table[(Int, String, Double)]("network_stats") {

  lazy val database = Database.forDataSource(DB.getDataSource())
  def iteration = column[Int]("iteration", O NotNull)
  def network = column[String]("network", O NotNull)
  def score = column[Double]("score", O NotNull)

  def * = iteration ~ network ~ score
  def findAll = database.withSession { implicit db: Session =>
    (for {
      t <- this
      _ <- Query orderBy(t.network)
     } yield t.iteration ~ t.network ~ t.score).list
  }
}

object MotifStats extends Table[(Int, String, Double)]("motif_stats") {

  lazy val database = Database.forDataSource(DB.getDataSource())
  def iteration = column[Int]("iteration", O NotNull)
  def seqtype = column[String]("seqtype", O NotNull)
  def pval = column[Double]("pval", O NotNull)

  def * = iteration ~ seqtype ~ pval

  def findAll = database.withSession { implicit db: Session =>
    (for {
      t <- this
      _ <- Query orderBy(t.seqtype)
     } yield t.iteration ~ t.seqtype ~ t.pval).list
  }
}

case class ClusterStat(numRows: Int, numColumns: Int, residual: Double)
case class IterationStat(clusters: Map[Int, ClusterStat], medianResidual: Double,
                         motifPValues: Map[String, Double],
                         networkScores: Map[String, Double],
                         fuzzyCoeff: Double)

class StatsReader {
  def readStats: HashMap[Int, IterationStat] = {

    // for I/O performance reasons, we use the findAll method to build up our stats
    val stats = new HashMap[Int, IterationStat]
    val iterationStats = IterationStats.findAll

    val clusterStats = ClusterStats.findAll
    val clusterIterationMap = new HashMap[Int, HashMap[Int, ClusterStat]]
    clusterStats.foreach { cstat =>
      if (!clusterIterationMap.contains(cstat._1)) {
        clusterIterationMap(cstat._1) = new HashMap[Int, ClusterStat]
      }
      clusterIterationMap(cstat._1)(cstat._2) = ClusterStat(cstat._3, cstat._4, cstat._5)
    }

    val motifStats = MotifStats.findAll
    val motifIterationMap = new HashMap[Int, HashMap[String, Double]]
    motifStats.foreach { mstat =>
      if (!motifIterationMap.contains(mstat._1)) {
        motifIterationMap(mstat._1) = new HashMap[String, Double]
      }
      motifIterationMap(mstat._1)(mstat._2) = mstat._3
    }

    val networkStats = NetworkStats.findAll
    val networkIterationMap = new HashMap[Int, HashMap[String, Double]]
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
