package controllers

import play.api.Play.current

import java.io._
import java.util.zip._
import java.util.regex._
import scala.collection.mutable.HashMap

import play.api.db._
import org.scalaquery.ql._
import org.scalaquery.ql.TypeMapper._
import org.scalaquery.ql.extended.{ExtendedTable => Table}
import org.scalaquery.ql.extended.SQLiteDriver.Implicit._
import org.scalaquery.session.{Database, Session}

object RowNames extends Table[(Int, String)]("row_names") {

  lazy val database = Database.forDataSource(DB.getDataSource())
  def orderNum = column[Int]("order_num", O NotNull)
  def name = column[String]("name", O NotNull)

  def * = orderNum ~ name
  def findAll = database.withSession { implicit db: Session =>
    (for {
      t <- this
      _ <- Query orderBy(t.name)
     } yield t.name).list.toArray
  }
}

object ColumnNames extends Table[(Int, String)]("column_names") {

  lazy val database = Database.forDataSource(DB.getDataSource())
  def orderNum = column[Int]("order_num", O NotNull)
  def name = column[String]("name", O NotNull)

  def * = orderNum ~ name
  def findAll = database.withSession { implicit db: Session =>
    (for {
      t <- this
      _ <- Query orderBy(t.name)
     } yield t.name).list.toArray
  }
}

object RowMembers extends Table[(Int, Int, Int)]("row_members") {

  lazy val database = Database.forDataSource(DB.getDataSource())
  def iteration = column[Int]("iteration", O NotNull)
  def cluster = column[Int]("cluster", O NotNull)
  def orderNum = column[Int]("order_num", O NotNull)

  def * = iteration ~ cluster ~ orderNum
  def findForIteration(iteration: Int) = database.withSession { implicit db: Session =>
    (for {
      t <- this if t.iteration === iteration
      _ <- Query orderBy(t.cluster)
     } yield t.cluster ~ t.orderNum).list
  }
}
object ColumnMembers extends Table[(Int, Int, Int)]("column_members") {

  lazy val database = Database.forDataSource(DB.getDataSource())
  def iteration = column[Int]("iteration", O NotNull)
  def cluster = column[Int]("cluster", O NotNull)
  def orderNum = column[Int]("order_num", O NotNull)

  def * = iteration ~ cluster ~ orderNum
  def findForIteration(iteration: Int) = database.withSession { implicit db: Session =>
    (for {
      t <- this if t.iteration === iteration
      _ <- Query orderBy(t.cluster)
     } yield t.cluster ~ t.orderNum).list
  }
}
object ClusterResiduals extends Table[(Int, Int, Double)]("cluster_residuals") {

  lazy val database = Database.forDataSource(DB.getDataSource())
  def iteration = column[Int]("iteration", O NotNull)
  def cluster = column[Int]("cluster", O NotNull)
  def residual = column[Double]("residual", O NotNull)

  def * = iteration ~ cluster ~ residual
  def findForIteration(iteration: Int) = database.withSession { implicit db: Session =>
    (for {
      t <- this if t.iteration === iteration
      _ <- Query orderBy(t.cluster)
     } yield t.cluster ~ t.residual).list
  }
}

case class Annotation(motifNum: Int, reverse: Boolean, position: Int, gene: String, pvalue: Double)
case class GeneAnnotations(gene: String, annotations: Seq[Annotation])
case class MotifInfo(motifNum: Int, evalue: Double, pssm: Array[Array[Float]], annotations: Array[Annotation])
case class Snapshot(rows: Map[Int, List[String]], columns: Map[Int, List[String]],
                    residuals: Map[Int, Double],
                    motifs: Map[String, Map[Int, Array[MotifInfo]]]) {
}

class SnapshotReader(Synonyms: SynonymsMap) {

  def readSnapshot(iteration: Int) : Option[Snapshot] = {
    val rowNames: Array[String] = RowNames.findAll
    val columnNames: Array[String] = ColumnNames.findAll
    printf("# row names: %d, # col names: %d\n", rowNames.length, columnNames.length)
    val rowMembers = RowMembers.findForIteration(iteration)
    val colMembers = ColumnMembers.findForIteration(iteration)
    val clusterResiduals = ClusterResiduals.findForIteration(iteration)

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

    printf("# row members: %d, # col members: %d, # residuals: %d\n",
           rowMembers.length, colMembers.length, clusterResiduals.length)
    Some(Snapshot(rows.toMap, columns.toMap, residuals.toMap,
                  Map[String, Map[Int, Array[MotifInfo]]]()))
  }
}
