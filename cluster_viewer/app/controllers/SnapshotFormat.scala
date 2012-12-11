package controllers

import play.api.Play.current

import java.io._
import java.util.zip._
import java.util.regex._
import scala.collection.Map
import scala.collection.mutable.{HashMap, ArrayBuffer}

import play.api.db._
import org.scalaquery.ql._
import org.scalaquery.ql.TypeMapper._
import org.scalaquery.ql.extended.{ExtendedTable => Table}
import org.scalaquery.ql.extended.SQLiteDriver.Implicit._
import org.scalaquery.session.{Database, Session}

// **********************************************************************
// **** Tables for results
// **********************************************************************
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

object MotifInfos extends Table[(Int, Int, Int, String, Int, Double)]("motif_infos") {

  lazy val database = Database.forDataSource(DB.getDataSource())
  def id = column[Int]("rowid", O NotNull)
  def iteration = column[Int]("iteration", O NotNull)
  def cluster = column[Int]("cluster", O NotNull)
  def seqtype = column[String]("seqtype", O NotNull)
  def motifNum = column[Int]("motif_num", O NotNull)
  def evalue = column[Double]("evalue", O NotNull)

  def * = id ~ iteration ~ cluster ~ seqtype ~ motifNum ~ evalue
  def findForIteration(iteration: Int) = database.withSession { implicit db: Session =>
    (for {
      t <- this if t.iteration === iteration
      _ <- Query orderBy(t.cluster)
     } yield t.id ~ t.cluster ~ t.seqtype ~ t.motifNum ~ t.evalue).list
  }
}

object MotifPSSMRows
extends Table[(Int, Int, Int, Float, Float, Float, Float)]("motif_pssm_rows") {

  lazy val database = Database.forDataSource(DB.getDataSource())
  def motifInfoId = column[Int]("motif_info_id", O NotNull)
  def iteration = column[Int]("iteration", O NotNull)
  def row = column[Int]("row", O NotNull)
  def a = column[Float]("a", O NotNull)
  def c = column[Float]("c", O NotNull)
  def g = column[Float]("g", O NotNull)
  def t = column[Float]("t", O NotNull)

  def * = motifInfoId ~ iteration ~ row ~ a ~ c ~ g ~ t
  def findForIteration(iteration: Int) = database.withSession { implicit db: Session =>
    (for {
      t <- this if t.iteration === iteration
      _ <- Query orderBy(t.row)
     } yield t.motifInfoId ~ t.a ~ t.c ~ t.g ~ t.t).list
  }
}

object MotifAnnotations
extends Table[(Int, Int, Int, Int, Boolean, Double)]("motif_annotations") {

  lazy val database = Database.forDataSource(DB.getDataSource())
  def motifInfoId = column[Int]("motif_info_id", O NotNull)
  def iteration = column[Int]("iteration", O NotNull)
  def geneNum = column[Int]("gene_num", O NotNull)
  def position = column[Int]("position", O NotNull)
  def reverse = column[Boolean]("reverse", O NotNull)
  def pvalue = column[Double]("pvalue", O NotNull)

  def * = motifInfoId ~ iteration ~ geneNum ~ position ~ reverse ~ pvalue
  def findForIteration(iteration: Int) = database.withSession { implicit db: Session =>
    (for {
      t <- this if t.iteration === iteration
      _ <- Query orderBy(t.geneNum)
     } yield t.motifInfoId ~ t.geneNum ~ t.position ~ t.reverse ~ t.pvalue).list
  }
}

object MotifQueries {
  
  lazy val database = Database.forDataSource(DB.getDataSource())
  def findForIteration(iteration: Int, rowNames: Array[String]) =
    // the result structure, a map from seqtype->cluster->motif info
    database.withSession { implicit db: Session =>

      // motif id -> pssm row
      val pssmMap = new HashMap[Int, ArrayBuffer[Array[Float]]]
      // motif id -> annotation
      val annotMap = new HashMap[Int, ArrayBuffer[(Int, Int, Boolean, Double)]]
      // sequence type -> list((motif id, motif num, evalue))
      val seqtypeMap = new HashMap[String, ArrayBuffer[(Int, Int, Double)]]
      // motif id -> cluster
      val clusterMap = new HashMap[Int, Int]

      val motifInfos = MotifInfos.findForIteration(iteration)
      val pssmRows = MotifPSSMRows.findForIteration(iteration)
      val annotations = MotifAnnotations.findForIteration(iteration)

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
// **** Result data
// **********************************************************************

case class Annotation(motifNum: Int, reverse: Boolean, position: Int, gene: String, pvalue: Double)
case class GeneAnnotations(gene: String, annotations: Seq[Annotation])
case class MotifInfo(motifNum: Int, evalue: Double, pssm: Array[Array[Float]], annotations: Array[Annotation])
case class Snapshot(rows: Map[Int, List[String]], columns: Map[Int, List[String]],
                    residuals: Map[Int, Double],
                    motifs: Map[String, Map[Int, Seq[MotifInfo]]]) {
}

// **********************************************************************
// **** Result data reader
// **********************************************************************

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
    val motifInfos = MotifQueries.findForIteration(iteration, rowNames)

    printf("# row members: %d, # col members: %d, # residuals: %d\n",
           rowMembers.length, colMembers.length, clusterResiduals.length)
    Some(Snapshot(rows.toMap, columns.toMap, residuals.toMap,
                  motifInfos.toMap))
  }
}
