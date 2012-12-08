package controllers

import play.api.Play.current
import java.io._
import java.sql.Timestamp

import play.api.db._
import org.scalaquery.ql._
import org.scalaquery.ql.TypeMapper._
import org.scalaquery.ql.extended.{ExtendedTable => Table}
import org.scalaquery.ql.extended.SQLiteDriver.Implicit._
import org.scalaquery.session.{Database, Session}

object RunInfos
extends Table[(Timestamp, Option[Timestamp], Int, Option[Int],
               String, String, Int, Int)]("run_infos") {

  lazy val database = Database.forDataSource(DB.getDataSource())
  def startTime = column[Timestamp]("start_time", O NotNull)
  def finishTime = column[Option[Timestamp]]("finish_time")
  def numIterations = column[Int]("num_iterations", O NotNull)
  def lastIteration = column[Option[Int]]("last_iteration")
  def organism = column[String]("organism", O NotNull)
  def species = column[String]("species", O NotNull)
  def numRows = column[Int]("num_rows", O NotNull)
  def numColumns = column[Int]("num_columns", O NotNull)

  def * = startTime ~ finishTime ~ numIterations ~ lastIteration ~ organism ~ species ~ numRows ~ numColumns
  def findAll = database.withSession { implicit db: Session =>
    (for {
      t <- this
     } yield t.startTime ~ t.finishTime ~ t.numIterations ~ t.lastIteration ~ t.organism ~ t.species ~
       t.numRows ~ t.numColumns).list
  }
}

case class RunStatus(startTime: Timestamp, finishTime: Option[Timestamp],
                     numIterations: Int, lastIteration: Option[Int],
                     organismCode: String,
                     species: String, numRows: Int, numColumns: Int) {
  def finished = finishTime != null
}

class RunStatusReader {

  def readRunStatus: Option[RunStatus] = {
    val runInfos = RunInfos.findAll
    if (runInfos.length > 0) {
      val info = runInfos(0)
      Some(RunStatus(info._1, info._2, info._3, info._4, info._5, info._6, info._7, info._8))
    } else None
  }
}
