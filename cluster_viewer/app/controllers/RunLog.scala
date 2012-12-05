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

object RunLogs extends Table[(Long, String)]("run_logs") {
  lazy val database = Database.forDataSource(DB.getDataSource())
  def id = column[Long]("rowid", O PrimaryKey,  O AutoInc, O NotNull)
  def name = column[String]("name", O NotNull)
  def * = id ~ name
  def findAll = database.withSession { implicit db: Session =>
    (for {
      t <- this
      _ <- Query orderBy(t.name)
     } yield t.id ~ t.name).list
  }
}

object LogEntries extends Table[(Long, Long, Long, Boolean, Float)]("log_entries") {
  lazy val database = Database.forDataSource(DB.getDataSource())
  def id = column[Long]("rowid", O PrimaryKey,  O AutoInc, O NotNull)
  def logId = column[Long]("log_id", O NotNull)
  def iteration = column[Long]("iteration", O NotNull)
  def wasActive = column[Boolean]("was_active", O NotNull)
  def scaling = column[Float]("scaling", O NotNull)

  def * = id ~ logId ~ iteration ~ wasActive ~ scaling

  def findAll = database.withSession { implicit db: Session =>
    (for {
      t <- this
      _ <- Query orderBy(t.iteration)
     } yield t.id ~ t.logId ~ t.iteration ~ t.wasActive ~ t.scaling).list
  }

  def entriesFor(runLogId: Long) = database.withSession { implicit db: Session =>
    (for {
      logEntry <- LogEntries if logEntry.logId === runLogId
      _ <- Query.orderBy(logEntry.iteration)
    } yield logEntry.wasActive ~ logEntry.scaling).list
  }
}

case class RunLog(functionName: String, scaling: Array[Float],
                  active: Array[Boolean]) {
  override def toString = {
    "RunLog('%s', %d scaling)".format(functionName, scaling.length)
  }
}

class RunLogReader {
  def readLogs(OutDirectory: File): Option[Array[RunLog]] = {
    val result = RunLogs.findAll.map { runLog =>
      val entries = LogEntries.entriesFor(runLog._1).unzip
      val active = entries._1.toArray
      val scaling = entries._2.toArray
      RunLog(runLog._2, scaling, active)
    }.toArray
    Some(result)
  }
}
