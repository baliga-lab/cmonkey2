package org.systemsbiology

import play.api.{Plugin, Application}
import javax.sql.DataSource
import java.sql.DriverManager
import java.io.PrintWriter

/**
 * First shot at a SQLite datasource plugin that can be configured at runtime,
 * after the application has been started.
 * This is my first Play 2 plugin.
 */
class SQLitePlugin(app: Application) extends Plugin {

  override def onStart {
    println("Starting SQLite plugin")
    Class.forName("org.sqlite.JDBC")
  }
}

class SQLiteDataSource(url: String) extends DataSource {
  def getConnection = {
    DriverManager.getConnection(url)
  }
  def getConnection(username: String, password: String) = {
    DriverManager.getConnection(url)
  }
  def getLoginTimeout = {
    println("getLoginTimeout()")
    100
  }
  def getLogWriter = {
    println("getLogWriter()")
    null
  }
  def getParentLogger = {
    println("getParentLogger()")
    null
  }
  def setLoginTimeout(seconds: Int) {
    println("setLoginTimeout()")
  }
  def setLogWriter(out: PrintWriter) {
    println("setLogWriter()")
  }
  def isWrapperFor(c: Class[_]) = {
    println("isWrapperFor()")
    false
  }
  def unwrap[T](x: Class[T]) = {
    println("unwrap()")
    x.newInstance
  }
}

object DatabaseConfig {
  var datasource: DataSource = null

  def createDataSource(url: String) {
    if (datasource == null) {
      datasource = new SQLiteDataSource(url)
    }
  }
}
