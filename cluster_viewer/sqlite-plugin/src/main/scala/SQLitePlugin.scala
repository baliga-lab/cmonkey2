package org.systemsbiology

import play.api.{Plugin, Application}
import javax.sql.DataSource
import java.sql.DriverManager
import java.io.{PrintWriter, File}

/**
 * First shot at a SQLite datasource plugin that can be configured at runtime,
 * after the application has been started.
 * This is my first Play 2 plugin.
 */
class SQLitePlugin(app: Application) extends Plugin {

  override def onStart {
    // automatic determination of output path
    println("Starting SQLite plugin, run from: " + System.getProperty("user.dir"))
    val dbfile = new File((new File(System.getProperty("user.dir"))).getParentFile,
                          "out/cmonkey_run.db")
    DatabaseConfig.url = "jdbc:sqlite:%s".format(dbfile.getAbsolutePath)
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

  var url: String = ""
  private var _datasource: DataSource = null

  def datasource: DataSource = {
    if (_datasource == null) _datasource = new SQLiteDataSource(url)
    _datasource
  }
}
