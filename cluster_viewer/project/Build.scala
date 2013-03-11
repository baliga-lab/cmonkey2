import sbt._
import Keys._
import play.Project._

object ApplicationBuild extends Build {

  val appName         = "clusters"
  val appVersion      = "1.0-SNAPSHOT"

  val appDependencies = Seq(
    // Add your project dependencies here,
    jdbc,
    anorm, 
    "com.typesafe.slick" %% "slick" % "1.0.0"
  )
  override lazy val settings = super.settings ++ clusterViewerSettings

  def clusterViewerSettings = Seq(
    organization := "org.systemsbiology",
    version := "1.0",
    scalaVersion := "2.10.0",
    resolvers += "official Maven mirror" at "http://repo.typesafe.com/typesafe/releases/"
  )

  def pluginDependencies = Seq(libraryDependencies ++= Seq("play" % "play_2.10" % "2.1.0"))

  val main = play.Project(appName, appVersion, appDependencies).settings(
    scalacOptions ++= Seq("-feature", "-deprecation", "-unchecked")
  ) dependsOn(sqlitePlugin)

  lazy val sqlitePlugin = Project("sqlite-plugin", file("sqlite-plugin")) settings(pluginDependencies :_*)
}
