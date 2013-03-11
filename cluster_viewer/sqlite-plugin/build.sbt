name         := "sqlite-plugin"

version      := "1.0"

organization := "org.systemsbiology"

scalaVersion := "2.10.0"

scalacOptions ++= Seq("-unchecked", "-deprecation")

resolvers += "official Maven mirror" at "http://mirrors.ibiblio.org/pub/mirrors/maven2/"

resolvers += "official Maven mirror" at "http://repo.typesafe.com/typesafe/releases/"

libraryDependencies ++= Seq("play" % "play_2.10" % "2.1.0",
                            "org.scalatest" % "scalatest_2.10" % "2.0.M5b",
                            "junit" % "junit" % "4.10")
