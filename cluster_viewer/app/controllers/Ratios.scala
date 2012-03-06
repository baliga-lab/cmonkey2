package controllers

import java.io._

case class RatioMatrix(rows: Array[String], columns: Array[String],
                       values: Array[Array[Float]])

class RatioMatrixFactory(ratiofile: File, synonyms: SynonymsMap) {

  var columnTitles: Array[String] = null
  var rowTitles: Array[String] = null
  val values = new java.util.ArrayList[Array[Float]]
  readFile

  def readFile {
    var in: BufferedReader = null
    try {
      in = new BufferedReader(new FileReader(ratiofile))
      // First line
      var line = in.readLine
      columnTitles = line.split("\t").tail.toArray
      val rows = new java.util.ArrayList[String]

      while (line != null) {
        line = in.readLine
        if (line != null) {
          val comps = line.split("\t")
          rows.add(synonyms(comps(0).replaceAll("\"", "").trim))
          values.add(comps.tail.map(str => str.toFloat).toArray)
        }
      }
      rowTitles = rows.toArray(new Array[String](0))
    } finally {
      if (in != null) in.close
    }
  }

  def readRatios(genes: Array[String]): RatioMatrix = {
    var i = 0
    val outValues = new java.util.ArrayList[Array[Float]]
    while (i < rowTitles.length) {
      if (genes.contains(rowTitles(i))) {
        outValues.add(values.get(i))
      }
      i += 1
    }
    RatioMatrix(genes, columnTitles,
                outValues.toArray(new Array[Array[Float]](0)))
  }
}
