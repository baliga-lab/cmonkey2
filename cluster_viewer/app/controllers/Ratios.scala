package controllers

import java.io._
import java.util.zip._

case class RatioMatrix(rows: Array[String], columns: Array[String],
                       values: Array[Array[Float]]) {
  def mean = {
    val numValues: Float = rows.length * columns.length
    var sum : Float = 0.0f
    for (i <- 0 until rows.length; j <- 0 until columns.length) sum += values(i)(j)
    sum / numValues
  }
/*
  def stdDev = {
    val numValues: Float = rows.length * columns.length
    var sum : Float = 0.0f
    val avg = mean
    val factor: Float = 1.0f / (numValues - 1)
    for (i <- 0 until rows.length; j <- 0 until columns.length) {
      val x = values(i)(j)
      sum += (x - avg) * (x - avg)
    }
    sum *= factor
    math.sqrt(sum)
  }
*/

  private def stdDev(vals: Seq[Float]) = {
    val numValues: Float = vals.length
    var sum : Float = 0.0f
    val mean = vals.sum / numValues
    val factor: Float = 1.0f / (numValues - 1)
    for (x <- vals) sum += (x - mean) * (x - mean)
    sum *= factor
    math.sqrt(sum)
  }

  def colStdDev(column: Int) = {
    val tmp = for (row <- 0 until rows.length) yield values(row)(column)
    stdDev(tmp)
  }
}

class RatioMatrixFactory(ratiofile: File, synonyms: SynonymsMap) {

  var columnTitles: Array[String] = null
  var rowTitles: Array[String] = null
  val values = new java.util.ArrayList[Array[Float]]
  readFile

  private def readFile {
    var in: BufferedReader = null
    try {
      if (ratiofile.getName.endsWith(".gz")) {
        in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(ratiofile))))
      } else {
        in = new BufferedReader(new FileReader(ratiofile))
      }
      read(in)
    } finally {
      if (in != null) in.close
    }
  }

  private def read(in: BufferedReader) {
    // First line
    var line = in.readLine
    columnTitles = line.split("\t").tail.toArray
    val rows = new java.util.ArrayList[String]

    while (line != null) {
      line = in.readLine
      if (line != null) {
        val comps = line.split("\t")
        rows.add(synonyms(comps(0).replaceAll("\"", "").trim))
        values.add(comps.tail.map(str => if (str != "nan") str.toFloat else Float.NaN).toArray)
      }
    }
    rowTitles = rows.toArray(new Array[String](0))
  }

  /**
   * Rearranges a row from the gene expression matrix according to
   * the given cluster column arrangement:
   * Columns inside the cluster go to the left, columns not in the cluster
   * go to the right.
   */
  private def rearrangeRow(row: Array[Float], inIndexes: Array[Int],
                           outIndexes: Array[Int]): Array[Float] = {
    val result = Array.ofDim[Float](row.length)
    var j = 0
    for (i <- 0 until inIndexes.length) {
      result(j) = row(inIndexes(i))
      j += 1
    }
    for (i <- 0 until outIndexes.length) {
      result(j) = row(outIndexes(i))
      j += 1
    }
    result
  }

  def readRatios(genes: Array[String],
                 conditions: Array[String]): RatioMatrix = {
    var i = 0
    val outValues = new java.util.ArrayList[Array[Float]]

    // partition columns according to cluster membership
    val columnToIndex: Map[String, Int] =
      (for (i <- 0 until columnTitles.length) yield (columnTitles(i) -> i)).toMap
    val inIndexes: Array[Int] = conditions.map(columnToIndex)
    val outIndexes:Array[Int] =
      (for (i <- 0 until columnTitles.length; if !inIndexes.contains(i)) yield i).toArray
    
    while (i < rowTitles.length) {
      if (genes.contains(rowTitles(i))) {
        // rearrange each according to the column partitioning
        outValues.add(rearrangeRow(values.get(i), inIndexes, outIndexes))
      }
      i += 1
    }
    RatioMatrix(genes, columnTitles,
                outValues.toArray(new Array[Array[Float]](0)))
  }
}
