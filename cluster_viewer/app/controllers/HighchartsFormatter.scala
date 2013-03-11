package controllers
import scala.collection.JavaConversions._
import scala.collection.mutable.ArrayBuffer

// Formatter object for creating common Javascript constructs for
// Highcharts graphs
object HighchartsFormatter {
  def formatStrings(strs: Seq[String]) = {
    "[%s]".format(strs.map(squote).mkString(", "))
  }
  def formatFloats(vals: Seq[Float]) = {
    "[%s]".format(vals.mkString(", "))
  }
  def formatInts(vals: Seq[Int]) = {
    "[%s]".format(vals.mkString(", "))
  }

  private def squote(value: Any): String = "'%s'".format(value.toString)

  def toHSSeries(ratios: RatioMatrix) = {
    val builder = new StringBuilder
    builder.append("[")
    for (row <- 0 until ratios.rows.length) {
      if (row > 0) builder.append(", ")
      builder.append(toHSSeriesEntry(ratios.rows(row), ratios.values(row)))
    }
    builder.append("]")
    builder.toString
  }

  def toHSSeriesEntry(name: String, values: Array[Float]) = {
    "{ name: '%s', data: %s }".format(name, formatFloats(values))
  }

  private def makeBoxPlotRow(tmp: IndexedSeq[Float]) = {
    val min = tmp(0)
    val max = tmp(tmp.length - 1)
    val median = tmp(tmp.length / 2 + tmp.length % 2)
    val quart = tmp.length / 4
    val lowerQuartile = tmp(quart)
    val upperQuartile = tmp(tmp.length - quart - 1)
    Array(min, lowerQuartile, median, upperQuartile, max)
  }

  private def getSortedHalfs(ratios: RatioMatrix, numColumnsIn: Int) = {
    val numColumns = ratios.columns.length
    val numRows = ratios.rows.length
    val outrowsIn = new ArrayBuffer[Array[Float]]
    val outrowsOut = new ArrayBuffer[Array[Float]]

    for (col <- 0 until numColumns) {
      val tmp = for (row <- 0 until numRows) yield ratios.values(row)(col)
      val outrow = makeBoxPlotRow(tmp.sortWith(_ < _))
      if (col < numColumnsIn) outrowsIn += outrow
      else outrowsOut += outrow
    }

    // split halfs
    val sortedIn = outrowsIn.sortWith((row1, row2) => row1(2) < row2(2))
    val sortedOut = outrowsOut.sortWith((row1, row2) => row1(2) < row2(2))
    (sortedIn, sortedOut)
  }
  def toHSSeriesBoxPlot(ratios: RatioMatrix, numColumnsIn: Int) = {
    val (sortedIn, sortedOut) = getSortedHalfs(ratios, numColumnsIn)
    val builder = new StringBuilder
    builder.append("[")

    for (row <- 0 until sortedIn.length) {
      if (row > 0) builder.append(",")
      builder.append(formatFloats(sortedIn(row)))
    }
    for (row <- 0 until sortedOut.length) {
      builder.append(",")
      builder.append(formatFloats(sortedOut(row)))
    }

    builder.append("]")
    builder.toString
  }
  // **********************************************************************
  private def stdDev(vals: Seq[Float]) = {
    val numValues: Float = vals.length
    var sum : Float = 0.0f
    val mean = vals.sum / numValues
    val factor: Float = 1.0f / (numValues - 1)
    for (x <- vals) sum += (x - mean) * (x - mean)
    sum *= factor
    math.sqrt(sum)
  }
/*
  private def getSortedHalfsStdDev(ratios: RatioMatrix, numColumnsIn: Int) = {
    val numColumns = ratios.columns.length
    val numRows = ratios.rows.length
    val outvalsIn = new ArrayBuffer[Float]
    val outvalsOut = new ArrayBuffer[Float]

    for (col <- 0 until numColumns) {
      val tmp = for (row <- 0 until numRows) yield ratios.values(row)(col)
      val outval = stdDev(tmp)
      if (col < numColumnsIn) outvalsIn += outval
      else outvalsOut += outval
    }

    // split halfs
    val sortedIn = outrowsIn.sortWith((row1, row2) => row1(2) < row2(2))
    val sortedOut = outrowsOut.sortWith((row1, row2) => row1(2) < row2(2))
    (sortedIn, sortedOut)
  }

  def toHSSeriesStdDev(ratios: RatioMatrix, numColumnsIn: Int) = {
  }*/
  // **********************************************************************

  def toHSSeries(values: Array[Double]) = {
    val builder = new StringBuilder
    builder.append("[ { name: 'mean resid', data: [")
    for (i <- 0 until values.length) {
      if (i > 0) builder.append(", ")
      builder.append("%f".format(values(i)))
    }
    builder.append("] }]")
    builder.toString
  }

  def toNRowNColHSSeries(stats: Map[Int, IterationStat]) = {
    val iterations = stats.keySet.toArray
    java.util.Arrays.sort(iterations)

    val meanRows = new Array[Double](iterations.length)
    val meanColumns = new Array[Double](iterations.length)

    for (i <- 0 until iterations.length) {
      var totalRows = 0.0
      var totalColumns = 0.0
      for (cluster <- stats(iterations(i)).clusters.keys) {
        totalRows += stats(iterations(i)).clusters(cluster).numRows
        totalColumns += stats(iterations(i)).clusters(cluster).numColumns
      }
      val numClusters = stats(iterations(i)).clusters.size
      meanRows(i) = totalRows / numClusters
      meanColumns(i) = totalColumns / numClusters
    }

    // columns
    val builder = new StringBuilder
    builder.append("[ { name: 'columns', data: [")
    for (i <- 0 until iterations.length) {
      if (i > 0) builder.append(", ")
      builder.append("%f".format(meanColumns(i)))
    }
    builder.append("] }, ")

    // rows
    builder.append("{ name: 'rows', data: [")
    for (i <- 0 until iterations.length) {
      if (i > 0) builder.append(", ")
      builder.append("%f".format(meanRows(i)))
    }
    builder.append("] }]")

    builder.toString
  }

  // motif-score-upstream, tricky: multiple sequence types
  def toMotifPValueSeries(stats: Map[Int, IterationStat],
                          runLogs: Array[RunLog]) = {
    if (!stats.isEmpty) {
      val iterations = stats.keySet.toArray
      java.util.Arrays.sort(iterations)
      val seqTypes = stats(iterations(0)).motifPValues.keys.toArray
      val builder = new StringBuilder

      builder.append("[")
      for (i <- 0 until seqTypes.length) {
        if (i > 0) builder.append(",")
        builder.append("{ name: '%s', data: [".format(seqTypes(i)))
        for (iteration <- iterations) {
          builder.append(stats(iteration).motifPValues(seqTypes(i)))
          builder.append(", ")
        }
        builder.append("] }")
      }
      builder.append("]\n")
      builder.toString
    } else "[]"
  }

  def toNetworkScoreSeries(stats: Map[Int, IterationStat],
                           runLogs: Array[RunLog]) = {
    val builder = new StringBuilder
    val scoreMap = new java.util.HashMap[String, java.util.ArrayList[Double]]    
    val iterations = stats.keySet.toArray
    java.util.Arrays.sort(iterations)
    for (iteration <- iterations) {
      val networkScores = stats(iteration).networkScores
      for (network <- networkScores.keys) {
        if (!scoreMap.containsKey(network)) {
          scoreMap(network) = new java.util.ArrayList[Double]
        }
        scoreMap(network).append(networkScores(network))
      }
    }
    builder.append("[")
    var started = false
    for (network <- scoreMap.keys) {
      if (started) builder.append(", ")
      else started = true
      builder.append("{ name: '%s', data: ".format(network))
      builder.append(scoreMap(network).mkString("[", ", ", "]"))
      builder.append("}")
    }
    builder.append("]")
    builder.toString
  }

  def toFuzzyCoeffSeries(stats: Map[Int, IterationStat]) = {
    val builder = new StringBuilder
    builder.append("[ { name: 'fuzzy coeff', data: [")
    val iterations = stats.keySet.toArray
    java.util.Arrays.sort(iterations)
    for (iteration <- iterations) {
      builder.append(stats(iteration).fuzzyCoeff)
      builder.append(", ")
    }
    builder.append("] } ]\n")
    builder.toString
  }

  def toHSSeriesEntry(runLog: RunLog): String = {
    val finalScaling = Array.ofDim[Float](runLog.active.length)
    for (i <- 0 until finalScaling.length) {
      finalScaling(i) = if (runLog.active(i)) runLog.scaling(i) else 0.0f
    }
    HighchartsFormatter.toHSSeriesEntry(runLog.functionName, finalScaling)
  }

}
