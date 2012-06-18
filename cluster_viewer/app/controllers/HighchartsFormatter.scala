package controllers
import scala.collection.JavaConversions._

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
    println("# rows: " + ratios.rows.length)
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

  def toNRowNColHSSeries(stats: Map[Int, IterationStats]) = {
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

  // motif-score-upstream
  def toMotifPValueSeries(stats: Map[Int, IterationStats],
                          runLogs: Array[RunLog]) = {
    val builder = new StringBuilder
    builder.append("[ { name: 'motif pvalue', data: [")
    val iterations = stats.keySet.toArray
    java.util.Arrays.sort(iterations)
    for (iteration <- iterations) {
      builder.append(stats(iteration).motifPValue)
      builder.append(", ")
    }
    builder.append("] } ]\n")
    builder.toString
  }

  def toNetworkScoreSeries(stats: Map[Int, IterationStats],
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

  def toFuzzyCoeffSeries(stats: Map[Int, IterationStats]) = {
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
