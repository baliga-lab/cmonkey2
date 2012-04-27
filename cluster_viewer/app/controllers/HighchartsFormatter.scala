package controllers

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
}
