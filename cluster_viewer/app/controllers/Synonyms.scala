package controllers

import java.io._
import java.util.zip._
import scala.collection.mutable.HashMap

/**
 * Providing a factory to read potentially limitless numbers of different
 * synonym formats
 */
trait SynonymsMap {
  def apply(str: String): String
}

/**
 * This synonyms map simply returns the name it was given as its alias
 */
object SynonymsMapIdentity extends SynonymsMap {
  def apply(str: String) = str
}

/**
 * Reads in a synonym file in the format
 * original,alternative1;alternative2;...
 */
class SynonymsMapCSV2(file: File) extends SynonymsMap {

  val entries = new HashMap[String, String]
  initEntries

  def apply(str: String): String = if (entries contains str) entries(str) else str

  private def initEntries {
    var in: BufferedReader = null
    try {
      in = new BufferedReader(new InputStreamReader(new GZIPInputStream(
        new FileInputStream(file))))

      var line = in.readLine
      while (line != null) {
        val comps = line.trim.split(",")
        val original = comps(0).toUpperCase
        entries(original) = original
        val alternatives = comps(1).split(";")
        for (alternative <- alternatives) {
          entries(alternative.toUpperCase) = original
        }
        line = in.readLine
      }
    } finally {
      if (in != null) in.close
    }
  }
}

/**
 * Use this factory to access the different formats which are specified
 * in the application.conf file.
 */
object SynonymsFactory {
  def getSynonyms(format: String, filepath: String): SynonymsMap = {
    if (format == null) SynonymsMapIdentity
    else {
      new SynonymsMapCSV2(new File(filepath))
    }
  }
}
