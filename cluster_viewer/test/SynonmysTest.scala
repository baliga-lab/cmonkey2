import org.specs2.mutable._

import play.api.test._
import play.api.test.Helpers._
import controllers._

class SynonymsMapCSV2Spec extends Specification {
  "SynonymsMapCSV2" should {
    "read a csv file" in {
      val file = new java.io.File("testdata/hal_synonyms.gz")
      val synonyms = new SynonymsMapCSV2(file)
      synonyms("NP_045946.1") must equalTo("VNG7001")
      synonyms("1446803") must equalTo("VNG7001")
      synonyms("foobar") must equalTo("foobar")
    }
  }
}
