package org.pharmgkb.pharmcat;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.stream.Collectors;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableMap;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;
import org.jsoup.select.Elements;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.pharmcat.reporter.MessageHelper;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.handlebars.ReportHelpers;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;
import org.pharmgkb.pharmcat.reporter.model.result.AnnotationReport;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.DrugReport;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.Genotype;
import org.pharmgkb.pharmcat.reporter.model.result.GuidelineReport;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.contains;
import static org.junit.jupiter.api.Assertions.*;


/**
 * This is a JUnit test for {@link Pipeline}.
 * This should the data generated from a full run of the PharmCAT matcher and reporter.
 *
 * @author Mark Woon
 */
class PipelineTest {
  private static final String sf_unknownDiplotype = Haplotype.UNKNOWN + TextConstants.GENOTYPE_DELIMITER + Haplotype.UNKNOWN;
  private static final List<String> sf_unknownCall = List.of(sf_unknownDiplotype);
  private static final List<String> sf_notCalled = List.of(TextConstants.UNCALLED);
  private static Path s_outsideCallFilePath;
  private static Path s_otherOutsideCallFilePath;


  @BeforeAll
  static void prepare() throws IOException {

    ReportHelpers.setDebugMode(true);
    s_outsideCallFilePath = TestUtils.createTestFile(PipelineTest.class, "outsideCall.tsv");
    try (BufferedWriter writer = Files.newBufferedWriter(s_outsideCallFilePath)) {
      writer.write("""
          ##Test Outside Call Data
          CYP2D6\tCYP2D6*1/CYP2D6*4\t\t\t0.6\t0.75\tp: 0.0\t\t\tv1.9-2017_02_09
              """);
    }

    s_otherOutsideCallFilePath = TestUtils.createTestFile(PipelineTest.class, "otherOutsideCall.tsv");
    try (BufferedWriter writer = Files.newBufferedWriter(s_otherOutsideCallFilePath)) {
      writer.write("""
          CYP2D6\t*3/*4
          G6PD\tB (wildtype)/B (wildtype)
          """);
    }
    //TestUtils.setSaveTestOutput(true);
  }

  @AfterEach
  void deleteDirectory(TestInfo testInfo) {
    TestUtils.deleteTestOutputDirectory(testInfo);
  }


  public static List<String> expectedCallsToRecommendedDiplotypes(List<String> expectedCalls) {
    Preconditions.checkArgument(expectedCalls.size() == 1);
    return expectedCalls.stream()
        .flatMap(s -> Arrays.stream(s.split(TextConstants.GENOTYPE_DELIMITER)))
        .toList();
  }


  /**
   * Checks for expected HTML output.
   *
   * @param expectedCalls - null if no data is provided, empty array if no call is expected
   */
  private void htmlChecks(Document document, String gene, @Nullable List<String> expectedCalls,
      @Nullable String drug, RecPresence cpicAnnPresence, RecPresence dpwgAnnPresence) {
    Map<String, List<String>> geneCallMap = new HashMap<>();
    geneCallMap.put(gene, expectedCalls);
    htmlChecks(document, geneCallMap, drug, cpicAnnPresence, dpwgAnnPresence);
  }

  /**
   * Checks for expected HTML output.
   *
   * @param expectedCalls - null if no data is provided, empty array if no call is expected
   */
  private void htmlChecks(Document document, Map<String, List<String>> expectedCalls,
      @Nullable String drug, RecPresence cpicAnnPresence, RecPresence dpwgAnnPresence) {

    for (String gene : expectedCalls.keySet()) {
      htmlCheckGene(document, gene, expectedCalls.get(gene));
    }

    if (drug != null) {
      htmlCheckDrug(document, expectedCalls, drug, cpicAnnPresence, dpwgAnnPresence);
    }
  }

  private void htmlCheckGene(Document document, String gene, List<String> expectedCalls) {
    if (expectedCalls == sf_unknownCall) {
      expectedCalls = null;
    }
    if (expectedCalls == null || expectedCalls.size() == 0) {
      // check section i
      assertEquals(0, document.select(".gs-" + gene + " .gs-dip").size());
      if (expectedCalls != null) {
        assertNotNull(document.getElementById("gs-uncallable-" + gene));
      }

      // check section iii
      Elements geneSection = document.select(".gene." + gene);
      assertEquals(1, geneSection.size());
      if (expectedCalls == null) {
        assertEquals(1, geneSection.get(0).getElementsByClass("no-data").size());
      } else {
        assertEquals("Not called", geneSection.select(".genotype-result").text());
      }

    } else {
      // check section i
      Elements gsDips = document.select(".gs-" + gene + " .gs-dip");
      assertEquals(expectedCalls.size(), gsDips.size());
      assertEquals(expectedCalls,
          gsDips.stream()
              .map(e -> e.child(0).text())
              .toList());

      // check section iii
      Elements geneSection = document.select(".gene." + gene);
      assertEquals(1, geneSection.size());
      assertEquals(0, geneSection.get(0).getElementsByClass("no-data").size());
    }
  }

  private void htmlCheckDrug(Document document, String gene, @Nullable List<String> expectedCalls,
      @Nullable String drug, RecPresence cpicAnnPresence, RecPresence dpwgAnnPresence) {
    Map<String, List<String>> geneCallMap = new HashMap<>();
    geneCallMap.put(gene, expectedCalls);
    htmlCheckDrug(document, geneCallMap, drug, cpicAnnPresence, dpwgAnnPresence);
  }

  private void htmlCheckDrug(Document document, Map<String, List<String>> expectedCalls, String drug,
      RecPresence cpicAnnPresence, RecPresence dpwgAnnPresence) {

    String sanitizedDrug = ReportHelpers.sanitizeCssSelector(drug);
    Elements drugSections = document.getElementsByClass(sanitizedDrug);

    if (cpicAnnPresence == RecPresence.NO && dpwgAnnPresence == RecPresence.NO) {
      assertEquals(0, drugSections.size());

    } else {
      assertEquals(1, drugSections.size());

      List<String> expectedRxCalls = new ArrayList<>();
      for (String gene : expectedCalls.keySet()) {
        List<String> calls = expectedCalls.get(gene);
        if (calls == null || calls.size() == 0) {
          Elements cpicDrugDips = drugSections.select(".cpic-" + sanitizedDrug + " .rx-dip");
          assertEquals(0, cpicDrugDips.size());

          Elements dpwgDrugDips = drugSections.select(".dpwg-" + sanitizedDrug + " .rx-dip");
          assertEquals(0, dpwgDrugDips.size());

          continue;
        }

        for (String call : calls) {
          expectedRxCalls.add(gene + ":" + call);
        }
      }

      htmlCheckDrugAnnotation(drugSections, "cpic", sanitizedDrug, cpicAnnPresence, expectedCalls.keySet(),
          expectedRxCalls);
      htmlCheckDrugAnnotation(drugSections, "dpwg", sanitizedDrug, dpwgAnnPresence, expectedCalls.keySet(),
          expectedRxCalls);
    }
  }

  private void htmlCheckDrugAnnotation(Elements drugSections, String src, String drug, RecPresence annPresence,
      Collection<String> genes, List<String> expectedRxCalls) {

    String baseSelector = "." + src + "-" + ReportHelpers.sanitizeCssSelector(drug);

    Elements drugDips = drugSections.select(baseSelector + " .rx-dip");
    if (annPresence == RecPresence.YES) {
      assertEquals(expectedRxCalls,
          drugDips.stream()
              .map(e -> cleanupRxDip(e, genes))
              .toList());
    } else {
      assertEquals(0, drugDips.size());

      if (drugSections.select(baseSelector).size() > 0) {
        Elements unmatchedDips = drugSections.select(baseSelector + " .rx-unmatched-dip");
        assertEquals(expectedRxCalls, unmatchedDips.stream()
            .map(e -> cleanupRxDip(e, genes))
            .toList());
      }
    }
  }

  public static String cleanupRxDip(Element rxDip, Collection<String> genes) {
    String dip = rxDip.text().replace("/ ", "/");
    for (String gene : genes) {
      dip = dip.replace(gene + ": ", gene + ":");
    }
    return dip;
  }


  /**
   * NOTE: if these assertions fail then new data may have been added from the DataManager because of an update to the
   * CPIC database. If that's true, then update these numbers to the current count. If the count changes with no known
   * change to the CPIC database then something may be wrong in code.
   */
  @Test
  void testCounts(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C9");
    testWrapper.execute(null);
    assertEquals(23, testWrapper.getContext().getGeneReports().keySet().stream()
        .flatMap((k) -> testWrapper.getContext().getGeneReports().get(k).values().stream()
            .map(GeneReport::getGeneDisplay))
        .collect(Collectors.toSet())
        .size()
    );
    assertEquals(128, testWrapper.getContext().getDrugReports().keySet().stream()
        .flatMap((k) -> testWrapper.getContext().getDrugReports().get(k).values().stream()
            .map(DrugReport::getName))
        .collect(Collectors.toSet())
        .size());
  }

  @Test
  void testAll(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println(
          """
              CYP2D6\t*3/*4
              HLA-A\t\t*31:01 positive
              HLA-B\t*15:02/*57:01"""
      );
    }
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("ABCG2")
        .reference("CACNA1S")
        .reference("CFTR")
        .reference("CYP2B6")
        .reference("CYP2C19")
        .variation("CYP2C19", "rs3758581", "G", "G") // to make it *1/*1
        .reference("CYP2C9")
        .reference("CYP3A4")
        .reference("CYP3A5")
        .reference("CYP4F2")
        .reference("DPYD")
        .reference("G6PD")
        .reference("IFNL3")
        .reference("NUDT15")
        .reference("RYR1")
        .reference("SLCO1B1")
        .reference("TPMT")
        .reference("UGT1A1")
        .reference("VKORC1");
    testWrapper.execute(outsideCallPath);

    testWrapper.testCalledByMatcher(
        "ABCG2",
        "CACNA1S",
        "CFTR",
        "CYP2B6",
        "CYP2C19",
        "CYP2C9",
        "CYP3A4",
        "CYP3A5",
        "CYP4F2",
        "DPYD",
        "G6PD",
        "IFNL3",
        "NUDT15",
        "RYR1",
        "SLCO1B1",
        "TPMT",
        "UGT1A1",
        "VKORC1"
    );
    testWrapper.testNotCalledByMatcher("CYP2D6", "HLA-A", "HLA-B");
  }


  @Test
  void testUndocumentedVariation(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .allowUnknownAllele()
        .variation("CYP2C19", "rs3758581", "G", "T");
    Path vcfFile = testWrapper.execute(null);

    testWrapper.testNotCalledByMatcher("CYP2C19");

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, "CYP2C19", Collections.emptyList(), null, RecPresence.YES, RecPresence.YES);
  }

  @Test
  void testUndocumentedVariationExtendedReport(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false)
        .extendedReport();
    testWrapper.getVcfBuilder()
        .allowUnknownAllele()
        .variation("CYP2C19", "rs3758581", "G", "T");
    Path vcfFile = testWrapper.execute(null);

    testWrapper.testNotCalledByMatcher("CYP2C19");

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, "CYP2C19", Collections.emptyList(), null, RecPresence.YES, RecPresence.NO);
  }

  @Test
  void testUndocumentedVariationsWithTreatAsReference(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .allowUnknownAllele()
        .variation("RYR1", "rs193922753", "G", "A")
        .variation("CYP2C19", "rs3758581", "G", "T");
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedRyr1Calls = List.of("No CPIC variants found");

    testWrapper.testNotCalledByMatcher("CYP2C19");
    testWrapper.testCalledByMatcher("RYR1");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "RYR1", expectedRyr1Calls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "RYR1", List.of("Reference", "Reference"));

    Document document = readHtmlReport(vcfFile);
    assertNotNull(document.getElementById("gs-undocVarAsRef-RYR1"));
    htmlChecks(document, "CYP2C19", Collections.emptyList(), null, RecPresence.YES, RecPresence.NO);
    htmlChecks(document, "RYR1", expectedRyr1Calls, null, RecPresence.YES, RecPresence.NO);
  }

  @Test
  void testUndocumentedVariationsWithTreatAsReferenceAndCombo(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, true, false, false);
    testWrapper.getVcfBuilder()
        .allowUnknownAllele()
        .variation("RYR1", "rs193922753", "G", "A");
    Path vcfFile = testWrapper.execute(null);

    // becomes Reference and custom snp because combo is enabled
    List<String> expectedCalls = List.of("g.38444212G>A (heterozygous)");

    testWrapper.testCalledByMatcher("RYR1");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "RYR1", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "RYR1", List.of("Reference", "g.38444212G>A"));

    Document document = readHtmlReport(vcfFile);
    assertNull(document.getElementById("gs-undocVarAsRef-RYR1"));
    htmlChecks(document, "RYR1", expectedCalls, null, RecPresence.YES, RecPresence.NO);
  }


  /**
   * This test illustrates when one gene in a two-gene guideline (amitriptyline) is not called that it should still be
   * able to come up with a matched annotation.
   */
  @Test
  void testCyp2c19(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs3758581", "G", "G");
    testWrapper.execute(s_otherOutsideCallFilePath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testPrintCpicCalls( "CYP2C19", "*1/*1");

    testWrapper.testMatchedAnnotations("amitriptyline", DataSource.CPIC, 1);
    testWrapper.testMatchedAnnotations("amitriptyline", DataSource.DPWG, 1);
    testWrapper.testMatchedAnnotations("citalopram", DataSource.CPIC, 1);
    testWrapper.testMatchedAnnotations("citalopram", DataSource.DPWG, 1);
    testWrapper.testMatchedAnnotations("ivacaftor", 0);
  }

  /**
   * This test case demos that an "ambiguity" {@link MessageAnnotation} which specifies a variant and a diplotype call
   * for a given drug report will be matched and added to the {@link DrugReport}
   */
  @Test
  void testCyp2c19_s1s2rs58973490het(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs12769205", "A", "G")
        .variation("CYP2C19", "rs58973490", "G", "A")
        .variation("CYP2C19", "rs4244285", "G", "A")
        .variation("CYP2C19", "rs3758581", "G", "G");
    testWrapper.execute(s_otherOutsideCallFilePath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testPrintCpicCalls( "CYP2C19", "*1/*2");
    testWrapper.testNotCalledByMatcher("CYP2D6");
    testWrapper.testPrintCpicCalls( "CYP2D6", "*3/*4");

    testWrapper.testMatchedAnnotations("amitriptyline", 2);
    testWrapper.testAnyMatchFromSource("amitriptyline", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("amitriptyline", DataSource.DPWG);
    testWrapper.testMatchedAnnotations("citalopram", 2);
    testWrapper.testAnyMatchFromSource("citalopram", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("citalopram", DataSource.DPWG);
    testWrapper.testMatchedAnnotations("clomipramine", 3);
    testWrapper.testAnyMatchFromSource("clomipramine", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("clomipramine", DataSource.DPWG);
    testWrapper.testMatchedAnnotations("ivacaftor", 0);

    GeneReport cyp2c19report = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2C19");
    assertNotNull(cyp2c19report);
    VariantReport vr = cyp2c19report.findVariantReport("rs58973490")
        .orElseThrow(() -> new RuntimeException("Variant missing from test data"));
    assertTrue(vr.isHetCall());

    // ambiguity message will not apply in this case because all variants are available for CYP2C19, but one message
    // should appear for the *1 call
    assertEquals(1, cyp2c19report.getMessages().stream()
        .filter(m -> m.getExceptionType().equals(MessageAnnotation.TYPE_AMBIGUITY) &&
            Objects.requireNonNull(m.getMatches().getVariant()).equals("rs58973490"))
        .count());

    testWrapper.testMessageCountForDrug(DataSource.CPIC, "amitriptyline", 1);
  }

  /**
   * This test case demos that an "ambiguity" {@link MessageAnnotation} which specifies a variant and a diplotype call
   * for a given drug report will not be matched when the variant in the message is homozygous
   */
  @Test
  void testCyp2c19_s1s2(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs12769205", "A", "G")
        .variation("CYP2C19", "rs4244285", "G", "A")
        .variation("CYP2C19", "rs3758581", "G", "G");
    testWrapper.execute(s_otherOutsideCallFilePath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testPrintCpicCalls( "CYP2C19", "*1/*2");

    testWrapper.testMatchedAnnotations("amitriptyline", 2);
    testWrapper.testAnyMatchFromSource("amitriptyline", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("amitriptyline", DataSource.DPWG);
    testWrapper.testMatchedAnnotations("citalopram", 2);
    testWrapper.testAnyMatchFromSource("citalopram", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("citalopram", DataSource.DPWG);
    testWrapper.testMatchedAnnotations("clomipramine", 3);
    testWrapper.testAnyMatchFromSource("clomipramine", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("clomipramine", DataSource.DPWG);
    testWrapper.testMatchedAnnotations("ivacaftor", 0);

    GeneReport cyp2c19report = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2C19");
    assertNotNull(cyp2c19report);

    // make sure the variant in question is not a het call
    VariantReport vr = cyp2c19report.findVariantReport("rs58973490")
        .orElseThrow(() -> new RuntimeException("Variant missing from test data"));
    assertFalse(vr.isHetCall());

    // the variant is hom so ambiguity message should not apply and, thus, no matching messages
    assertEquals(0, cyp2c19report.getMessages().stream()
        .filter(m -> m.getExceptionType().equals(MessageAnnotation.TYPE_AMBIGUITY) &&
            Objects.requireNonNull(m.getMatches().getVariant()).equals("rs58973490"))
        .count());

    testWrapper.testAnyMatchFromSource("amitriptyline", DataSource.CPIC);
    // the variant is hom so ambiguity message should not match
    testWrapper.testMessageCountForDrug(DataSource.CPIC, "amitriptyline", 0);

    // CYP2C19 reference is *38, not *1, so should not have reference message
    testWrapper.testMessageCountForGene(DataSource.CPIC, "CYP2C19", 0);
  }

  @Test
  void testClomipramineCall(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs12769205", "G", "G")
        .variation("CYP2C19", "rs4244285", "A", "A")
        .variation("CYP2C19", "rs3758581", "G", "G");
    testWrapper.execute(s_otherOutsideCallFilePath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testPrintCpicCalls( "CYP2C19", "*2/*2");

    testWrapper.testMatchedAnnotations("amitriptyline", 2);
    testWrapper.testAnyMatchFromSource("amitriptyline", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("amitriptyline", DataSource.DPWG);
    testWrapper.testMatchedAnnotations("clomipramine", 3);
    testWrapper.testAnyMatchFromSource("clomipramine", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("clomipramine", DataSource.DPWG);
    testWrapper.testMatchedAnnotations("desipramine", 1);
    testWrapper.testAnyMatchFromSource("desipramine", DataSource.CPIC);
    testWrapper.testMatchedAnnotations("doxepin", 2);
    testWrapper.testAnyMatchFromSource("doxepin", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("doxepin", DataSource.DPWG);
    testWrapper.testMatchedAnnotations("imipramine", 3);
    testWrapper.testAnyMatchFromSource("imipramine", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("imipramine", DataSource.DPWG);
    testWrapper.testMatchedAnnotations("nortriptyline", 2);
    testWrapper.testAnyMatchFromSource("nortriptyline", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("nortriptyline", DataSource.DPWG);
    testWrapper.testMatchedAnnotations("trimipramine", 1);

    testWrapper.testMatchedAnnotations("clopidogrel", 4);
    testWrapper.testAnyMatchFromSource("clopidogrel", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("clopidogrel", DataSource.DPWG);

    testWrapper.testMatchedAnnotations("lansoprazole", 2);
    testWrapper.testAnyMatchFromSource("lansoprazole", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("lansoprazole", DataSource.DPWG);

    // voriconazole has 2 populations with recommendations so should have 2 matching annotations from CPIC
    // and 1 from DPWG
    testWrapper.testMatchedAnnotations("voriconazole", DataSource.CPIC, 2);
    testWrapper.testMatchedAnnotations("voriconazole", DataSource.DPWG, 1);
  }

  @Test
  void testCyp2c19noCall(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs12769205", "A", "G")
        .variation("CYP2C19", "rs4244285", "A", "A");
    testWrapper.execute(s_otherOutsideCallFilePath);

    testWrapper.testNotCalledByMatcher("CYP2C19");

    testWrapper.testNoMatchFromSource("citalopram", DataSource.CPIC);
    testWrapper.testNoMatchFromSource("citalopram", DataSource.DPWG);
    testWrapper.testNoMatchFromSource("ivacaftor", DataSource.CPIC);
    testWrapper.testNoMatchFromSource("ivacaftor", DataSource.DPWG);
  }

  @Test
  void testCyp2c19s4bs17rs28399504missing(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs12248560", "T", "T")
        .missing("CYP2C19", "rs28399504")
        .variation("CYP2C19", "rs3758581", "G", "G");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testPrintCpicCalls("CYP2C19", "*4/*4", "*4/*17", "*17/*17");

    testWrapper.testMatchedAnnotations("citalopram", 6);
    testWrapper.testMatchedAnnotations("citalopram", DataSource.CPIC, 3);
    testWrapper.testMatchedAnnotations("citalopram", DataSource.DPWG, 3);
    testWrapper.testAnyMatchFromSource("citalopram", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("citalopram", DataSource.DPWG);
  }

  @Test
  void testCyp2c19s1s4het(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs12248560", "T", "T")
        .variation("CYP2C19", "rs28399504", "A", "G")
        .variation("CYP2C19", "rs3758581", "G", "G");
    testWrapper.execute(s_outsideCallFilePath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testReportable("CYP2D6");

    testWrapper.testPrintCpicCalls("CYP2D6", "*1/*4");
    testWrapper.testPrintCpicCalls("CYP2C19", "*4/*17");

    GeneReport geneReport = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2D6");
    assertNotNull(geneReport);
    assertTrue(geneReport.isOutsideCall());
  }

  @Test
  void testCyp2c19s1s4missingS1(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs12248560", "C", "T")
        .variation("CYP2C19", "rs28399504", "A", "G")
        .missing("CYP2C19", "rs3758581");
    testWrapper.execute(s_outsideCallFilePath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testReportable("CYP2D6");

    testWrapper.testPrintCpicCalls("CYP2D6", "*1/*4");
    testWrapper.testPrintCpicCalls("CYP2C19",  "*1/*4", "*4/*38");

    GeneReport cyp2d6Report = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2D6");
    assertNotNull(cyp2d6Report);
    assertTrue(cyp2d6Report.isOutsideCall());

    GeneReport cyp2c19report = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2C19");
    assertNotNull(cyp2c19report);
    assertTrue(cyp2c19report.isMissingVariants());

    assertFalse(cyp2c19report.isPhased());
    assertTrue(cyp2c19report.findVariantReport("rs12248560").map(VariantReport::isHetCall).orElse(false));
    assertTrue(cyp2c19report.findVariantReport("rs3758581").map(VariantReport::isMissing).orElse(false));

    assertTrue(cyp2c19report.hasHaplotype("*38"));

    // message is for *1/*4 being ambiguous with unphased data
    assertEquals(2, cyp2c19report.getMessages().stream()
        .filter(m -> m.getExceptionType().equals(MessageAnnotation.TYPE_AMBIGUITY))
        .count());

    testWrapper.testAnyMatchFromSource("amitriptyline", DataSource.CPIC);
    testWrapper.testMessageCountForDrug(DataSource.CPIC, "amitriptyline", 2);
  }

  @Test
  void testCyp2c19s4s17(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs12248560", "C", "T")
        .variation("CYP2C19", "rs28399504", "A", "G")
        .variation("CYP2C19", "rs3758581", "G", "G");
    testWrapper.execute(s_outsideCallFilePath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testReportable("CYP2D6");

    testWrapper.testPrintCpicCalls("CYP2D6", "*1/*4");
    testWrapper.testPrintCpicCalls("CYP2C19", "*1/*4");

    GeneReport cyp2d6Report = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2D6");
    assertNotNull(cyp2d6Report);
    assertTrue(cyp2d6Report.isOutsideCall());
  }

  @Test
  void testCyp2c19SingleGeneMatch(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19")
        .variation("CYP2C19", "rs3758581", "A", "G")
        .missing("CYP2C19", "rs56337013");
    testWrapper.execute(s_outsideCallFilePath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testReportable("CYP2D6");

    testWrapper.testPrintCpicCalls("CYP2D6", "*1/*4");
    testWrapper.testPrintCpicCalls("CYP2C19", "*1/*38");

    GeneReport cyp2d6Report = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2D6");
    assertNotNull(cyp2d6Report);
    assertTrue(cyp2d6Report.isOutsideCall());
  }


  @Test
  void testCftrRefRef(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CFTR");
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("No CPIC variants found");

    testWrapper.testCalledByMatcher("CFTR");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "CFTR", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "CFTR",
        List.of("ivacaftor non-responsive CFTR sequence", "ivacaftor non-responsive CFTR sequence"));
    testWrapper.testPrintCalls(DataSource.CPIC, "CFTR", expectedCalls);

    testWrapper.testMatchedAnnotations("ivacaftor", 1);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, "CFTR", expectedCalls, "ivacaftor", RecPresence.YES, RecPresence.NO);
  }

  @Test
  void testCftrD1270NHet(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CFTR", "rs11971167", "G", "A");
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("D1270N (heterozygous)");

    testWrapper.testCalledByMatcher("CFTR");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "CFTR", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "CFTR", List.of("ivacaftor non-responsive CFTR sequence", "D1270N"));
    testWrapper.testPrintCalls(DataSource.CPIC, "CFTR", expectedCalls);

    testWrapper.testMatchedAnnotations("ivacaftor", 1);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, "CFTR", expectedCalls, "ivacaftor", RecPresence.YES, RecPresence.NO);
  }

  @Test
  void testCftrD1270NG551D(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CFTR", "rs11971167", "G", "A")
        .variation("CFTR", "rs75527207", "G", "A");
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("D1270N/G551D");

    testWrapper.testCalledByMatcher("CFTR");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "CFTR", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "CFTR", expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testPrintCalls(DataSource.CPIC, "CFTR", expectedCalls);

    testWrapper.testMatchedAnnotations("ivacaftor", 1);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, "CFTR", expectedCalls, "ivacaftor", RecPresence.YES, RecPresence.NO);
  }


  @Test
  void testRyr1HomoRef(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("RYR1");
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("No CPIC variants found");

    testWrapper.testCalledByMatcher("RYR1");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "RYR1", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "RYR1", List.of("Reference", "Reference"));
    testWrapper.testPrintCalls(DataSource.CPIC, "RYR1", expectedCalls);

    testWrapper.testMatchedAnnotations("desflurane", 0);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, "RYR1", expectedCalls, null, RecPresence.YES, RecPresence.NO);
  }

  @Test
  void testRyr1HomoVar(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("RYR1", "rs193922747", "T", "C");
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("c.103T>C (heterozygous)");

    testWrapper.testCalledByMatcher("RYR1");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "RYR1", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "RYR1", List.of("Reference", "c.103T>C"));
    testWrapper.testPrintCalls(DataSource.CPIC, "RYR1", expectedCalls);

    testWrapper.testMatchedAnnotations("desflurane", 0);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, "RYR1", expectedCalls, null, RecPresence.YES, RecPresence.NO);
  }

  @Test
  void testRyr1Het(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("RYR1", "rs193922747", "T", "C")
        .variation("RYR1", "rs193922748", "C", "T");
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("c.103T>C/c.130C>T");

    testWrapper.testCalledByMatcher("RYR1");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "RYR1", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "RYR1", expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testPrintCalls(DataSource.CPIC, "RYR1", expectedCalls);

    testWrapper.testMatchedAnnotations("desflurane", 0);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, "RYR1", expectedCalls, null, RecPresence.YES, RecPresence.NO);
  }


  @Test
  void testRyr1HomoRefCacna1sRef(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CACNA1S")
        .reference("RYR1");
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("No CPIC variants found");

    testWrapper.testCalledByMatcher("CACNA1S", "RYR1");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "CACNA1S", expectedCalls);
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "RYR1", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "CACNA1S", List.of("Reference", "Reference"));
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "RYR1", List.of("Reference", "Reference"));
    testWrapper.testPrintCalls(DataSource.CPIC, "CACNA1S", expectedCalls);
    testWrapper.testPrintCalls(DataSource.CPIC, "RYR1", expectedCalls);

    testWrapper.testMatchedAnnotations("desflurane", 1);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, new ImmutableMap.Builder<String, List<String>>()
            .put("CACNA1S", expectedCalls)
            .put("RYR1", expectedCalls)
            .build(),
        "desflurane", RecPresence.YES, RecPresence.NO);
  }


  @Test
  void testCacna1sHomoRef(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CACNA1S");
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("No CPIC variants found");

    testWrapper.testCalledByMatcher("CACNA1S");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "CACNA1S", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "CACNA1S", List.of("Reference", "Reference"));
    testWrapper.testPrintCalls(DataSource.CPIC, "CACNA1S", expectedCalls);

    testWrapper.testMatchedAnnotations("desflurane", 0);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, new ImmutableMap.Builder<String, List<String>>()
            .put("CACNA1S", expectedCalls)
            .put("RYR1", sf_unknownCall)
            .build(),
        "desflurane", RecPresence.YES_NO_MATCH, RecPresence.NO);
  }

  @Test
  void testCacna1sHomoVar(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CACNA1S", "rs1800559", "C", "T");
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("c.3257G>A (heterozygous)");

    testWrapper.testCalledByMatcher("CACNA1S");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "CACNA1S", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "CACNA1S", List.of("Reference", "c.3257G>A"));
    testWrapper.testPrintCalls(DataSource.CPIC, "CACNA1S", expectedCalls);

    testWrapper.testMatchedAnnotations("desflurane", 0);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, "CACNA1S", expectedCalls, null, RecPresence.YES, RecPresence.NO);
  }

  @Test
  void testCacna1sHet(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CACNA1S", "rs1800559", "C", "T")
        .variation("CACNA1S", "rs772226819", "G", "A");
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("c.520C>T/c.3257G>A");

    testWrapper.testCalledByMatcher("CACNA1S");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "CACNA1S", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "CACNA1S", expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testPrintCalls(DataSource.CPIC, "CACNA1S", expectedCalls);

    testWrapper.testMatchedAnnotations("desflurane", 0);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, "CACNA1S", expectedCalls, null, RecPresence.YES, RecPresence.NO);
  }


  @Test
  void testRosuvastatin(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("ABCG2", "rs2231142", "G", "T")
        .variation("SLCO1B1", "rs56101265", "T", "C");
    Path vcfFile = testWrapper.execute(null);

    testWrapper.testCalledByMatcher("ABCG2", "SLCO1B1");
    testWrapper.testPrintCpicCalls("SLCO1B1", "*1/*2");

    testWrapper.testMatchedAnnotations("rosuvastatin", 1);

    // no dpyd - should not have DPYD warning
    Document document = readHtmlReport(vcfFile);
    Elements capecitabineSection = document.getElementsByClass("capecitabine");
    assertEquals(0, capecitabineSection.size());

    Elements dpydSection = document.select(".gene.dpyd");
    assertEquals(1, dpydSection.size());
    assertEquals(1, dpydSection.get(0).getElementsByClass("no-data").size());
  }

  @Test
  void testAmitryptylineCallWoCyp2c19(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("DPYD");
    testWrapper.execute(s_outsideCallFilePath);

    testWrapper.testReportable("CYP2D6");
    testWrapper.testPrintCpicCalls("CYP2D6", "*1/*4");
    testWrapper.testRecommendedDiplotypes("CYP2D6", "*1", "*4");

    testWrapper.testMatchedAnnotations("amitriptyline", 2);
    testWrapper.testAnyMatchFromSource("amitriptyline", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("amitriptyline", DataSource.DPWG);
  }


  @Test
  void testSlco1b1HomWild(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("SLCO1B1");
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("*1/*1");

    testWrapper.testCalledByMatcher("SLCO1B1");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "SLCO1B1", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "SLCO1B1", expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testPrintCalls(DataSource.CPIC, "SLCO1B1", expectedCalls);

    testWrapper.testMatchedAnnotations("simvastatin", 1);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, "SLCO1B1", expectedCalls, "simvastatin", RecPresence.YES, RecPresence.NO);
  }

  @Test
  void testSlco1b1HomVar(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("SLCO1B1", "rs2306283", "A", "G")
        .variation("SLCO1B1", "rs4149056", "C", "C");
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("*5/*15");

    testWrapper.testCalledByMatcher("SLCO1B1");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "SLCO1B1", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "SLCO1B1", expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testPrintCalls(DataSource.CPIC, "SLCO1B1", expectedCalls);

    testWrapper.testMatchedAnnotations("simvastatin", 1);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, "SLCO1B1", expectedCalls, "simvastatin", RecPresence.YES, RecPresence.NO);
  }

  @Test
  void testSlco1b1Test5(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("SLCO1B1", "rs2306283", "A", "G")
        .variation("SLCO1B1", "rs11045852", "A", "G")
        .variation("SLCO1B1", "rs74064213", "A", "G");
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("*1/*44");

    testWrapper.testCalledByMatcher("SLCO1B1");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "SLCO1B1", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "SLCO1B1", expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testPrintCalls(DataSource.CPIC, "SLCO1B1", expectedCalls);

    testWrapper.testMatchedAnnotations("simvastatin", 1);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, "SLCO1B1", expectedCalls, "simvastatin", RecPresence.YES, RecPresence.NO);
  }

  @Test
  void testSlco1b1Test3(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("SLCO1B1", "rs2306283", "A", "G")
        .variation("SLCO1B1", "rs4149056", "T", "C");
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("*1/*15");

    testWrapper.testCalledByMatcher("SLCO1B1");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "SLCO1B1", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "SLCO1B1", expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testPrintCalls(DataSource.CPIC, "SLCO1B1", expectedCalls);

    testWrapper.testMatchedAnnotations("simvastatin", 1);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, "SLCO1B1", expectedCalls, "simvastatin", RecPresence.YES, RecPresence.NO);
  }

  @Test
  void testSlco1b1Test4(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("SLCO1B1", "rs4149056", "T", "C")
        .variation("SLCO1B1", "rs71581941", "C", "T");
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("*5/*45");

    testWrapper.testCalledByMatcher("SLCO1B1");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "SLCO1B1", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "SLCO1B1", expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testPrintCalls(DataSource.CPIC, "SLCO1B1", expectedCalls);

    testWrapper.testMatchedAnnotations("simvastatin", 1);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, "SLCO1B1", expectedCalls, "simvastatin", RecPresence.YES, RecPresence.NO);
  }

  /**
   * This tests a special case of SLCO1B1. The gene in this scenario is "uncalled" by the matcher due the sample VCF
   * data. However, SLCO1B1 has an override that will display the rs4149056 diplotype regardless of call status. That
   * same override will assign alleles to use for recommendation lookup
   */
  @Test
  void testSlco1b1UncalledOverride(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("SLCO1B1", "rs2306283", "G", "G")
        .variation("SLCO1B1", "rs4149056", "T", "C")
        .variation("SLCO1B1", "rs11045853", "A", "A")
        .variation("SLCO1B1", "rs72559748", "G", "G");
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("rs4149056 C/rs4149056 T");

    testWrapper.testNotCalledByMatcher("SLCO1B1");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "SLCO1B1", sf_unknownCall);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "SLCO1B1", List.of("*1", "*5"));
    testWrapper.testPrintCalls(DataSource.CPIC, "SLCO1B1", expectedCalls);

    testWrapper.testMatchedAnnotations("simvastatin", 1);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, "SLCO1B1", expectedCalls, "simvastatin", RecPresence.YES, RecPresence.NO);
  }


  @Test
  void testUgt1a1Phased(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("UGT1A1", "rs887829", "C", "T");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*1/*80");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*1", "*80");
  }

  @Test
  void testUgt1a1Unphased(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("UGT1A1", "rs887829", "C", "T");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*1/*80");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*1", "*80");
  }

  @Test
  void testUgt1a1s1s1(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("UGT1A1");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*1/*1");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*1", "*1");
  }

  @Test
  void testUgt1a1S1S80S28(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("UGT1A1", "rs887829", "C", "T")
        .variation("UGT1A1", "rs3064744", "TA(7)", "TA(8)");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*1/*80+*28");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*1", "*80+*28");
  }

  @Test
  void testUgt1a1S28S37(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(9)");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*28/*37");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*37", "*28");
  }

  @Test
  void testUgt1a1s28s80phased(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("UGT1A1", "rs887829", "C", "T")
        .variation("UGT1A1", "rs3064744", "TA(7)", "TA(8)");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*1/*80+*28");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*1", "*80+*28");

    // the guideline should not have an ambiguity message
    testWrapper.testMatchedAnnotations("atazanavir", 1);
    testWrapper.testMessageCountForDrug(DataSource.CPIC, "atazanavir", 0);

    testWrapper.testMessageCountForGene(DataSource.CPIC, "UGT1A1", 2);
    testWrapper.testGeneHasMessage(DataSource.CPIC, "UGT1A1", "reference-allele");
  }

  @Test
  void testUgt1a1s28s80s6s60phased(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("UGT1A1", "rs887829", "T", "C")
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(7)")
        .variation("UGT1A1", "rs4148323", "G", "A");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*6/*80+*28");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*6", "*80+*28");
  }

  @Test
  void testUgt1a1s28s80s6s60unphased(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("UGT1A1", "rs887829", "T", "C")
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(7)")
        .variation("UGT1A1", "rs4148323", "G", "A");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*6/*80+*28");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*6", "*80+*28");
  }

  @Test
  void testUgt1a1s6s6(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("UGT1A1", "rs4148323", "A", "A");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*6/*6");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*6", "*6");
  }

  @Test
  void testUgt1a1s6s60s80s28MissingPhased(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .missing("UGT1A1", "rs3064744")
        .variation("UGT1A1", "rs887829", "C", "T")
        .variation("UGT1A1", "rs4148323", "A", "G");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*6/*80", "*6/*80+*28", "*6/*80+*37");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*6", "*80");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*6", "*80+*28");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*6", "*80+*37");
  }

  @Test
  void testUgt1a1s6s60s80s28MissingUnphased(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .missing("UGT1A1", "rs3064744")
        .variation("UGT1A1", "rs887829", "C", "T")
        .variation("UGT1A1", "rs4148323", "A", "G");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*6/*80", "*6/*80+*28", "*6/*80+*37");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*6", "*80");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*6", "*80+*28");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*6", "*80+*37");
  }

  @Test
  void testUgt1a1s80s28missing(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .missing("UGT1A1", "rs3064744")
        .variation("UGT1A1", "rs887829", "C", "T");
    Path vcfFile = testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*1", "*80");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*1", "*80+*28");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*1", "*80+*37");
    testWrapper.testPrintCpicCalls("UGT1A1", "*1/*80", "*1/*80+*28", "*1/*80+*37");

    Document document = readHtmlReport(vcfFile);
    Elements drugSection = document.select(".cpic-atazanavir");
    assertEquals(2, drugSection.size());

    Elements d1 = drugSection.get(0).select(".rx-dip");
    assertEquals(1, d1.size());
    assertEquals("UGT1A1:*1/*80", d1.get(0).text());

    Elements d2 = drugSection.get(1).select(".rx-dip");
    assertEquals(2, d2.size());
    assertEquals(List.of("UGT1A1:*1/*80+*28", "UGT1A1:*1/*80+*37"), d2.stream().map(Element::text).toList());
  }

  @Test
  void testUgt1a1na12717(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("UGT1A1", "rs887829", "T", "T")
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(7)");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*80/*80+*28");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*80", "*80+*28");
  }

  @Test
  void testUgt1a1s28homMissing(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .missing("UGT1A1", "rs887829")
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(8)");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*28/*28", "*28/*80+*28", "*80+*28/*80+*28");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*28", "*28");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*28", "*80+*28");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*80+*28", "*80+*28");
  }

  @Test
  void testUgt1a1s28s60Hom(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(7)");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*1/*28");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*1", "*28");
  }

  @Test
  void testUgt1a1s27s28unphaseds80s60missing(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .missing("UGT1A1", "rs887829")
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(7)")
        .variation("UGT1A1", "rs35350960", "C", "A");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*27/*28", "*27/*80+*28");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*27", "*28");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*27", "*80+*28");
  }

  @Test
  void testUgt1a1HG00436(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("UGT1A1", "rs887829", "T", "C")
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(7)")
        .variation("UGT1A1", "rs35350960", "A", "C");
    testWrapper.execute(null);

    testWrapper.testNotCalledByMatcher("UGT1A1");
  }

  @Test
  void testUgt1a1s1s80s27s60s28missingphased(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .missing("UGT1A1", "rs3064744")
        .variation("UGT1A1", "rs887829", "T", "C")
        .variation("UGT1A1", "rs35350960", "A", "C");
    testWrapper.execute(null);

    testWrapper.testNotCalledByMatcher("UGT1A1");
    GeneReport geneReport = testWrapper.getContext().getGeneReport(DataSource.CPIC, "UGT1A1");
    assertNotNull(geneReport);
    assertTrue(geneReport.isPhased());
  }

  @Test
  void testUgt1a1s1s60s80s6phased(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("UGT1A1", "rs887829", "T", "C")
        .variation("UGT1A1", "rs35350960", "A", "C");
    testWrapper.execute(null);

    testWrapper.testNotCalledByMatcher("UGT1A1");
    GeneReport geneReport = testWrapper.getContext().getGeneReport(DataSource.CPIC, "UGT1A1");
    assertNotNull(geneReport);
    assertTrue(geneReport.isPhased());
  }

  @Test
  void testUgt1a1s1s60s80s28s6phased(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("UGT1A1", "rs887829", "T", "C")
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(7)")
        .variation("UGT1A1", "rs35350960", "A", "C");
    testWrapper.execute(null);

    testWrapper.testNotCalledByMatcher("UGT1A1");
    GeneReport geneReport = testWrapper.getContext().getGeneReport(DataSource.CPIC, "UGT1A1");
    assertNotNull(geneReport);
    assertTrue(geneReport.isPhased());
  }

  @Test
  void testUgt1a1s1s37s80s60phased(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("UGT1A1", "rs887829", "T", "C")
        .variation("UGT1A1", "rs3064744", "TA(9)", "TA(7)");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testReportable("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*1/*80+*37");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*1", "*80+*37");
    GeneReport geneReport = testWrapper.getContext().getGeneReport(DataSource.CPIC, "UGT1A1");
    assertNotNull(geneReport);
    assertTrue(geneReport.isPhased());
  }

  @Test
  void testCyp3a5Missing3Message(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .missing("CYP3A5", "rs776746");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP3A5");
    testWrapper.testReportable("CYP3A5");
    testWrapper.testPrintCpicCalls("CYP3A5", "*1/*1");
    testWrapper.testRecommendedDiplotypes("CYP3A5", "*1", "*1");

    GeneReport gene = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP3A5");
    assertNotNull(gene);
    // rs776746 should be missing from this report
    assertNotNull(gene.getVariantReports());
    assertTrue(gene.getVariantReports().stream().anyMatch(v -> v.isMissing() && v.getDbSnpId().equals("rs776746")));

    // the guideline should have a matching message
    assertTrue(testWrapper.getContext().getDrugReports().get(DataSource.CPIC).values().stream()
        .filter(r -> r.getName().equals("tacrolimus"))
        .allMatch(r -> r.getMessages().size() > 0));

    assertFalse(gene.isPhased());
  }

  @Test
  void testCyp3a5v1(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP3A5", "rs776746", "T", "C");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP3A5");
    testWrapper.testReportable("CYP3A5");
    testWrapper.testPrintCpicCalls("CYP3A5", "*1/*3");
    testWrapper.testRecommendedDiplotypes("CYP3A5", "*1", "*3");
  }

  @Test
  void testCyp3a5v2(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP3A5", "rs28383479", "C", "T")
        .variation("CYP3A5", "rs776746", "C", "T")
    ;
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP3A5");
    testWrapper.testPrintCpicCalls("CYP3A5", "*3/*9");
    testWrapper.testRecommendedDiplotypes("CYP3A5", "*3", "*9");
  }

  @Test
  void testCyp3a5v3(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP3A5", "rs776746", "C", "C")
    ;
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP3A5");
    testWrapper.testPrintCpicCalls("CYP3A5", "*3/*3");
    testWrapper.testRecommendedDiplotypes("CYP3A5", "*3", "*3");
  }

  @Test
  void testCyp3a5v4(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP3A5", "rs776746", "T", "C")
    ;
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP3A5");
    testWrapper.testPrintCpicCalls("CYP3A5", "*1/*3");
    testWrapper.testRecommendedDiplotypes("CYP3A5", "*1", "*3");
  }

  @Test
  void testCyp3a5v5(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP3A5", "rs28383479", "T", "C")
        .variation("CYP3A5", "rs776746", "T", "C")
    ;
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP3A5");
    testWrapper.testPrintCpicCalls("CYP3A5", "*3/*9");
    testWrapper.testRecommendedDiplotypes("CYP3A5", "*3", "*9");
  }

  @Test
  void testHlab(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("HLA-B\t*15:02/*57:01");
    }

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C9");
    testWrapper.execute(outsideCallPath);

    testWrapper.testCalledByMatcher("CYP2C9");
    testWrapper.testReportable("CYP2C9");
    testWrapper.testPrintCalls(DataSource.CPIC, "CYP2C9", "*1/*1");
    testWrapper.testPrintCalls(DataSource.DPWG, "CYP2C9", "*1/*1");

    testWrapper.testNotCalledByMatcher("HLA-B");
    testWrapper.testReportable("HLA-B");
    testWrapper.testSourcePhenotype(DataSource.CPIC, "HLA-B", "*57:01 positive");
    testWrapper.testSourcePhenotype(DataSource.CPIC, "HLA-B", "*58:01 negative");
    testWrapper.testSourcePhenotype(DataSource.CPIC, "HLA-B", "*15:02 positive");
    testWrapper.testSourcePhenotype(DataSource.DPWG, "HLA-B", "*57:01 positive");
    testWrapper.testSourcePhenotype(DataSource.DPWG, "HLA-B", "*58:01 negative");
    testWrapper.testSourcePhenotype(DataSource.DPWG, "HLA-B", "*15:02 positive");

    // *57:01 guideline
    testWrapper.testMatchedAnnotations("abacavir", DataSource.CPIC, 1);
    // TODO: revert when DPWG HLA's are supported again
    //testWrapper.testMatchedAnnotations("abacavir", DataSource.DPWG, 1);
    // *58:01 guideline
    testWrapper.testMatchedAnnotations("allopurinol", DataSource.CPIC, 1);
    // TODO: revert when DPWG HLA's are supported again
    //testWrapper.testMatchedAnnotations("allopurinol", DataSource.DPWG, 1);
    // *15:02 guideline (along with CYP2C9)
    // TODO: revert when DPWG HLA's are supported again
    //testWrapper.testMatchedAnnotations("phenytoin", 4);
    testWrapper.testMatchedAnnotations("phenytoin", 3);
    testWrapper.testAnyMatchFromSource("phenytoin", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("phenytoin", DataSource.DPWG);
  }

  @Test
  void testSingleHlabAllele(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("HLA-B\t*15:02");
    }

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C9");
    testWrapper.execute(outsideCallPath);

    testWrapper.testNotCalledByMatcher("HLA-B");
    testWrapper.testReportable("HLA-B");
    testWrapper.testSourcePhenotype(DataSource.CPIC, "HLA-B", "*57:01 negative");
    testWrapper.testSourcePhenotype(DataSource.CPIC, "HLA-B", "*58:01 negative");
    testWrapper.testSourcePhenotype(DataSource.CPIC, "HLA-B", "*15:02 positive");
    testWrapper.testSourcePhenotype(DataSource.DPWG, "HLA-B", "*57:01 negative");
    testWrapper.testSourcePhenotype(DataSource.DPWG, "HLA-B", "*58:01 negative");
    testWrapper.testSourcePhenotype(DataSource.DPWG, "HLA-B", "*15:02 positive");

    testWrapper.testMatchedAnnotations("abacavir", DataSource.CPIC, 1);
    testWrapper.testMatchedAnnotations("allopurinol", DataSource.CPIC, 1);
    testWrapper.testMatchedAnnotations("phenytoin", 3);
    testWrapper.testAnyMatchFromSource("phenytoin", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("phenytoin", DataSource.DPWG);
  }

  @Test
  void testHlabPhenotype(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("HLA-B\t\t*57:01 positive");
    }

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C9");
    testWrapper.execute(outsideCallPath);

    testWrapper.testCalledByMatcher("CYP2C9");
    testWrapper.testNotCalledByMatcher("HLA-B");
    testWrapper.testReportable("CYP2C9");
    testWrapper.testReportable("HLA-B");
    // TODO: revert when DPWG HLA's are supported again
    //testWrapper.testMatchedAnnotations("abacavir", 2);
    testWrapper.testMatchedAnnotations("abacavir", 1);
    testWrapper.testMatchedAnnotations("abacavir", DataSource.CPIC, 1);
    // TODO: revert when DPWG HLA's are supported again
    //testWrapper.testMatchedAnnotations("abacavir", DataSource.DPWG, 1);
    // allopurinol relies on a different allele for recs so no matches
    testWrapper.testMatchedAnnotations("allopurinol", 0);
    // phenytoin also relies on a different allele but there will be a match for DPWG since the recommendations are
    // split between the two genes on that side
    testWrapper.testMatchedAnnotations("phenytoin", 1);
    testWrapper.testNoMatchFromSource("phenytoin", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("phenytoin", DataSource.DPWG);
  }

  /**
   * An example report that shows a few different types of recommendation scenarios all in one report. The examples
   * shown are:
   * <ul>
   *   <li>celecoxib = 1 CPIC recommenation</li>
   *   <li>citalopram = 2 recommendations: 1 CPIC, 1 DPWG, 1 gene and it's called</li>
   *   <li>clomipramine = 2 recommendations: 1 CPIC, 1 DPWG, 2 gene but only 1 called</li>
   *   <li>carbamezepine = 3 CPIC recommendations on different populations</li>
   *   <li>clopidogrel = 4 recommendations: 3 CPIC on different pops, 1 DPWG</li>
   *   <li>flucloxacillin = 0 recommendations but the gene is reportable</li>
   *   <li>fluvoxamine = 0 recommendations, no gene reportable</li>
   *   <li>siponimod = 1 DPWG recommendation</li>
   * </ul>
   */
  @Test
  void testRecommendationExamples(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("HLA-A\t\t*31:01 positive");
      writer.println("HLA-B\t*57:01/*58:01\t");
    }

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C9")
        .variation("CYP2C19", "rs12769205", "G", "G")
        .variation("CYP2C19", "rs4244285", "A", "A")
        .variation("CYP2C19", "rs3758581", "G", "G");
    testWrapper.execute(outsideCallPath);

    testWrapper.testRecommendedDiplotypes("CYP2C19", "*2", "*2");
    testWrapper.testPrintCpicCalls("CYP2C19", "*2/*2");

    GeneReport cyp2c9 = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2C9");
    assertNotNull(cyp2c9);
    assertEquals(1, cyp2c9.getRecommendationDiplotypes().size());
    assertTrue(cyp2c9.getRecommendationDiplotypes().stream().allMatch(d -> d.getActivityScore().equals("2.0")));

    testWrapper.testReportable("CYP2C19", "CYP2C9", "HLA-A", "HLA-B");
    testWrapper.testMatchedAnnotations("celecoxib", 1);
    testWrapper.testAnyMatchFromSource("celecoxib", DataSource.CPIC);
    testWrapper.testMatchedAnnotations("citalopram", 2);
    testWrapper.testMatchedAnnotations("clomipramine", 2);
    testWrapper.testMatchedAnnotations("clopidogrel", 4);
    testWrapper.testMatchedAnnotations("clopidogrel", DataSource.CPIC, 3);
    testWrapper.testMatchedAnnotations("clopidogrel", DataSource.DPWG, 1);
    testWrapper.testNoMatchFromSource("flucloxacillin", DataSource.CPIC);
    // TODO: revert when DPWG HLA's are supported again
    //testWrapper.testMatchedAnnotations("flucloxacillin", DataSource.DPWG, 1);
    testWrapper.testNoMatchFromSource("fluvoxamine", DataSource.CPIC);
    testWrapper.testNoMatchFromSource("fluvoxamine", DataSource.DPWG);
    testWrapper.testMatchedAnnotations("siponimod", 1);
    testWrapper.testAnyMatchFromSource("siponimod", DataSource.DPWG);

    // TODO: revert when DPWG HLA's are supported again
    //testWrapper.testMatchedAnnotations("carbamazepine", 5);
    testWrapper.testMatchedAnnotations("carbamazepine", 3);
  }

  @Test
  void testTpmtStar1s(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("TPMT", "rs1800460", "C", "T")
        .variation("TPMT", "rs1142345", "T", "C");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("TPMT");
    testWrapper.testPrintCpicCalls("TPMT", "*1/*3A");
    testWrapper.testRecommendedDiplotypes("TPMT", "*1", "*3A");

    GeneReport tpmtReport = testWrapper.getContext().getGeneReport(DataSource.CPIC, "TPMT");
    assertNotNull(tpmtReport);
    assertEquals(43, tpmtReport.getVariantReports().size());
  }


  @Test
  void testCyp2c9star61(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C9", "rs1799853", "C", "T")
        .variation("CYP2C9", "rs202201137", "A", "G");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP2C9");
    testWrapper.testPrintCpicCalls("CYP2C9", "*1/*61");
    testWrapper.testRecommendedDiplotypes("CYP2C9", "*1", "*61");
  }

  @Test
  void testCyp2c9star1Hom(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C9");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP2C9");
    testWrapper.testPrintCpicCalls("CYP2C9", "*1/*1");
    testWrapper.testRecommendedDiplotypes("CYP2C9", "*1", "*1");
    testWrapper.testMatchedAnnotations("celecoxib", 1);
    testWrapper.testMatchedAnnotations("ibuprofen", 1);
    testWrapper.testMatchedAnnotations("lornoxicam", 1);
  }


  /**
   * Test CYP2B6 for a het *34 sample file. When doing the "top match" scenario this will only match to a 1/34 and,
   * thus, only match to a single recommendation. This test will have a different outcome when run in "all matches" mode
   * and should be compared with {@link #testCyp2b6star1star34AllMatch(TestInfo)}.
   */
  @Test
  void testCyp2b6star1star34(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP2B6", "rs34223104", "T", "C")
        .variation("CYP2B6", "rs3211371", "C", "A")
        .variation("CYP2B6", "rs3745274", "G", "T")
        .variation("CYP2B6", "rs2279343", "A", "G")
    ;
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP2B6");
    testWrapper.testPrintCpicCalls("CYP2B6", "*1/*34");
    testWrapper.testRecommendedDiplotypes("CYP2B6", "*1", "*34");
    testWrapper.testMatchedAnnotations("efavirenz", 1);
  }

  /**
   * This test is just like {@link #testCyp2b6star1star34(TestInfo)} but run in "all matches" mode. This should result in 2
   * possible different calls coming from the matcher. These two have different phenotypes and, thus, match to different
   * recommendations.
   */
  @Test
  void testCyp2b6star1star34AllMatch(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, true);
    testWrapper.getVcfBuilder()
        .variation("CYP2B6", "rs34223104", "T", "C")
        .variation("CYP2B6", "rs3211371", "C", "A")
        .variation("CYP2B6", "rs3745274", "G", "T")
        .variation("CYP2B6", "rs2279343", "A", "G")
    ;
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP2B6");
    testWrapper.testPrintCpicCalls("CYP2B6", "*1/*34", "*33/*36");
    testWrapper.testRecommendedDiplotypes("CYP2B6", "*1", "*34");
    testWrapper.testRecommendedDiplotypes("CYP2B6", "*33", "*36");
    testWrapper.testMatchedAnnotations("efavirenz", 2);
  }


  /* NUDT15 */
  @Test
  void testNudt15Ref(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("NUDT15");
    testWrapper.execute(null);

    testWrapper.testPrintCpicCalls("NUDT15", "*1/*1");
    testWrapper.testRecommendedDiplotypes("NUDT15", "*1", "*1");

    testWrapper.testMatchedAnnotations("azathioprine", 2);
    testWrapper.testMatchedAnnotations("mercaptopurine", 2);
    testWrapper.testAnyMatchFromSource("mercaptopurine", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("mercaptopurine", DataSource.DPWG);
    testWrapper.testMatchedAnnotations("thioguanine", 2);
    testWrapper.testAnyMatchFromSource("thioguanine", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("thioguanine", DataSource.DPWG);
  }

  @Test
  void testNudt15S2(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("NUDT15", "rs746071566", "GAGTCG(3)", "GAGTCG(4)")
        .variation("NUDT15", "rs116855232", "C", "T")
    ;
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("NUDT15");
    testWrapper.testPrintCpicCalls("NUDT15", "*1/*2");
    testWrapper.testRecommendedDiplotypes("NUDT15", "*1", "*2");

    DrugReport azaReport = testWrapper.getContext().getDrugReport(DataSource.CPIC, "azathioprine");
    assertNotNull(azaReport);
    GuidelineReport azaCpicGuideline = azaReport.getGuidelines().iterator().next();
    List<Genotype> genotypes = Genotype.makeGenotypes(azaCpicGuideline.getGeneReports());
    assertEquals(1, genotypes.size());

    testWrapper.testMatchedAnnotations("azathioprine", 2);
  }

  @Test
  void testNudt15S3(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("NUDT15", "rs116855232", "C", "T")
    ;
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("NUDT15");
    testWrapper.testPrintCpicCalls("NUDT15", "*1/*3");
    testWrapper.testRecommendedDiplotypes("NUDT15", "*1", "*3");

    testWrapper.testMatchedAnnotations("azathioprine", 2);
    testWrapper.testMatchedAnnotations("mercaptopurine", 2);
    testWrapper.testMatchedAnnotations("thioguanine", 2);
  }


  /* MT-RNR1 */
  @Test
  void testMtrnr1(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("MT-RNR1\t1555A>G");
    }
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19")
        .reference("CYP2C9")
    ;
    testWrapper.execute(outsideCallPath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testCalledByMatcher("CYP2C9");
    testWrapper.testReportable("MT-RNR1");
    testWrapper.testMatchedAnnotations("amikacin", 1);
  }


  /* IFNL3/4 */
  @Test
  void testIfnl3(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("IFNL3")
    ;
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("IFNL3");
    testWrapper.testReportable("IFNL3");
    testWrapper.testPrintCpicCalls("IFNL3", "rs12979860 reference (C)/rs12979860 reference (C)");
    testWrapper.testMatchedAnnotations("peginterferon alfa-2a", 0);
    testWrapper.testMatchedAnnotations("peginterferon alfa-2b", 0);
  }


  /**
   * Tests whether an allele that is unknown to PharmCAT/CPIC will still go through the system without throwing an
   * exception and will be reported properly.
   */
  @Test
  void testCyp2d6AlleleWithNoFunction(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo,".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("CYP2D6\t*1/*XXX");
    }

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19");
    testWrapper.execute(outsideCallPath);

    testWrapper.testPrintCalls(DataSource.CPIC, "CYP2D6", "*1/*XXX");
    testWrapper.testPrintCalls(DataSource.DPWG, "CYP2D6", "*1/*XXX");
    testWrapper.testSourcePhenotype(DataSource.CPIC, "CYP2D6", "n/a");

    // this nonsense allele will still match to "Indeterminate" phenotypes in guidelines for CYP2D6
    testWrapper.testMatchedAnnotations("atomoxetine", DataSource.CPIC, 2);
    DrugReport atoReport = testWrapper.getContext().getDrugReport(DataSource.CPIC, "atomoxetine");
    assertNotNull(atoReport);
    assertNotNull(atoReport.getGuidelines());
    assertEquals(1, atoReport.getGuidelines().size());
    assertTrue(atoReport.getGuidelines().stream()
        .flatMap((g) -> g.getAnnotations().stream())
        .allMatch((a) -> a.getPhenotypes().containsKey("CYP2D6") && a.getPhenotypes().containsValue("Indeterminate")));

    GeneReport geneReport = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2D6");
    assertNotNull(geneReport);
    assertEquals(1, geneReport.getRecommendationDiplotypes().size());
    Diplotype diplotype = geneReport.getRecommendationDiplotypes().first();
    assertEquals("One Normal function allele and one Unassigned function allele", ReportHelpers.gsFunction(diplotype));
  }


  /**
   * Should have call multimatch message.
   */
  @Test
  void testCyp2d6EquivalentDoubleCall(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo,".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("CYP2D6\t*1/*1");
      writer.println("CYP2D6\t*1/*2");
    }

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19");
    Path vcfFile = testWrapper.execute(outsideCallPath);

    GeneReport geneReport = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2D6");
    assertNotNull(geneReport);
    assertEquals(2, geneReport.getRecommendationDiplotypes().size());

    Diplotype diplotype = geneReport.getRecommendationDiplotypes().first();
    assertThat(diplotype.getPhenotypes(), contains("Normal Metabolizer"));
    assertEquals("Two Normal function alleles", ReportHelpers.gsFunction(diplotype));

    Document document = readHtmlReport(vcfFile);
    Elements clomipramineSection = document.select(".guideline.clomipramine");
    assertEquals(1, clomipramineSection.size());
    assertEquals(1, clomipramineSection.get(0).getElementsByClass(MessageHelper.MSG_MUlTI_CALL).size());
  }

  /**
   * Should not have call multimatch message.
   * Should have inferred CYP2D6 copy number.
   */
  @Test
  void testCyp2d6DoubleCall(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo,".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("CYP2D6\t*1/*1");
      writer.println("CYP2D6\t*1/*1x7");
    }

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19");
    Path vcfFile = testWrapper.execute(outsideCallPath);

    GeneReport geneReport = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2D6");
    assertNotNull(geneReport);
    assertEquals(2, geneReport.getRecommendationDiplotypes().size());

    Iterator<Diplotype> recommendationDipIt = geneReport.getRecommendationDiplotypes().iterator();
    Diplotype diplotype = recommendationDipIt.next();
    assertThat(diplotype.getPhenotypes(), contains("Normal Metabolizer"));
    assertEquals("Two Normal function alleles", ReportHelpers.gsFunction(diplotype));

    diplotype = recommendationDipIt.next();
    assertNotNull(diplotype.getAllele2());
    assertEquals("*1x" + TextConstants.GTE + "3", diplotype.getAllele2().getName());
    assertEquals("One Increased function allele and one Normal function allele", ReportHelpers.gsFunction(diplotype));

    Document document = readHtmlReport(vcfFile);
    Elements clomipramineSection = document.select(".guideline.clomipramine");
    assertEquals(1, clomipramineSection.size());
    assertEquals(0, clomipramineSection.get(0).getElementsByClass(MessageHelper.MSG_MUlTI_CALL).size());
  }

  /**
   * Check that a single CYP2D6 phenotype has all the activity scores that could possibly be mapped to it.
   * Specifically, "Intermediate Metabolizer" maps to both "0.25", "0.5", "0.75" and "1.0" activity scores for CYP2D6.
   */
  @Test
  void testCyp2d6PhenotypeHasMultipleActivityScore(TestInfo testInfo) throws Exception {
    Path outDir = TestUtils.getTestOutputDir(testInfo, false);
    Path outsideCallPath = outDir.resolve(TestUtils.getTestName(testInfo) + ".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("CYP2D6\t\tIntermediate Metabolizer");
    }

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19");
    Path vcfFile = testWrapper.execute(outsideCallPath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testPrintCalls(DataSource.CPIC, "CYP2C19", "*38/*38");

    testWrapper.testNotCalledByMatcher("CYP2D6");
    testWrapper.testPrintCalls(DataSource.CPIC, "CYP2D6", "Intermediate Metabolizer");

    testWrapper.testMessageCountForGene(DataSource.CPIC, "CYP2C19", 1);
    testWrapper.testGeneHasMessage(DataSource.CPIC, "CYP2C19", "reference-allele");
    testWrapper.testGeneHasMessage(DataSource.CPIC, "CYP2D6", MessageHelper.MSG_OUTSIDE_CALL);

    // Check that a single CYP2D6 phenotype has all the activity scores that could possibly be mapped to it.
    // Specifically, "Intermediate Metabolizer" maps to both "0.25", "0.5", "0.75" and "1.0" activity scores for CYP2D6.
    GeneReport geneReport = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2D6");
    assertNotNull(geneReport);
    // we expect one diplotype to exist with multiple lookup keys
    assertEquals(1, geneReport.getRecommendationDiplotypes().size());
    Diplotype diplotype = geneReport.getRecommendationDiplotypes().iterator().next();
    // there should be no single activity score specified since this phenotype maps to more than one
    assertTrue(TextConstants.isUnspecified(diplotype.getActivityScore()));
    // there should be two and only two lookup keys, one for each activity score
    assertEquals(4, diplotype.getLookupKeys().size());
    // the two lookup keys should be the two activity scores that correspond to Intermediate Metabolizer
    assertThat(diplotype.getLookupKeys(), contains("0.25", "0.5", "0.75", "1.0"));

    Document document = readHtmlReport(vcfFile);
    assertNotNull(document.getElementById("CYP2D6"));
    Elements cyp2d6Section = document.select(".gene.CYP2D6");
    assertEquals(1, cyp2d6Section.size());
    assertEquals(1, cyp2d6Section.get(0).getElementsByClass(MessageHelper.MSG_OUTSIDE_CALL).size());
    assertEquals(0, cyp2d6Section.get(0).getElementsByClass(MessageHelper.MSG_CYP2D6_MODE).size());

    Elements clomipramineSection = document.select(".guideline.clomipramine");
    assertEquals(1, clomipramineSection.size());
    assertEquals(1, clomipramineSection.get(0).getElementsByClass(MessageHelper.MSG_MULTI_SCORE).size());
  }


  @Test
  void testCyp2d6CpicVsDpwg(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo,".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("CYP2D6\t*1x2/*9");
      writer.println("CYP2D6\t*1x2/*10");
      writer.println("CYP2D6\t*1x2/*17");
      writer.println("CYP2D6\t*1x3/*1");
      writer.println("CYP2D6\t*1/*1");
      writer.println("CYP2D6\t*4/*10");
      writer.println("CYP2D6\t*4/*4");
    }

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19");
    testWrapper.execute(outsideCallPath);

    List<String> expectedCyp2d6Calls = List.of("*1/*1", "*1x2/*9", "*1x2/*10", "*1x2/*17", "*1/*1x3", "*4/*4", "*4/*10");

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testNotCalledByMatcher("CYP2D6");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "CYP2D6", expectedCyp2d6Calls);
    testWrapper.testPrintCalls(DataSource.CPIC, "CYP2D6", expectedCyp2d6Calls);
    testWrapper.testPrintCalls(DataSource.DPWG, "CYP2D6", expectedCyp2d6Calls);

    // TODO: finish this!
  }


  @Test
  void testCyp3a4(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP3A4", "rs72552799", "T", "T")
    ;
    testWrapper.execute(null);
    testWrapper.testCalledByMatcher("CYP3A4");
    testWrapper.testReportable("CYP3A4");
    testWrapper.testPrintCalls(DataSource.DPWG, "CYP3A4", "*8/*8");
    testWrapper.testMatchedAnnotations("quetiapine", 1);
  }

  /**
   * Added to check the output of a partial match for CYP2C19 and make sure messages are applied
   */
  @Test
  void testPartialCall(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, true, true, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19")
        .variation("CYP2C19", "rs367543002", "C", "T")
        .variation("CYP2C19", "rs3758581", "G", "G")
        .missing("CYP2C19", "rs367543003");
    testWrapper.execute(null);
    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testReportable("CYP2C19");
  }


  /**
   * Added to have an example of running in CYP2D6-matching mode and make sure messages are applied
   */
  @Test
  void testCyp2d6Matcher(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false, true, true);
    testWrapper.getVcfBuilder()
        .reference("CYP2D6")
        .reference("CYP2C19")
        .variation("CYP2C19", "rs3758581", "G", "G");
    Path vcfFile = testWrapper.execute(null);
    testWrapper.testCalledByMatcher("CYP2C19", "CYP2D6");
    testWrapper.testReportable("CYP2C19", "CYP2D6");

    Document document = readHtmlReport(vcfFile);
    assertNotNull(document.getElementById("CYP2D6"));
    Elements cyp2d6Section = document.select(".gene.CYP2D6");
    assertEquals(1, cyp2d6Section.size());
    assertEquals(0, cyp2d6Section.get(0).getElementsByClass(MessageHelper.MSG_OUTSIDE_CALL).size());
    assertEquals(1, cyp2d6Section.get(0).getElementsByClass(MessageHelper.MSG_CYP2D6_MODE).size());
  }


  /**
   * Tests how PharmCAT handles that state when sample VCF data exists for a gene and an outside call also exists for
   * that gene. Currently, this should execute successfully by ignoring VCF data and using the outside call
   */
  @Test
  void testOutsideCallCollision(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("CYP2C19\t*2/*2");
    }

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19");
    testWrapper.execute(outsideCallPath);

    testWrapper.testNotCalledByMatcher("CYP2C19");
    // this is the diplotype indicated in the outside call, not the one matched
    testWrapper.testPrintCpicCalls( "CYP2C19", "*2/*2");

    testWrapper.testMessageCountForGene(DataSource.CPIC, "CYP2C19", 2);
    testWrapper.testGeneHasMessage(DataSource.CPIC, "CYP2C19", "prefer-sample-data",
        MessageHelper.MSG_OUTSIDE_CALL);
  }


  /**
   * Tests that an "unordered" diplotype should normalize to the ordered version then it can be used for matching
   */
  @Test
  void testOutsideCallDiplotypeNormalization(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      // dipltoype in backwards order
      writer.println("CYP2C19\t*2/*1");
    }

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C9");
    testWrapper.execute(outsideCallPath);

    testWrapper.testNotCalledByMatcher("CYP2C19");
    // this should be a normalized version of hte given diplotype
    testWrapper.testPrintCpicCalls( "CYP2C19", "*1/*2");
  }

  @Test
  void testOutsideCallPhenotypeOverridesDiplotype(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo,".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("CYP2D6\t*1/*1\tPM\t" + TextConstants.GTE + "4.0");
    }

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19");
    testWrapper.execute(outsideCallPath);

    GeneReport geneReport = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2D6");
    assertNotNull(geneReport);
    assertEquals(1, geneReport.getRecommendationDiplotypes().size());

    Diplotype diplotype = geneReport.getRecommendationDiplotypes().first();
    assertThat(diplotype.getPhenotypes(), contains("Poor Metabolizer"));

    DrugReport drugReport = testWrapper.getContext().getDrugReport(DataSource.CPIC, "clomipramine");
    assertNotNull(drugReport);
    assertEquals(1, drugReport.getGuidelines().size());
    GuidelineReport guidelineReport = drugReport.getGuidelines().first();
    assertEquals(1, guidelineReport.getAnnotations().size());
    AnnotationReport annotationReport = guidelineReport.getAnnotations().first();
    assertEquals("Poor Metabolizer", annotationReport.getPhenotypes().get("CYP2D6"));
  }

  /**
   * Can we use activity scores in outside call files? It should be specified in the column for "phenotype"
   */
  @Test
  void testOutsideCallActivityScore(TestInfo testInfo) throws Exception {
    Path outDir = TestUtils.getTestOutputDir(testInfo, false);
    Path outsideCallPath = outDir.resolve(TestUtils.getTestName(testInfo) + ".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("CYP2D6\t\t\t1.25");
    }

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19");
    testWrapper.execute(outsideCallPath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testPrintCalls(DataSource.CPIC, "CYP2C19", "*38/*38");

    testWrapper.testNotCalledByMatcher("CYP2D6");
    testWrapper.testPrintCalls(DataSource.CPIC, "CYP2D6", "Normal Metabolizer (1.25)");

    testWrapper.testMessageCountForGene(DataSource.CPIC, "CYP2C19", 1);
    testWrapper.testGeneHasMessage(DataSource.CPIC, "CYP2C19", "reference-allele");
    testWrapper.testGeneHasMessage(DataSource.CPIC, "CYP2D6", MessageHelper.MSG_OUTSIDE_CALL);
  }

  /**
   * This test ensures that a user can specify both a diplotype AND a phenotype from an outside call. This also tests
   * to make sure the user can override the internally-known phenotype with their own phenotype assignment. *2/*10 would
   * normally be a Normal Metabolizer but this outside call overrides it as an Intermediate Metabolizer.
   */
  @Test
  void testOutsideCallActivityScoreAndPhenotype(TestInfo testInfo) throws Exception {
    Path outDir = TestUtils.getTestOutputDir(testInfo, false);
    Path outsideCallPath = outDir.resolve(TestUtils.getTestName(testInfo) + ".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("CYP2D6\t*2/*10\tIntermediate Metabolizer\t1.25");
    }

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19");
    testWrapper.execute(outsideCallPath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testPrintCalls(DataSource.CPIC, "CYP2C19", "*38/*38");

    testWrapper.testNotCalledByMatcher("CYP2D6");
    testWrapper.testPrintCalls(DataSource.CPIC, "CYP2D6", "*2/*10");
    GeneReport cyp2d6Report = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2D6");
    assertNotNull(cyp2d6Report);
    assertEquals(1, cyp2d6Report.getSourceDiplotypes().size());
    assertTrue(cyp2d6Report.getSourceDiplotypes().stream().allMatch((d) -> d.getPhenotypes().contains("Intermediate Metabolizer")));

    testWrapper.testMessageCountForGene(DataSource.CPIC, "CYP2C19", 1);
    testWrapper.testGeneHasMessage(DataSource.CPIC, "CYP2C19", "reference-allele");
    testWrapper.testGeneHasMessage(DataSource.CPIC, "CYP2D6", MessageHelper.MSG_OUTSIDE_CALL);
  }


  @Test
  void testF5(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    String gene = "F5";
    testWrapper.getVcfBuilder()
        .reference(gene);
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("rs6025 C/rs6025 C");

    testWrapper.testCalledByMatcher(gene);
    testWrapper.testSourceDiplotypes(DataSource.DPWG, gene, expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.DPWG, gene, expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testPrintCalls(DataSource.DPWG, gene, expectedCalls);

    testWrapper.testMatchedAnnotations("hormonal contraceptives for systemic use", 1);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, gene, expectedCalls, "hormonal contraceptives for systemic use", RecPresence.NO, RecPresence.YES);
  }


  /**
   * Can we use phenotype for F5 DPWG matching?
   */
  @Test
  void testF5OutsideCallPhenotype(TestInfo testInfo) throws Exception {
    Path outDir = TestUtils.getTestOutputDir(testInfo, false);
    Path outsideCallFile = outDir.resolve(TestUtils.getTestName(testInfo) + ".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallFile))) {
      writer.println("F5\t\tFactor V Leiden absent");
    }

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.execute(outsideCallFile);

    String gene = "F5";
    List<String> expectedCalls = List.of("Factor V Leiden absent");

    testWrapper.testNotCalledByMatcher(gene);
    testWrapper.testSourceDiplotypes(DataSource.DPWG, gene, expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.DPWG, gene, expectedCalls);
    testWrapper.testReportable(gene);
    testWrapper.testPrintCalls(DataSource.DPWG, gene, expectedCalls);

    testWrapper.testMatchedAnnotations("hormonal contraceptives for systemic use", DataSource.CPIC, 0);
    testWrapper.testMatchedAnnotations("hormonal contraceptives for systemic use", DataSource.DPWG, 1);

    Document document = readHtmlReport(outsideCallFile);
    htmlCheckGene(document, gene, sf_notCalled);
    htmlCheckDrug(document, gene, expectedCalls, "hormonal contraceptives for systemic use",
        RecPresence.NO, RecPresence.YES);
  }


  @Test
  void testF5OutsideCallDiplotype(TestInfo testInfo) throws Exception {
    Path outDir = TestUtils.getTestOutputDir(testInfo, false);
    Path outsideCallFile = outDir.resolve(TestUtils.getTestName(testInfo) + ".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallFile))) {
      writer.println("F5\trs6025 C/rs6025 T (Factor V Leiden)");
    }

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.execute(outsideCallFile);

    String gene = "F5";
    List<String> expectedCalls = List.of("rs6025 C/rs6025 T (Factor V Leiden)");

    testWrapper.testNotCalledByMatcher(gene);
    testWrapper.testSourceDiplotypes(DataSource.DPWG, gene, expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.DPWG, gene, expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testSourcePhenotype(DataSource.DPWG, gene, "Factor V Leiden heterozygous");
    testWrapper.testRecommendedPhenotype(DataSource.DPWG, gene, "Factor V Leiden heterozygous");
    testWrapper.testReportable(gene);
    testWrapper.testPrintCalls(DataSource.DPWG, gene, expectedCalls);

    testWrapper.testMatchedAnnotations("hormonal contraceptives for systemic use", DataSource.CPIC, 0);
    testWrapper.testMatchedAnnotations("hormonal contraceptives for systemic use", DataSource.DPWG, 1);

    Document document = readHtmlReport(outsideCallFile);
    htmlChecks(document, gene, expectedCalls, "hormonal contraceptives for systemic use", RecPresence.NO, RecPresence.YES);
  }


  @Test
  void testWarfarinMissingRs12777823(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C9")
        .reference("CYP4F2")
        .reference("VKORC1")
        .missingExtraPosition("CYP2C9", "rs12777823")
    ;
    Path vcfFile = testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP2C9");
    testWrapper.testReportable("CYP2C9");

    Document document = readHtmlReport(vcfFile);
    assertNotNull(document.getElementById("CYP2C9"));

    Elements warfarinCpicDips = document.select(".cpic-warfarin .rx-dip");
    assertEquals(3, warfarinCpicDips.size());
    assertEquals(warfarinCpicDips.stream()
            .map(Element::text)
            .toList(),
        List.of("CYP2C9:*1/*1", "CYP4F2:*1/*1", "VKORC1: rs9923231 reference (C)/ rs9923231 reference (C)"));

    Elements cpicWarfarinHighlightedVars = document.select(".cpic-warfarin .rx-hl-var");
    assertEquals(1, cpicWarfarinHighlightedVars.size());
    assertEquals("rs12777823:Unknown", cpicWarfarinHighlightedVars.get(0).text());
  }


  @Test
  void testG6pdRef_male(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .male()
        .reference("G6PD");
    Path vcfFile = testWrapper.execute(null);

    testWrapper.testCalledByMatcher("G6PD");
    testWrapper.testReportable("G6PD");

    Document document = readHtmlReport(vcfFile);
    Elements g6pdSections = document.select(".gene.G6PD");
    assertEquals(1, g6pdSections.size());
    Elements g6pdCallElems = g6pdSections.get(0).getElementsByClass("genotype-result");
    assertEquals(1, g6pdCallElems.size());
    assertEquals("B (reference)", g6pdCallElems.text());
  }

  @Test
  void testG6pd_Ref_female(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .female()
        .reference("G6PD");
    Path vcfFile = testWrapper.execute(null);

    testWrapper.testCalledByMatcher("G6PD");
    testWrapper.testReportable("G6PD");

    Document document = readHtmlReport(vcfFile);
    Elements g6pdSections = document.select(".gene.G6PD");
    assertEquals(1, g6pdSections.size());
    Elements g6pdCallElems = g6pdSections.get(0).getElementsByClass("genotype-result");
    assertEquals(1, g6pdCallElems.size());
    assertEquals("B (reference)/B (reference)", g6pdCallElems.text());
  }

  @Test
  void testG6pd_MDPSCB_male(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .male()
        .reference("G6PD")
        .variation("G6PD", "rs5030868", "A");
    Path vcfFile = testWrapper.execute(null);

    String gene = "G6PD";
    List<String> expectedCalls = List.of("Mediterranean, Dallas, Panama, Sassari, Cagliari, Birmingham");

    testWrapper.testCalledByMatcher(gene);
    testWrapper.testSourceDiplotypes(DataSource.CPIC, gene, expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, gene, expectedCalls);
    testWrapper.testReportable(gene);
    testWrapper.testPrintCalls(DataSource.CPIC, gene, expectedCalls);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, gene, expectedCalls, "aspirin", RecPresence.YES, RecPresence.NO);
  }

  @Test
  void testG6pd_MDPSCB_female_homo(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("G6PD")
        .variation("G6PD", "rs5030868", "A", "A");
    Path vcfFile = testWrapper.execute(null);

    String gene = "G6PD";
    List<String> expectedCalls = List.of("Mediterranean, Dallas, Panama, Sassari, Cagliari, Birmingham/Mediterranean, Dallas, Panama, Sassari, Cagliari, Birmingham");

    testWrapper.testCalledByMatcher(gene);
    testWrapper.testSourceDiplotypes(DataSource.CPIC, gene, expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, gene, expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testReportable(gene);
    testWrapper.testPrintCalls(DataSource.CPIC, gene, expectedCalls);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, gene, expectedCalls, "aspirin", RecPresence.YES, RecPresence.NO);
  }

  @Test
  void testG6pd_MDPSCB_female_het(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("G6PD")
        .variation("G6PD", "rs5030868", "G", "A");
    Path vcfFile = testWrapper.execute(null);

    String gene = "G6PD";
    List<String> expectedCalls = List.of("B (reference)/Mediterranean, Dallas, Panama, Sassari, Cagliari, Birmingham");

    testWrapper.testCalledByMatcher(gene);
    testWrapper.testSourceDiplotypes(DataSource.CPIC, gene, expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, gene, expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testReportable(gene);
    testWrapper.testPrintCalls(DataSource.CPIC, gene, expectedCalls);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, gene, expectedCalls, "aspirin", RecPresence.YES, RecPresence.NO);
  }


  @Test
  void testG6pd_chatham_male(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .male()
        .reference("G6PD")
        .variation("G6PD", "rs5030869", "T");
    Path vcfFile = testWrapper.execute(null);

    String gene = "G6PD";
    List<String> expectedCalls = List.of("Chatham");

    testWrapper.testCalledByMatcher(gene);
    testWrapper.testSourceDiplotypes(DataSource.CPIC, gene, expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, gene, expectedCalls);
    testWrapper.testReportable(gene);
    testWrapper.testPrintCalls(DataSource.CPIC, gene, expectedCalls);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, gene, expectedCalls, "aspirin", RecPresence.YES, RecPresence.NO);
  }

  @Test
  void testG6pd_chatham_female_homo(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("G6PD")
        .variation("G6PD", "rs5030869", "T", "T");
    Path vcfFile = testWrapper.execute(null);

    String gene = "G6PD";
    List<String> expectedCalls = List.of("Chatham/Chatham");

    testWrapper.testCalledByMatcher(gene);
    testWrapper.testSourceDiplotypes(DataSource.CPIC, gene, expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, gene, expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testReportable(gene);
    testWrapper.testPrintCalls(DataSource.CPIC, gene, expectedCalls);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, gene, expectedCalls, "aspirin", RecPresence.YES, RecPresence.NO);
  }

  @Test
  void testG6pd_chatham_female_het(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("G6PD")
        .variation("G6PD", "rs5030869", "C", "T");
    Path vcfFile = testWrapper.execute(null);

    String gene = "G6PD";
    List<String> expectedCalls = List.of("B (reference)/Chatham");

    testWrapper.testCalledByMatcher(gene);
    testWrapper.testSourceDiplotypes(DataSource.CPIC, gene, expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, gene, expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testReportable(gene);
    testWrapper.testPrintCalls(DataSource.CPIC, gene, expectedCalls);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, gene, expectedCalls, "aspirin", RecPresence.YES, RecPresence.NO);
  }

  @Test
  void testG6pd_Arakawa_male(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .male()
        .reference("G6PD")
        .variation("G6PD", "chrX", 154532082, "A");
    Path vcfFile = testWrapper.execute(null);

    testWrapper.testCalledByMatcher("G6PD");
    testWrapper.testReportable("G6PD");

    Document document = readHtmlReport(vcfFile);
    Elements g6pdSections = document.select(".gene.G6PD");
    assertEquals(1, g6pdSections.size());
    Elements g6pdCallElems = g6pdSections.get(0).getElementsByClass("genotype-result");
    assertEquals(1, g6pdCallElems.size());
    assertEquals("Arakawa", g6pdCallElems.text());
  }

  @Test
  void testG6pd_Arakawa_female_het(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .female()
        .reference("G6PD")
        .variation("G6PD", "chrX", 154532082, "G", "A");
    Path vcfFile = testWrapper.execute(null);

    testWrapper.testCalledByMatcher("G6PD");
    testWrapper.testReportable("G6PD");

    Document document = readHtmlReport(vcfFile);
    Elements g6pdSections = document.select(".gene.G6PD");
    assertEquals(1, g6pdSections.size());
    Elements g6pdCallElems = g6pdSections.get(0).getElementsByClass("genotype-result");
    assertEquals(1, g6pdCallElems.size());
    assertEquals("Arakawa/B (reference)", g6pdCallElems.text());
  }

  @Test
  void testG6pd_Arakawa_female_homo(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .female()
        .reference("G6PD")
        .variation("G6PD", "chrX", 154532082, "A", "A");
    Path vcfFile = testWrapper.execute(null);

    testWrapper.testCalledByMatcher("G6PD");
    testWrapper.testReportable("G6PD");

    Document document = readHtmlReport(vcfFile);
    Elements g6pdSections = document.select(".gene.G6PD");
    assertEquals(1, g6pdSections.size());
    Elements g6pdCallElems = g6pdSections.get(0).getElementsByClass("genotype-result");
    assertEquals(1, g6pdCallElems.size());
    assertEquals("Arakawa/Arakawa", g6pdCallElems.text());
  }
  
  public static Document readHtmlReport(@Nullable Path file) throws IOException {
    if (file == null) {
      throw new IOException("No such file: " + file);
    }
    Path reporterOutput = file.getParent().resolve(BaseConfig.getBaseFilename(file) +
        BaseConfig.REPORTER_SUFFIX + ".html");
    return Jsoup.parse(reporterOutput.toFile());
  }
}
