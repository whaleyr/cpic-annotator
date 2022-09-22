package org.pharmgkb.pharmcat.reporter;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.cpic.Drug;
import org.pharmgkb.pharmcat.reporter.model.result.DrugReport;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.util.CliUtils;


/**
 * This class acts as a central context for all data needed to generate the final report.
 *
 * @author greytwist
 * @author Ryan Whaley
 */
public class ReportContext {
  @Expose
  @SerializedName("title")
  private final String f_title;
  @Expose
  @SerializedName("generatedOn")
  private final Date f_generatedOn = new Date();
  @Expose
  @SerializedName("pharmcatVersion")
  private final String f_pharmcatVersion = CliUtils.getVersion();
  @Expose
  @SerializedName("cpicVersion")
  private String m_cpicVersion;
  @Expose
  @SerializedName("genes")
  private final Map<DataSource, SortedMap<String, GeneReport>> m_geneReports;
  @Expose
  @SerializedName("drugs")
  private final Map<DataSource, SortedMap<String, DrugReport>> m_drugReports = new HashMap<>();
  @Expose
  @SerializedName("messages")
  private final List<MessageAnnotation> f_messages = new ArrayList<>();

  /**
   * Public constructor. Compiles all the incoming data into useful objects to be held for later reporting.
   *
   * @param geneReports {@link GeneReport} objects, non-null but can be empty
   * @param title the optional text to show as a user-friendly title or identifier for this report
   */
  public ReportContext(Env env, Map<DataSource, SortedMap<String, GeneReport>> geneReports, String title) throws IOException {
    f_title = title;
    m_geneReports = geneReports;

    // get CPIC drug data
    DrugCollection drugCollection = new DrugCollection();
    m_cpicVersion = validateCpicVersions(drugCollection);
    Map<String, DrugReport> cpicReports = m_drugReports.computeIfAbsent(DataSource.CPIC, (s) -> new TreeMap<>());
    for (Drug drug : drugCollection.listReportable()) {
      cpicReports.put(drug.getDrugName().toLowerCase(), new DrugReport(drug, this));
    }

    // get DPWG/PharmGKB drug data
    PgkbGuidelineCollection pgkbGuidelineCollection = new PgkbGuidelineCollection();
    Map<String, DrugReport> dpwgReports = m_drugReports.computeIfAbsent(DataSource.DPWG, (s) -> new TreeMap<>());
    // go through all DPWG-PharmGKB drugs, we iterate this way because one guideline may have multiple chemicals/drugs
    for (String drugName : pgkbGuidelineCollection.getChemicals()) {
      dpwgReports.put(drugName.toLowerCase(),
          new DrugReport(drugName, pgkbGuidelineCollection.findGuidelinePackages(drugName), this));
    }

    // now that all reports are generated, apply applicable messages
    MessageHelper messageHelper = env.getMessageHelper();
    // to gene reports
    geneReports.values().stream()
        .flatMap((m) -> m.values().stream())
        .forEach(messageHelper::addMatchingMessagesTo);
    // to drug reports
    for (DataSource source : m_drugReports.keySet()) {
      for (DrugReport drugReport : m_drugReports.get(source).values()) {
        messageHelper.addMatchingMessagesTo(drugReport, this, source);

        // add a message for any gene that has missing data
        drugReport.getRelatedGeneSymbols().stream()
            .map((s) -> getGeneReport(source, s))
            .filter((gr) -> gr != null && !gr.isOutsideCall() && gr.isMissingVariants())
            .forEach((gr) -> drugReport.addMessage(new MessageAnnotation(MessageAnnotation.TYPE_NOTE,
                "missing-variants",
                "Some position data used to define " + gr.getGeneDisplay() +
                    " alleles is missing which may change the matched genotype. See <a href=\"" +
                    gr.getGeneDisplay() + "\">" + gr.getGeneDisplay() +
                    "</a> in Section III for for more information.")));
      }
    }
  }

  private String validateCpicVersions(DrugCollection drugCollection) {
    Set<String> cpicVersions = new HashSet<>();
    // check GeneReports from the Phenotyper
    for (GeneReport geneReport : m_geneReports.get(DataSource.CPIC).values()) {
      if (geneReport.getPhenotypeVersion() != null) {
        cpicVersions.add(geneReport.getPhenotypeVersion());
      }
      if (geneReport.getAlleleDefinitionVersion() != null) {
        cpicVersions.add(geneReport.getAlleleDefinitionVersion());
      }
    }

    // check CPIC drug data
    for (Drug drug : drugCollection.list()) {
      if (drug.getCpicVersion() == null || drug.getCpicVersion().equals("n/a")) {
        continue;
      }
      cpicVersions.add(drug.getCpicVersion());
    }

    if (cpicVersions.size() == 0) {
      return TextConstants.NA;
    } else {
      if (cpicVersions.size() > 1) {
        addMessage(new MessageAnnotation(MessageAnnotation.TYPE_NOTE, "multiple-cpic-versions",
            "Multiple CPIC versions used to generate gene and drug reports: " + cpicVersions + "."
        ));
      }
      return String.join(", ", cpicVersions);
    }
  }



  /**
   * Gets the set of all {@link DrugReport} objects that hold drug information and thier recommendations
   * @return a set of {@link DrugReport} objects
   */
  public Map<DataSource, SortedMap<String, DrugReport>> getDrugReports() {
    return m_drugReports;
  }

  /**
   * Gets the set of all {@link GeneReport} objects that are reported in this context
   */
  public Map<DataSource, SortedMap<String, GeneReport>> getGeneReports() {
    return m_geneReports;
  }

  public List<DrugReport> getDrugReports(String drug) {
    return m_drugReports.keySet().stream()
        .map((k) -> m_drugReports.get(k).get(drug))
        .filter(Objects::nonNull)
        .toList();
  }

  public @Nullable DrugReport getDrugReport(DataSource source, String drug) {
    return m_drugReports.get(source).get(drug);
  }

  public List<GeneReport> getGeneReports(String gene) {
    return m_geneReports.keySet().stream()
        .map((k) -> m_geneReports.get(k).get(gene))
        .filter(Objects::nonNull)
        .toList();
  }

  public GeneReport getGeneReport(DataSource source, String gene) {
    GeneReport geneReport = m_geneReports.get(source).get(gene);
    if (geneReport == null) {
      throw new IllegalStateException("No gene report for " + gene);
    }
    return geneReport;
  }

  /**
   * The user-freindly title for the report
   * @return the title string
   */
  public String getTitle() {
    return f_title;
  }

  /**
   * Gets the timestamp this context was compiled
   * @return the timestamp this context was compiled
   */
  public Date getGeneratedOn() {
    return f_generatedOn;
  }

  /**
   * Gets the PharmCAT version tag this context was created with
   * @return a verstion tag string in the form vX.Y
   */
  public String getPharmcatVersion() {
    return f_pharmcatVersion;
  }

  public String getCpicVersion() {
    return m_cpicVersion;
  }

  public List<MessageAnnotation> getMessages() {
    return f_messages;
  }

  public void addMessage(MessageAnnotation message) {
    f_messages.add(message);
  }
}
