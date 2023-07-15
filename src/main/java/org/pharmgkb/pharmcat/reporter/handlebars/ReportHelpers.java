package org.pharmgkb.pharmcat.reporter.handlebars;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import com.google.common.html.HtmlEscapers;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.format.html.Report;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;
import org.pharmgkb.pharmcat.reporter.model.cpic.Publication;
import org.pharmgkb.pharmcat.reporter.model.result.AnnotationReport;
import org.pharmgkb.pharmcat.reporter.model.result.CallSource;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.Genotype;
import org.pharmgkb.pharmcat.reporter.model.result.GuidelineReport;

import static org.pharmgkb.pharmcat.reporter.TextConstants.*;
import static org.pharmgkb.pharmcat.reporter.caller.DpydCaller.isDpyd;


/**
 * Class to hold methods that can be used as helper methods in the Handlebars templating system.
 *
 * @author Ryan Whaley
 */
@SuppressWarnings("unused")
public class ReportHelpers {
  private static final String sf_variantAlleleTemplate = "<td class=\"%s\">%s%s</td>";
  private static boolean s_debugMode;

  /**
   * Private constructor for utility class.
   */
  private ReportHelpers() {
  }


  public static void setDebugMode(boolean debugMode) {
    s_debugMode = debugMode;
  }



  public static String listGenes(Collection<String> genes) {
    StringBuilder builder = new StringBuilder();
    for (String gene : genes) {
      if (builder.length() > 0) {
        builder.append(", ");
      }
      builder.append("<a href=\"#")
          .append(gene)
          .append("\">")
          .append(gene)
          .append("</a>");
    }
    return builder.toString();
  }

  public static String listSources(Collection<String> urls) {
    StringBuilder builder = new StringBuilder();
    int x = 0;
    for (String url : urls) {
      x += 1;
      if (x > 1) {
        builder.append(", ");
      }
      builder.append(externalHref(url, Integer.toString(x)));
    }
    return builder.toString();
  }


  public static String printRecMap(Map<String, String> data) {
    if (data.size() == 1) {
      return "<p>" + capitalizeNA(data.values().iterator().next()) + "</p>";
    }

    StringBuilder builder = new StringBuilder()
        .append("<dl class=\"compact mt-0\">");
    for (String key : data.keySet()) {
      builder.append("<dt>")
          .append(key)
          .append(":</dt><dd>")
          .append(capitalizeNA(data.get(key)))
          .append("</dd>");
    }
    builder.append("</dl>");
    return builder.toString();
  }


  public static String pluralize(String word, Object col) {
    if (moreThanOne(col)) {
      return word + "s";
    }
    return word;
  }


  private static final Pattern sf_4chPattern = Pattern.compile(".{9}|.+$");

  public static String variantAlleles(VariantReport variantReport) {
    String cellStyle = variantReport.isNonWildType() ? "nonwild" : "";
    String mismatch = variantReport.isHasUndocumentedVariations() ? "<div class=\"callMessage\">Undocumented variation</div>" : "";
    if (variantReport.isHasUndocumentedVariations()) {
      cellStyle = StringUtils.strip(cellStyle + " mismatch");
    }
    String call = formatCall(variantReport.getCall());
    return String.format(sf_variantAlleleTemplate, cellStyle, call, mismatch);
  }

  public static String wildtypeAllele(VariantReport variantReport) {
    return formatCall(variantReport.getWildTypeAllele());
  }

  private static String formatCall(String call) {
    if (call.length() <= 9) {
      return call;
    }
    String[] alleles = call.split("/");
    StringBuilder builder = new StringBuilder();
    for (String allele : alleles) {
      if (builder.length() > 0) {
        builder.append("/<br />");
      }
      if (allele.length() > 8) {
        Matcher m = sf_4chPattern.matcher(allele);
        List<String> parts = new ArrayList<>();
        while (m.find()) {
          parts.add(m.group());
        }
        builder.append(String.join("<br />", parts));
      } else {
        builder.append(allele);
      }
    }
    return builder.toString();
  }


  public static boolean gsShowDrug(String drug, Collection<String> drugs) {
    return drugs.contains(drug);
  }


  public static String gsCall(Diplotype diplotype) {
    if (diplotype.isUnknownAlleles()) {
      return "<span class=\"gs-uncalled-" + HtmlEscapers.htmlEscaper().escape(diplotype.getGene()) + "\">" + UNCALLED + "</span>";
    }
    if (diplotype.getInferredSourceDiplotypes() != null) {
      if (diplotype.getInferredSourceDiplotypes().size() > 1) {
        throw new IllegalStateException("Cannot determine genotype call (inferred from " +
            diplotype.getInferredSourceDiplotypes().size() + " diplotypes!)");
      }
      diplotype = diplotype.getInferredSourceDiplotypes().first();
    }
    return HtmlEscapers.htmlEscaper().escape(diplotype.getLabel());
  }

  public static String gsFunction(Diplotype diplotype) {
    if (diplotype.isCombination() && isDpyd(diplotype.getGene())) {
      return TextConstants.SEE_DRUG;
    }

    String f1 = diplotype.getAllele1() != null && diplotype.getAllele1().getFunction() != null ?
        diplotype.getAllele1().getFunction() : null;
    String f2 = diplotype.getAllele2() != null && diplotype.getAllele2().getFunction() != null ?
        diplotype.getAllele2().getFunction() : null;
    boolean isSinglePloidy = diplotype.getAllele1() != null && diplotype.getAllele2() == null;
    if (isSinglePloidy) {
      if (StringUtils.isNotBlank(f1)) {
        return f1;
      }
    } else {
      if (StringUtils.isNotBlank(f1) && StringUtils.isNotBlank(f2)) {
        if (f1.equals(f2)) {
          return String.format("Two %s alleles", f1);
        }
        else {
          String[] functions = new String[]{f1, f2};
          Arrays.sort(functions);
          return String.format("One %s allele and one %s allele", functions[0], functions[1]);
        }
      }
    }
    return TextConstants.NA.toUpperCase();
  }

  public static String gsPhenotype(Diplotype diplotype) {
    if (isDpyd(diplotype.getGene())) {
      return TextConstants.SEE_DRUG;
    }
    if (diplotype.getPhenotypes().size() == 0) {
      return TextConstants.NA.toUpperCase();
    }
    return String.join("; ", diplotype.getPhenotypes());
  }


  public static String rxAnnotationClass(DataSource source, String drug) {
    StringBuilder builder = new StringBuilder();
    if (source == DataSource.CPIC) {
      builder.append("cpic-");
    } else {
      builder.append("dpwg-");
    }
    builder.append(sanitizeCssSelector(drug));
    return builder.toString();
  }

  public static boolean rxIsCpicWarfarin(String drug, DataSource source) {
    return drug.equals("warfarin") && source == DataSource.CPIC;
  }

  public static boolean rxDpydInferred(Genotype genotype) {
    return genotype.isInferred() && genotype.getDiplotypes().stream()
        .map(Diplotype::getGene)
        .anyMatch(g -> g.equals("DPYD"));
  }

  public static boolean rxInferred(Genotype genotype) {
    return genotype.isInferred() && genotype.getDiplotypes().stream()
        .map(Diplotype::getGene)
        .noneMatch(g -> g.equals("DPYD"));
  }

  public static String rxGenotype(Genotype genotype, AnnotationReport annotationReport,
      GuidelineReport guidelineReport) {
    if (genotype.getDiplotypes().size() == 0) {
      return TextConstants.UNKNOWN_GENOTYPE;
    }
    StringBuilder builder = new StringBuilder()
      .append(renderRxDiplotypes(genotype.getDiplotypes(), false));
    if (annotationReport.getHighlightedVariants().size() > 0) {
      for (String var : annotationReport.getHighlightedVariants()) {
        if (builder.length() > 0) {
          builder.append(";<br />");
        }
        builder.append("<span class=\"rx-hl-var\">")
            .append(var)
            .append("</span>");
      }
    }
    return builder.toString();
  }

  public static String rxGenotypeDebug(Genotype genotype, GuidelineReport guidelineReport) {
    if (!s_debugMode) {
      return "";
    }

    boolean hasInferred = genotype.getDiplotypes().stream()
        .anyMatch(d -> d.getInferredSourceDiplotypes() != null && d.getInferredSourceDiplotypes().size() > 0);
    if (hasInferred) {
      return "<div class=\"alert alert-debug\">" +
          "<div class=\"hint\">Inferred:</div>" +
          "<span class=\"nowrap\">" +
          renderRxDiplotypes(genotype.getDiplotypes(), true) +
          "</span></div>";
    }
    return "";
  }

  public static String rxUnmatchedDiplotypes(SortedSet<Diplotype> diplotypes) {
    return renderRxDiplotypes(diplotypes, true, false, "rx-unmatched-dip");
  }

  public static boolean rxUnmatchedDiplotypesInferred(Report report) {
    // checking this to add top padding to account for superscript dagger
    return report.isUnmatchedInferred() || report.isUnmatchedDpydInferred();
  }


  private static String renderRxDiplotypes(Collection<Diplotype> diplotypes, boolean forDebug) {
    return renderRxDiplotypes(diplotypes, false, forDebug, "rx-dip");
  }

  private static String renderRxDiplotypes(Collection<Diplotype> diplotypes, boolean noLengthLimit, boolean forDebug,
      String dipClass) {
    SortedSet<Diplotype> displayDiplotypes = new TreeSet<>();
    for (Diplotype diplotype : diplotypes) {
      if (!forDebug && diplotype.getInferredSourceDiplotypes() != null) {
        displayDiplotypes.addAll(diplotype.getInferredSourceDiplotypes());
      } else {
        displayDiplotypes.add(diplotype);
      }
    }

    StringBuilder builder = new StringBuilder();
    for (Diplotype diplotype : displayDiplotypes) {
      if (builder.length() > 0) {
        builder.append(";<br />");
      }
      builder.append("<span");
      if (!forDebug) {
        builder.append(" class=\"")
            .append(dipClass)
            .append("\"");
      }
      builder.append("><a href=\"#")
          .append(diplotype.getGene())
          .append("\">")
          .append(diplotype.getGene())
          .append("</a>:");
      String call = diplotype.getLabel();
      if (noLengthLimit || call.length() <= 15) {
        builder.append(call);
      } else {
        int idx = call.indexOf("/");
        if (idx == -1) {
          builder.append("<br />")
              .append(call);
        } else {
          String a = call.substring(0, idx + 1);
          String b = call.substring(idx + 1);
          if (a.length() > 15) {
            builder.append("<br />");
          }
          builder.append(a)
              .append("<br />")
              .append(b);
        }
      }
      builder.append("</span>");
    }
    return builder.toString();
  }



  public static String rxImplications(SortedMap<String, String> implications) {
    if (implications.size() == 0) {
      return "";
    }
    if (implications.size() == 1) {
      Map.Entry<String, String> entry = implications.entrySet().iterator().next();
      return entry.getKey() + ": " + capitalizeNA(entry.getValue());
    }
    StringBuilder builder = new StringBuilder()
        .append("<ul class=\"noPadding mt-0\">");
    for (Map.Entry<String, String> entry : implications.entrySet()) {
      builder.append("<li>")
          .append(entry.getKey())
          .append(": ")
          .append(capitalizeNA(entry.getValue()))
          .append("</li>");
    }
    builder.append("</ul>");
    return builder.toString();
  }

  public static List<MessageAnnotation> rxAnnotationMessages(AnnotationReport annotationReport) {
    return annotationReport.getMessages().stream()
        .filter(MessageAnnotation.isMessage)
        .toList();
  }


  public static String amdSubtitle(GeneReport geneReport) {
    StringBuilder builder = new StringBuilder();

    if (isDpyd(geneReport.getGene()) && geneReport.getMatcherComponentHaplotypes().size() == 0) {
      builder.append("Haplotype");
    } else {
      builder.append("Genotype");
    }
    if (geneReport.getSourceDiplotypes().size() > 1) {
      builder.append("s");
    }
    builder.append(" ");
    if (geneReport.isOutsideCall()) {
      builder.append("Reported");
    } else {
      builder.append("Matched");
    }
    return builder.toString();
  }

  public static boolean amdNoCall(GeneReport report) {
    if (report.getCallSource() == CallSource.NONE) {
      return true;
    } else if (report.getCallSource() == CallSource.MATCHER) {
      return report.getVariantReports().size() == 0 ||
          report.getVariantReports().stream().allMatch(VariantReport::isMissing);
    }
    return false;
  }

  public static String amdAlleleFunction(Map<String, Map<String, String>> functionMap, String gene, String allele) {
    Map<String, String> functions = functionMap.get(gene);
    if (functions != null) {
      String function = functions.get(allele);
      if (function != null) {
        return allele + " - " + function;
      }
    }
    return allele;
  }

  public static String amdNoDataMessage(Collection<String> compactNoDataGenes) {
    return "<p class=\"noGeneData\">No data provided for " +
        compactNoDataGenes.stream()
            .map((g) -> "<span class=\"gene " + g.toLowerCase() + "\"><span class=\"no-data\">" + g + "</span></span>")
            .collect(Collectors.joining(", ")) +
        ".</p>";
  }

  public static boolean amdIsSingleCall(GeneReport report) {
    return amdGeneCalls(report).size() == 1;
  }

  public static List<String> amdGeneCalls(GeneReport geneReport) {
    if (geneReport.getCallSource() == CallSource.NONE) {
      return LIST_UNCALLED_NO_DATA;
    }
    if (geneReport.getCallSource() == CallSource.MATCHER) {
      if (geneReport.isNoData()) {
        return LIST_UNCALLED_NO_DATA;
      } else if (!geneReport.isReportable()) {
        return LIST_UNCALLED;
      }
    }
    return geneReport.getRecommendationDiplotypes().stream()
        .flatMap(d -> d.getInferredSourceDiplotypes() == null ? Stream.of(d)
            : d.getInferredSourceDiplotypes().stream())
        .map(Diplotype::getLabel)
        .toList();
  }

  public static String amdGeneCall(GeneReport report) {
    return amdGeneCalls(report).get(0);
  }


  public static String amdPhaseStatus(GeneReport geneReport) {
    if (geneReport.isOutsideCall()) {
      return "Unavailable for calls made outside PharmCAT";
    }
    return geneReport.isPhased() ? "Phased" : "Unphased";
  }

  public static boolean amdShowUnphasedNote(GeneReport geneReport) {
    return !geneReport.isPhased() && !isDpyd(geneReport.getGeneDisplay());
  }

  public static boolean amdHasUncalledHaps(GeneReport geneReport) {
    return geneReport.getUncalledHaplotypes() != null &&
        geneReport.getUncalledHaplotypes().size() > 0;
  }

  public static String amdUncalledHaps(GeneReport geneReport) {
    return String.join(", ", geneReport.getUncalledHaplotypes());
  }


  public static long amdTotalMissingVariants(GeneReport geneReport) {
    return geneReport.getVariantReports().stream()
        .filter(VariantReport::isMissing)
        .count();
  }

  public static int amdTotalVariants(GeneReport geneReport) {
    return geneReport.getVariantReports().size();
  }

  public static List<MessageAnnotation> amdMessages(GeneReport geneReport) {
    return geneReport.getMessages().stream()
        .filter(MessageAnnotation.isMessage)
        .toList();
  }

  public static List<String> amdExtraPositionNotes(GeneReport geneReport) {
    return geneReport.getMessages().stream()
        .filter(MessageAnnotation.isExtraPositionNote)
        .map(MessageAnnotation::getMessage)
        .toList();
  }


  public static boolean moreThanOne(Object obj) {
    if (obj == null) {
      return false;
    }
    //noinspection rawtypes
    if (obj instanceof Collection col) {
      return col.size() > 1;
    }
    //noinspection rawtypes
    if (obj instanceof Map map) {
      return map.size() > 1;
    }
    return false;
  }

  public static boolean thereAre(Object obj) {
    if (obj == null) {
      return false;
    }
    //noinspection rawtypes
    if (obj instanceof Collection col) {
      return col.size() > 0;
    }
    //noinspection rawtypes
    if (obj instanceof Map map) {
      return map.size() > 0;
    }
    return false;
  }


  public static boolean contains(Collection col, Object item) {
    return col.contains(item);
  }

  public static int add(int x, int y) {
    return x + y;
  }

  public static String externalHref(String href, String text) {
    return "<a href=\"" + href + "\" target=\"_blank\" rel=\"noopener noreferrer\">" + text + "</a>";
  }

  public static String sanitizeCssSelector(String value) {
    return value.replaceAll("\\W+", "_")
        .replaceAll("_+", "_")
        .replaceAll("^_+", "")
        .replaceAll("_+$", "");
  }

  private static final Pattern sf_endsWithPuncuationPattern = Pattern.compile(".*?\\p{Punct}$");

  public static String printCitation(Publication pub) {
    String url = "https://www.ncbi.nlm.nih.gov/pubmed/" + pub.getPmid();
    String period = "";
    if (!sf_endsWithPuncuationPattern.matcher(pub.getTitle()).matches()) {
      period = ".";
    }
    return externalHref(url, pub.getTitle()) + period + " <i>" + pub.getJournal() + "</i>. " + pub.getYear() +
        ". PMID:" + pub.getPmid();
  }

  public static String capitalizeNA(String text) {
    if (TextConstants.NA.equalsIgnoreCase(text)) {
      return "N/A";
    }
    return text;
  }

  public static String messageClass(MessageAnnotation msg) {
    return msg.getName()
        .replaceAll("\\p{Punct}", "-")
        .replaceAll("\\s+", "-")
        .replaceAll("-+", "-");
  }
}
