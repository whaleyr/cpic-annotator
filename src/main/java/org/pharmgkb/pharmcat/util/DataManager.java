package org.pharmgkb.pharmcat.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.lang.invoke.MethodHandles;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Properties;
import java.util.zip.ZipFile;
import javax.annotation.Nonnull;
import com.google.common.base.Preconditions;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.http.HttpEntity;
import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;
import org.pharmgkb.common.io.util.CliHelper;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.definition.CuratedDefinitionParser;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;


/**
 * This class manages external resources (e.g. allele definition files, dosing guideline annotations)
 *
 * @author Mark Woon
 */
public class DataManager {
  public static final Path DEFAULT_DEFINITION_DIR = PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/alleles");
  public static final Path DEFAULT_REPORTER_DIR = PathUtils.getPathToResource("org/pharmgkb/pharmcat/reporter");
  public static final Path DEFAULT_GUIDELINE_DIR = DEFAULT_REPORTER_DIR.resolve("guidelines");
  public static final String DOSING_GUIDELINE_URL = "https://api.pharmgkb.org/v1/download/file/data/dosingGuidelines.json.zip?ref=pharmcat";
  public static final String EXEMPTIONS_JSON_FILE_NAME = "exemptions.json";
  public static final String MESSAGES_JSON_FILE_NAME = "messages.json";
  private String m_googleUser;
  private String m_googleKey;
  private DataSerializer m_dataSerializer = new DataSerializer();
  private boolean m_verbose;



  private DataManager(Path propertyFile, boolean verbose) throws IOException {

    Properties properties = new Properties();
    try (BufferedReader reader = Files.newBufferedReader(propertyFile)) {
      properties.load(reader);
    }
    m_googleUser = StringUtils.stripToNull((String)properties.get("google.user"));
    Preconditions.checkState(m_googleUser != null, "Missing property: 'google.user");
    m_googleKey = StringUtils.stripToNull((String)properties.get("google.key"));
    Preconditions.checkState(m_googleKey != null, "Missing property: 'google.key");
    m_verbose = verbose;
  }


  public static void main(String[] args) {

    try {
      CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
          .addOption("p", "properties-file", "PharmCAT properties file", false, "p")
          .addOption("dl", "download-dir", "directory to save downloaded files", false, "dl")
          .addOption("a", "alleles-dir", "directory to save generated allele definition files", false, "a")
          .addOption("m", "messages-dir", "directory to write messages to", false, "m")
          .addOption("g", "guidelines-dir", "directory to save guideline annotations to", false, "g")
          .addOption("sd", "skip-download", "skip downloading")
          .addOption("sa", "skip-alleles", "skip alleles")
          .addOption("sm", "skip-messages", "skip messages")
          .addOption("sg", "skip-guidelines", "skip guidelines");


      if (!cliHelper.parse(args)) {
        System.exit(1);
      }

      Path propsFile = CliUtils.getPropsFile(cliHelper, "p");
      Path downloadDir;
      if (cliHelper.hasOption("dl")) {
        downloadDir = cliHelper.getValidDirectory("dl", true);
      } else {
        downloadDir = Files.createTempDirectory("pharmcat");
        downloadDir.toFile().deleteOnExit();
      }

      try {
        Path allelesDir;
        if (cliHelper.hasOption("a")) {
          allelesDir = cliHelper.getValidDirectory("a", true);
        } else {
          allelesDir = DEFAULT_DEFINITION_DIR;
        }
        Path messageDir;
        if (cliHelper.hasOption("m")) {
          messageDir = cliHelper.getValidDirectory("m", true);
        } else {
          messageDir = DEFAULT_REPORTER_DIR;
        }
        Path guidelinesDir;
        if (cliHelper.hasOption("g")) {
          guidelinesDir = cliHelper.getValidDirectory("g", true);
        } else {
          guidelinesDir = DEFAULT_GUIDELINE_DIR;
          if (Files.exists(guidelinesDir)) {
            if (!Files.isDirectory(guidelinesDir)) {
              System.out.println(guidelinesDir + " is not a directory");
              return;
            }
          } else {
            Files.createDirectories(guidelinesDir);
          }
        }

        Path exemptionsTsv = downloadDir.resolve("exemptions.tsv");
        Path exemptionsJson = allelesDir.resolve(EXEMPTIONS_JSON_FILE_NAME);
        Path messagesTsv = downloadDir.resolve("messages.tsv");
        Path messagesJson = messageDir.resolve(MESSAGES_JSON_FILE_NAME);
        Path guidelinesZip = downloadDir.resolve("guidelines.zip");

        DataManager manager = new DataManager(propsFile, cliHelper.isVerbose());
        if (!cliHelper.hasOption("sd")) {
          manager.download(downloadDir, exemptionsTsv, messagesTsv);
          if (!cliHelper.hasOption("sg")) {
            manager.downloadGuidelines(guidelinesZip);
          }
        }

        if (!cliHelper.hasOption("sa")) {
          manager.transformAlleleDefinitions(downloadDir, allelesDir);
          manager.transformExemptions(exemptionsTsv, exemptionsJson);
        }
        if (!cliHelper.hasOption("sm")) {
          manager.transformMessages(messagesTsv, messagesJson);
        }
        if (!cliHelper.hasOption("sg")) {
          manager.transformGuidelines(guidelinesZip, guidelinesDir);
        }

      } finally {
        if (!cliHelper.hasOption("dl")) {
          FileUtils.deleteQuietly(downloadDir.toFile());
        }
      }

    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }


  private void download(@Nonnull Path downloadDir, @Nonnull Path exemptionsFile, @Nonnull Path messagesFile) throws Exception {

    SheetsHelper sh = new SheetsHelper(m_googleUser, m_googleKey);
    sh.downloadAlleleDefinitions(downloadDir);
    sh.downloadAlleleExemptionsFile(exemptionsFile);
    sh.downloadMessagesFile(messagesFile);
  }


  /**
   * Does the work for stepping through the files and applying the format.
   */
  private void transformAlleleDefinitions(@Nonnull Path downloadDir, @Nonnull Path outDir) throws IOException {

    System.out.println();
    System.out.println("Saving allele definitions in " + outDir.toString());
    try (DirectoryStream<Path> files = Files.newDirectoryStream(downloadDir, f -> f.toString().endsWith("_translation.tsv"))) {
      for (Path file : files) {
        if (m_verbose) {
          System.out.println("Parsing " + file);
        }
        CuratedDefinitionParser parser = new CuratedDefinitionParser(file);

        DefinitionFile definitionFile = parser.parse();
        if (!parser.getWarnings().isEmpty()) {
          System.out.println("Warnings for " + file);
          parser.getWarnings()
              .forEach(System.out::println);
        }

        Path jsonFile = outDir.resolve(PathUtils.getBaseFilename(file) + ".json");
        m_dataSerializer.serializeToJson(definitionFile, jsonFile);
        if (m_verbose) {
          System.out.println("Wrote " + jsonFile);
        }
      }
    }
  }

  private void transformExemptions(@Nonnull Path tsvFile, @Nonnull Path jsonFile) throws IOException {

    System.out.println();
    System.out.println("Saving exemptions to " + jsonFile.toString());
    m_dataSerializer.serializeToJson(m_dataSerializer.deserializeExemptionsFromTsv(tsvFile), jsonFile);
  }

  private void transformMessages(@Nonnull Path tsvFile, @Nonnull Path jsonFile) throws IOException {

    System.out.println();
    System.out.println("Saving messages to " + jsonFile.toString());
    m_dataSerializer.serializeToJson(m_dataSerializer.deserializeMessagesFromTsv(tsvFile), jsonFile);
  }


  private void downloadGuidelines(@Nonnull Path guidelinesZip) throws Exception {

    try (CloseableHttpClient client = HttpClients.createDefault()) {
      try (CloseableHttpResponse response = client.execute(new HttpGet(DOSING_GUIDELINE_URL))) {
        HttpEntity entity = response.getEntity();
        if (entity != null) {
          try (InputStream in = entity.getContent();
               OutputStream out = Files.newOutputStream(guidelinesZip)) {
            IOUtils.copy(in, out);
          }
        }
      }
    }
  }

  private void transformGuidelines(@Nonnull Path zipFile, @Nonnull Path guidelinesDir) throws IOException {

    System.out.println();
    System.out.println("Saving guidelines to " + guidelinesDir.toString());
    try (ZipFile zip = new ZipFile(zipFile.toFile())) {
      zip.stream()
          .filter(ze -> {
            String name = ze.getName().toLowerCase();
            return name.matches("cpic_.*\\.json") || name.matches("created_.*\\.txt");
          })
          .forEachOrdered(ze -> {
            try (InputStream inputStream = zip.getInputStream(ze);
                 OutputStream out = Files.newOutputStream(guidelinesDir.resolve(ze.getName()))) {
              if (m_verbose) {
                System.out.println("Extracting " + ze.getName());
              }
              IOUtils.copy(inputStream, out);
            } catch (IOException ex) {
              throw new RuntimeException("Error extracting " + ze.getName(), ex);
            }
          });
    }
  }
}
