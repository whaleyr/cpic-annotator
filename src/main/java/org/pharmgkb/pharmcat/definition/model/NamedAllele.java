package org.pharmgkb.pharmcat.definition.model;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.SortedSet;
import java.util.regex.Pattern;
import com.google.common.base.Objects;
import com.google.common.base.Preconditions;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.apache.commons.lang3.ObjectUtils;
import org.checkerframework.checker.nullness.qual.NonNull;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.pharmcat.haplotype.Iupac;
import org.pharmgkb.pharmcat.haplotype.MatchData;
import org.pharmgkb.pharmcat.util.HaplotypeNameComparator;


/**
 * A named allele (aka Haplotype, Star Allele, etc.).
 *
 * @author Ryan Whaley
 */
public class NamedAllele implements Comparable<NamedAllele> {
  @Expose
  @SerializedName("name")
  private final String m_name;
  @Expose
  @SerializedName("id")
  private final String m_id;
  @Expose
  @SerializedName("alleles")
  private final String[] m_alleles;
  @Expose
  @SerializedName("cpicAlleles")
  private final String[] m_cpicAlleles;
  @Expose
  @SerializedName(value = "reference", alternate = {"matchesreferencesequence"})
  private final boolean m_isReference;
  /**
   * This should always be false. The only time it might be true is on initial data ingestion by
   * {@link org.pharmgkb.pharmcat.util.DataManager}, after which they are immediately deleted.
   */
  @Expose
  @SerializedName("structuralVariant")
  private final boolean m_isStructuralVariant;
  //-- variables after this point are used by NamedAlleleMatcher --//
  /** The set of positions that are missing from this copy of the NamedAllele **/
  private final SortedSet<VariantLocus> m_missingPositions;
  // generated by initialize()
  // m_alleleMap and m_cpicAlleleMap must be HashMap so that it doesn't use VariantLocus.compareTo()
  // because VariantLocus.position can change
  private boolean m_isInitialized;
  private HashMap<VariantLocus, String> m_alleleMap;
  private HashMap<VariantLocus, String> m_cpicAlleleMap;
  private List<Integer> m_wobblePositions;
  private int m_score;
  private transient Pattern m_permutations;
  // generated by combinations code
  @Expose
  @SerializedName("numCombinations")
  private int m_numCombinations;
  @Expose
  @SerializedName("numPartials")
  private int m_numPartials;


  /**
   * Primary constructor.
   * Use this when reading in allele definitions.
   */
  public NamedAllele(String id, String name, String[] alleles, String[] cpicAlleles, boolean isReference) {
    this(id, name, alleles, cpicAlleles, Collections.emptySortedSet(), isReference, 0, 0);
  }

  public NamedAllele(String id, String name, String[] alleles, String[] cpicAlleles,
      SortedSet<VariantLocus> missingPositions, boolean isReference) {
    this(id, name, alleles, cpicAlleles, missingPositions, isReference, 0, 0);
  }

  /**
   * Constructor for duplicating/modifying a {@link NamedAllele}.
   */
  public NamedAllele(String id, String name, String[] alleles, String[] cpicAlleles,
      SortedSet<VariantLocus> missingPositions, boolean isReference, int numCombinations, int numPartials) {
    Preconditions.checkNotNull(id);
    Preconditions.checkNotNull(name);
    Preconditions.checkNotNull(cpicAlleles);
    Preconditions.checkNotNull(missingPositions);
    m_id = id;
    m_name = name;
    m_alleles = alleles;
    m_cpicAlleles = cpicAlleles;
    m_missingPositions = missingPositions;
    m_isReference = isReference;
    m_isStructuralVariant = false;
    m_numCombinations = numCombinations;
    m_numPartials = numPartials;

    // m_alleles will be null after initial import when removing ignored positions
    if (isReference && m_alleles != null && Arrays.stream(m_alleles).anyMatch(java.util.Objects::isNull)) {
      try {
        throw new IllegalStateException("check: " + name + " (" + id + ")\n" + String.join(", ", m_alleles));
      } catch (Exception ex) {
        ex.printStackTrace(System.out);
        System.exit(1);
      }
    }
  }


  /**
   * Call this to initialize {@link NamedAllele} for use after initial import.
   */
   void initializeForImport(VariantLocus[] refVariants) {
    Preconditions.checkNotNull(refVariants);
    Preconditions.checkNotNull(m_cpicAlleles);
    Preconditions.checkState(refVariants.length == m_cpicAlleles.length, "Mismatch of variants for " + this +
        ", " + refVariants.length + " reference variants does not match " + m_cpicAlleles.length);

    m_cpicAlleleMap = new HashMap<>();
    if (m_isReference) {
      for (int x = 0; x < refVariants.length; x += 1) {
        Preconditions.checkNotNull(refVariants[x]);
        Preconditions.checkNotNull(m_cpicAlleles[x]);
        m_cpicAlleleMap.put(refVariants[x], m_cpicAlleles[x]);
      }
    } else {
      for (int x = 0; x < refVariants.length; x += 1) {
        Preconditions.checkNotNull(refVariants[x]);
        m_cpicAlleleMap.put(refVariants[x], m_cpicAlleles[x]);
      }
    }
    m_isInitialized = true;
  }


  /**
   * Call this to initialize this object for use.
   */
  public void initialize(VariantLocus[] refVariants) {
    Preconditions.checkNotNull(refVariants);
    Preconditions.checkNotNull(m_alleles);
    Preconditions.checkState(refVariants.length == m_alleles.length, "Mismatch of variants for " +
        this + ", " + refVariants.length + " reference variants does not match " + m_alleles.length);

    if (m_isInitialized) {
      return;
    }
    m_alleleMap = new HashMap<>();
    m_cpicAlleleMap = new HashMap<>();
    m_wobblePositions = new ArrayList<>();
    for (int x = 0; x < refVariants.length; x += 1) {
      m_alleleMap.put(refVariants[x], m_alleles[x]);
      m_cpicAlleleMap.put(refVariants[x], m_cpicAlleles[x]);
      if (m_alleles[x] != null) {
        if (Iupac.isWobble(m_alleles[x])) {
          m_wobblePositions.add(x);
        }
        m_score++;
      }
    }
    m_score -= m_numPartials;
    calculatePermutations(refVariants);
    m_isInitialized = true;
  }

  /**
   * Call this to initialize this object for use.
   * This variant of {@link #initialize} allows an arbitrary {@code score} to be set.
   *
   * @param score the score for this {@link NamedAllele}
   */
  public void initialize(VariantLocus[] refVariants, int score) {

    initialize(refVariants);
    m_score = score;
  }


  /**
   * Take sample alleles into consideration for scoring.
   * This takes wobbles into consideration.
   */
  public int scoreForSample(MatchData matchData, Collection<String> sequences) {
    if (m_wobblePositions == null || m_wobblePositions.isEmpty()) {
      return m_score;
    }
    int score = m_score;
    for (Integer idx : m_wobblePositions) {
      VariantLocus vl = matchData.getPositions()[idx];
      int numRefs = 0;
      for (String seq : sequences) {
        String allele = matchData.getAllele(seq, idx);
        if (allele.equals(vl.getRef())) {
          numRefs += 1;
        }
      }
      // if all alleles at position is ref, don't score this wobble
      if (numRefs == sequences.size()) {
        score -= 1;
      }
    }
    return score;
  }


  /**
   * The name of this named allele (e.g. *1, Foo123Bar)
   */
  public String getName() {
    return m_name;
  }

  /**
   * The CPIC identifier for this named allele (e.g. CA10000.1)
   */
  public String getId() {
    return m_id;
  }

  public boolean isReference() {
    return m_isReference;
  }


  /**
   * Gets if this is a structural variant.
   * This should always be false. The only time it might be true is on initial data ingestion by
   * {@link org.pharmgkb.pharmcat.util.DataManager}, after which they are immediately deleted.
   */
  boolean isStructuralVariant() {
    return m_isStructuralVariant;
  }


  public boolean isCombination() {
    return m_numCombinations > 1;
  }

  public int getNumCombinations() {
    return m_numCombinations;
  }

  public boolean isPartial() {
    return m_numPartials > 0;
  }

  public int getNumPartials() {
    return m_numPartials;
  }


  /**
   * The array of alleles that define this allele.
   * <p>
   * <em>Note:</em> use this in conjunction with {@link DefinitionFile#getVariants()} to get the name of the variant
   */
  public String[] getAlleles() {
    return m_alleles;
  }

  public String getAllele(int idx) {
    return m_alleles[idx];
  }

  public @Nullable String getAllele(VariantLocus variantLocus) {
    Preconditions.checkState(m_isInitialized, "This NamedAllele has not been initialized()");
    return m_alleleMap.get(variantLocus);
  }

  public String[] getCpicAlleles() {
    return m_cpicAlleles;
  }

  public String getCpicAllele(int x) {
    return m_cpicAlleles[x];
  }

  public @Nullable String getCpicAllele(VariantLocus variantLocus) {
    Preconditions.checkState(m_isInitialized, "This NamedAllele has not been initialized()");
    return m_cpicAlleleMap.get(variantLocus);
  }


  /**
   * Gets the score (the number of alleles that matched) for this allele if it is matched.
   * It is usually the same as the number of non-null alleles, but can be set to anything via
   * {@link #initialize(VariantLocus[], int)}.
   */
  public int getScore() {
    return m_score;
  }


  /**
   * Gets the positions that are missing from this NamedAllele.
   */
  public SortedSet<VariantLocus> getMissingPositions() {
    //noinspection ReplaceNullCheck
    if (m_missingPositions == null) {
      // this is possible if marshalled via GSON
      return Collections.emptySortedSet();
    }
    return m_missingPositions;
  }


  @Override
  public String toString() {
    return m_name + " [" + m_id + "]";
  }

  @Override
  public int compareTo(@NonNull NamedAllele o) {
    if (m_isReference && !o.isReference()) {
      return -1;
    } else if (o.isReference() && !m_isReference) {
      return 1;
    }
    int rez = HaplotypeNameComparator.getComparator().compare(m_name, o.m_name);
    if (rez != 0) {
      return rez;
    }
    return ObjectUtils.compare(m_id, o.m_id);
  }


  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (!(o instanceof NamedAllele that)) {
      return false;
    }
    return Objects.equal(m_name, that.getName()) &&
        Objects.equal(m_id, that.getId()) &&
        Arrays.equals(m_alleles, that.getAlleles());
  }

  @Override
  public int hashCode() {
    return Objects.hashCode(m_name, m_id, m_alleles);
  }



  //-- permutation code --//

  public Pattern getPermutations() {
    return m_permutations;
  }


  private void calculatePermutations(VariantLocus[] refVariants) {

    List<VariantLocus> sortedRefVariants = Arrays.stream(refVariants).sorted().toList();
    StringBuilder builder = new StringBuilder();
    for (VariantLocus variant : sortedRefVariants) {
      builder.append(variant.getPosition())
          .append(":");
      String allele = m_alleleMap.get(variant);
      if (allele != null) {
        if (allele.length() == 1) {
          allele = Iupac.lookup(allele).getRegex();
        }
        builder.append(allele);
      } else {
        builder.append(".*?");
      }
      builder.append(";");
    }
    m_permutations = Pattern.compile(builder.toString());
  }
}
