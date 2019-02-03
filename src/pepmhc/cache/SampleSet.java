
package pepmhc.cache;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import jam.app.JamLogger;
import jam.app.JamProperties;
import jam.peptide.Peptide;
import jam.util.ListUtil;

/**
 * Generates and maintains a global set of peptides sampled from a
 * larger pool.
 */
public final class SampleSet {
    private final int sampleCount;
    private final int sampleSize;
    private final String peptideFile;
    private final boolean withReplacement;

    private final List<Set<Peptide>> sampleList;

    private static SampleSet global = null;

    /**
     * Name of the system property that specifies the full path
     * name of the flat file containing the peptides to sample.
     */
    public static final String PEPTIDE_FILE_PROPERTY = "pepmhc.sample.peptideFile";

    /**
     * Name of the system property that specifies the number of
     * distinct peptide samples to generate.
     */
    public static final String SAMPLE_COUNT_PROPERTY = "pepmhc.sample.sampleCount";

    /**
     * Name of the system property that specifies the number of
     * peptides to include in each sample.
     */
    public static final String SAMPLE_SIZE_PROPERTY = "pepmhc.sample.sampleSize";

    /**
     * Name of the system property that specifies whether to sample
     * with replacement or require all sampled peptides to be unique.
     */
    public static final String WITH_REPLACEMENT_PROPERTY = "pepmhc.sample.withReplacement";

    private SampleSet() {
        this.sampleSize = resolveSampleSize();
        this.sampleCount = resolveSampleCount();
        this.peptideFile = resolvePeptideFile();
        this.withReplacement = resolveWithReplacement();

        this.sampleList = new ArrayList<Set<Peptide>>(sampleCount);
        samplePeptides();
    }

    private static String resolvePeptideFile() {
        return JamProperties.getRequired(PEPTIDE_FILE_PROPERTY);
    }

    private static int resolveSampleCount() {
        return JamProperties.getRequiredInt(SAMPLE_COUNT_PROPERTY);
    }

    private static int resolveSampleSize() {
        return JamProperties.getRequiredInt(SAMPLE_SIZE_PROPERTY);
    }

    private static boolean resolveWithReplacement() {
        return JamProperties.getOptionalBoolean(WITH_REPLACEMENT_PROPERTY, true);
    }

    private void samplePeptides() {
        List<Peptide> allPeptides = Peptide.loadFlatFile(peptideFile);

        if (withReplacement)
            sampleWithReplacement(allPeptides);
        else
            sampleUnique(allPeptides);
    }

    private void sampleWithReplacement(List<Peptide> allPeptides) {
        while (sampleList.size() < sampleCount)
            sampleList.add(sampleUnique(allPeptides, sampleSize));
    }

    private static Set<Peptide> sampleUnique(List<Peptide> allPeptides, int uniqueSize) {
        JamLogger.info("Sampling [%d] of [%d] peptides...", uniqueSize, allPeptides.size());

        Set<Peptide> sampledPeptides = new LinkedHashSet<Peptide>(uniqueSize);

        while (sampledPeptides.size() < uniqueSize)
            sampledPeptides.add(ListUtil.select(allPeptides));

        return sampledPeptides;
    }

    private void sampleUnique(List<Peptide> allPeptides) {
        //
        // Ensure that the peptides are unique across all samples by
        // sampling all peptides first...
        //
        Set<Peptide> allSampled = sampleUnique(allPeptides, sampleCount * sampleSize);

        // Now assign each peptide to a sample...
        int sampleIndex = 0;
        int peptideIndex = 0;

        for (Peptide peptide : allSampled) {
            if (peptideIndex == 0)
                sampleList.add(new LinkedHashSet<Peptide>(sampleSize));

            sampleList.get(sampleIndex).add(peptide);
            ++peptideIndex;

            if (peptideIndex == sampleSize) {
                peptideIndex = 0;
                ++sampleIndex;
            }
        }
    }

    /**
     * Returns the global sample set.
     *
     * @return the global sample set.
     */
    public static SampleSet global() {
        if (global == null)
            global = new SampleSet();

        return global;
    }

    /**
     * Returns a read-only view of a peptide sample.
     *
     * @param sampleIndex the index of the sample to view.
     *
     * @return a read-only view of the specified sample.
     *
     * @throws IndexOutOfBoundsException if the sample index is not
     * valid.
     */
    public Set<Peptide> viewSample(int sampleIndex) {
        return Collections.unmodifiableSet(sampleList.get(sampleIndex));
    }
}
