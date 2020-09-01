
package pepmhc.calc;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import jam.app.JamLogger;
import jam.app.JamProperties;
import jam.math.IntUtil;

import jene.peptide.Peptide;

/**
 * Generates and maintains sets of peptides sampled from a larger pool.
 */
public final class PeptideSample {
    private final List<Set<Peptide>> sampleList;

    private static PeptideSample global = null;

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

    private PeptideSample(List<Set<Peptide>> sampleList) {
        this.sampleList = sampleList;
    }

    /**
     * Generates new non-overlapping samples from a flat file of
     * peptides.
     *
     * @param peptideFile the full path name of the flat file
     * containing the peptides to sample.
     *
     * @param sampleCount the number of distinct sample sets.
     *
     * @param sampleSize the number of peptides to include in each
     * sample.
     *
     * @return the new sample set.
     *
     * @throws RuntimeException if the total number of sampled
     * peptides ({@code sampleCount * sampleSize}) exceeds the
     * number in the flat file.
     */
    public static PeptideSample create(String peptideFile, int sampleCount, int sampleSize) {
        return create(Peptide.load(peptideFile), sampleCount, sampleSize);
    }

    /**
     * Generates new non-overlapping samples from a pool of peptides.
     *
     * @param peptidePool the peptides to sample from.
     *
     * @param sampleCount the number of distinct sample sets.
     *
     * @param sampleSize the number of peptides to include in each
     * sample.
     *
     * @return the new sample set.
     *
     * @throws RuntimeException if the total number of sampled
     * peptides ({@code sampleCount * sampleSize}) exceeds the
     * number in the input list.
     */
    public static PeptideSample create(List<Peptide> peptidePool, int sampleCount, int sampleSize) {
        int poolCount = peptidePool.size();
        int totalCount = sampleCount * sampleSize;

        if (totalCount > poolCount)
            throw new IllegalArgumentException("Number of peptide samples exceeds the pool size.");

        JamLogger.info("Sampling [%d] of [%d] peptides...", totalCount, poolCount);
        
        int[] poolIndexes = IntUtil.sample(poolCount, totalCount);
        List<Set<Peptide>> sampleList = new ArrayList<Set<Peptide>>(sampleCount);

        for (int sampleIndex = 0; sampleIndex < sampleCount; ++sampleIndex) {
            Set<Peptide> sampleSet = new HashSet<Peptide>(sampleSize);

            for (int peptideIndex = 0; peptideIndex < sampleSize; ++peptideIndex)
                sampleSet.add(peptidePool.get(poolIndexes[peptideIndex + sampleSize * sampleIndex]));

            sampleList.add(Collections.unmodifiableSet(sampleSet));
        }

        return new PeptideSample(sampleList);
    }

    /**
     * Returns the global sample set.
     *
     * @return the global sample set.
     */
    public static PeptideSample global() {
        if (global == null)
            global = createGlobal();

        return global;
    }

    private static PeptideSample createGlobal() {
        int sampleSize = resolveSampleSize();
        int sampleCount = resolveSampleCount();
        String peptideFile = resolvePeptideFile();

        return create(peptideFile, sampleCount, sampleSize);
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

    /**
     * Returns the number of independent samples.
     *
     * @return the number of independent samples.
     */
    public int countSamples() {
        return sampleList.size();
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
        //
        // The sample sets are already unmodifiable when this instance
        // is constructed...
        //
        return sampleList.get(sampleIndex);
    }
}
