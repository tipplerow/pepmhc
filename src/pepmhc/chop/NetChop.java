
package pepmhc.chop;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import jam.app.JamEnv;
import jam.app.JamProperties;
import jam.math.DoubleRange;
import jam.math.IntRange;
import jam.math.IntUtil;
import jam.peptide.Peptide;
import jam.util.RegexUtil;

/**
 * Simulates proteasomal processing of a peptide: Predicts cleavage
 * sites using {@code netchop} and assembles cleaved fragments with
 * prescribed lengths.
 */
public final class NetChop {
    private final Peptide peptide;
    private final int[]   lengths;
    private final double  threshold;

    private List<Double> scores;
    private List<Peptide> fragments;

    private NetChop(Peptide peptide, int[] lengths, double threshold) {
        this.peptide = peptide;
        this.lengths = lengths;
        this.threshold = threshold;

        Arrays.sort(this.lengths);
        THRESHOLD_PROBABILITY_RANGE.validate("Threshold probability", threshold);
    }

    /**
     * Name of the environment variable that defines the absolute path
     * of the {@code netchop} executable file.  If the system property
     * {@code pepmhc.chop.netchop} is also defined, it will override
     * the environment variable.
     */
    public static final String EXECUTABLE_PATH_ENV = "NETCHOP_EXE";

    /**
     * Name of the system property that defines the absolute path of
     * the {@code netchop} executable file.
     */
    public static final String EXECUTABLE_PATH_PROPERTY = "pepmhc.chop.netchop";

    /**
     * Default value for the threshold probability for assigned
     * cleavage sites.
     */
    public static final double THRESHOLD_PROBABILITY_DEFAULT = 0.5;

    /**
     * Valid range of threshold probabilities.
     */
    public static final DoubleRange THRESHOLD_PROBABILITY_RANGE = DoubleRange.FRACTIONAL;

    /**
     * Default value for the lengths of the cleaved peptides.
     */
    public static final int[] PEPTIDE_LENGTHS_DEFAULT = new int[] { 9, 10 };

    /**
     * Simulates proteasomal processing of a peptide.
     *
     * <p>Predicts cleavage sites using the default threshold
     * probability and assembles cleaved fragments having the
     * prescribed lengths.
     *
     * @param peptide the peptide to chop.
     *
     * @param lengths the lengths of the cleaved peptides to
     * generate.
     *
     * @return a list containing the cleaved peptide fragments.
     */
    public static List<Peptide> chop(Peptide peptide, int[] lengths) {
        return chop(peptide, lengths, THRESHOLD_PROBABILITY_DEFAULT);
    }

    /**
     * Simulates proteasomal processing of a peptide.
     *
     * <p>Predicts cleavage sites using the assigned threshold
     * probability and assembles cleaved fragments having the
     * prescribed lengths.
     *
     * @param peptide the peptide to chop.
     *
     * @param lengths the lengths of the cleaved peptides to
     * generate.
     *
     * @param threshold the threshold probability for assigned
     * cleavage sites.
     *
     * @return a list containing the cleaved peptide fragments.
     */
    public static List<Peptide> chop(Peptide peptide, int[] lengths, double threshold) {
        NetChop chopper = new NetChop(peptide, lengths, threshold);
        return chopper.chop();
    }

    /**
     * Determines whether the {@code netchop} executable is installed.
     *
     * @return {@code true} iff the {@code netchop} executable is
     * installed at the location specified by the required system
     * property or environment variable.
     */
    public static boolean isInstalled() {
        try {
            return resolveExecutableFile().canExecute();
        }
        catch (Exception ex) {
            return false;
        }
    }

    /**
     * Resolves the absolute path to the {@code netchop} executable file.
     *
     * @return the path specified by the {@code EXECUTABLE_PATH_PROPERTY},
     * if set, or the {@code EXECUTABLE_PATH_ENV} environment variable.
     *
     * @throws RuntimeException unless a path is specified.
     */
    public static File resolveExecutableFile() {
        return new File(resolveExecutableName());
    }

    /**
     * Resolves the absolute path to the {@code netchop} executable file.
     *
     * @return the path specified by the {@code EXECUTABLE_PATH_PROPERTY},
     * if set, or the {@code EXECUTABLE_PATH_ENV} environment variable.
     *
     * @throws RuntimeException unless a path is specified.
     */
    public static String resolveExecutableName() {
        if (JamProperties.isSet(EXECUTABLE_PATH_PROPERTY))
            return JamProperties.getRequired(EXECUTABLE_PATH_PROPERTY);
        else
            return JamEnv.getRequired(EXECUTABLE_PATH_ENV);
    }

    private List<Peptide> chop() {
        computeScores();
        assembleFragments();

        return fragments;
    }

    private void computeScores() {
        scores = NetChopRunner.chop(peptide);
    }

    private void assembleFragments() {
        fragments = new ArrayList<Peptide>();

        for (int length : lengths)
            assembleFragments(length);
    }

    private void assembleFragments(int fragLen) {
        //
        // No cleavage fragments unless the peptide is as long as the
        // fragment length...
        //
        if (peptide.length() < fragLen)
            return;
        
        // A cleavage fragment of length "L" is generated with
        // C-terminus at index "L - 1" iff the residue at index
        // "L - 1" is a cleavage site.
        int cterm = fragLen - 1;

        if (isCleavageSite(cterm))
            fragments.add(cleavageFragment(cterm, fragLen));

        for (cterm = fragLen; cterm < peptide.length(); ++cterm) {
            //
            // A cleavage fragment of length "L" is generated with
            // C-terminus at index "i" iff:
            //
            // (1) The residue at index "i" is a cleavage site, and
            // (2) The residue at index "i - L" is a cleavage site.
            //
            if (isCleavageSite(cterm) && isCleavageSite(cterm - fragLen))
                fragments.add(cleavageFragment(cterm, fragLen));
        }
    }

    private boolean isCleavageSite(int index) {
        return scores.get(index) >= threshold;
    }

    private Peptide cleavageFragment(int cterm, int fragLen) {
        return peptide.fragment(new IntRange(cterm - fragLen + 1, cterm));
    }
}
