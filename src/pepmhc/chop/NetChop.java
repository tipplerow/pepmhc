
package pepmhc.chop;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import jam.app.JamEnv;
import jam.app.JamProperties;
import jam.math.DoubleRange;
import jam.math.IntRange;
import jam.peptide.Peptide;

/**
 * Simulates proteasomal processing of a peptide:  Predicts cleavage
 * sites using {@code netchop} and assembles cleaved fragments having
 * a desired length.
 */
public final class NetChop {
    private final Peptide peptide;
    private final int[]   lengths;
    private final double  threshold;

    private List<Double> scores;
    private List<Peptide> fragments;

    private NetChop(Peptide peptide, int[] lengths) {
        this.peptide   = peptide;
        this.lengths   = lengths;
        this.threshold = resolveThreshold();
    }

    /**
     * Name of the environment variable that defines the absolute path
     * of the {@code netchop} executable file.  If the system property
     * {@code pepmhc.chop.netchop} is also defined, it will take
     * precedence.
     */
    public static final String EXECUTABLE_PATH_ENV = "NETCHOP_EXE";

    /**
     * Name of the system property that defines the absolute path of
     * the {@code netchop} executable file.
     */
    public static final String EXECUTABLE_PATH_PROPERTY = "pepmhc.chop.netchop";

    /**
     * Name of the system property that defines the threshold
     * probability for assigned cleavage sites.
     */
    public static final String THRESHOLD_PROBABILITY_PROPERTY = "pepmhc.chop.thresholdProbability";

    /**
     * Default value for the threshold probability for assigned
     * cleavage sites.
     */
    public static final double THRESHOLD_PROBABILITY_DEFAULT = 0.5;

    private static double resolveThreshold() {
        return JamProperties.getOptionalDouble(THRESHOLD_PROBABILITY_PROPERTY,
                                               DoubleRange.FRACTIONAL,
                                               THRESHOLD_PROBABILITY_DEFAULT);
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

    /**
     * Simulates proteasomal processing of a peptide: Predicts
     * cleavage sites and assembles cleaved fragments having a
     * desired length.
     *
     * @param peptide the peptide to chop.
     *
     * @param lengths the fragment lengths.
     *
     * @return a set containing the cleaved peptide fragments.
     */
    public static List<Peptide> chop(Peptide peptide, int... lengths) {
        NetChop chopper = new NetChop(peptide, lengths);
        return chopper.chop();
    }

    private List<Peptide> chop() {
        scores = NetChopRunner.chop(peptide);
        fragments = new ArrayList<Peptide>();

        for (int cterm = 0; cterm < peptide.length(); ++cterm) {
            //
            // A cleavage fragment of length "N" is generated with
            // C-terminus at index "i" iff:
            //
            // (1) The residue at index "i" is a cleavage site, and
            // (2) The residue at index "i - N" is a cleavage site.
            //
            if (isCleavageSite(cterm))
                for (int length : lengths)
                    if (isCleavageSite(cterm - length))
                        fragments.add(cleavageFragment(cterm, length));
        }

        return fragments;
    }

    private boolean isCleavageSite(int index) {
        return (index == -1)
            || (index == peptide.length() - 1)
            || (index >= 0 && scores.get(index) >= threshold);
    }

    private Peptide cleavageFragment(int cterm, int length) {
        return peptide.fragment(new IntRange(cterm - length + 1, cterm));
    }
}
