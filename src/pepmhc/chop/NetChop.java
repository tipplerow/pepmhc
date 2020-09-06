
package pepmhc.chop;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import jam.app.JamEnv;
import jam.app.JamProperties;
import jam.math.Probability;
import jam.math.UnitIndex;
import jam.math.UnitIndexRange;
import jam.util.RegexUtil;

import jene.peptide.Peptide;

/**
 * Simulates proteasomal processing of a peptide: Predicts cleavage
 * sites using {@code netchop} and assembles cleaved fragments with
 * prescribed lengths.
 */
public final class NetChop {
    private final int[] lengths;
    private final Peptide peptide;
    private final Probability threshold;

    private List<Peptide> fragments;
    private List<Probability> scores;

    private NetChop(Peptide peptide, int[] lengths, Probability threshold) {
        this.peptide = peptide;
        this.lengths = lengths;
        this.threshold = threshold;

        Arrays.sort(this.lengths);
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
    public static final Probability THRESHOLD_PROBABILITY_DEFAULT = Probability.ONE_HALF;

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
    public static List<Peptide> chop(Peptide peptide, int[] lengths, Probability threshold) {
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
        scores = NetChopRunner.score(peptide);
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

        // Special logic for the first fragment...
        addFirstFragment(fragLen);

        // For all other fragments, beginning with N-terminus residue
        // at position 2 in the main peptide:
        UnitIndex nterm = UnitIndex.instance(2);
        UnitIndex cterm = UnitIndex.instance(fragLen + 1);
        UnitIndex lastC = UnitIndex.instance(peptide.length());

        while (cterm.LE(lastC)) {
            //
            // The peptide with N-terminus at position "nterm" and
            // C-terminus at position "cterm" is a cleavage fragment
            // iff the residue at position "cterm" is a cleavage site
            // and the residue at the position to the left of "nterm"
            // is a cleavage site...
            //
            if (isCleavageSite(cterm) && isCleavageSite(nterm.prev()))
                fragments.add(cleavageFragment(nterm, cterm));

            nterm = nterm.next();
            cterm = cterm.next();
        }
    }

    private void addFirstFragment(int fragLen) {
        //
        // The peptide with N-terminus at position 1 and C-terminus
        // at position "L" is a cleavage fragment iff the residue at
        // position "L" is a cleavage site...
        //
        UnitIndex nterm = UnitIndex.instance(1);
        UnitIndex cterm = UnitIndex.instance(fragLen);

        if (isCleavageSite(cterm))
            fragments.add(cleavageFragment(nterm, cterm));
    }

    private boolean isCleavageSite(UnitIndex position) {
        return position.get(scores).GE(threshold);
    }

    private Peptide cleavageFragment(UnitIndex nterm, UnitIndex cterm) {
        return peptide.fragment(UnitIndexRange.instance(nterm, cterm));
    }
}
