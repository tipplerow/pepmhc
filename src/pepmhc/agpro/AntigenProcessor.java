
package pepmhc.agpro;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import jam.app.JamProperties;
import jam.app.PropertyList;
import jam.math.IntUtil;
import jam.peptide.Peptide;
import jam.util.RegexUtil;

import pepmhc.chop.NetChop;
import pepmhc.tap.TAP;

/**
 * Encapsulates the components of antigen processing prior to
 * presentation by MHC class I molecules: peptide cleavage by
 * proteasomes and TAP-assisted transport into the endoplasmic
 * reticulum (prior to MHC loading).
 */
public final class AntigenProcessor {
    private final int[] cleavageLength;

    private final boolean useNetChop;
    private final boolean useTAPConsensus;

    private final double netChopThreshold;
    private final double tapAlphaShrinkage;
    private final double tapScoreThreshold;

    private final TAP tap;
    private final NetChop netChop;

    private static AntigenProcessor global = null;

    private AntigenProcessor(int[]   cleavageLength,
                             boolean useNetChop,
                             boolean useTAPConsensus,
                             double  netChopThreshold,
                             double  tapAlphaShrinkage,
                             double  tapScoreThreshold) {
        this.cleavageLength = cleavageLength;

        this.useNetChop = useNetChop;
        this.useTAPConsensus = useTAPConsensus;

        this.netChopThreshold = netChopThreshold;
        this.tapAlphaShrinkage = tapAlphaShrinkage;
        this.tapScoreThreshold = tapScoreThreshold;

        if (useNetChop)
            this.netChop = new NetChop(cleavageLength, netChopThreshold);
        else
            this.netChop = null;

        if (useTAPConsensus)
            this.tap = new TAP(tapAlphaShrinkage, tapScoreThreshold);
        else
            this.tap = null;
    }

    /**
     * Name of the system property that defines the lengths of the
     * peptides produced by proteasomal peptide cleavage.
     */
    public static final String CLEAVAGE_LENGTH_PROPERTY = "pepmhc.agpro.cleavageLength";

    /**
     * Name of the system property that specifies whether to use the
     * {@code netchop} program to predict peptide cleavage products.
     */
    public static final String USE_NETCHOP_PROPERTY = "pepmhc.agpro.useNetChop";

    /**
     * Name of the system property that specifies whether to use the
     * TAP consensus matrix to compute binding affinity and to filter
     * cleavage products by affinity score.
     */
    public static final String USE_TAP_PROPERTY = "pepmhc.agpro.useTAP";

    /**
     * Name of the system property that defines the threshold cleavage
     * probability passed to the {@code netchop} predictor.  (Required
     * if using {@code netchop}.)
     */
    public static final String NETCHOP_THRESHOLD_PROPERTY = "pepmhc.agpro.netChopThreshold";

    /**
     * Name of the system property that defines the alpha shrinkage
     * factor passed to the TAP predictor.  (Required if using TAP).
     */
    public static final String TAP_ALPHA_PROPERTY = "pepmhc.agpro.tapAlpha";

    /**
     * Name of the system property that defines the score selection
     * threshold passed to the TAP predictor. (Required if using TAP).
     */
    public static final String TAP_SCORE_THRESHOLD_PROPERTY = "pepmhc.agpro.tapScoreThreshold";

    /**
     * Default value for the lengths of the peptides produced by
     * proteasomal peptide cleavage.
     */
    public static final int[] CLEAVAGE_LENGTH_DEFAULT = new int[] { 9, 10 };

    /**
     * Default value specifying whether to use the {@code netchop}
     * program to predict peptide cleavage products.
     */
    public static final boolean USE_NETCHOP_DEFAULT = true;

    /**
     * Default value specifying whether to use the TAP consensus
     * matrix to compute binding affinity and to filter cleavage
     * products by affinity score.
     */
    public static final boolean USE_TAP_DEFAULT = true;

    /**
     * Default value for the threshold cleavage probability passed to
     * the {@code netchop} predictor.
     */
    public static final double NETCHOP_THRESHOLD_DEFAULT = 0.5;

    /**
     * Default value for the alpha shrinkage factor passed to the TAP
     * predictor.
     */
    public static final double TAP_ALPHA_DEFAULT = 0.2;

    /**
     * Default value for the score selection threshold passed to the
     * TAP predictor.
     */
    public static final double TAP_SCORE_THRESHOLD_DEFAULT = 1.0;

    /**
     * Returns a new antigen processor with the default configuration.
     *
     * @return a new antigen processor with the default configuration.
     */
    public static AntigenProcessor defaultProcessor() {
        return new AntigenProcessor(CLEAVAGE_LENGTH_DEFAULT,
                                    USE_NETCHOP_DEFAULT,
                                    USE_TAP_DEFAULT,
                                    NETCHOP_THRESHOLD_DEFAULT,
                                    TAP_ALPHA_DEFAULT,
                                    TAP_SCORE_THRESHOLD_DEFAULT);
    }

    /**
     * Creates a new <em>exhaustive</em> antigen processor: all
     * possible native peptides are generated with no filtering
     * for cleavage site prediction or TAP transport efficiency.
     *
     * @param cleavageLength the lengths of the peptides produced
     * by proteasomal cleavage.
     *
     * @return an exhaustive antigen processor with the specified
     * cleavage lengths.
     */
    public static AntigenProcessor exhaustive(int... cleavageLength) {
        return new AntigenProcessor(cleavageLength, false, false,
                                    Double.NaN, Double.NaN, Double.NaN);
    }

    /**
     * Returns the global antigen processor defined through system
     * properties.
     *
     * @return the global antigen processor defined through system
     * properties.
     *
     * @throws RuntimeException unless the system properties define
     * a valid antigen processor.
     */
    public static AntigenProcessor global() {
        if (global == null)
            global = createGlobal();

        return global;
    }

    private static AntigenProcessor createGlobal() {
        int[] cleavageLength = resolveCleavageLength();

        boolean useNetChop = resolveUseNetChop();
        boolean useTAPConsensus = resolveUseTAPConsensus();

        double netChopThreshold = useNetChop ? resolveNetChopThreshold() : Double.NaN;

        double tapAlphaShrinkage = useTAPConsensus ? resolveTAPAlphaShrinkage() : Double.NaN;
        double tapScoreThreshold = useTAPConsensus ? resolveTAPScoreThreshold() : Double.NaN;

        return new AntigenProcessor(cleavageLength, useNetChop, useTAPConsensus,
                                    netChopThreshold, tapAlphaShrinkage, tapScoreThreshold);
    }

    private static int[] resolveCleavageLength() {
        return IntUtil.parseIntArray(JamProperties.getRequired(CLEAVAGE_LENGTH_PROPERTY), RegexUtil.COMMA);
    }

    private static boolean resolveUseNetChop() {
        return JamProperties.getRequiredBoolean(USE_NETCHOP_PROPERTY);
    }

    private static boolean resolveUseTAPConsensus() {
        return JamProperties.getRequiredBoolean(USE_TAP_PROPERTY);
    }

    private static double resolveNetChopThreshold() {
        return JamProperties.getRequiredDouble(NETCHOP_THRESHOLD_PROPERTY, NetChop.THRESHOLD_PROBABILITY_RANGE);
    }

    private static double resolveTAPAlphaShrinkage() {
        return JamProperties.getRequiredDouble(TAP_ALPHA_PROPERTY, TAP.ALPHA_RANGE);
    }

    private static double resolveTAPScoreThreshold() {
        return JamProperties.getRequiredDouble(TAP_SCORE_THRESHOLD_PROPERTY, TAP.THRESHOLD_SCORE_RANGE);
    }

    /**
     * Creates a new antigen processor by resolving the properties
     * defined in a file.
     *
     * @param file the property file containing the antigen processor
     * definitition.
     *
     * @return a new antigen processor defined by the properties in
     * the specified file.
     *
     * @throws RuntimeException unless the file can be opened for
     * reading and contains properties that define a valid antigen
     * processor.
     */
    public static AntigenProcessor resolve(File file) {
        PropertyList properties = JamProperties.parseFile(file);

        int[] cleavageLength = IntUtil.parseIntArray(properties.require(CLEAVAGE_LENGTH_PROPERTY), RegexUtil.COMMA);

        boolean useNetChop      = Boolean.valueOf(properties.require(USE_NETCHOP_PROPERTY));
        boolean useTAPConsensus = Boolean.valueOf(properties.require(USE_TAP_PROPERTY));

        double netChopThreshold  = Double.NaN;
        double tapAlphaShrinkage = Double.NaN;
        double tapScoreThreshold = Double.NaN;

        if (useNetChop)
            netChopThreshold = Double.valueOf(properties.require(NETCHOP_THRESHOLD_PROPERTY));

        if (useTAPConsensus) {
            tapAlphaShrinkage = Double.valueOf(properties.require(TAP_ALPHA_PROPERTY));
            tapScoreThreshold = Double.valueOf(properties.require(TAP_SCORE_THRESHOLD_PROPERTY));
        }

        return new AntigenProcessor(cleavageLength, useNetChop, useTAPConsensus,
                                    netChopThreshold, tapAlphaShrinkage, tapScoreThreshold);
    }

    /**
     * Creates a new antigen processor by resolving the properties
     * defined in a file.
     *
     * @param fileName the name of the property file containing the
     * antigen processor definitition.
     *
     * @return a new antigen processor defined by the properties in
     * the specified file.
     *
     * @throws RuntimeException unless the file can be opened for
     * reading and contains properties that define a valid antigen
     * processor.
     */
    public static AntigenProcessor resolve(String fileName) {
        return resolve(new File(fileName));
    }

    /**
     * Returns the active {@code netchop} predictor.
     *
     * @return the active {@code netchop} predictor.
     */
    public NetChop getNetChop() {
        return netChop;
    }

    /**
     * Returns the active TAP predictor.
     *
     * @return the active TAP predictor.
     */
    public TAP getTAP() {
        return tap;
    }

    /**
     * Simulates the antigen processing for a single protein.
     *
     * @param protein the protein to process.
     *
     * @return the peptide fragments to be presented to MHC molecules.
     */
    public List<Peptide> process(Peptide protein) {
        return transport(cleave(protein));
    }

    private List<Peptide> cleave(Peptide protein) {
        if (useNetChop)
            return netChop.chop(protein);
        else
            return cleaveNative(protein);
    }

    private List<Peptide> cleaveNative(Peptide protein) {
        List<Peptide> fragments = new ArrayList<Peptide>();

        for (int L : cleavageLength)
            fragments.addAll(protein.nativeFragments(L));

        return fragments;
    }

    private List<Peptide> transport(List<Peptide> fragments) {
        if (useTAPConsensus)
            return tap.transport(fragments);
        else
            return fragments;
    }
}
