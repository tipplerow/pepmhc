
package pepmhc.engine;

import jam.app.JamProperties;

/**
 * Enumerates peptide-MHC binding affinity prediction methods.
 */
public enum PredictionMethod {
    NET_MHC, NET_MHC_PAN, NET_MHC_STAB_PAN, SMM, SMM_PMBEC;

    private static PredictionMethod global = null;

    /**
     * Name of the system property that specifies the global
     * prediction method.
     */
    public static final String METHOD_PROPERTY = "pepmhc.engine.predictionMethod";

    /**
     * Returns the global prediction method specified through system
     * properties.
     *
     * @return the global prediction method.
     */
    public static PredictionMethod global() {
        if (global == null)
            global = JamProperties.getRequiredEnum(METHOD_PROPERTY, PredictionMethod.class);

        return global;
    }
}
