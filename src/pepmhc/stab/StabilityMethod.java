
package pepmhc.stab;

import jam.app.JamProperties;

/**
 * Enumerates peptide-MHC stability prediction methods.
 */
public enum StabilityMethod {
    /**
     * Direct prediction of half-life by the {@code netMHCstabpan}
     * engine.
     */
    NET_MHC_STAB_PAN,

    /**
     * Inferred prediction of half-life by computing binding affinity
     * with {@code netMHCpan} and using a regression model to convert
     * from affinity to half-life.
     */
    NET_MHC_PAN_AFFINITY_PROXY;

    private static StabilityMethod global = null;

    /**
     * Name of the system property that specifies the global stability
     * prediction method.
     */
    public static final String METHOD_PROPERTY = "pepmhc.stab.method";

    /**
     * Returns the global prediction method specified through system
     * properties.
     *
     * @return the global prediction method.
     */
    public static StabilityMethod global() {
        if (global == null)
            global = JamProperties.getRequiredEnum(METHOD_PROPERTY, StabilityMethod.class);

        return global;
    }
}
