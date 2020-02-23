
package pepmhc.stab;

import jam.app.JamProperties;

import pepmhc.stab.net.NetStabPredictor;
import pepmhc.stab.proxy.AffinityProxyPredictor;

/**
 * Enumerates peptide-MHC stability prediction methods.
 */
public enum StabilityMethod {
    /**
     * Inferred prediction of half-life by computing binding affinity
     * with {@code netMHCpan} and using a regression model to convert
     * from affinity to half-life.
     */
    AFFINITY_PROXY {
        @Override public StabilityPredictor getPredictor() {
            return AffinityProxyPredictor.INSTANCE;
        }
    },

    /**
     * Direct prediction of half-life by {@code netMHCstabpan}.
     */
    NET_MHC_STAB_PAN {
        @Override public StabilityPredictor getPredictor() {
            return NetStabPredictor.INSTANCE;
        }
    };

    private static StabilityMethod global = null;

    /**
     * System property that specifies the global stability prediction
     * method.
     */
    public static final String METHOD_PROPERTY = "pepmhc.stabilityMethod";

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

    /**
     * Returns the prediction engine for this method.
     *
     * @return the prediction engine for this method.
     */
    public abstract StabilityPredictor getPredictor();
}
