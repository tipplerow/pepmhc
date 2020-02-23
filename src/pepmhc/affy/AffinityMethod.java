
package pepmhc.affy;

import jam.app.JamProperties;

import pepmhc.affy.net.NetMHCPredictor;
import pepmhc.affy.net.NetMHCPanPredictor;
import pepmhc.affy.smm.SMMPredictor;

/**
 * Enumerates peptide-MHC binding affinity prediction methods.
 */
public enum AffinityMethod {
    /**
     * Affinity prediction of half-life by {@code netMHC}.
     */
    NET_MHC {
        @Override public AffinityPredictor getPredictor() {
            return NetMHCPredictor.INSTANCE;
        }
    },

    /**
     * Affinity prediction of half-life by {@code netMHCpan}.
     */
    NET_MHC_PAN {
        @Override public AffinityPredictor getPredictor() {
            return NetMHCPanPredictor.INSTANCE;
        }
    },

    /**
     * Affinity prediction by the original stabilized matrix method.
     */
    SMM {
        @Override public AffinityPredictor getPredictor() {
            return SMMPredictor.INSTANCE;
        }
    },

    /**
     * Affinity prediction by the stabilized matrix method with a
     * peptide-MHC binding energy covariance (PMBEC) matrix prior.
     */
    SMM_PMBEC {
        @Override public AffinityPredictor getPredictor() {
            return SMMPredictor.INSTANCE;
        }
    };

    private static AffinityMethod global = null;

    /**
     * System property that specifies the global affinity prediction
     * method.
     */
    public static final String METHOD_PROPERTY = "pepmhc.affinityMethod";

    /**
     * Returns the global prediction method specified through system
     * properties.
     *
     * @return the global prediction method.
     */
    public static AffinityMethod global() {
        if (global == null)
            global = JamProperties.getRequiredEnum(METHOD_PROPERTY, AffinityMethod.class);

        return global;
    }

    /**
     * Returns the prediction engine for this method.
     *
     * @return the prediction engine for this method.
     */
    public abstract AffinityPredictor getPredictor();
}
