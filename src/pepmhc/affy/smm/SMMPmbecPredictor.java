
package pepmhc.affy.smm;

import pepmhc.affy.AffinityMethod;

public final class SMMPmbecPredictor extends MatrixPredictor {
    private SMMPmbecPredictor() {}

    /**
     * The single instance.
     */
    public static final SMMPmbecPredictor INSTANCE = new SMMPmbecPredictor();

    @Override public AffinityMethod getMethod() {
        return AffinityMethod.SMM_PMBEC;
    }
}
