
package pepmhc.affy.smm;

import pepmhc.affy.AffinityMethod;

public final class SMMPredictor extends MatrixPredictor {
    private SMMPredictor() {}

    /**
     * The single instance.
     */
    public static final SMMPredictor INSTANCE = new SMMPredictor();

    @Override public AffinityMethod getMethod() {
        return AffinityMethod.SMM;
    }
}
