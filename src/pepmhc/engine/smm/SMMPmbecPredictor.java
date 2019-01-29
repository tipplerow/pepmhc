
package pepmhc.engine.smm;

import pepmhc.engine.PredictionMethod;

public final class SMMPmbecPredictor extends MatrixPredictor {
    private SMMPmbecPredictor() {}

    /**
     * The single instance.
     */
    public static final SMMPmbecPredictor INSTANCE = new SMMPmbecPredictor();

    @Override public PredictionMethod getMethod() {
        return PredictionMethod.SMM_PMBEC;
    }
}
