
package pepmhc.engine.smm;

import pepmhc.engine.PredictionMethod;

public final class SMMPredictor extends MatrixPredictor {
    private SMMPredictor() {}

    /**
     * The single instance.
     */
    public static final SMMPredictor INSTANCE = new SMMPredictor();

    @Override public PredictionMethod getMethod() {
        return PredictionMethod.SMM;
    }
}
