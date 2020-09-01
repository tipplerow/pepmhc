
package pepmhc.stab;

import jam.util.PairKeyTable;

import jene.hla.Allele;

import pepmhc.bind.BindStore;

/**
 * Provides a compute-on-demand service and persistent storage for
 * peptide-MHC stability records.
 */
public final class StabilityStore extends BindStore<StabilityRecord> {
    private static final PairKeyTable<StabilityMethod, Allele, StabilityStore> instances = PairKeyTable.hash();

    private StabilityStore(StabilityTable table, StabilityPredictor predictor, Allele allele) {
        super(table, predictor, allele);
    }

    /**
     * Returns the stability store for a given allele and prediction
     * method.
     *
     * @param method the stability prediction method.
     *
     * @param allele the allele of the binding MHC molecule.
     *
     * @return the stability store for the specified allele and
     * prediction method.
     */
    public static synchronized StabilityStore instance(StabilityMethod method, Allele allele) {
        StabilityStore instance = instances.get(method, allele);

        if (instance == null)
            instance = newInstance(method, allele);

        return instance;
    }

    private static StabilityStore newInstance(StabilityMethod method, Allele allele) {
        StabilityStore instance =
            new StabilityStore(StabilityTable.create(method, allele), method.getPredictor(), allele);

        instances.put(method, allele, instance);
        return instance;
    }

    @Override public StabilityMethod getMethod() {
        return getPredictor().getMethod();
    }

    @Override public StabilityPredictor getPredictor() {
        return (StabilityPredictor) super.predictor;
    }
}
