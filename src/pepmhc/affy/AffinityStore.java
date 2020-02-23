
package pepmhc.affy;

import jam.hla.Allele;
import jam.util.PairKeyTable;

import pepmhc.bind.BindStore;

/**
 * Provides a compute-on-demand service and persistent storage for
 * peptide-MHC affinity records.
 */
public final class AffinityStore extends BindStore<AffinityRecord> {
    private static final PairKeyTable<AffinityMethod, Allele, AffinityStore> instances = PairKeyTable.hash();

    private AffinityStore(AffinityTable table, AffinityPredictor predictor, Allele allele) {
        super(table, predictor, allele);
    }

    /**
     * Returns the affinity store for a given allele and prediction
     * method.
     *
     * @param method the affinity prediction method.
     *
     * @param allele the allele of the binding MHC molecule.
     *
     * @return the affinity store for the specified allele and
     * prediction method.
     */
    public static synchronized AffinityStore instance(AffinityMethod method, Allele allele) {
        AffinityStore instance = instances.get(method, allele);

        if (instance == null)
            instance = newInstance(method, allele);

        return instance;
    }

    private static AffinityStore newInstance(AffinityMethod method, Allele allele) {
        AffinityStore instance =
            new AffinityStore(AffinityTable.create(method, allele), method.getPredictor(), allele);

        instances.put(method, allele, instance);
        return instance;
    }

    @Override public AffinityMethod getMethod() {
        return getPredictor().getMethod();
    }

    @Override public AffinityPredictor getPredictor() {
        return (AffinityPredictor) super.predictor;
    }
}
