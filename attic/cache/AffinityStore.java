
package pepmhc.cache;

import java.util.Collection;
import java.util.List;

import jam.app.JamLogger;
import jam.hla.Allele;
import jam.peptide.Peptide;
import jam.sql.SQLStore;
import jam.util.PairKeyTable;

import pepmhc.binder.BindingRecord;
import pepmhc.engine.Predictor;
import pepmhc.engine.PredictionMethod;

/**
 * Provides a compute-on-demand service and persistent storage for
 * peptide-MHC binding affinity records.
 */
public final class AffinityStore extends SQLStore<Peptide, BindingRecord> {
    private final Allele allele;
    private final Predictor predictor;
    private final PredictionMethod method;

    private static final PairKeyTable<PredictionMethod, Allele, AffinityStore> instances = PairKeyTable.hash();

    private AffinityStore(AffinityTable table, PredictionMethod method, Allele allele) {
        super(table);

        this.method = method;
        this.allele = allele;
        this.predictor = Predictor.instance(method);
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
    public static AffinityStore instance(PredictionMethod method, Allele allele) {
        JamLogger.info("Requesting affinity store: [%s, %s]...", method, allele);
        return instanceSync(method, allele);
    }

    private static synchronized AffinityStore instanceSync(PredictionMethod method, Allele allele) {
        AffinityStore instance = instances.get(method, allele);

        if (instance == null)
            instance = newInstance(method, allele);

        JamLogger.info("Returning affinity store: [%s, %s].", method, allele);
        return instance;
    }

    private static AffinityStore newInstance(PredictionMethod method, Allele allele) {
        AffinityStore instance =
            new AffinityStore(AffinityTable.create(method, allele), method, allele);

        instances.put(method, allele, instance);
        return instance;
    }

    /**
     * Returns the HLA allele served by this store.
     *
     * @return the HLA allele served by this store.
     */
    public Allele getAllele() {
        return allele;
    }

    /**
     * Returns the prediction method used by this store.
     *
     * @return the prediction method used by this store.
     */
    public PredictionMethod getMethod() {
        return method;
    }

    @Override protected BindingRecord compute(Peptide peptide) {
        return compute(List.of(peptide)).get(0);
    }

    @Override protected List<BindingRecord> compute(Collection<Peptide> peptides) {
        return predictor.predict(allele, peptides);
    }
}
