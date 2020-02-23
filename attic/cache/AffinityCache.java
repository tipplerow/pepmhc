
package pepmhc.cache;

import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import jam.app.JamLogger;
import jam.hla.Allele;
import jam.peptide.Peptide;
import jam.util.MapUtil;
import jam.util.PairKeyTable;
import jam.util.SetUtil;

import pepmhc.binder.BindingRecord;
import pepmhc.engine.Predictor;
import pepmhc.engine.PredictionMethod;

/**
 * Maintains an in-memory cache of peptide-MHC binding records backed
 * by a persistent database store.
 */
public final class AffinityCache {
    private final AffinityStore store;
    private final Map<Peptide, BindingRecord> cache;

    private static final PairKeyTable<PredictionMethod, Allele, AffinityCache> instances = PairKeyTable.hash();

    private AffinityCache(AffinityStore store) {
        this.store = store;
        this.cache = new HashMap<Peptide, BindingRecord>();
    }

    /**
     * Returns the affinity cache for a given allele and prediction
     * method.
     *
     * @param method the affinity prediction method.
     *
     * @param allele the allele of the binding MHC molecule.
     *
     * @return the affinity cache for the specified allele and
     * prediction method.
     */
    public static AffinityCache instance(PredictionMethod method, Allele allele) {
        JamLogger.info("Requesting affinity cache: [%s, %s]...", method, allele);
        return instanceSync(method, allele);
    }
    
    private static synchronized AffinityCache instanceSync(PredictionMethod method, Allele allele) {
        AffinityCache instance = instances.get(method, allele);

        if (instance == null)
            instance = newInstance(method, allele);

        JamLogger.info("Returning affinity cache: [%s, %s].", method, allele);
        return instance;
    }

    private static AffinityCache newInstance(PredictionMethod method, Allele allele) {
        AffinityCache instance =
            new AffinityCache(AffinityStore.instance(method, allele));

        instances.put(method, allele, instance);
        return instance;
    }

    /**
     * Clears this cache: removes all cached records from memory and
     * removes this cache from the underlying instance map.
     */
    public void clear() {
        JamLogger.info("Clearing affinity cache [%s, %s]...", getMethod(), getAllele());
        cache.clear();

        synchronized (instances) {
            instances.remove(getMethod(), getAllele());
        }
    }

    /**
     * Returns the affinity record for a given peptide.
     *
     * @param peptide the peptide of interest.
     *
     * @return the affinity record for the specified peptide.
     */
    public BindingRecord get(Peptide peptide) {
        return get(List.of(peptide)).get(0);
    }

    /**
     * Returns the affinity records for a collection of peptides.
     *
     * @param peptides the peptides of interest.
     *
     * @return the affinity records for the specified peptides (in the
     * order returned by the collection iterator).
     */
    public List<BindingRecord> get(Collection<Peptide> peptides) {
        JamLogger.info("Requesting [%d] records from affinity cache [%s, %s]...",
                       peptides.size(), getMethod(), getAllele());

        return getSync(peptides);
    }

    private List<BindingRecord> getSync(Collection<Peptide> peptides) {
        //
        // Identify peptides from the input collection that are not
        // present in the cache...
        //
        Set<Peptide> missing = SetUtil.missing(cache.keySet(), peptides);

        if (!missing.isEmpty()) {
            //
            // Query (compute on-demand) affinity for each missing
            // peptide...
            //
            List<BindingRecord> records = store.get(missing);

            // Add the new records to the in-memory cache (the
            // persistent store is updated automatically)...
            updateCache(records);
        }

        JamLogger.info("Returning [%d] records from affinity cache [%s, %s]...",
                       peptides.size(), getMethod(), getAllele());

        // Pull the records from the cache in the order requested...
        return MapUtil.get(cache, peptides);
    }

    private void updateCache(List<BindingRecord> records) {
        JamLogger.info("Adding [%d] records to the affinity cache [%s, %s]...",
                       records.size(), getMethod(), getAllele());

        for (BindingRecord record : records)
            cache.put(record.getPeptide(), record);
    }

    /**
     * Returns the HLA allele served by this cache.
     *
     * @return the HLA allele served by this cache.
     */
    public Allele getAllele() {
        return store.getAllele();
    }

    /**
     * Returns the prediction method used by this cache.
     *
     * @return the prediction method used by this cache.
     */
    public PredictionMethod getMethod() {
        return store.getMethod();
    }
}
