
package pepmhc.affy;

import jam.app.JamProperties;
import jam.hla.Allele;
import jam.util.PairKeyTable;

import pepmhc.bind.BindCache;

/**
 * Maintains an in-memory cache of affinity records backed by a
 * persistent database store.
 */
public final class AffinityCache extends BindCache<AffinityRecord> {
    private static final PairKeyTable<AffinityMethod, Allele, AffinityCache> instances = PairKeyTable.hash();

    private AffinityCache(AffinityTable table, AffinityPredictor predictor, Allele allele) {
        super(table, predictor, allele);
    }

    /**
     * Name of the environment variable that specifies the directory
     * containing the persistent database store. The system property
     * {@code pepmhc.stab.cacheDir} will take precedence if both are
     * specified.
     */
    public static final String CACHE_DIRECTORY_ENV = "PEPMHC_AFFINITY_CACHE";

    /**
     * Name of the system property that specifies the directory
     * containing the persistent database store.
     */
    public static final String CACHE_DIRECTORY_PROPERTY = "pepmhc.affinityCache";

    /**
     * Returns the name of the directory containing the persistent
     * affinity database.
     *
     * @return the name of the directory containing the persistent
     * affinity database.
     */
    public static String cacheDir() {
        return JamProperties.resolve(CACHE_DIRECTORY_PROPERTY, CACHE_DIRECTORY_ENV, null);
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
    public static synchronized AffinityCache instance(AffinityMethod method, Allele allele) {
        AffinityCache instance = instances.get(method, allele);

        if (instance == null)
            instance = newInstance(method, allele);

        return instance;
    }

    private static AffinityCache newInstance(AffinityMethod method, Allele allele) {
        AffinityCache instance =
            new AffinityCache(AffinityTable.create(method, allele), method.getPredictor(), allele);

        instances.put(method, allele, instance);
        return instance;
    }

    /**
     * Clears this cache: removes all cached records from memory and
     * removes this cache from the underlying instance map.
     */
    @Override public void clear() {
        super.clear();

        synchronized (instances) {
            instances.remove(getMethod(), getAllele());
        }
    }

    @Override public AffinityMethod getMethod() {
        return (AffinityMethod) super.getMethod();
    }
}
