
package pepmhc.stab;

import jam.app.JamProperties;
import jam.hla.Allele;
import jam.util.PairKeyTable;

import pepmhc.bind.BindCache;

/**
 * Maintains an in-memory cache of stability records backed by a
 * persistent database store.
 */
public final class StabilityCache extends BindCache<StabilityRecord> {
    private static final PairKeyTable<StabilityMethod, Allele, StabilityCache> instances = PairKeyTable.hash();

    private StabilityCache(StabilityStore store) {
        super(store);
    }

    /**
     * Name of the environment variable that specifies the directory
     * containing the persistent database store. The system property
     * {@code pepmhc.stab.cacheDir} will take precedence if both are
     * specified.
     */
    public static final String CACHE_DIRECTORY_ENV = "PEPMHC_STABILITY_CACHE";

    /**
     * Name of the system property that specifies the directory
     * containing the persistent database store.
     */
    public static final String CACHE_DIRECTORY_PROPERTY = "pepmhc.stabilityCache";

    /**
     * Returns the name of the directory containing the persistent
     * stability database.
     *
     * @return the name of the directory containing the persistent
     * stability database.
     */
    public static String cacheDir() {
        return JamProperties.resolve(CACHE_DIRECTORY_PROPERTY, CACHE_DIRECTORY_ENV, null);
    }

    /**
     * Returns the stability cache for a given allele and prediction
     * method.
     *
     * @param method the stability prediction method.
     *
     * @param allele the allele of the binding MHC molecule.
     *
     * @return the stability cache for the specified allele and
     * prediction method.
     */
    public static synchronized StabilityCache instance(StabilityMethod method, Allele allele) {
        StabilityCache instance = instances.get(method, allele);

        if (instance == null)
            instance = newInstance(method, allele);

        return instance;
    }

    private static StabilityCache newInstance(StabilityMethod method, Allele allele) {
        StabilityCache instance =
            new StabilityCache(StabilityStore.instance(method, allele));

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

    @Override public StabilityMethod getMethod() {
        return (StabilityMethod) super.getMethod();
    }
}
