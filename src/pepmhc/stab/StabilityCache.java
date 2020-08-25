
package pepmhc.stab;

import jam.app.JamProperties;
import jam.util.PairKeyTable;

import jean.hla.Allele;

import pepmhc.bind.BindCache;

/**
 * Provides a compute-on-demand service, in-memory caching, and
 * persistent storage for peptide-MHC stability records.
 */
public final class StabilityCache extends BindCache<StabilityRecord> {
    private static final PairKeyTable<StabilityMethod, Allele, StabilityCache> instances = PairKeyTable.hash();

    private StabilityCache(StabilityTable table, StabilityPredictor predictor, Allele allele) {
        super(table, predictor, allele);
    }

    /**
     * Name of the environment variable that specifies the directory
     * containing the persistent database store. The system property
     * {@code pepmhc.stabilityCache} will take precedence if both are
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
     * Clears any active caches for a given HLA allele (e.g., after
     * the allele has been processed).
     *
     * @param allele the allele that has been processed.
     */
    public static void clear(Allele allele) {
        for (StabilityMethod method : StabilityMethod.values())
            clear(method, allele);
    }

    /**
     * Removes a stability cache from memory (but retains all
     * persistent stability records).  method.
     *
     * @param method the stability prediction method.
     *
     * @param allele the allele of the binding MHC molecule.
     */
    public static synchronized void clear(StabilityMethod method, Allele allele) {
        StabilityCache instance = instances.get(method, allele);

        if (instance != null)
            instance.clear();

        instances.remove(method, allele);
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
            new StabilityCache(StabilityTable.create(method, allele), method.getPredictor(), allele);

        instances.put(method, allele, instance);
        return instance;
    }

    @Override public StabilityMethod getMethod() {
        return getPredictor().getMethod();
    }

    @Override public StabilityPredictor getPredictor() {
        return (StabilityPredictor) super.predictor;
    }

    @Override public Class getRecordClass() {
        return StabilityRecord.class;
    }
}
