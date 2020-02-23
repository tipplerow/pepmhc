
package pepmhc.bind;

import jam.hla.Allele;
import jam.peptide.Peptide;
import jam.sql.SQLCache;

/**
 * Maintains an in-memory cache of peptide-MHC binding records backed
 * by a persistent database store.
 */
public abstract class BindCache<R extends BindRecord> extends SQLCache<Peptide, R> {
    /**
     * Creates a new cache with a backing database store to provide
     * persistent storage and to compute records on demand.
     *
     * @param store the backing database store.
     */
    protected BindCache(BindStore<R> store) {
        super(store);
    }

    protected BindStore<R> getStore() {
        return (BindStore<R>) store;
    }

    /**
     * Returns the HLA allele served by this cache.
     *
     * @return the HLA allele served by this cache.
     */
    public Allele getAllele() {
        return getStore().getAllele();
    }

    /**
     * Returns the prediction method used to compute new records.
     *
     * @return the prediction method used to compute new records.
     */
    public Enum getMethod() {
        return getStore().getPredictor().getMethod();
    }

    private String methodName() {
        return getMethod().name();
    }

    /**
     * Returns the predictor used to compute new records.
     *
     * @return the predictor used to compute new records.
     */
    public BindPredictor<R> getPredictor() {
        return getStore().getPredictor();
    }

    @Override public String getName() {
        return getAllele().longKey();
    }
}
