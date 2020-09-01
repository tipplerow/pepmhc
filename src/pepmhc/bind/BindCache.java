
package pepmhc.bind;

import java.util.Collection;
import java.util.List;

import jam.sql.SQLCache;
import jam.sql.SQLKeyTable;

import jene.hla.Allele;
import jene.peptide.Peptide;

/**
 * Provides a compute-on-demand service, in-memory caching, and
 * persistent storage for peptide-MHC binding records.
 *
 * @param <R> the type of binding record (affinity or stability)
 * to cache.
 */
public abstract class BindCache<R extends BindRecord> extends SQLCache<Peptide, R> {
    /**
     * The allele served by this cache.
     */
    protected final Allele allele;

    /**
     * The predictor serving this cache.
     */
    protected final BindPredictor<R> predictor;

    /**
     * Creates a new record cache for a given allele and predictor.
     *
     * @param table the database table providing persistent storage.
     *
     * @param predictor the predictor serving this cache.
     *
     * @param allele the allele served by this cache.
     */
    protected BindCache(SQLKeyTable<Peptide, R> table, BindPredictor<R> predictor, Allele allele) {
        super(table);

        this.allele = allele;
        this.predictor = predictor;
    }

    /**
     * Returns the HLA allele served by this cache.
     *
     * @return the HLA allele served by this cache.
     */
    public Allele getAllele() {
        return allele;
    }

    /**
     * Returns the prediction method used to compute new records.
     *
     * @return the predictor method used to compute new records.
     */
    public Enum getMethod() {
        return predictor.getMethod();
    }

    /**
     * Returns the predictor used to compute new records.
     *
     * @return the predictor used to compute new records.
     */
    public BindPredictor<R> getPredictor() {
        return predictor;
    }

    @Override public String getName() {
        return getAllele() + ":" + getMethod();
    }

    @Override protected R compute(Peptide peptide) {
        return compute(List.of(peptide)).get(0);
    }

    @Override protected List<R> compute(Collection<Peptide> peptides) {
        return predictor.predict(allele, peptides);
    }

    @Override public Class getKeyClass() {
        return Peptide.class;
    }
}
