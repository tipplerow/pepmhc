
package pepmhc.bind;

import java.util.Collection;
import java.util.List;

import jam.hla.Allele;
import jam.peptide.Peptide;
import jam.sql.SQLStore;
import jam.sql.SQLTable;

/**
 * Provides a compute-on-demand service and persistent storage for
 * peptide-MHC binding records.
 *
 * @param <R> the type of binding record (affinity or stability)
 * produced by the predictor.
 */
public abstract class BindStore<R extends BindRecord> extends SQLStore<Peptide, R> {
    /**
     * The allele served by this store.
     */
    protected final Allele allele;

    /**
     * The predictor serving this store.
     */
    protected final BindPredictor<R> predictor;

    /**
     * Creates a new record store for a given allele and predictor.
     *
     * @param table the database table providing persistent storage.
     *
     * @param predictor the predictor serving this store.
     *
     * @param allele the allele served by this store.
     */
    protected BindStore(SQLTable<Peptide, R> table, BindPredictor<R> predictor, Allele allele) {
        super(table);

        this.allele = allele;
        this.predictor = predictor;
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

    @Override protected R compute(Peptide peptide) {
        return compute(List.of(peptide)).get(0);
    }

    @Override protected List<R> compute(Collection<Peptide> peptides) {
        return predictor.predict(allele, peptides);
    }
}
