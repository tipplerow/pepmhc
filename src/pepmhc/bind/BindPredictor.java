
package pepmhc.bind;

import java.util.Collection;
import java.util.List;

import jene.hla.Allele;
import jene.peptide.Peptide;

/**
 * Defines an interface for computational engines or algorithms that
 * predict binding affinity or stability for peptide-MHC complexes.
 *
 * @param <R> the type of binding record (affinity or stability)
 * produced by the predictor.
 */
public abstract class BindPredictor<R extends BindRecord> {
    /**
     * Returns the enumerated prediction method.
     *
     * @return the enumerated prediction method.
     */
    public abstract Enum getMethod();

    /**
     * Determines whether the underlying engine for this predictor
     * is installed and available for the JVM.
     *
     * @return {@code true} iff the underlying engine for this
     * predictor is installed and available for the JVM.
     */
    public abstract boolean isInstalled();

    /**
     * Computes the binding record for a single peptide presented to
     * an HLA allele.
     *
     * @param allele the receiving allele.
     *
     * @param peptide the presented peptide.
     *
     * @return the binding record for the presented peptide.
     */
    public R predict(Allele allele, Peptide peptide) {
        return predict(allele, List.of(peptide)).get(0);
     }

    /**
     * Computes binding records for a collection of peptides presented
     * to an HLA allele.
     *
     * @param allele the receiving allele.
     *
     * @param peptides the presented peptides.
     *
     * @return the binding records indexed by the presented peptides.
     */
    public BindRecordMap map(Allele allele, Collection<Peptide> peptides) {
        return BindRecordMap.hash(predict(allele, peptides));
    }

    /**
     * Computes binding records for a collection of peptides presented
     * to an HLA allele.
     *
     * @param allele the receiving allele.
     *
     * @param peptides the presented peptides.
     *
     * @return the binding records for the presented peptides (in the
     * order returned by the collection iterator)
     */
    public abstract List<R> predict(Allele allele, Collection<Peptide> peptides);
}
