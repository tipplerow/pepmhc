
package pepmhc.neo;

import java.util.Collection;

import jam.hla.PeptideSource;
import jam.hugo.HugoPeptideTable;
import jam.hugo.HugoSymbol;
import jam.lang.JamException;
import jam.peptide.Peptide;

/**
 * Stores neo-peptides generated by antigen processing of mutated
 * genes along with self-peptides generated from the same genes.
 */
public final class PeptideSourceTable {
    private final HugoPeptideTable neoPeptideTable;
    private final HugoPeptideTable selfPeptideTable;

    private PeptideSourceTable(HugoPeptideTable neoPeptideTable,
                               HugoPeptideTable selfPeptideTable) {
        this.neoPeptideTable = neoPeptideTable;
        this.selfPeptideTable = selfPeptideTable;
    }

    /**
     * The single empty peptide-source table.
     */
    public static PeptideSourceTable EMPTY =
        create(HugoPeptideTable.EMPTY, HugoPeptideTable.EMPTY);

    /**
     * Creates a new peptide-source table from its constituent data
     * tables.
     *
     * @param neoPeptideTable the neo-peptides generated by antigen
     * processing of mutated genes.
     *
     * @param selfPeptideTable the self-peptides generated by antigen
     * processing of mutated genes.
     *
     * @return the new peptide source table.
     */
    public static PeptideSourceTable create(HugoPeptideTable neoPeptideTable,
                                            HugoPeptideTable selfPeptideTable) {
        return new PeptideSourceTable(neoPeptideTable, selfPeptideTable);
    }

    /**
     * Returns the peptides of a given type generated by antigen
     * processing of mutated genes.
     *
     * @param source the peptide type: "neo" or "self".
     *
     * @return the peptides of the specified type.
     */
    public HugoPeptideTable getPeptideTable(PeptideSource source) {
        switch (source) {
        case NEO:
            return neoPeptideTable;

        case SELF:
            return selfPeptideTable;

        default:
            throw JamException.runtime("Unknown peptide source: [%s].", source);
        }
    }

    /**
     * Returns the neo-peptides generated by antigen processing of
     * mutated genes.
     *
     * @return the neo-peptides generated by antigen processing of
     * mutated genes.
     */
    public HugoPeptideTable getNeoPeptideTable() {
        return neoPeptideTable;
    }

    /**
     * Returns the neo-peptides generated by antigen processing for
     * a given gene.
     *
     * @param hugoSymbol the HUGO symbol for the gene of interest.
     *
     * @return all neo-peptides generated from antigen processing of
     * the specified gene (will be empty for unmutated genes).
     */
    public Collection<Peptide> getNeoPeptides(HugoSymbol hugoSymbol) {
        return neoPeptideTable.get(hugoSymbol);
    }

    /**
     * Returns the self-peptides generated by antigen processing of
     * mutated genes.
     *
     * @return the self-peptides generated by antigen processing of
     * mutated genes.
     */
    public HugoPeptideTable getSelfPeptideTable() {
        return selfPeptideTable;
    }

    /**
     * Returns the self-peptides generated by antigen processing for
     * a given gene.
     *
     * @param hugoSymbol the HUGO symbol for the gene of interest.
     *
     * @param selfReference the reference self-peptidome to be used
     * for unmutated genes.
     *
     * @return all self-peptides generated from antigen processing of
     * the specified gene (the reference self-peptidome for unmutated
     * genes).
     */
    public Collection<Peptide> getSelfPeptides(HugoSymbol hugoSymbol, HugoPeptideTable selfReference) {
        if (selfPeptideTable.contains(hugoSymbol))
            return selfPeptideTable.get(hugoSymbol);
        else
            return selfReference.get(hugoSymbol);
    }
}