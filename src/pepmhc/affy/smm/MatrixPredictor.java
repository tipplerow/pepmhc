
package pepmhc.affy.smm;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import jam.hla.Allele;
import jam.peptide.Peptide;
import jam.util.PairKeyTable;

import pepmhc.affy.AffinityPredictor;
import pepmhc.affy.AffinityRecord;

public abstract class MatrixPredictor extends AffinityPredictor {
    private final PairKeyTable<Allele, Integer, StabilizedMatrix> table;

    protected MatrixPredictor() {
        this.table = PairKeyTable.hash();
    }

    @Override public boolean isInstalled() {
        return true;
    }

    @Override public AffinityRecord predict(Allele allele, Peptide peptide) {
        StabilizedMatrix matrix = getMatrix(allele, peptide.length());
        return new AffinityRecord(peptide, matrix.computeIC50(peptide));
    }

    @Override public List<AffinityRecord> predict(Allele allele, Collection<Peptide> peptides) {
        List<AffinityRecord> records = new ArrayList<AffinityRecord>(peptides.size());

        for (Peptide peptide : peptides)
            records.add(predict(allele, peptide));

        return records;
    }

    private StabilizedMatrix getMatrix(Allele allele, Integer length) {
        StabilizedMatrix matrix = table.get(allele, length);

        if (matrix == null) {
            matrix = StabilizedMatrix.instance(getMethod(), allele, length);
            table.put(allele, length, matrix);
        }

        return matrix;
    }
}
