
package pepmhc.engine.smm;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jam.hla.Allele;
import jam.peptide.Peptide;

import pepmhc.binder.BindingRecord;
import pepmhc.engine.Predictor;

public abstract class MatrixPredictor implements Predictor {
    private final Map<Allele, Map<Integer, StabilizedMatrix>> outerMap;

    protected MatrixPredictor() {
        this.outerMap = new HashMap<Allele, Map<Integer, StabilizedMatrix>>();
    }

    @Override public boolean isInstalled() {
        return true;
    }

    @Override public BindingRecord predict(Allele allele, Peptide peptide) {
        return new BindingRecord(peptide, getMatrix(allele, peptide.length()).computeIC50(peptide));
    }

    @Override public List<BindingRecord> predict(Allele allele, Collection<Peptide> peptides) {
        List<BindingRecord> records = new ArrayList<BindingRecord>(peptides.size());

        for (Peptide peptide : peptides)
            records.add(predict(allele, peptide));

        return records;
    }

    private StabilizedMatrix getMatrix(Allele allele, int length) {
        StabilizedMatrix matrix = getInnerMap(allele).get(length);

        if (matrix == null) {
            matrix = StabilizedMatrix.instance(getMethod(), allele, length);
            getInnerMap(allele).put(length, matrix);
        }

        return matrix;
    }

    private Map<Integer, StabilizedMatrix> getInnerMap(Allele allele) {
        Map<Integer, StabilizedMatrix> innerMap = outerMap.get(allele);

        if (innerMap == null) {
            innerMap = new HashMap<Integer, StabilizedMatrix>();
            outerMap.put(allele, innerMap);
        }

        return innerMap;
    }
}
