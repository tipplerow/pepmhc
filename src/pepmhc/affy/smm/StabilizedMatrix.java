
package pepmhc.affy.smm;

import java.util.List;
import java.util.Map;

import jam.app.JamEnv;
import jam.hla.Allele;
import jam.io.FileUtil;
import jam.lang.JamException;
import jam.peptide.Peptide;
import jam.peptide.Residue;

import pepmhc.affy.AffinityMethod;

/**
 * Implements a <em>stabilized matrix</em> for predicting peptide-MHC
 * binding affinities.
 */
public final class StabilizedMatrix {
    private final double intercept;
    private final List<Map<Residue, Double>> elements;

    StabilizedMatrix(List<Map<Residue, Double>> elements, double intercept) {
        this.elements = elements;
        this.intercept = intercept;
    }

    /**
     * Returns the stabilized matrix for a given method, allele, and
     * peptide length.
     *
     * @param method the desired prediction method.
     *
     * @param allele the desired MHC allele.
     *
     * @param length the desired peptide length.
     *
     * @return the stabilized matrix for the specified method and
     * allele.
     */
    public static StabilizedMatrix instance(AffinityMethod method, Allele allele, int length) {
        return load(resolveFileName(method, allele, length));
    }

    /**
     * Reads a stabilized matrix from a formatted data file.
     *
     * @param fileName the name of the data file.
     *
     * @return the stabilized matrix with parameters contained in the
     * specified data file.
     */
    public static StabilizedMatrix load(String fileName) {
        return MatrixReader.load(fileName);
    }

    private static String resolveFileName(AffinityMethod method, Allele allele, int length) {
        return FileUtil.join(dirName(method), baseName(allele, length));
    }

    private static String dirName(AffinityMethod method) {
        String homeDir = JamEnv.getRequired("PEPMHC_HOME");

        switch (method) {
        case SMM:
            return FileUtil.join(homeDir, "data", "smm");

        case SMM_PMBEC:
            return FileUtil.join(homeDir, "data", "smm_pmbec");

        default:
            throw JamException.runtime("Unsupported prediction method: [%s].", method);
        }
    }

    private static String baseName(Allele allele, int length) {
        return String.format("%s-%d.txt", allele.longKey().replace('*', '-'), length);
    }

    /**
     * Computes the predicted binding affinity for a given peptide.
     *
     * @param peptide the peptide of interest.
     *
     * @return the predicted binding affinity for the given peptide.
     */
    public double computeIC50(Peptide peptide) {
        if (peptide.length() != elements.size())
            throw new IllegalArgumentException("Invalid peptide length.");

        double logsum = intercept;

        for (int index = 0; index < peptide.length(); ++index)
            logsum += elements.get(index).get(peptide.at(index));

        return Math.pow(10.0, logsum);
    }

    /**
     * Returns the matrix element (log10 contribution) for a given
     * residue and binding position.
     *
     * @param residue the peptide residue of interest.
     *
     * @param position the position where the residue binds.
     *
     * @return the matrix element (log10 contribution) for the given
     * residue and binding position.
     */
    public double getElement(Residue residue, int position) {
        return elements.get(position).get(residue);
    }

    /**
     * Returns the value of the intercept for this stabilized matrix.
     *
     * @return the value of the intercept for this stabilized matrix.
     */
    public double getIntercept() {
        return intercept;
    }

    /**
     * Returns the peptide length covered by this stabilized matrix.
     *
     * @return the peptide length covered by this stabilized matrix.
     */
    public int getPeptideLength() {
        return elements.size();
    }
}
