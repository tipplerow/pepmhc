
package pepmhc.tap;

import java.io.File;
import java.util.List;

import jam.app.JamEnv;
import jam.data.DenseDataMatrix;
import jam.data.DenseDataMatrixLoader;
import jam.io.FileUtil;
import jam.matrix.JamMatrix;

import jene.peptide.Residue;

/**
 * Stores the scoring matrices used in the TAP transport predictions.
 */
public final class TAPMatrix extends DenseDataMatrix<Residue, TAPPosition> {
    private TAPMatrix(List<Residue> residues, List<TAPPosition> positions, JamMatrix scores) {
        super(residues, positions, scores, false, false);
    }

    private static TAPMatrix consensus = null;

    /**
     * Returns the consensus scoring matrix from Peters et al.,
     * J. Immunol. 171, 1741--1749 (2003).
     *
     * @return the consensus scoring matrix.
     */
    public static TAPMatrix consensus() {
        if (consensus == null)
            consensus = load(resolveConsensusFile());

        return consensus;
    }

    private static String resolveConsensusFile() {
        return FileUtil.join(JamEnv.getRequired("PEPMHC_HOME"), "data", "tap", "consensus.tsv");
    }

    /**
     * Loads a TAP scoring matrix from a data file.
     *
     * @param file the file to load.
     *
     * @return the TAP scoring matrix contained in the specified
     * file.
     *
     * @throws RuntimeException unless the file can be opened for
     * reading and contains a valid scoring matrix.
     */
    public static TAPMatrix load(File file) {
        TAPMatrixLoader loader = new TAPMatrixLoader(file);
        return loader.load();
    }

    /**
     * Loads a TAP scoring matrix from a data file.
     *
     * @param fileName the name of the file to load.
     *
     * @return the TAP scoring matrix contained in the specified
     * file.
     *
     * @throws RuntimeException unless the file can be opened for
     * reading and contains a valid scoring matrix.
     */
    public static TAPMatrix load(String fileName) {
        return load(new File(fileName));
    }

    private static final class TAPMatrixLoader extends DenseDataMatrixLoader<Residue, TAPPosition> {
        TAPMatrixLoader(File file) {
            super(file);
        }

        @Override public TAPMatrix load() {
            return (TAPMatrix) super.load();
        }

        @Override public TAPMatrix newMatrix(List<Residue> residues, List<TAPPosition> positions, JamMatrix scores) {
            return new TAPMatrix(residues, positions, scores);
        }

        @Override public TAPPosition parseColKey(String colKey) {
            return TAPPosition.valueOf(colKey);
        }

        @Override public Residue parseRowKey(String rowKey) {
            return Residue.valueOfCode1(rowKey);
        }
    }
}
