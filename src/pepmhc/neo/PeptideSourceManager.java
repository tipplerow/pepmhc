
package pepmhc.neo;

import java.io.File;

import jam.hugo.HugoPeptideTable;
import jam.tcga.TumorBarcode;

/**
 * Manages...
 */
public final class PeptideSourceManager {
    private final String dirName;

    // The reference self-peptidome (for genes that are not mutated)...
    private final HugoPeptideTable selfReference;

    private PeptideSourceManager(String dirName, HugoPeptideTable selfReference) {
        this.dirName = dirName;
        this.selfReference = selfReference;
    }

    private static final String NEO_SUFFIX = "_neo_peptide.tsv.gz";
    private static final String SELF_SUFFIX = "_self_peptide.tsv.gz";

    /**
     * Creates a new peptide-source data manager.
     *
     * @param dirName the directory containing individual tumor
     * data files.
     *
     * @param selfRefName the name of the file containing the
     * reference self-peptidome (for genes that are not mutated).
     *
     * @return a new peptide-source data manager for the specified
     * directory.
     */
    public static PeptideSourceManager create(String dirName, String selfRefName) {
        return create(dirName, HugoPeptideTable.load(selfRefName));
    }

    /**
     * Creates a new peptide-source data manager.
     *
     * @param dirName the directory containing individual tumor
     * data files.
     *
     * @param selfRefTable a table containing the reference
     * self-peptidome (for genes that are not mutated).
     *
     * @return a new peptide-source data manager for the specified
     * directory.
     */
    public static PeptideSourceManager create(String dirName, HugoPeptideTable selfRefTable) {
        return new PeptideSourceManager(dirName, selfRefTable);
    }

    /**
     * Determines whether the peptide-source data for a given tumor
     * sample exists.
     *
     * @param barcode the barcode of the tumor sample.
     *
     * @return {@code true} iff the peptide-source data for the
     * specified tumor sample exists in the data directory.
     */
    public boolean exists(TumorBarcode barcode) {
        return neoPeptideFile(barcode).exists();
    }

    /**
     * Loads the peptide-source table for a given tumor sample.
     *
     * @param barcode the barcode of the desired tumor sample.
     *
     * @return the peptide-source table for the specified tumor sample
     * (an empty table if the data is not present).
     */
    public PeptideSourceTable load(TumorBarcode barcode) {
        File neoFile = neoPeptideFile(barcode);
        File selfFile = selfPeptideFile(barcode);

        if (neoFile.canRead() && selfFile.canRead())
            return load(neoFile, selfFile);
        else
            return PeptideSourceTable.EMPTY;
    }

    private File neoPeptideFile(TumorBarcode barcode) {
        return new File(dirName, neoBaseName(barcode));
    }

    private static String neoBaseName(TumorBarcode barcode) {
        return barcode.getKey() + NEO_SUFFIX;
    }

    private File selfPeptideFile(TumorBarcode barcode) {
        return new File(dirName, selfBaseName(barcode));
    }

    private static String selfBaseName(TumorBarcode barcode) {
        return barcode.getKey() + SELF_SUFFIX;
    }

    private static PeptideSourceTable load(File neoFile, File selfFile) {
        HugoPeptideTable neoTable = HugoPeptideTable.load(neoFile);
        HugoPeptideTable selfTable = HugoPeptideTable.load(selfFile);

        return PeptideSourceTable.create(neoTable, selfTable);
    }

    /**
     * Stores the peptide-source table for a given tumor sample.
     *
     * @param barcode the barcode for the tumor sample.
     *
     * @param source the peptide-source table for the specified tumor
     * sample.
     */
    public void store(TumorBarcode barcode, PeptideSourceTable source) {
        source.getNeoPeptideTable().store(neoPeptideFile(barcode));
        source.getSelfPeptideTable().store(selfPeptideFile(barcode));
    }
}
