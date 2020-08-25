
package pepmhc.neo;

import java.util.List;

import jam.app.JamApp;
import jam.app.JamLogger;

import jean.hugo.HugoPeptideTable;
import jean.hugo.HugoSymbol;
import jean.maf.MAFFastaList;
import jean.missense.MissenseManager;
import jean.peptide.Peptide;
import jean.tcga.TumorBarcode;

/**
 * Generates neo-peptides and self-peptides derived from mutated
 * protein structures by passing the proteins through the default
 * antigen processor.
 */
public final class PeptideSourceProcessor extends JamApp {
    private final String missenseDir;
    private final String selfPepFile;
    private final String barcodeFile;
    private final String pepSourceDir;
    private final String agproPropFile;

    private int processed = 0;
    private List<TumorBarcode> barcodes;
    private MissenseManager missenseManager;
    private HugoPeptideTable selfPeptideRef;
    private PeptideSourceManager pepSourceManager;

    private PeptideSourceProcessor(String missenseDir,
                                   String selfPepFile,
                                   String barcodeFile,
                                   String pepSourceDir,
                                   String agproPropFile) {
        this.missenseDir = missenseDir;
        this.selfPepFile = selfPepFile;
        this.barcodeFile = barcodeFile;
        this.pepSourceDir = pepSourceDir;
        this.agproPropFile = agproPropFile;
    }

    /**
     * Runs the peptide-source processor.
     *
     * @param missenseDir the directory containing the missense FASTA
     * files.
     *
     * @param selfPepFile the name of the file containing the reference
     * self-peptidome.
     *
     * @param barcodeFile the name of the file containing tumor sample
     * barcodes to process.
     *
     * @param pepSourceDir the directory where the peptide-source files
     * will be written.
     *
     * @param agproPropFile optional property file used to customized
     * the antigen processor settings ({@code null} to use the default
     * parameters).
     */
    public static void run(String missenseDir,
                           String selfPepFile,
                           String barcodeFile,
                           String pepSourceDir,
                           String agproPropFile) {
        PeptideSourceProcessor processor =
            new PeptideSourceProcessor(missenseDir,
                                       selfPepFile,
                                       barcodeFile,
                                       pepSourceDir,
                                       agproPropFile);
        processor.run();
    }

    private void run() {
        barcodes = TumorBarcode.load(barcodeFile);
        selfPeptideRef = HugoPeptideTable.load(selfPepFile);
        missenseManager = MissenseManager.create(missenseDir);
        pepSourceManager = PeptideSourceManager.create(pepSourceDir, selfPeptideRef);

        writeRuntimeEnv("JAM_", "JEAN_", "PEPMHC_");
        writeRuntimeProperties("jam.", "jean.", "pepmhc.");

        processBarcodes();
        JamLogger.info("DONE!");
    }

    private void processBarcodes() {
        barcodes.parallelStream().forEach(barcode -> processBarcode(barcode));
    }

    private void processBarcode(TumorBarcode barcode) {
        ++processed;
        JamLogger.info("Processing barcode [%s] (%d of %d)...", barcode.getKey(), processed, barcodes.size());

        if (pepSourceManager.exists(barcode)) {
            JamLogger.info("Peptide source already exists for barcode [%s]; skipping...", barcode.getKey());
            return;
        }

        MAFFastaList fastaList =
            missenseManager.load(barcode).sort();

        JamLogger.info("Found [%d] mutated peptides for barcode [%s]...", fastaList.size(), barcode.getKey());

        if (fastaList.isEmpty())
            return;

        PeptideSourceView sourceView =
            PeptideSourceEngine.process(barcode, fastaList, selfPeptideRef, agproPropFile);

        pepSourceManager.store(barcode, sourceView);
    }

    private static void usage() {
        System.err.print("Usage: pepmhc.neo.PeptideSourceProcessor");
        System.err.print(" MISSENSE_DIR");
        System.err.print(" SELF_PEP_FILE");
        System.err.print(" BARCODE_FILE");
        System.err.print(" PEP_SOURCE_DIR");
        System.err.print(" [AGPRO_PROP_FILE]");
        System.err.println();
        System.exit(1);
    }

    public static void main(String[] args) {
        if (args.length < 4 || args.length > 5)
            usage();

        String missenseDir = args[0];
        String selfPepFile = args[1];
        String barcodeFile = args[2];
        String pepSourceDir = args[3];
        String agproPropFile = null;

        if (args.length == 5)
            agproPropFile = args[4];

        run(missenseDir, selfPepFile, barcodeFile, pepSourceDir, agproPropFile);
    }
}
