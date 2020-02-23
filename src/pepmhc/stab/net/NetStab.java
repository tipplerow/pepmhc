
package pepmhc.stab.net;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import jam.app.JamProperties;
import jam.hla.Allele;
import jam.lang.JamException;
import jam.peptide.Peptide;

import pepmhc.stab.StabilityRecord;

/**
 * Predicts the stability of peptide-MHC complexes using the
 * {@code netMHCstabpan} engine.
 */ 
public final class NetStab {
    private final Allele allele;
    private final Collection<Peptide> peptides;
    private final List<StabilityRecord> records;

    // Peptides passed to the netMHCstabpan executable must have
    // the same length, so we need to separate the peptide input
    // by length and then recompose the results in the original
    // peptide order...
    private List<Peptide> pep8  = new ArrayList<Peptide>();
    private List<Peptide> pep9  = new ArrayList<Peptide>();
    private List<Peptide> pep10 = new ArrayList<Peptide>();
    private List<Peptide> pep11 = new ArrayList<Peptide>();

    private List<StabilityRecord> stab8;
    private List<StabilityRecord> stab9;
    private List<StabilityRecord> stab10;
    private List<StabilityRecord> stab11;

    private NetStab(Allele allele, Collection<Peptide> peptides) {
        this.allele   = allele;
        this.peptides = peptides;
        this.records  = new ArrayList<StabilityRecord>(peptides.size());
    }

    /**
     * Name of the environment variable that defines the full path to
     * the {@code netMHCstsabpan} executable.  If the system property
     * {@code pepmhc.netMHCstabpan} is also defined, it will take
     * precedence.
     */
    public static final String EXECUTABLE_PATH_ENV = "NET_MHC_STAB_PAN_EXE";

    /**
     * Name of the system property that defines the full path to the
     * {@code netMHCpan} executable file.
     */
    public static final String EXECUTABLE_PATH_PROPERTY = "pepmhc.netMHCstabpan";

    /**
     * Determines whether the {@code netMHCstabpan} command-line
     * program is installed and is executable.
     *
     * @return {@code true} iff the {@code netMHCstabpan} executable
     * is installed at the resolved executable file location.
     */
    public static boolean isInstalled() {
        return resolveExecutableFile().canExecute();
    }

    /**
     * Resolves the full path to the {@code netMHCstabpan} executable file.
     *
     * @return the path specified by the {@code EXECUTABLE_PATH_PROPERTY}
     * (if set), or the path specified by the {@code EXECUTABLE_PATH_ENV}
     * environment variable (otherwise).
     *
     * @throws RuntimeException if neither property nor environment
     * variable is set.
     */
    public static File resolveExecutableFile() {
        return new File(resolveExecutableName());
    }

    /**
     * Resolves the full path to the {@code netMHCstabpan} executable file.
     *
     * @return the path specified by the {@code EXECUTABLE_PATH_PROPERTY}
     * (if set), or the path specified by the {@code EXECUTABLE_PATH_ENV}
     * environment variable (otherwise).
     *
     * @throws RuntimeException if neither property nor environment
     * variable is set.
     */
    public static String resolveExecutableName() {
        return JamProperties.resolve(EXECUTABLE_PATH_PROPERTY, EXECUTABLE_PATH_ENV, null);
    }

    /**
     * Predicts the MHC-peptide complex stability for a given allele
     * and target peptide.
     *
     * @param allele the binding MHC allele.
     *
     * @param peptide the peptide target.
     *
     * @return the stability record for the specified allele and
     * peptide.
     */
    public static StabilityRecord run(Allele allele, Peptide peptide) {
        return run(allele, List.of(peptide)).get(0);
    }

    /**
     * Predicts the MHC-peptide complex stability for a given allele
     * and a collection of target peptides.
     *
     * @param allele the binding MHC allele.
     *
     * @param peptides the peptide targets.
     *
     * @return the stability records for the specified allele and
     * peptides.
     */
    public static List<StabilityRecord> run(Allele allele, Collection<Peptide> peptides) {
        NetStab netStab = new NetStab(allele, peptides);
        return netStab.run();
    }

    private List<StabilityRecord> run() {
        //
        // Peptides passed to the netMHCstabpan executable must have
        // the same length, so we need to separate the peptide input
        // by length and then recompose the results in the original
        // peptide order...
        //
        mapPeptides();
        predictUniform();
        reduceRecords();

        return records;
    }

    private void mapPeptides() {
        for (Peptide peptide : peptides) {
            switch (peptide.length()) {
            case 8:
                pep8.add(peptide);
                break;

            case 9:
                pep9.add(peptide);
                break;

            case 10:
                pep10.add(peptide);
                break;

            case 11:
                pep11.add(peptide);
                break;

            default:
                throw JamException.runtime("Invalid peptide length: [%d].", peptide.length());
            }
        }
    }

    private void predictUniform() {
        stab8  = predictUniform(pep8);
        stab9  = predictUniform(pep9);
        stab10 = predictUniform(pep10);
        stab11 = predictUniform(pep11);
    }

    private List<StabilityRecord> predictUniform(List<Peptide> uniform) {
        if (uniform.isEmpty())
            return Collections.emptyList();
        else
            return NetStabBatch.run(allele, uniform);
    }

    private void reduceRecords() {
        Iterator<StabilityRecord> iter8  = stab8.iterator();
        Iterator<StabilityRecord> iter9  = stab9.iterator();
        Iterator<StabilityRecord> iter10 = stab10.iterator();
        Iterator<StabilityRecord> iter11 = stab11.iterator();

        for (Peptide peptide : peptides) {
            switch (peptide.length()) {
            case 8:
                records.add(iter8.next());
                break;

            case 9:
                records.add(iter9.next());
                break;

            case 10:
                records.add(iter10.next());
                break;

            case 11:
                records.add(iter11.next());
                break;

            default:
                throw new IllegalStateException("Invalid peptide passed the length screen.");
            }
        }

        if (iter8.hasNext())
            throw new IllegalStateException("Not all 8-mers were processed.");

        if (iter9.hasNext())
            throw new IllegalStateException("Not all 9-mers were processed.");

        if (iter10.hasNext())
            throw new IllegalStateException("Not all 10-mers were processed.");

        if (iter11.hasNext())
            throw new IllegalStateException("Not all 11-mers were processed.");

        if (records.size() != peptides.size())
            throw new IllegalStateException("Not all peptides were processed.");
    }
}
