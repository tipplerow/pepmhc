
package pepmhc.agpro;

import java.io.File;
import java.util.Collection;
import java.util.Set;

import jam.app.JamProperties;
import jam.hugo.HugoSymbol;
import jam.hugo.HugoPeptideTable;
import jam.peptide.Peptide;

/**
 * Provides a database of self-peptides derived from germline protein
 * structures processed by the antigen processing machinery.
 */
public final class SelfPeptideDb {
    private final AntigenProcessor agProcessor;
    private final HugoPeptideTable peptideTable;

    private static SelfPeptideDb global = null;

    private SelfPeptideDb(AntigenProcessor agProcessor,
                          HugoPeptideTable peptideTable) {
        this.agProcessor  = agProcessor;
        this.peptideTable = peptideTable;
    }

    /**
     * Name of the system property that specifies the configuration
     * file for the global antigen processor.
     */
    public static final String CONFIG_FILE_PROPERTY = "pepmhc.agpro.selfDbConfigFile";

    /**
     * Name of the system property that specifies the file containing
     * self-peptides mapped to HUGO symbol.
     */
    public static final String PEPTIDE_FILE_PROPERTY = "pepmhc.agpro.selfDbPeptideFile";

    /**
     * Returns the global database defined by system properties.
     *
     * @return the global database defined by system properties.
     */
    public static SelfPeptideDb global() {
        if (global == null)
            global = createGlobal();

        return global;
    }

    private static SelfPeptideDb createGlobal() {
        return load(resolveConfigFile(), resolvePeptideFile());
    }

    private static String resolveConfigFile() {
        return JamProperties.getRequired(CONFIG_FILE_PROPERTY);
    }

    private static String resolvePeptideFile() {
        return JamProperties.getRequired(PEPTIDE_FILE_PROPERTY);
    }

    /**
     * Loads a database from configuration and peptide files.
     *
     * @param configFile the antigen processor configuration file.
     *
     * @param peptideFile the file containing self-peptides mapped to
     * HUGO symbol.
     *
     * @return a new database containing the antigen processor defined
     * by the configuration file and gene-peptide mappings contained
     * in the peptide file.
     *
     * @throws RuntimeException unless the input files can be opened
     * for reading and contain valid data.
     */
    public static SelfPeptideDb load(File configFile, File peptideFile) {
        AntigenProcessor agProcessor = AntigenProcessor.resolve(configFile);
        HugoPeptideTable peptideTable = HugoPeptideTable.load(peptideFile);

        return new SelfPeptideDb(agProcessor, peptideTable);
    }

    /**
     * Loads a database from configuration and peptide files.
     *
     * @param configFile the name of the antigen processor configuration
     * file.
     *
     * @param peptideFile the name of the file containing self-peptides
     * mapped to HUGO symbol.
     *
     * @return a new database containing the antigen processor defined
     * by the configuration file and gene-peptide mappings contained
     * in the peptide file.
     *
     * @throws RuntimeException unless the input files can be opened
     * for reading and contain valid data.
     */
    public static SelfPeptideDb load(String configFile, String peptideFile) {
        return load(new File(configFile), new File(peptideFile));
    }

    /**
     * Returns the antigen processor that generated the self-peptides
     * in this database.
     *
     * @return the antigen processor that generated the self-peptides
     * in this database.
     */
    public AntigenProcessor antigenProcessor() {
        return agProcessor;
    }

    /**
     * Identifies HUGO symbols contained in this database.
     *
     * @param symbol the HUGO symbol of interest.
     *
     * @return {@code true} iff this database contains one or more
     * peptides mapped to the specified symbol.
     */
    public boolean contains(HugoSymbol symbol) {
        return peptideTable.contains(symbol);
    }

    /**
     * Identifies peptides contained in this database.
     *
     * @param peptide the peptide of interest.
     *
     * @return {@code true} iff the peptide is contained in this
     * database (is derived from one or more germline proteins).
     */
    public boolean contains(Peptide peptide) {
        return peptideTable.contains(peptide);
    }

    /**
     * Returns all HUGO symbols contained in this database.
     *
     * @return an unmodifiable set containing all HUGO symbols
     * contained in this database.
     */
    public Set<HugoSymbol> geneSet() {
        return peptideTable.viewSymbols();
    }

    /**
     * Returns all peptides mapped to a given HUGO symbol.
     *
     * @param symbol the HUGO symbol of interest.
     *
     * @return an unmodifiable collection containing all peptides
     * mapped to the specified symbol (an empty list if this database
     * does not contain the symbol).
     */
    public Collection<Peptide> lookup(HugoSymbol symbol) {
        return peptideTable.get(symbol);
    }

    /**
     * Returns the number of gene-peptide mappings in this database.
     *
     * @return the number of gene-peptide mappings in this database.
     */
    public int size() {
        return peptideTable.size();
    }
}
