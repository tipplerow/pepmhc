
package pepmhc.cache;

import jam.app.JamProperties;
import jam.hla.Allele;
import jam.io.FileUtil;

import pepmhc.engine.PredictionMethod;

/**
 * Determines unique database files for HLA alleles and prediction
 * methods.
 */
public final class AffinityDbFile {
    private final Allele allele;
    private final PredictionMethod method;

    private AffinityDbFile(PredictionMethod method, Allele allele) {
        this.allele = allele;
        this.method = method;
    }

    /**
     * Name of the environment variable that specifies the directory
     * containing the persistent database store.  The system property
     * {@code pepmhc.cache.directory} will take precedence if both are
     * specified.
     */
    public static final String CACHE_DIRECTORY_ENV = "PEPMHC_AFFINITY_CACHE";

    /**
     * Name of the system property that specifies the directory
     * containing the persistent database store.
     */
    public static final String CACHE_DIRECTORY_PROPERTY = "pepmhc.affinity.cacheDir";

    /**
     * Determines the unique database file name for a given
     * prediction method and allele.
     *
     * @param method the affinity prediction method.
     *
     * @param allele the allele of the binding MHC molecule.
     *
     * @return the unique database file name for the specified
     * prediction method and allele.
     */
    public static String resolve(PredictionMethod method, Allele allele) {
        AffinityDbFile dbFile = new AffinityDbFile(method, allele);
        return dbFile.resolve();
    }

    private String resolve() {
        return FileUtil.join(resolveDbDir(), formatDbFile());
    }

    private static String resolveDbDir() {
        String dirName = JamProperties.resolve(CACHE_DIRECTORY_PROPERTY, CACHE_DIRECTORY_ENV, null);
        FileUtil.ensureDir(dirName);
        return dirName;
    }

    private String formatDbFile() {
        return String.format("%s_%s.db", method, formatAllele());
    }

    private String formatAllele() {
        return allele.longKey().replace("*", "-").replace(":", "-");
    }
}
