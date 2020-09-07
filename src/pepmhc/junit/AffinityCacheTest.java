
package pepmhc.junit;

import java.io.File;

import jene.hla.Allele;
import jene.peptide.Peptide;

import pepmhc.affy.AffinityCache;
import pepmhc.affy.AffinityMethod;
import pepmhc.affy.AffinityPredictor;
import pepmhc.affy.AffinityRecord;
import pepmhc.affy.AffinityTable;

import org.junit.*;
import static org.junit.Assert.*;

public class AffinityCacheTest {
    static {
        System.setProperty(AffinityCache.CACHE_DIRECTORY_PROPERTY, "data/cache");
    }

    private static final Allele allele = Allele.instance("HLA-A*02:01");
    private static final AffinityMethod method = AffinityMethod.NET_MHC_PAN;
    private static final AffinityPredictor predictor = method.getPredictor();

    @Before public void setUp() {
        deleteDbFile();
    }

    @AfterClass public static void tearDownClass() {
        deleteDbFile();
    }

    private static void deleteDbFile() {
        File file = new File(AffinityTable.dbFile(method, allele));

        if (file.exists())
            file.delete();
    }

    @Test public void testNetMHC() {
        if (!predictor.isInstalled())
            return;

        AffinityRecord record = AffinityCache.instance(method, allele).require(Peptide.instance("AEFGPWQTV"));

        AffinityCache.instance(method, allele).require(Peptide.instance("AEFGPWQTV"));
        AffinityCache.instance(method, allele).require(Peptide.instance("AEFGPWQTV"));
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.AffinityCacheTest");
    }
}
