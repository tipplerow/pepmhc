
package pepmhc.junit;

import jam.hla.Allele;
import jam.peptide.Peptide;

import pepmhc.binder.BindingRecord;
import pepmhc.cache.AffinityCache;
import pepmhc.engine.PredictionMethod;
import pepmhc.engine.Predictor;
import pepmhc.engine.net.NetMHCPredictor;

import org.junit.*;
import static org.junit.Assert.*;

public class AffinityCacheTest {
    static {
        System.setProperty(AffinityCache.CACHE_DIRECTORY_PROPERTY, "data/cache");
        System.setProperty(NetMHCPredictor.EXECUTABLE_PATH_PROPERTY, "/Users/scott/local/netMHC-4.0/netMHC");
    }

    private final Allele allele = Allele.instance("HLA-A*02:01");
    private final PredictionMethod method = PredictionMethod.NET_MHC;

    @Test public void testNetMHC() {
        if (!Predictor.isInstalled(method))
            return;

        BindingRecord record = AffinityCache.get(method, allele, Peptide.parse("AEFGPWQTV"));
        System.out.println(record);

        AffinityCache.get(method, allele, Peptide.parse("AEFGPWQTV"));
        AffinityCache.get(method, allele, Peptide.parse("AEFGPWQTV"));
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.AffinityCacheTest");
    }
}
