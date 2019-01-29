
package pepmhc.junit;

import java.util.List;

import jam.peptide.Peptide;

import pepmhc.cache.AffinityCache;
import pepmhc.engine.BindingRecord;
import pepmhc.engine.PredictionMethod;
import pepmhc.engine.Predictor;
import pepmhc.engine.net.NetMHCPredictor;

import org.junit.*;
import static org.junit.Assert.*;

public class AffinityCacheTest {
    static {
        System.setProperty(AffinityCache.CACHE_DIRECTORY_PROPERTY, "data/test");
        System.setProperty(NetMHCPredictor.EXECUTABLE_PATH_PROPERTY, "/Users/scott/local/netMHC-4.0/netMHC");
    }

    private final PredictionMethod method = PredictionMethod.NET_MHC;

    @Test public void testNetMHC() {
        if (!Predictor.isInstalled(method))
            return;

        BindingRecord record = AffinityCache.get(method, "HLA-A*02:01", Peptide.parse("AEFGPWQTV"));
        System.out.println(record);

        AffinityCache.get(method, "HLA-A*02:01", Peptide.parse("AEFGPWQTV"));
        AffinityCache.get(method, "HLA-A*02:01", Peptide.parse("AEFGPWQTV"));
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.AffinityCacheTest");
    }
}
