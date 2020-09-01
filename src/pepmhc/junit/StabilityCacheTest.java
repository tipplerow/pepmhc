
package pepmhc.junit;

import jene.hla.Allele;
import jene.peptide.Peptide;

import pepmhc.stab.net.NetStab;
import pepmhc.stab.StabilityCache;
import pepmhc.stab.StabilityMethod;
import pepmhc.stab.StabilityRecord;

import org.junit.*;
import static org.junit.Assert.*;

public class StabilityCacheTest {
    static {
        System.setProperty(StabilityCache.CACHE_DIRECTORY_PROPERTY, "data/cache");
    }

    private final Allele allele = Allele.instance("HLA-A*02:01");
    private final StabilityMethod method = StabilityMethod.NET_MHC_STAB_PAN;

    @Test public void testNetMHC() {
        if (!NetStab.isInstalled())
            return;

        StabilityCache cache = StabilityCache.instance(method, allele);
        StabilityRecord record = cache.require(Peptide.instance("AEFGPWQTV"));
        System.out.println(record);

        cache.require(Peptide.instance("AEFGPWQTV"));
        cache.require(Peptide.instance("AEFGPWQTV"));
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.StabilityCacheTest");
    }
}
