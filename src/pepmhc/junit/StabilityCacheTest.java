
package pepmhc.junit;

import jam.hla.Allele;
import jam.peptide.Peptide;

import jam.stab.NetStab;
import jam.stab.StabilityCache;
import jam.stab.StabilityRecord;

import org.junit.*;
import static org.junit.Assert.*;

public class StabilityCacheTest {
    static {
        System.setProperty(StabilityCache.CACHE_DIRECTORY_PROPERTY, "data/cache");
    }

    private final Allele allele = Allele.instance("HLA-A*02:01");

    @Test public void testNetMHC() {
        if (!NetStab.isInstalled())
            return;

        StabilityRecord record = StabilityCache.get(allele, Peptide.parse("AEFGPWQTV"));
        System.out.println(record);

        StabilityCache.get(allele, Peptide.parse("AEFGPWQTV"));
        StabilityCache.get(allele, Peptide.parse("AEFGPWQTV"));
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.StabilityCacheTest");
    }
}
