
package pepmhc.junit;

import jam.hla.Allele;
import jam.peptide.Peptide;
import jam.stab.NetStab;

import org.junit.*;
import static org.junit.Assert.*;

public class NetStabTest {
    private final Allele A0101 = Allele.instance("HLA-A*01:01");
    private final Allele A0201 = Allele.instance("HLA-A*02:01");

    @Test public void testExecutable() {
        System.out.print("NetStab ");
        System.out.print(NetStab.isInstalled() ? "IS" : "is NOT");
        System.out.print(" installed.");
        System.out.println();
    }

    @Test public void testPredict1() {
        if (!NetStab.isInstalled())
            return;

        assertEquals(0.43, NetStab.run(A0101, Peptide.parse("AAAWYLWEV")).getHalfLife(), 0.01);
        assertEquals(0.11, NetStab.run(A0101, Peptide.parse("AEFGPWQTV")).getHalfLife(), 0.01);

        assertEquals(7.61, NetStab.run(A0201, Peptide.parse("AAAWYLWEV")).getHalfLife(), 0.01);
        assertEquals(0.50, NetStab.run(A0201, Peptide.parse("AEFGPWQTV")).getHalfLife(), 0.01);
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.NetStabTest");
    }
}
