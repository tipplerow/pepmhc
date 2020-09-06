
package pepmhc.junit;

import java.util.List;

import jam.math.Probability;

import jene.peptide.Peptide;

import pepmhc.chop.NetChop;
import pepmhc.chop.NetChopRunner;

import org.junit.*;
import static org.junit.Assert.*;

public class NetChopRunnerTest {
    //0123456789 0123456789 0123456789 0123456789
    //    |        ||     |
    //MAGRSGDNDE ELLKAVRIIK ILYKSNPYPE PKGSRQARKN

    private static final Peptide PEPTIDE =
        Peptide.instance("MAGRSGDNDEELLKAVRIIKILYKSNPYPEPKGSRQARKN");

    @Test public void testChop() {
        if (!NetChop.isInstalled())
            return;

        List<Probability> scores = NetChopRunner.score(PEPTIDE);

        System.out.println();
        for (int index = 0; index < scores.size(); ++index) {
            double score = scores.get(index).doubleValue();
            String symbol = score > 0.5 ? "*" : "";

            System.out.println(String.format("%d: %6.4f %s", index, score, symbol));
        }

        assertEquals(PEPTIDE.length(), scores.size());

        assertTrue(scores.get(3).equals(0.7831, 0.0001));
        assertTrue(scores.get(4).equals(0.0255, 0.0001));
        assertTrue(scores.get(11).equals(0.5506, 0.0001));
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.NetChopRunnerTest");
    }
}
