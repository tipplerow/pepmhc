
package pepmhc.junit;

import java.util.List;

import jam.math.Probability;

import jene.peptide.Peptide;

import pepmhc.chop.NetChop;
import pepmhc.chop.NetChopRunner;

import org.junit.*;
import static org.junit.Assert.*;

public class NetChopRunnerTest {
    private static final Peptide PEPTIDE =
        Peptide.instance("MAGRSGDNDEELLKAVRIIKILYKSNPYPEPKGSRQARKN");

    @Test public void testChop() {
        if (!NetChop.isInstalled())
            return;

        List<Probability> scores = NetChopRunner.score(PEPTIDE);

        assertEquals(PEPTIDE.length(), scores.size());

        assertTrue(scores.get(3).equals(0.7831, 0.0001));
        assertTrue(scores.get(4).equals(0.0255, 0.0001));
        assertTrue(scores.get(11).equals(0.5506, 0.0001));
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.NetChopRunnerTest");
    }
}
