
package pepmhc.junit;

import java.util.List;

import jam.math.Probability;
import jam.math.UnitIndex;
import jam.math.UnitIndexRange;

import jene.peptide.Peptide;

import pepmhc.chop.NetChop;
import pepmhc.chop.NetChopEngine;

import org.junit.*;
import static org.junit.Assert.*;

public class NetChopEngineTest {
    // 123456789012345678901234
    // *  *       **      * ***
    // MAGRSGDNDEELLKAVRIIKILYK

    private static final Peptide PEPTIDE =
        Peptide.instance("MAGRSGDNDEELLKAVRIIKILYK");

    private static final List<Peptide> FRAGMENTS =
        List.of(Peptide.instance("SGDNDEELL"),
                Peptide.instance("KAVRIIKIL"));

    private static final double TOLERANCE = 0.0001;

    @Test public void testChop() {
        if (!NetChop.isInstalled())
            return;

        NetChopEngine engine = NetChopEngine.run(PEPTIDE);
        assertEquals(FRAGMENTS, engine.chop(9, Probability.ONE_HALF));
    }

    @Test public void testCleavageProb() {
        if (!NetChop.isInstalled())
            return;

        NetChopEngine engine = NetChopEngine.run(PEPTIDE);

        assertTrue(engine.getNTerminusCleavageProb(UnitIndex.instance(1)).equals(1.0,    TOLERANCE));
        assertTrue(engine.getNTerminusCleavageProb(UnitIndex.instance(2)).equals(0.7606, TOLERANCE));
        assertTrue(engine.getNTerminusCleavageProb(UnitIndex.instance(3)).equals(0.4834, TOLERANCE));

        assertTrue(engine.getCTerminusCleavageProb(UnitIndex.instance(1)).equals(0.7606, TOLERANCE));
        assertTrue(engine.getCTerminusCleavageProb(UnitIndex.instance(2)).equals(0.4834, TOLERANCE));
        assertTrue(engine.getCTerminusCleavageProb(UnitIndex.instance(3)).equals(0.0885, TOLERANCE));
        /*
        for (UnitIndex index = UnitIndex.first(); PEPTIDE.contains(index); index = index.next()) {
            double score = engine.getCTerminusCleavageProb(index).doubleValue();
            String symbol = score > 0.5 ? "*" : "";

            System.out.println(String.format("%s: %6.4f %s", index, score, symbol));
        }
        */
        assertTrue(engine.computeCleavageProb(UnitIndexRange.instance(1, 9)).equals(0.0229, TOLERANCE));
        assertTrue(engine.computeCleavageProb(UnitIndexRange.instance(2, 10)).equals(0.7606 * 0.0249, TOLERANCE));
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.NetChopEngineTest");
    }
}
