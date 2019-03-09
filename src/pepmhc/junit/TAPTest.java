
package pepmhc.junit;

import java.util.ArrayList;
import java.util.List;

import jam.io.LineReader;
import jam.io.ObjectReader;
import jam.math.StatUtil;
import jam.peptide.Peptide;
import jam.util.RegexUtil;
import jam.vector.VectorView;

import pepmhc.tap.TAP;

import org.junit.*;
import static org.junit.Assert.*;

public class TAPTest {
    private static final String DIEZ_S1_FILE = "data/test/DiezRiveroS1.csv";
    private static final String IEDB_FILE    = "data/test/IEDB.csv";

    @Test public void testDiezCorr() {
        List<Double> measured = new ArrayList<Double>();
        List<Double> predicted = new ArrayList<Double>();

        LineReader reader = LineReader.open(DIEZ_S1_FILE);
        reader.next(); // Skip header line

        for (String line : reader) {
            String[] fields = RegexUtil.split(RegexUtil.COMMA, line, 2);
            Peptide peptide = Peptide.parse(fields[0]);

            measured.add(Double.parseDouble(fields[1]));
            predicted.add(TAP.INSTANCE.score(peptide));
        }

        reader.close();
        assertTrue(StatUtil.cor(VectorView.wrap(measured), VectorView.wrap(predicted)) > 0.75);
    }

    @Test public void testIEDBCorr() {
        List<Double> IEDB = new ArrayList<Double>();
        List<Double> ours = new ArrayList<Double>();

        LineReader reader = LineReader.open(IEDB_FILE);
        reader.next(); // Skip header line

        for (String line : reader) {
            String[] fields = RegexUtil.split(RegexUtil.COMMA, line, 2);
            Peptide peptide = Peptide.parse(fields[0]);

            IEDB.add(Double.parseDouble(fields[1]));
            ours.add(TAP.INSTANCE.score(peptide));
        }

        reader.close();
        assertTrue(StatUtil.cor(VectorView.wrap(IEDB), VectorView.wrap(ours)) < -0.95);
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.TAPTest");
    }
}
