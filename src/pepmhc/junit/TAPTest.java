
package pepmhc.junit;

import java.util.ArrayList;
import java.util.List;

import jam.io.TableReader;
import jam.math.StatUtil;
import jam.peptide.Peptide;
import jam.vector.VectorView;

import pepmhc.tap.TAP;

import org.junit.*;
import static org.junit.Assert.*;

public class TAPTest {
    private static final String DIEZ_S1_FILE = "data/test/DiezRiveroS1.csv";
    private static final String IEDB_FILE    = "data/test/IEDB.csv";

    @Test public void testDiezCorr() {
        List<Double> measuredList = new ArrayList<Double>();
        List<Double> predictedList = new ArrayList<Double>();

        TableReader reader = TableReader.open(DIEZ_S1_FILE);

        for (List<String> fields : reader) {
            Peptide peptide = Peptide.instance(fields.get(0));

            double measured  = Double.parseDouble(fields.get(1));
            double predicted = TAP.consensus().score(peptide);

            measuredList.add(measured);
            predictedList.add(predicted);
        }

        reader.close();
        assertTrue(StatUtil.cor(VectorView.wrap(measuredList), VectorView.wrap(predictedList)) > 0.75);
    }

    @Test public void testIEDBCorr() {
        List<Double> IEDB = new ArrayList<Double>();
        List<Double> ours = new ArrayList<Double>();

        TableReader reader = TableReader.open(IEDB_FILE);

        for (List<String> fields : reader) {
            Peptide peptide = Peptide.instance(fields.get(0));

            IEDB.add(Double.parseDouble(fields.get(1)));
            ours.add(TAP.consensus().score(peptide));
        }

        reader.close();
        // IEDB flips the sign...
        assertTrue(StatUtil.cor(VectorView.wrap(IEDB), VectorView.wrap(ours)) < -0.95);
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("pepmhc.junit.TAPTest");
    }
}
