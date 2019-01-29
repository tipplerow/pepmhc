
package pepmhc.engine.smm;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

import jam.io.IOUtil;
import jam.lang.JamException;
import jam.peptide.Residue;
import jam.util.RegexUtil;

final class MatrixReader {
    private final String fileName;
    private final BufferedReader buffReader;

    private int pepLength;
    private double intercept;
    private List<Map<Residue, Double>> elements;

    private final static Pattern DELIM = RegexUtil.TAB;

    private MatrixReader(String fileName) {
        this.fileName = fileName;
        this.buffReader = IOUtil.openReader(fileName);
    }

    static StabilizedMatrix load(String fileName) {
        MatrixReader matrixReader = new MatrixReader(fileName);
        return matrixReader.load();
    }

    private StabilizedMatrix load() {
        try {
            parsePeptideLength();
            initMatrixElements();
            readMatrixElements();
            parseIntercept();
        }
        catch (Exception ex) {
            throw JamException.runtime(ex);
        }
        finally {
            IOUtil.close(buffReader);
        }

        return new StabilizedMatrix(elements, intercept);
    }

    private void parsePeptideLength() throws IOException {
        String   line   = buffReader.readLine();
        String[] fields = RegexUtil.split(DELIM, line, 2);

        pepLength = Integer.parseInt(fields[1]);
    }

    private void initMatrixElements() {
        elements = new ArrayList<Map<Residue, Double>>(pepLength);

        while (elements.size() < pepLength)
            elements.add(new EnumMap<Residue, Double>(Residue.class));
    }

    private void readMatrixElements() throws IOException {
        for (int index = 0; index < Residue.countNative(); ++index)
            readMatrixRow();
    }

    private void readMatrixRow() throws IOException {
        String   line   = buffReader.readLine();
        String[] fields = RegexUtil.split(DELIM, line, pepLength + 1);

        Residue residue = Residue.valueOfCode1(fields[0]);

        for (int position = 0; position < pepLength; ++position)
            elements.get(position).put(residue, Double.parseDouble(fields[position + 1]));
    }

    private void parseIntercept() throws IOException {
        String   line   = buffReader.readLine();
        String[] fields = RegexUtil.split(DELIM, line, 2);

        intercept = Double.parseDouble(fields[1]);
    }
}
