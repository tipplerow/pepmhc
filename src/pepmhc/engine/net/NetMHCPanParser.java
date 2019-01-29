
package pepmhc.engine.net;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.regex.Pattern;

import jam.io.IOUtil;
import jam.lang.JamException;
import jam.peptide.Peptide;
import jam.util.RegexUtil;

final class NetMHCPanParser {
    private final BufferedReader reader;
    private final Map<Peptide, Double> affinities;

    private static final String DASHED_LINE_MATCH = "----------";

    private static final Pattern DATA_LINE_DELIM = RegexUtil.MULTI_WHITE_SPACE;

    private static final int PEPTIDE_FIELD_INDEX = 2;
    private static final int AFFINITY_FIELD_INDEX = 12;

    private NetMHCPanParser(BufferedReader reader) {
        this.reader = reader;
        this.affinities = new LinkedHashMap<Peptide, Double>();
    }

    static Map<Peptide, Double> parse(BufferedReader reader) {
        NetMHCPanParser parser = new NetMHCPanParser(reader);
        return parser.parse();
    }

    private Map<Peptide, Double> parse() {
        try {
            skipHeader();
            parseData();
        }
        catch (IOException ioex) {
            throw JamException.runtime(ioex);
        }
        finally {
            IOUtil.close(reader);
        }

        return affinities;
    }

    private void skipHeader() throws IOException {
        //
        // The header contains two dashed lines; the data begins
        // immediately after the second dashed line...
        //
        readToDashedLine();
        readToDashedLine();
    }

    private void readToDashedLine() throws IOException {
        while (true) {
            String line = reader.readLine();

            if (line == null || isDashedLine(line))
                return;
        }
    }

    private static boolean isDashedLine(String line) {
        return line != null && line.startsWith(DASHED_LINE_MATCH);
    }

    private void parseData() throws IOException {
        //
        // The data section ends with a dashed line...
        //
        while (true) {
            String line = reader.readLine();

            if (line == null || isDashedLine(line))
                return;
            else
                parseLine(line);
        }
    }

    private void parseLine(String line) {
        String[] fields = DATA_LINE_DELIM.split(line.trim());

        if (fields.length <= AFFINITY_FIELD_INDEX)
            throw JamException.runtime("Invalid data line [%s].", line);

        Peptide peptide  = Peptide.parse(fields[PEPTIDE_FIELD_INDEX]);
        double  affinity = Double.valueOf(fields[AFFINITY_FIELD_INDEX]);

        affinities.put(peptide, affinity);
    }
}
