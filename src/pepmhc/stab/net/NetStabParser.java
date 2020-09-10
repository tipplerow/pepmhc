
package pepmhc.stab.net;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

import jam.io.IOUtil;
import jam.lang.JamException;
import jam.math.Percentile;
import jam.util.RegexUtil;

import jene.chem.HalfLife;
import jene.peptide.Peptide;

import pepmhc.stab.StabilityRecord;

/**
 * Parses output written by the {@code netMHCstabpan} stability
 * predictor program.
 */
public final class NetStabParser {
    private final BufferedReader reader;
    private final List<StabilityRecord> records = new ArrayList<StabilityRecord>();

    private static final String DASHED_LINE_MATCH = "----------";

    private static final Pattern DATA_LINE_DELIM = RegexUtil.MULTI_WHITE_SPACE;

    private static final int PEPTIDE_FIELD_INDEX = 2;
    private static final int HALF_LIFE_FIELD_INDEX = 5;
    private static final int PERCENTILE_FIELD_INDEX = 6;

    protected NetStabParser(BufferedReader reader) {
        this.reader = reader;
    }

    /**
     * Parses an output file written by {@code netMHCstabpan}.
     *
     * @param file the path to the output file to parse.
     *
     * @return a list containing every binding record in the output file.
     *
     * @throws RuntimeException if any I/O errors occur.
     */
    public static List<StabilityRecord> parse(File file) {
        return parse(IOUtil.openReader(file));
    }

    /**
     * Parses an output file written by {@code netMHCstabpan}.
     *
     * @param fileName the name to the output file to parse.
     *
     * @return a list containing every binding record in the output file.
     *
     * @throws RuntimeException if any I/O errors occur.
     */
    public static List<StabilityRecord> parse(String fileName) {
        return parse(IOUtil.openReader(fileName));
    }

    /**
     * Parses an output file written by {@code netMHCstabpan}.
     *
     * @param reader an open reader for the output file.
     *
     * @return a list containing every binding record in the output file.
     *
     * @throws RuntimeException if any I/O errors occur.
     */
    public static List<StabilityRecord> parse(BufferedReader reader) {
        NetStabParser parser = new NetStabParser(reader);
        return parser.parse();
    }

    private List<StabilityRecord> parse() {
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

        return records;
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

        if (fields.length <= PERCENTILE_FIELD_INDEX)
            throw JamException.runtime("Invalid data line [%s].", line);

        Peptide peptide = Peptide.instance(fields[PEPTIDE_FIELD_INDEX]);
        HalfLife halfLife = HalfLife.parse(fields[HALF_LIFE_FIELD_INDEX]);
        Percentile percentile = Percentile.parse(fields[PERCENTILE_FIELD_INDEX]);

        records.add(new StabilityRecord(peptide, halfLife, percentile));
    }
}
