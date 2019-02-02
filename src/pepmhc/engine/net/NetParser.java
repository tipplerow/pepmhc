
package pepmhc.engine.net;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

import jam.io.IOUtil;
import jam.lang.JamException;
import jam.peptide.Peptide;
import jam.util.RegexUtil;

import pepmhc.binder.BindingRecord;

/**
 * Parses output written by the {@code netMHC} and {@code netMHCpan}
 * binding predictor programs.
 */
public final class NetParser {
    private final BufferedReader reader;
    private final List<BindingRecord> records = new ArrayList<BindingRecord>();

    private static final String DASHED_LINE_MATCH = "----------";

    private static final Pattern DATA_LINE_DELIM = RegexUtil.MULTI_WHITE_SPACE;

    private static final int PEPTIDE_FIELD_INDEX = 2;
    private static final int AFFINITY_FIELD_INDEX = 12;
    private static final int PERCENTILE_FIELD_INDEX = 13;

    private NetParser(BufferedReader reader) {
        this.reader = reader;
    }

    /**
     * Parses an output file written by {@code netMHC} or {@code netMHCpan}.
     *
     * @param file the path to the output file to parse.
     *
     * @return a list containing every binding record in the output file.
     *
     * @throws RuntimeException if any I/O errors occur.
     */
    public static List<BindingRecord> parse(File file) {
        return parse(IOUtil.openReader(file));
    }

    /**
     * Parses an output file written by {@code netMHC} or {@code netMHCpan}.
     *
     * @param fileName the name to the output file to parse.
     *
     * @return a list containing every binding record in the output file.
     *
     * @throws RuntimeException if any I/O errors occur.
     */
    public static List<BindingRecord> parse(String fileName) {
        return parse(IOUtil.openReader(fileName));
    }

    /**
     * Parses an output file written by {@code netMHC} or {@code netMHCpan}.
     *
     * @param reader a buffered reader opened to the beginning of the output
     * file; the reader is closed before returning from this method.
     *
     * @return a list containing every binding record in the output file.
     *
     * @throws RuntimeException if any I/O errors occur.
     */
    public static List<BindingRecord> parse(BufferedReader reader) {
        NetParser parser = new NetParser(reader);
        return parser.parse();
    }

    private List<BindingRecord> parse() {
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

        Peptide peptide    = Peptide.parse(fields[PEPTIDE_FIELD_INDEX]);
        double  affinity   = Double.valueOf(fields[AFFINITY_FIELD_INDEX]);
        double  percentile = Double.valueOf(fields[PERCENTILE_FIELD_INDEX]);

        records.add(new BindingRecord(peptide, affinity, percentile));
    }
}
