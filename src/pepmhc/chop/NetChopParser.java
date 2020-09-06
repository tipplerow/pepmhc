
package pepmhc.chop;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

import jam.io.IOUtil;
import jam.lang.JamException;
import jam.math.Probability;
import jam.util.RegexUtil;

/**
 * Parses output written by the {@code netchop} program.
 */
public final class NetChopParser {
    private final BufferedReader reader;
    private final List<Probability> scores = new ArrayList<Probability>();

    private static final String DASHED_LINE_MATCH = "--------------------";
    private static final Pattern DATA_LINE_DELIM = RegexUtil.MULTI_WHITE_SPACE;

    private static final int SCORE_INDEX = 3;
    private static final int FIELD_COUNT = 5;

    private NetChopParser(BufferedReader reader) {
        this.reader = reader;
    }

    /**
     * Parses an output file written by {@code netchop}.
     *
     * @param file the path to the output file to parse.
     *
     * @return a list containing every cleavage score in the
     * output file.
     *
     * @throws RuntimeException if any I/O errors occur.
     */
    public static List<Probability> parse(File file) {
        return parse(IOUtil.openReader(file));
    }

    /**
     * Parses an output file written by {@code netchop}.
     *
     * @param fileName the name to the output file to parse.
     *
     * @return a list containing every cleavage score in the
     * output file.
     *
     * @throws RuntimeException if any I/O errors occur.
     */
    public static List<Probability> parse(String fileName) {
        return parse(IOUtil.openReader(fileName));
    }

    /**
     * Parses an output file written by {@code netchop}.
     *
     * @param reader a buffered reader opened to the beginning of
     * the output file; the reader is closed before returning from
     * this method.
     *
     * @return a list containing every cleavage score in the output
     * file.
     *
     * @throws RuntimeException if any I/O errors occur.
     */
    public static List<Probability> parse(BufferedReader reader) {
        NetChopParser parser = new NetChopParser(reader);
        return parser.parse();
    }

    private List<Probability> parse() {
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

        return scores;
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

        if (fields.length != FIELD_COUNT)
            throw JamException.runtime("Invalid data line [%s].", line);

        scores.add(Probability.parse(fields[SCORE_INDEX]));
    }
}
