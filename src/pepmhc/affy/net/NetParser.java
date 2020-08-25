
package pepmhc.affy.net;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

import jam.io.IOUtil;
import jam.lang.JamException;
import jam.util.RegexUtil;

import jean.peptide.Peptide;

import pepmhc.affy.AffinityRecord;

/**
 * Parses output written by the {@code netMHC}, {@code netMHCpan}, and
 * {@code netMHCstabpan} programs.
 */
public abstract class NetParser {
    private final BufferedReader reader;
    private final List<AffinityRecord> records = new ArrayList<AffinityRecord>();

    private static final String DASHED_LINE_MATCH = "----------";

    private static final Pattern DATA_LINE_DELIM = RegexUtil.MULTI_WHITE_SPACE;

    /**
     * Wraps a parser around an open stream reader.
     *
     * @param reader an open stream reader.
     *
     * @throws RuntimeException if the file cannot be opened.
     */
    protected NetParser(BufferedReader reader) {
        this.reader = reader;
    }

    /**
     * Returns the zero-offset index of the field containing the
     * peptide structure.
     *
     * @return the zero-offset index of the field containing the
     * peptide structure.
     */
    public abstract int getPeptideFieldIndex();

    /**
     * Returns the zero-offset index of the field containing the
     * binding affinity.
     *
     * @return the zero-offset index of the field containing the
     * binding affinity.
     */
    public abstract int getAffinityFieldIndex();

    /**
     * Returns the zero-offset index of the field containing the
     * affinity percentile.
     *
     * @return the zero-offset index of the field containing the
     * affinity percentile.
     */
    public abstract int getPercentileFieldIndex();

    /**
     * Parses the file opened by the constructor.
     *
     * @return the affinity records present in the output file.
     *
     * @throws RuntimeException if any I/O errors occur.
     */
    public List<AffinityRecord> parse() {
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

        if (fields.length <= getPercentileFieldIndex())
            throw JamException.runtime("Invalid data line [%s].", line);

        Peptide peptide    = Peptide.instance(fields[getPeptideFieldIndex()]);
        double  affinity   = Double.valueOf(fields[getAffinityFieldIndex()]);
        double  percentile = Double.valueOf(fields[getPercentileFieldIndex()]);

        records.add(new AffinityRecord(peptide, affinity, percentile));
    }
}
