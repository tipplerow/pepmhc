
package pepmhc.bind;

import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;

import jam.sql.SQLColumn;
import jam.sql.SQLDb;
import jam.sql.SQLKeyTable;

import jene.hla.Allele;
import jene.peptide.Peptide;

/**
 * Provides a base class for persistent database tables of binding
 * records indexed by peptide.
 *
 * @param <R> the type of binding record (affinity or stability)
 * stored in the table.
 */
public abstract class BindTable<R extends BindRecord> extends SQLKeyTable<Peptide, R> {
    /**
     * Creates a new binding record table with a fixed database
     * manager.
     *
     * @param db the database manager.
     */
    protected BindTable(SQLDb db) {
        super(db);
    }

    /**
     * The name of the {@code peptide} column.
     */
    public static final String PEPTIDE_NAME = "peptide";

    /**
     * The name of the {@code percentile} column.
     */
    public static final String PERCENTILE_NAME = "percentile";

    /**
     * Meta-data for the {@code peptide} column.
     */
    public static final SQLColumn PEPTIDE_COLUMN =
        SQLColumn.create(PEPTIDE_NAME, "string")
        .primaryKey();

    /**
     * Meta-data for the {@code percentile} column.
     */
    public static final SQLColumn PERCENTILE_COLUMN =
        SQLColumn.create(PERCENTILE_NAME, "double");

    /**
     * Replaces special characters in allele names to allow them to be
     * used in file names.
     *
     * @param allele the allele to format.
     *
     * @return a properly formatted string unique to the allele.
     */
    public static String formatAllele(Allele allele) {
        return allele.longKey().replace("*", "-").replace(":", "-");
    }

    @Override public Peptide getKey(R record) {
        return record.getPeptide();
    }

    @Override public Peptide getKey(ResultSet resultSet, String columnLabel) throws SQLException {
        return Peptide.instance(resultSet.getString(columnLabel));
    }

    @Override public void prepareKey(PreparedStatement statement, int index, Peptide peptide) throws SQLException {
        statement.setString(index, peptide.formatString());
    }
}
