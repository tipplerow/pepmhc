
package pepmhc.neo;

import java.util.List;

import jam.app.JamLogger;
import jam.hugo.HugoPeptideTable;
import jam.hugo.HugoSymbol;
import jam.lang.JamException;
import jam.maf.MAFFastaList;
import jam.maf.MAFFastaRecord;
import jam.peptide.Peptide;
import jam.tcga.TumorBarcode;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

import pepmhc.agpro.AntigenProcessor;

final class PeptideSourceEngine {
    private final TumorBarcode barcode;
    private final MAFFastaList fastaList;
    private final HugoPeptideTable selfPepReference;

    private final AntigenProcessor antigenProcessor =
        AntigenProcessor.defaultProcessor();

    private final Multimap<HugoSymbol, Peptide> neoPeptideMap = HashMultimap.create();
    private final Multimap<HugoSymbol, Peptide> selfPeptideMap = HashMultimap.create();

    private PeptideSourceEngine(TumorBarcode barcode,
                                MAFFastaList fastaList,
                                HugoPeptideTable selfPepReference) {
        this.barcode = barcode;
        this.fastaList = fastaList;
        this.selfPepReference = selfPepReference;
    }

    static PeptideSourceView process(TumorBarcode barcode,
                                     MAFFastaList fastaList,
                                     HugoPeptideTable selfPepReference) {
        PeptideSourceEngine engine =
            new PeptideSourceEngine(barcode, fastaList, selfPepReference);

        return engine.process();
    }

    private PeptideSourceView process() {
        processRecords();
        
        HugoPeptideTable neoPeptideTable = HugoPeptideTable.create(neoPeptideMap);
        HugoPeptideTable selfPeptideTable = HugoPeptideTable.create(selfPeptideMap);

        return PeptideSourceView.create(neoPeptideTable, selfPeptideTable, selfPepReference);
    }

    private void processRecords() {
        for (MAFFastaRecord record : fastaList)
            processRecord(record);
    }

    private void processRecord(MAFFastaRecord record) {
        HugoSymbol symbol = record.getHugoSymbol();
        JamLogger.info("Processing gene [%s:%s]...", barcode.getKey(), symbol.getKey());
        
        List<Peptide> fragments =
            antigenProcessor.process(record.getPeptide());

        for (Peptide fragment : fragments)
            if (isSelfPeptide(fragment))
                selfPeptideMap.put(symbol, fragment);
            else
                neoPeptideMap.put(symbol, fragment);
    }

    private boolean isSelfPeptide(Peptide peptide) {
        return selfPepReference.contains(peptide);
    }
}
