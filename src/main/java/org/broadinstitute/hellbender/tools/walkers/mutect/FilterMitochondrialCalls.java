package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.io.File;
import java.util.ArrayList;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * <p>Filter variants in a Mutect2 Mitochondrial VCF.</p>
 * <p>
 * <h3>Usage example</h3>
 * <pre>
 * gatk FilterMitochondrialCalls \
 *   -V mitochondria.vcf.gz \
 *   -O filtered.vcf.gz
 * </pre>
 */
@CommandLineProgramProperties(
        summary = "Filter mitochondrial SNVs and indels called by Mutect2",
        oneLineSummary = "Filter mitochondrial SNVs and indels called by Mutect2",
        programGroup = VariantFilteringProgramGroup.class
)
@DocumentedFeature
public final class FilterMitochondrialCalls extends VariantWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "The output filtered VCF file", optional = false)
    private final String outputVcf = null;

    @ArgumentCollection
    protected MitochondrialFiltersArgumentCollection MTFAC = new MitochondrialFiltersArgumentCollection();

    private VariantContextWriter vcfWriter;

    private Mutect2FilteringEngine filteringEngine;

    @Override
    public void onTraversalStart() {
        final VCFHeader inputHeader = getHeaderForVariants();
        final Set<VCFHeaderLine> headerLines = inputHeader.getMetaDataInSortedOrder().stream()
                .filter(line -> !line.getKey().equals(Mutect2FilteringEngine.FILTERING_STATUS_VCF_KEY)) //remove header line from Mutect2 stating that calls are unfiltered.
                .collect(Collectors.toSet());
        headerLines.add(new VCFHeaderLine(Mutect2FilteringEngine.FILTERING_STATUS_VCF_KEY, "These calls have been filtered by " + FilterMitochondrialCalls.class.getSimpleName() + " to label false positives with a list of failed filters and true positives with PASS."));

        GATKVCFConstants.MITOCHONDRIAL_FILTER_NAMES.stream().map(GATKVCFHeaderLines::getFilterLine).forEach(headerLines::add);

        headerLines.addAll(getDefaultToolVCFHeaderLines());

        final VCFHeader vcfHeader = new VCFHeader(headerLines, inputHeader.getGenotypeSamples());
        vcfWriter = createVCFWriter(new File(outputVcf));
        vcfWriter.writeHeader(vcfHeader);

        ArrayList<String> allSampleNames = vcfHeader.getSampleNamesInOrder();
        if(allSampleNames.size() != 1) {
            throw new UserException(String.format("Expected single sample VCF to filter, but there were %s samples.", allSampleNames.size()));
        }

        final String sampleName = allSampleNames.get(0);
        filteringEngine = new Mutect2FilteringEngine(MTFAC, sampleName, Optional.empty(), GATKVCFConstants.LOD_KEY);
    }

    @Override
    public void apply(final VariantContext vc, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        final FilterResult filterResult = filteringEngine.calculateMitochondrialFilters(MTFAC, vc);

        final VariantContextBuilder vcb = new VariantContextBuilder(vc);

        vcb.filters(filterResult.getFilters());
        filterResult.getAttributes().entrySet().forEach(e -> vcb.attribute(e.getKey(), e.getValue()));

        vcfWriter.add(vcb.make());
    }

    @Override
    public void closeTool() {
        if (vcfWriter != null) {
            vcfWriter.close();
        }
    }
}
