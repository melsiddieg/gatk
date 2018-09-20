package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.downsampling.MutectDownsampler;
import org.broadinstitute.hellbender.utils.downsampling.ReadsDownsampler;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.io.File;
import java.util.*;

/**
 * A wrapper for Mutect2 which calls SNPs and indels at various allele fractions in mitochondrial data. See {@Mutect2}
 * for more details.
 *
 *  <pre>
 *  gatk MitochondrialCaller \
 *   -R reference.fa \
 *   -I sample.bam \
 *   --median-autosomal-coverage 30 \
 *   -O unfiltered.vcf
 * </pre>
 *
 * The filtering tool {@link FilterMitochondrialCalls} should be run on the output vcf from this tool to generate
 * the final vcf.
 */
@CommandLineProgramProperties(
        summary = "Call mitochondrial SNVs and indels via local assembly of haplotypes",
        oneLineSummary = "Call mitochondrial SNVs and indels via local assembly of haplotypes",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)
@DocumentedFeature
public class MitochondrialCaller extends AssemblyRegionWalker {

    @ArgumentCollection
    protected MitochondrialCallerArgumentCollection MTAC = new MitochondrialCallerArgumentCollection();

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "File to which variants should be written")
    public File outputVCF;

    private VariantContextWriter vcfWriter;

    private Mutect2Engine m2Engine;

    @Override
    protected int defaultMinAssemblyRegionSize() { return 50; }

    @Override
    protected int defaultMaxAssemblyRegionSize() { return 300; }

    @Override
    protected int defaultAssemblyRegionPadding() { return 100; }

    @Override
    protected int defaultMaxReadsPerAlignmentStart() { return 50; }

    @Override
    protected double defaultActiveProbThreshold() { return 0.002; }

    @Override
    protected int defaultMaxProbPropagationDistance() { return 50; }

    @Override
    protected boolean includeReadsWithDeletionsInIsActivePileups() { return true; }

    @Override
    public boolean useVariantAnnotations() { return true;}

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        List<ReadFilter> readFilters = new ArrayList<>(Mutect2Engine.makeStandardMutect2ReadFilters());
        readFilters.add(ReadFilterLibrary.NON_CHIMERIC_OA_FILTER);
        return readFilters;
    }

    @Override
    public ReadTransformer makePostReadFilterTransformer() {
        return super.makePostReadFilterTransformer().andThen(Mutect2Engine.makeStandardMutect2PostFilterReadTransformer(referenceArguments.getReferencePath(), !MTAC.dontClipITRArtifacts));
    }

    @Override
    public List<Class<? extends Annotation>> getDefaultVariantAnnotationGroups() {
        List<Class<? extends Annotation>> annotations = new ArrayList<>(Mutect2Engine.getStandardMutect2AnnotationGroups());
        annotations.add(StandardMitochondrialAnnotation.class);
        return new ArrayList<>(annotations);
    }

    @Override
    protected ReadsDownsampler createDownsampler() {
        return new MutectDownsampler(maxReadsPerAlignmentStart, 0, 1);
    }

    @Override
    public AssemblyRegionEvaluator assemblyRegionEvaluator() { return m2Engine; }

    @Override
    public void onTraversalStart() {
        Set<String> samples =  ReadUtils.getSamplesFromHeader(getHeaderForReads());
        if (samples.size() != 1 ){
            throw new UserException(String.format("The input bam has more than one sample: %s", Arrays.toString(samples.toArray())));
        }
        String sampleName = samples.iterator().next();
        M2ArgumentCollection m2Args =  new M2ArgumentCollection(MTAC, sampleName);
        VariantAnnotatorEngine annotatorEngine = new VariantAnnotatorEngine(makeVariantAnnotations(), null, Collections.emptyList(), false);
        m2Engine = new Mutect2Engine(m2Args, createOutputBamIndex, createOutputBamMD5, getHeaderForReads(), referenceArguments.getReferenceFileName(), annotatorEngine, GATKVCFConstants.LOD_KEY);
        vcfWriter = createVCFWriter(outputVCF);
        m2Engine.writeHeader(vcfWriter, getMitochondrialCallerVCFHeaderLines());
    }

    @Override
    public Collection<Annotation> makeVariantAnnotations(){
        final Collection<Annotation> annotations = super.makeVariantAnnotations();
        if (MTAC.autosomalCoverage > 0) {
            annotations.add(new PolymorphicNuMT(MTAC.autosomalCoverage));
        }
        return annotations;
    }

    @Override
    public Object onTraversalSuccess() {
        return "SUCCESS";
    }

    @Override
    public void apply(final AssemblyRegion region, final ReferenceContext referenceContext, final FeatureContext featureContext ) {
        m2Engine.callRegion(region, referenceContext, featureContext).forEach(vcfWriter::add);
    }

    @Override
    public void closeTool() {
        if (vcfWriter != null) {
            vcfWriter.close();
        }
        if (m2Engine != null) {
            m2Engine.shutdown();
        }
    }

    private Set<VCFHeaderLine> getMitochondrialCallerVCFHeaderLines() {
        final Set<VCFHeaderLine> headerInfo = new HashSet<>();
        headerInfo.add(new VCFHeaderLine("Mutect Version", m2Engine.getVersion()));
        headerInfo.add(new VCFHeaderLine(Mutect2FilteringEngine.FILTERING_STATUS_VCF_KEY, "Warning: unfiltered Mutect 2 calls.  Please run " + FilterMitochondrialCalls.class.getSimpleName() + " to remove false positives."));
        headerInfo.addAll(m2Engine.getAnnotationEngine().getVCFAnnotationDescriptions(false));
        headerInfo.addAll(getDefaultToolVCFHeaderLines());
        GATKVCFConstants.STANDARD_MITO_INFO_FIELDS.stream().map(GATKVCFHeaderLines::getInfoLine).forEach(headerInfo::add);
        return headerInfo;
    }
}
