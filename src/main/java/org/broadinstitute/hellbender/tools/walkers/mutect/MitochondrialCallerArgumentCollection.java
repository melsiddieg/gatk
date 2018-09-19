package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerArgumentCollection;

public class MitochondrialCallerArgumentCollection extends AssemblyBasedCallerArgumentCollection {
    private static final long serialVersionUID = 9341L;
    public static final String EMISSION_LOD_LONG_NAME = "lod-to-emit";
    public static final String EMISSION_LOG_SHORT_NAME = "emit-lod";
    public static final String INITIAL_LOD_LONG_NAME = "initial-lod";
    public static final String INITIAL_LOD_SHORT_NAME = "init-lod";
    public static final String MEDIAN_AUTOSOMAL_COVERAGE_LONG_NAME = "median-autosomal-coverage";

    /**
     * Only variants with LODs exceeding this threshold will be written to the VCF, regardless of filter status.
     * Set to less than or equal to lod. Increase argument value to reduce false positives in the callset.
     * Default setting of 0 is optimized for high sensitivity and will emit negative training data that
     * {@link FilterMitochondrialCalls} should then filter.
     */
    @Argument(fullName = EMISSION_LOD_LONG_NAME, shortName = EMISSION_LOG_SHORT_NAME, optional = true, doc = "LOD threshold to emit variant to VCF.")
    public double emissionLod = 0;

    /**
     * Only variants with estimated LODs exceeding this threshold will be considered active. Default of 0 set to
     * maximize sensitivity. Increase this argument to decrease runtime at the cost of sensitivity.
     */
    @Argument(fullName = INITIAL_LOD_LONG_NAME, shortName = INITIAL_LOD_SHORT_NAME, optional = true, doc = "LOD threshold to consider pileup active.")
    public double initialLod = 0;

    /**
     * Two or more phased substitutions separated by this distance or less are merged into MNPs.
     */
    @Advanced
    @Argument(fullName = M2ArgumentCollection.MAX_MNP_DISTANCE_LONG_NAME, shortName = M2ArgumentCollection.MAX_MNP_DISTANCE_SHORT_NAME,
            doc = "Two or more phased substitutions separated by this distance or less are merged into MNPs.", optional = true)
    public int maxMnpDistance = 1;

    /**
     * When opposite ends of a fragment are inverted tandem repeats of each other, the sequence past one end may be copied onto the other
     * during library prep.  These artifacts are especially prevalent when DNA is damaged as in the case of FFPE samples and ancient DNA.
     * Artifact turned off by default in MitochondrialCaller.
     */
    @Argument(fullName= M2ArgumentCollection.IGNORE_ITR_ARTIFACTS_LONG_NAME, doc="Turn on read transformer that clips artifacts associated with end repair insertions near inverted tandem repeats.", optional = true)
    public boolean dontClipITRArtifacts = true;

    @Advanced
    @Argument(fullName = M2ArgumentCollection.GET_AF_FROM_AD_LONG_NAME, doc="Use allelic depth to calculate allele fraction", optional = true)
    public boolean calculateAFfromAD = true;

    /**
     * Used to model autosomal coverage when calling mitochondria. The median tends to be a more robust center statistic.
     */
    @Advanced
    @Argument(fullName = MEDIAN_AUTOSOMAL_COVERAGE_LONG_NAME, doc="Annotate possible polymorphic NuMT based on Poisson distribution given median autosomal coverage", optional = true)
    public double autosomalCoverage;
}
