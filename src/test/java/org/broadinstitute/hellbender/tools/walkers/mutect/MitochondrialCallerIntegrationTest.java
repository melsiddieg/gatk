package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

@Test(groups = {"variantcalling"})
public class MitochondrialCallerIntegrationTest extends CommandLineProgramTest {
    private static final File NA12878_MITO_BAM = new File(toolsTestDir, "mutect/mito/NA12878.bam");
    private static final File MITO_REF = new File(toolsTestDir, "mutect/mito/Homo_sapiens_assembly38.mt_only.fasta");

    @Test
    public void testMitochondria() throws Exception {
        Utils.resetRandomGenerator();
        final File unfilteredVcf = createTempFile("unfiltered", ".vcf");

        final List<String> args = Arrays.asList("-I", NA12878_MITO_BAM.getAbsolutePath(),
                "-R", MITO_REF.getAbsolutePath(),
                "-L", "chrM:1-1000",
                "-min-pruning", "5",
                "--" + MitochondrialCallerArgumentCollection.MEDIAN_AUTOSOMAL_COVERAGE_LONG_NAME, "1556", //arbitrary "autosomal" mean coverage used only for testing
                "-O", unfilteredVcf.getAbsolutePath());
        runCommandLine(args);


        final List<VariantContext> variants = VariantContextTestUtils.streamVcf(unfilteredVcf).collect(Collectors.toList());
        final Set<String> variantKeys = variants.stream().map(vc -> Mutect2IntegrationTest.keyForVariant(vc)).collect(Collectors.toSet());

        final List<String> expectedKeys = Arrays.asList(
                "chrM:152-152 [T*, C]",
                "chrM:263-263 [A*, G]",
                "chrM:301-301 [A*, AC]",
                "chrM:302-302 [A*, AC, C, ACC]",
                "chrM:310-310 [T*, TC]",
                "chrM:750-750 [A*, G]");
        Assert.assertTrue(expectedKeys.stream().allMatch(variantKeys::contains));

        Assert.assertEquals(variants.get(0).getGenotype("NA12878").getAnyAttribute(GATKVCFConstants.ORIGINAL_CONTIG_MISMATCH_KEY), "0");
        Assert.assertEquals(variants.get(0).getGenotype("NA12878").getAnyAttribute(GATKVCFConstants.POTENTIAL_POLYMORPHIC_NUMT_KEY), "true");
    }
}
