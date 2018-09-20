package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.annotations.Test;

public final class FilterMitochondrialCallsIntegrationTest extends CommandLineProgramTest {
    private File localTestData = new File(getTestDataDir(), "mitochondria/unfiltered.vcf");
    private static final List<Set<String>> expectedFilters = new ArrayList<>(Arrays.<Set<String>>asList(Collections.EMPTY_SET,
            new HashSet<>(Arrays.asList(GATKVCFConstants.CHIMERIC_ORIGINAL_ALIGNMENT_FILTER_NAME)),
            new HashSet<>(Arrays.asList(GATKVCFConstants.LOW_LOD_FILTER_NAME, GATKVCFConstants.LOW_AVG_ALT_QUALITY_FILTER_NAME)),
            Collections.EMPTY_SET, Collections.EMPTY_SET, Collections.EMPTY_SET));

    @Test()
    public void testFilter() {
        final File filteredVcf = createTempFile("filtered", ".vcf");

        final List<String> args = Arrays.asList("-V", localTestData.getPath(),
                "-O", filteredVcf.getPath());
        runCommandLine(args);

        final List<VariantContext> variants = VariantContextTestUtils.streamVcf(filteredVcf).collect(Collectors.toList());
        final Iterator<Set<String>> expected = expectedFilters.iterator();

        for(VariantContext v : variants){
            assertEqualsSet(v.getFilters(), expected.next(), "filters don't match expected");
        }
    }
}
