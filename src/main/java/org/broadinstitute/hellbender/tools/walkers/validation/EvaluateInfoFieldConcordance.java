package org.broadinstitute.hellbender.tools.walkers.validation;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.EnumMap;
import java.io.IOException;

import org.apache.commons.collections4.Predicate;
import org.apache.commons.lang.mutable.MutableLong;

import htsjdk.variant.variantcontext.VariantContext;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;

import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.engine.AbstractConcordanceWalker;

import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

@CommandLineProgramProperties(
        summary=EvaluateInfoFieldConcordance.USAGE_SUMMARY,
        oneLineSummary=EvaluateInfoFieldConcordance.USAGE_ONE_LINE_SUMMARY,
        programGroup=VariantEvaluationProgramGroup.class)
@DocumentedFeature
@BetaFeature
public class EvaluateInfoFieldConcordance extends AbstractConcordanceWalker {
    static final String USAGE_ONE_LINE_SUMMARY = "Evaluate concordance of info fields in an input VCF against a validated truth VCF";
    static final String USAGE_SUMMARY = "This tool evaluates info fields from an input VCF against a VCF that has been validated and is considered to represent ground truth.\n";
    public static final String SUMMARY_LONG_NAME = "summary";
    public static final String SUMMARY_SHORT_NAME = "S";
    @Argument(doc="A table of summary statistics (true positives, sensitivity, etc.)", fullName="summary", shortName=SUMMARY_SHORT_NAME)
    protected String summary;
    @Argument(fullName="eval-info-key", shortName="eval-info-key", doc="Info key from eval vcf", optional=true)
    protected String evalInfoKey = GATKVCFConstants.CNN_2D_KEY;
    @Argument(fullName="truth-info-key", shortName="truth-info-key", doc="Info key from truth vcf", optional=true)
    protected String truthInfoKey = GATKVCFConstants.CNN_2D_KEY;
    @Argument(fullName="epsilon", shortName="epsilon", doc="Difference tolerance", optional=true)
    protected double epsilon = 0.1;
    private final EnumMap<ConcordanceState, MutableLong> snpCounts = new EnumMap<>(ConcordanceState.class);
    private final EnumMap<ConcordanceState, MutableLong> indelCounts = new EnumMap<>(ConcordanceState.class);

    private double snpSumDelta = 0.0;
    private double snpSumDeltaSquared = 0.0;
    private double indelSumDelta = 0.0;
    private double indelSumDeltaSquared = 0.0;

    public void onTraversalStart() {
        for (ConcordanceState state : ConcordanceState.values()) {
            this.snpCounts.put(state, new MutableLong(0));
            this.indelCounts.put(state, new MutableLong(0));
        }
    }

    protected void apply(AbstractConcordanceWalker.TruthVersusEval truthVersusEval, ReadsContext readsContext, ReferenceContext refContext) {
        ConcordanceState concordanceState = truthVersusEval.getConcordance();
        if (truthVersusEval.getTruthIfPresentElseEval().isSNP()) {
            this.snpCounts.get(concordanceState).increment();
        } else {
            this.indelCounts.get(concordanceState).increment();
        }
        switch (concordanceState) {
            case TRUE_POSITIVE: {
                this.infoDifference(truthVersusEval.getEval(), truthVersusEval.getTruth());
                break;
            }
            case FALSE_POSITIVE:
            case FALSE_NEGATIVE:
            case FILTERED_TRUE_NEGATIVE:
            case FILTERED_FALSE_NEGATIVE: {
                break;
            }
            default: {
                throw new IllegalStateException("Unexpected ConcordanceState: " + concordanceState.toString());
            }
        }
    }

    private void infoDifference(VariantContext eval, VariantContext truth) {
        double evalVal = Double.valueOf((String)eval.getAttribute(this.evalInfoKey));
        double truthVal = Double.valueOf((String)truth.getAttribute(this.truthInfoKey));
        double delta = evalVal - truthVal;
        double deltaSquared = delta * delta;
        if (eval.isSNP()) {
            this.snpSumDelta += Math.sqrt(deltaSquared);
            this.snpSumDeltaSquared += deltaSquared;
        } else if (eval.isIndel()) {
            this.indelSumDelta += Math.sqrt(deltaSquared);
            this.indelSumDeltaSquared += deltaSquared;
        }
        if (Math.abs(delta) > this.epsilon) {
            this.logger.warn(String.format("Difference (%f) greater than epsilon (%f) at %s:%d %s:", delta, this.epsilon, eval.getContig(), eval.getStart(), eval.getAlleles().toString()));
            this.logger.warn(String.format("\t\tTruth info: " + truth.getAttributes().toString()));
            this.logger.warn(String.format("\t\t Eval info: " + eval.getAttributes().toString()));
        }
    }

    public Object onTraversalSuccess() {
        double snpN = this.snpCounts.get(ConcordanceState.TRUE_POSITIVE).doubleValue();
        double snpMean = this.snpSumDelta / snpN;
        double snpVariance = (this.snpSumDeltaSquared - this.snpSumDelta * this.snpSumDelta / snpN) / snpN;
        double snpStd = Math.sqrt(snpVariance);
        double indelN = this.indelCounts.get(ConcordanceState.TRUE_POSITIVE).doubleValue();
        double indelMean = this.indelSumDelta / indelN;
        double indelVariance = (this.indelSumDeltaSquared - this.indelSumDelta * this.indelSumDelta / indelN) / indelN;
        double indelStd = Math.sqrt(indelVariance);

        this.logger.info(String.format("SNP average delta %f and standard deviation: %f", snpMean, snpStd));
        this.logger.info(String.format("INDEL average delta %f and standard deviation: %f", indelMean, indelStd));

        try {
            InfoConcordanceRecord.InfoConcordanceWriter concordanceWriter = InfoConcordanceRecord.getWriter(Paths.get(this.summary));
            concordanceWriter.writeRecord(new InfoConcordanceRecord(VariantContext.Type.SNP, this.evalInfoKey, this.truthInfoKey, snpMean, snpStd));
            concordanceWriter.writeRecord(new InfoConcordanceRecord(VariantContext.Type.INDEL, this.evalInfoKey, this.truthInfoKey, indelMean, indelStd));
            concordanceWriter.close();
        }
        catch (IOException e) {
            throw new UserException("Encountered an IO exception writing the concordance summary table", e);
        }

        return "SUCCESS";
    }

    protected boolean areVariantsAtSameLocusConcordant(VariantContext truth, VariantContext eval) {
        boolean sameRefAllele = truth.getReference().equals(eval.getReference());
        boolean containsAltAllele = eval.getAlternateAlleles().contains(truth.getAlternateAllele(0));
        return sameRefAllele && containsAltAllele;
    }

    protected Predicate<VariantContext> makeTruthVariantFilter() {
        return vc -> !vc.isFiltered() && !vc.isSymbolicOrSV();
    }

}
