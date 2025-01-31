# DDRi_RNAseq
RNAseq analysis repository for analysis
# Experimental Design

## Samples and Conditions
Each sample has two key attributes:
- Gene Target (e.g. GAL4, PRDX1, NTC, etc.)
- DNAPKi Treatment (TRUE/FALSE)

We combine these into a single factor called condition. For example:
- GAL4_FALSE = "GAL4 with no DNAPKi treatment"
- GAL4_TRUE = "GAL4 with DNAPKi treatment"
- PRDX1_FALSE = "PRDX1 with no DNAPKi treatment"
- PRDX1_TRUE = "PRDX1 with DNAPKi treatment"
- NTC_FALSE = "Negative control (no DNAPKi)"

## Reference Level
We make NTC_FALSE the reference (baseline) level. This ensures that any coefficient for "conditionXYZ" in the linear model is interpreted as:

(XYZ) - (NTC_FALSE)

## Design Matrix
We set up a linear model with the formula `~ replicate + condition`, where:
- replicate accounts for batch or biological replicate effects
- condition captures all combinations of gene target and DNAPKi status

Concretely, our design matrix columns include:
- Intercept (which corresponds to NTC_FALSE)
- Replicate terms
- conditionGAL4_FALSE, conditionGAL4_TRUE, conditionPRDX1_FALSE, conditionPRDX1_TRUE, etc.

## Single Comparisons
For any given conditionX, we do a single comparison X vs NTC_FALSE by extracting that specific coefficient in the model. This tells us whether X is significantly different from the negative control (no DNAPKi).

## Double Comparisons (DNAPKi Effects)
We define contrasts to compare (GeneTarget_TRUE - GeneTarget_FALSE).

For example:
- GAL4_DNAPKi_effect = conditionGAL4_TRUE - conditionGAL4_FALSE
- PRDX1_DNAPKi_effect = conditionPRDX1_TRUE - conditionPRDX1_FALSE

Each such contrast tests whether DNAPKi treatment (for a specific gene target) significantly changes expression compared to the same gene target without DNAPKi.

## Interpretation

### Single Comparison vs. NTC_FALSE:
- logFC > 0 ⇒ Genes are up-regulated relative to the negative control
- logFC < 0 ⇒ Genes are down-regulated relative to the negative control

### Double Comparison (DNAPKi effect):
- logFC > 0 ⇒ Genes are further increased by DNAPKi treatment (for that gene target)
- logFC < 0 ⇒ Genes are decreased by DNAPKi for that gene target
