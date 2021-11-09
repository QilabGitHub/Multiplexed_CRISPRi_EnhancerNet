# Multiplexed_CRISPRi_EnhancerNet
# 1. Calculate the depletion score for multiplexed CRISPRi screen.
Firstly, to obtain sufficient representative reads for the downstream calculation, counts of sgRNA pairs were first filtered by requiring all paired constructors with at least 30 reads in D0, then added with a pseudocount of 10. The log2 enrichment scores of all sgRNA pairs were then calculated in D30 relative to D0 by using the relative frequencies in D0 and D30. The mean and the standard deviation of enrichment scores for paired non-targeting (neg-neg) sgRNAs were used to normalize the log2 enrichment scores. The normalized log2 enrichment scores were represented as the depletion score (see below).
![This is an image](https://myoctocat.com/assets/images/base-octocat.svg)
# 2. Calculate the epistasis score for sgRNA pairs targeting enhancers
## The defination of epistasis types with interaction score
## An additive model to calculate the interaction score to measure epistasis 
