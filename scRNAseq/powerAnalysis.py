# Zeal Jinwala
# May 30, 2022

# Problem Statement:
# Given a scRNAseq dataset: 
#   Assume that the total count for every sample is 5000 (i.e., sum of column = 5000).
#   Imagine there was a row (gene G) in this dataset for which the count is expected to be 1 in 10% of samples and 0 in the remaining 90% of samples.  We are doing an experiment where we would like to know if the expression of gene G changes in experimental vs control conditions, and we will measure n samples (single cells) from each condition.
#   Plot the statistical power to detect a 10% increase in the expression of G in experimental vs control at Bonferroni-corrected p < 0.05 as a function of n, assuming that we will be performing a similar test for significance on 1000 genes total.  How many samples from each condition do we need to measure to achieve a power of 95%?
#   (Make the simplifying assumption that the counts for this gene follow a Poisson distribution. There are many potential statistical tests you may use; for simplicity feel free to use a t-test)

from statsmodels.stats.power import TTestIndPower
import math
import numpy as np
from matplotlib import pyplot

def get_num_samples(effect: float, alpha: float, power: float) -> int:
    '''
    this function solves for any one parameter of the power of a two sample t-test
    :effect: expected difference between two groups
    :alpha: Bonferroni adjusted alpha value
    :power: desired power 
    '''
    analysis = TTestIndPower()
    results = analysis.solve_power(effect, power=power, nobs1=None, ratio=1.0, alpha=alpha)
    return math.floor(results)

def generate_power_plot() -> None:
    '''
    Plot the statistical power to detect a 10% increase in the expression of G in 
    experimental vs control at Bonferroni-corrected p < 0.05 as a function 
    of n (number of samples from each condition) 
    '''
    # range was decided based on number of samples obtained in the prev function (get_num_samples()) 
    sample_sizes = np.array(range(100, 8000, 100))
    effect_sizes = np.array([0.1])
    analysis = TTestIndPower()
    # dep_var: which var is used for x-axis
    # nobs: values of number of observations
    # effect_sizes: values of effect size in the plot
    analysis.plot_power(dep_var="nobs", nobs=sample_sizes, effect_size=effect_sizes)
    pyplot.savefig("question2.png")

if __name__ == "__main__":
    # 10% increase
    EFFECT = 0.1 
    # desired power
    POWER = 0.95 
    # Bonferroni adjusted p<0.05
    # assuming that we will be performing a similar test for significance on 1000 genes total
    ALPHA = 0.05/1000 
    num_samples = get_num_samples(effect=EFFECT, alpha=ALPHA, power=POWER)
    print("Num samples:", num_samples)
    print("Generating a plot now...")
    _ = generate_power_plot()
