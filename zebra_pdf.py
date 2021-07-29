# Author: Daniel Hakim
# July 28 2021
#  We estimate confidence intervals for expected
#  genome coverage given
#       -genome length,
#       -read length,
#       -and number of reads
#  Estimates are created by random trial,
#  simplified random trial, and
#  closed form algebraic solution

import random
from scipy.stats import binom
from math import ceil, sqrt
import copy


if __name__ == "__main__":
    # Modify to plot coverage for different scenarios
    GENOME_LENGTH = 3878815
    READ_LENGTH = 149
    NUM_READS = 173
    # Modify number of trials done in empirical simulations
    NUM_TRIALS = 1000


def compute_probability_coverage_greater_than(coverage, genome_length, mean_read_length, num_reads):
    mean_read_length = ceil(mean_read_length)

    # Builds the pdf with the dynamic programming algorithm
    # then tells the fraction of the pdf to the right of the specified coverage
    full_pdf = dyn_prog_algebra(genome_length, mean_read_length, num_reads)
    num_buckets = ceil(genome_length / mean_read_length)

    # Because of the bucketing assumption,
    # the fact mean read length might not be an integer
    # because mean read length might not evenly divide genome length
    # and just due to floating point errors,
    # we can hit a situation where coverage is
    # sitting in a bucket that takes up
    # 99.99% of the pdf, and accidentally miss it.
    # This is why it would be really nice to have a
    # closed form formula for the pdf
    # We just add a small fudge factor to coverage
    # threshold that we report p value for to account for it.
    coverage = coverage + 0.00001

    p_left = 0
    for k in range(num_buckets+1):
        if k / num_buckets <= coverage:
            p_left += full_pdf[k]
        else:
            if genome_length == 5452778:
                print(full_pdf)
                print(num_buckets)
                print(k / num_buckets)
                print(coverage)
                print((k-1) / num_buckets)
            break
    return p_left


def random_trial(num_trials, genome_length, read_length, num_reads):
    # In this trial, we assume reads start at a uniformly
    # random position in the genome
    # and mimic the exact actions of the zebra filter.
    trial_results = []

    for trial in range(num_trials):
        coverset = set()
        for read in range(num_reads):
            start_index = random.randint(0, genome_length-read_length)
            for i in range(start_index, start_index+read_length):
                coverset.add(i)
        coverage = len(coverset)/genome_length
        trial_results.append(coverage)

    return trial_results


def simplified_trial(num_trials, genome_length, read_length, num_reads):
    # In this trial, we simplify the zebra filter by assuming coverage
    # is calculated against buckets, each bucket of the genome is
    # the same size as read length, and reads are uniformly distributed across
    # buckets.
    # While not exactly the same as random trial above,
    # this will give a close approximation for
    # the zebra filter's coverage outputs,
    # and is easily solved algebraically.
    num_buckets = ceil(genome_length / read_length)
    trial_results = []

    for trial in range(num_trials):
        coverset = set()
        for read in range(num_reads):
            bucket = random.randint(0, num_buckets - 1)
            coverset.add(bucket)
        coverage = len(coverset)/num_buckets
        trial_results.append(coverage)

    return trial_results


# This exactly computes the pdf for the simplified
# trial with dynamic programming
def dyn_prog_algebra(genome_length, read_length, num_reads):
    num_buckets = ceil(genome_length / read_length)
    # p_arr[i] represents the probability that i unique buckets are picked
    # by the reads so far.
    p_arr = [1]
    for r in range(num_buckets):
        p_arr.append(0)

    # Our scratch buffer swaps back and forth between two arrays
    p_arr_swap = copy.copy(p_arr)

    for r in range(num_reads):
        for i in range(len(p_arr)):
            # Each read must either go into a new bucket or a used bucket
            if i > 0:
                p_new = p_arr[i-1] * ((num_buckets - (i-1)) / num_buckets)
            else:
                p_new = 0
            p_old = p_arr[i] * (i / num_buckets)
            p_arr_swap[i] = p_old + p_new

        p_temp = p_arr
        p_arr = p_arr_swap
        p_arr_swap = p_temp

    return p_arr


# This makes some false assumptions to create a
# binomial distribution that is very close to the pdf
# but is a little bit too wide.  (Mean coverage is correct)
def simplified_algebra_wrong(genome_length, read_length, num_reads):
    # This wrongly assumes a read can be in multiple buckets
    # This model gets the right mean, but slightly overestimates
    # the bounds on the confidence interval.  Its very easy to calculate.
    num_buckets = ceil(genome_length / read_length)
    # Odds of a read being in a bucket is 1/num_buckets.
    # Odds of a read not being in a bucket is (num_buckets-1)/num_buckets
    # Odds of a bucket not being covered by any of n reads is:
    #   ((num_buckets-1)/num_buckets)**n

    # Expected coverage is the fraction of buckets
    # that are covered by at least 1 read
    # coverage = 1 - ((num_buckets-1)/num_buckets)**n
    coverage = 1 - ((num_buckets-1)/num_buckets)**num_reads

    # We (falsely) assume all buckets have the same
    # independent chance of being covered,
    # we can calculate the confidence interval for coverage
    # using binomial distribution where num_buckets
    # is the number of trials and coverage is the
    # weighting of the coin
    n, p = num_buckets, coverage
    mean, var, skew, kurt = binom.stats(n, p, moments='mvsk')

    # Now we build a confidence interval.  scipy can only do two sided,
    # but we can upgrade this with a one sided one later if we so choose.
    # we just have to do it ourselves.
    min_bound, max_bound = binom.interval(0.90, n, p)
    return mean/num_buckets, min_bound/num_buckets, max_bound/num_buckets


if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt

    num_buckets = ceil(GENOME_LENGTH / READ_LENGTH)

    t1 = random_trial(NUM_TRIALS, GENOME_LENGTH, READ_LENGTH, NUM_READS)
    t2 = simplified_trial(NUM_TRIALS, GENOME_LENGTH, READ_LENGTH, NUM_READS)
    full_pdf = dyn_prog_algebra(GENOME_LENGTH, READ_LENGTH, NUM_READS)
    mean_coverage, min_bound_binom, max_bound_binom = \
        simplified_algebra_wrong(GENOME_LENGTH, READ_LENGTH, NUM_READS)

    print("Algebraic predictions")
    print("Mean coverage:")
    print(mean_coverage)
    print("90% confidence interval for coverage (binom approximation)")
    print([min_bound_binom, max_bound_binom])
    print("Full PDF")
    print(full_pdf)

    total_p = 0
    min_bound_95_exact = 0
    for i in range(len(full_pdf)):
        if total_p + full_pdf[i] > 0.05:
            min_bound_95_exact = i/num_buckets
            break
        total_p += full_pdf[i]

    print("Total better be 1: ", sum(full_pdf))

    print("Exact min_bound:" + str(1-total_p))
    print(min_bound_95_exact)

    print("Empirical Means")
    print("Randomized Trial")
    print(np.mean(t1))
    print("Simplified Randomized Trial")
    print(np.mean(t2))

    in_confidence_t1 = [x for x in t1 if min_bound_binom <= x <= max_bound_binom]
    in_confidence_t2 = [x for x in t2 if min_bound_binom <= x <= max_bound_binom]

    in_confidence_t1_exact = [x for x in t1 if min_bound_95_exact <= x]
    in_confidence_t2_exact = [x for x in t2 if min_bound_95_exact <= x]

    print("Empirical %Trials in Confidence Interval:")
    print("Randomized Trial")
    print(len(in_confidence_t1) / len(t1))
    print(len(in_confidence_t1_exact) / len(t1))
    print("Simplified Randomized Trial")
    print(len(in_confidence_t2) / len(t2))
    print(len(in_confidence_t1_exact) / len(t1))


    plt.hist(t1, bins=20)
    plt.title("Randomized Trial")
    plt.xlabel("Coverage")
    plt.ylabel("Trials")
    plt.axvline(min_bound_95_exact, label="p=.05", color="black", linestyle="-")
    plt.show()
    plt.hist(t2, bins=20)
    plt.title("Randomized Trial- Simplified Binning")
    plt.xlabel("Coverage")
    plt.ylabel("Trials")
    plt.axvline(min_bound_95_exact, label="p=.05", color="black", linestyle="-")
    plt.show()

    left = 0
    right = 0
    for i in range(len(full_pdf)):
        if full_pdf[i] > 0.0001:
            left = i
            break
    for i in reversed(range(len(full_pdf))):
        if full_pdf[i] > 0.0001:
            right = i
            break

    x = range(left, right)
    print([b/num_buckets for b in x])
    plt.bar([b/num_buckets for b in x], full_pdf[left:right], width=0.8/num_buckets, color='red')
    plt.title("Simplified Binning Full PDF")
    plt.xlabel("Coverage")
    plt.ylabel("Probability")
    plt.axvline(min_bound_95_exact, label="p=.05", color="black", linestyle="-")
    plt.show()
