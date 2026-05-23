#pragma once

#include <math.h>
#include <stddef.h>

/**
 * Values ​​of the Student's t-distribution for a two-tailed 95% CI for df = 1..30.
 * Values ​​correspond to the 0.975 quantile (α/2 = 0.025 in each tail).
 */
static const double T_TABLE[] = {
    12.706, 4.303, 3.182, 2.776, 2.571, // df 1-5
    2.447,  2.365, 2.306, 2.262, 2.228, // df 6-10
    2.201,  2.179, 2.160, 2.145, 2.131, // df 11-15
    2.120,  2.110, 2.101, 2.093, 2.086, // df 16-20
    2.080,  2.074, 2.069, 2.064, 2.060, // df 21-25
    2.056,  2.052, 2.048, 2.045, 2.042, // df 26-30
};

/**
 * @brief Returns the value of the t-distribution for the 0.975 level (95% CI).
 *
 * Implementation:
 * - For df <= 30: exact values ​​from the precalculated table.
 * - For df > 30: Cornish–Fisher expansion, refining the normal distribution.
 *
 * @param df Degrees of freedom, typically n - 1.
 * @return The t-quantile value.
 */
static inline double t_ppf_975(int df) {
    if (df <= 0) {
        return 0.0;
    }
    if (df <= 30) {
        return T_TABLE[df - 1];
    }
    // Cornish-Fisher approximation for large df
    double z = 1.959964; // normal quantile for 0.975
    return z + (z * z * z + z) / (4.0 * df) + (5 * z * z * z * z * z + 16 * z * z * z + 3 * z) / (96.0 * df * df);
}

/**
 * Result of statistics calculation.
 */
typedef struct {
    double mean;  // Sample mean
    double error; // 95% CI half-width
} StatResult;

/**
 * @brief Calculates the mean and half-width of the 95% confidence interval.
 *
 * Algorithm:
 * 1. Calculate the arithmetic mean.
 * 2. Calculate the unbiased standard deviation (divided by n-1).
 * 3. Standard error of the mean (SEM) = std_dev / sqrt(n).
 * 4. Half-width of the CI = t_critical * SEM.
 *
 * @param data Array with sample data.
 * @param n Number of elements in the sample.
 * @return Structure with mean and error.
 */
static inline StatResult calc_stats(const double *data, size_t n) {
    if (n == 0) {
        return (StatResult){0.0, 0.0};
    }

    double sum = 0.0;
    for (size_t i = 0; i < n; i++) {
        sum += data[i];
    }
    double mean = sum / (double)n;

    if (n == 1) {
        return (StatResult){mean, 0.0};
    }

    double sq_sum = 0.0;
    for (size_t i = 0; i < n; i++) {
        double d = data[i] - mean;
        sq_sum += d * d;
    }
    double std_dev = sqrt(sq_sum / (double)(n - 1)); // sample std
    double sem = std_dev / sqrt((double)n);          // standard error
    double ci = t_ppf_975((int)(n - 1)) * sem;

    return (StatResult){mean, ci};
}
