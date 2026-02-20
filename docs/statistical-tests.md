# Statistical Tests

To properly examine the relationships between the data collected, the following is a rough description of the logic behind the analysis.

1. Assess if data is from a normal distribution:
   1. [Kologorov-Smirnov](https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test) (KS) will assess if the data comes from a referbce distribution (normal)
   2. If this is invalidated, then parametric tests are no longer appilcable
2. The test for unequal variance must be conducted as well. Due to the KS test failing to show normality, a Flinger-Klinger (FK) test was used to examine variances between groups
   1. The FK test is a nonparametric test to show if groups contain equal variance
   2. Failing this tests reinforces the use of non-paramateric tests
3. Assess if data across groups comes from the same distrubtion, regardless of shape
   1. [Kruskal-Wallis](https://en.wikipedia.org/wiki/Kruskal%E2%80%93Wallis_one-way_analysis_of_variance) (KW) test is a nonparametric test that will compare two or more independent samples (sample size does not matter) to see if they come from the same distribution. *A significant result will indicate taht at least one sample in the test stochastically dominates the other sample(s), though the dominant sample is not clear from the test*. A post-hoc test such as Dunn's or Conover-Imam highlights which sample is dominant.
4. If applicable, which sample(s) are statistically different from each other
   1. The [Conover-Imam](https://stats.stackexchange.com/questions/141856/what-is-the-difference-between-various-kruskal-wallis-post-hoc-tests) test provides a method to uncover specific significant relationships between groups of samples
   2. 