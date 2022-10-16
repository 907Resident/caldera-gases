# Statistical Tests

To properly examine the relationships between the data collected, the following is a rough description of the logic behind the analysis.

1. Assess if data is from a normal distribution:
   1. [Kologorov-Smirnov](https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test) (KS) will assess if the data comes from a referbce distribution (normal)
   2. If this is invalidated, then parametric tests are no longer appilcable
2. Assess if data across groups comes from the same distrubtion, regardless of shape
   1. [Kruskal-Wallis](https://en.wikipedia.org/wiki/Kruskal%E2%80%93Wallis_one-way_analysis_of_variance) (KW) 
3. Determine if variance is equal across groups