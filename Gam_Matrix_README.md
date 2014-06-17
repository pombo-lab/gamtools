Description of GAM Matrices:


When we do a GAM experiment, we cut the nucleus into very thin slices and then 
ask how often two genomic bins are present in the same slice. For two loci A 
and B, the result of a GAM experiment is therefore a 2x2 table:

          | A Absent | A Present
----------|----------|-----------
B Absent  | 240      | 20
B Present | 13       | 50

In a perfect experiment, the number of times A and B are seen together should 
be inversely proportional to their distance (since the further apart they are, 
the less likely they are to appear in the same thin slice). Conversely, the 
number of times they appear together should correlate positively with Hi-C 
ligation frequency, since the closer they are the more frequently they ought to 
be ligated. In the case of the above example, A and B are together 50 times out 
of 323 experiments.

However, we cannot simply use the number of times both are present, since we 
may be better at detecting the presence of one locus than the other. In the 
above example, we detect locus A in 70 tubes, but we detect locus B in only 63 
tubes. Therefore, we may fail to detect the presence of locus B in 10% of the 
"true" cases, and the true number of times that A and B are present together 
may be higher than what we observe.

To normalise for differences in the *marginal* probability of detection for A 
or B, we use the normalised linkage disequilibrium or D' (sometimes written as 
Dprime). For a detailed description, see 
http://en.wikipedia.org/wiki/Linkage_disequilibrium

In short, we first calculate the linkage disequilibrium. We calculate the 
probability of seeing A and B together if they were distributed independently:

Prob seeing A = 70 / 323
              = 0.216

Prob seeing B = 63 / 323
              = 0.195

Prob seeing AB at random = 0.216 * 0.195
                         = 0.0421

Then we simply subtract the random expected probability from the observed 
probability to get the linkage:

Observed AB = 50 / 323
            = 0.155

Linkage = 0.155 - 0.0421
        = 0.113

Note: If we observe A and B together *less* frequently than expected by chance, 
we get a negative linkage.

Now to correct for the marginal probabilities of A and B, we divide by the 
maximum possible linkage we could observe.

In the example above, we detect B 63 times, and thus the maximum possible 
observations of AB is 63 (i.e. when B is present A is always present). If this 
were the case, the linkage would be:

Observed AB = 63 / 323
            = 0.195

Linkage = 0.195 - 0.0421
        = 0.153

And the normalized linkage is the observed linkage divided by the maximum 
linkage:

Dprime = 0.113 / 0.153
       = 0.739


This means that a Dprime of 1.0 means A and B were always together, 0.0 means 
they were together exactly as often as expected if independent and -1.0 means 
they were never seen together.
