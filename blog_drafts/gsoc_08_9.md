# A few brief life updates

In a moment typically reserved for 'As seen on TV ads...', I knocked over an entire glass of water onto my [laptop](https://www.newegg.com/icicle-silver-asus-zenbook-um433da-nh74-mainstream/p/N82E16834235285?Item=N82E16834235285) a few days ago. While the laptop sports one of the new Ryzen 7 8-thread CPUs, the technology remains quite sensitive to water. A nearly-full 8 ounce glass of water came tumbling down on the keyboard. Like most times of brief crisis, I was stunned at what was happening before me. Although I tried to resuscitate the laptop with rice, vigorous shaking, and some prayers, it eventually went black after about an hour or so. Unbeknownst to me ASUS has a limited warranty, free replacement program for accidental damage. The process will take about two months, but I can't complain as I caused the flood! Luckily, all of my work files are in cloud storage so I can continue on my desktop. The laptop was the primary vehicle for nearly all of the code completed during GSOC 2020, so it feels a bit strange heading into the final leg of the race without your normal shoes. 

On a separate note, I recently received an email about upcoming talks at the Open Geospatial Consortium (OGC). The upcoming health-related GIS talks are on the intersection of ethics and GIS, hosted by Ajay K. Gupta and Dr. Ram Peruvemba (Health Solutions Research, Inc.), as well as Ed Parsons (Google). These topics are of immense interest to me and I hope to integrate them into my future work. Preview of the talks by Gupta and Peruvemba [here](https://portal.ogc.org/files/?artifact_id=94352), and Parsons [here](https://www.edparsons.com/2020/06/the-ethics-of-geospatial-the-four-es/). 

# GSOC updates

As stated on a previous blog post, the majority of the GSOC work is complete and the following objectives are additional tasks I hope to contribute to PySAL. 

## Updating LOSH inference to new the `__crand()` engine

This task has proven a bit tricky. On a GSOC ESDA call with Levi and Serge, they highlighted how the actual known behavior of conditional randomization as it relates to the LOSH equation is still unknown. So we are in a bit of uncharted territory. We have decided to try and implement the conditional randomization in pure python before attempting to migrate it to `numba`. 

The key part of the LOSH equation is getting the sum of the lagged neighbors. We are hoping that the following snippet is an effective way to do so: 

```
(zrand - (zrand @ weights_i)) @ weights_i)
```

where `zrand` is a subset of observations around an area `i`, and `zi` is the held-out `i` value. 

While I am still working on this, it has been interesting to compare to the conditional randomization procedure available in the R `spdep` package. Specifically, the authors use the bootstrap procedure proposed by [Xu et al (2014)](https://link.springer.com/article/10.1007%2Fs00168-014-0605-5), or a slightly modified version of it. I wonder if they already tried to do a conditional randomization procedure as we are doing and determined it to be unsuitable? Only time will tell!

## Draft univariate and multivariate local Geary estimator

Another set of estimators that were identified as relevant to PySAL are the univariate and multivariate local Geary estimator. These are incredibly interesting statistics that I largely misunderstood on my first read. There is a subtle 'switch' (to avoid using the word transformation) that happens between the univariate and multivariate version of the statistic. 

In the univariate case, we are interested in:

$c_i = \sum_j w_{ij} (x_i - x_j)^2$

However, when we extend the equation to the multivariate case,

$c_{k,i} = \sum_{v=1}^k c_{v,i}$

and Anselin notes that

>The expression in equation (10) [or the $c_{k,i}$ statistic above] can be generalized in a number of ways. As mentioned, other distance measures can be applied, such as a Manhattan distance (absolute differences), or, in general, any Minkowski distance metric. [*Anselin page 9, 2017*](https://geodacenter.github.io/docs/LA_multivariateGeary1.pdf)

Whereas the univariate implementation can, more or less, but fixed in a single form, the multivariate implementation must be implemented in such a way that is flexible to different types of distance metrics. My current implementation of the multivariate local Geary can correctly local Geary values based on a squared differences, although I have not yet figured out inference.

The inference for the multivariate local Geary statistic is challenging for a few reasons. Most of the templates for conditional randomization in PySAL are univariate or bivariate at most. Indeed, the new `_crand()` engine currently *blocks* multivariate inputs. The multivariate local join counts estimator was successful because I could simplify the inputs into a single vector (due to the nature of the statistic), and then pass it through existing conditional randomization engines. However, the multivariate local Geary requires that

>A conditional permutation approach consists of holding the *tuple* of values observed at *i* fixed, and computing the statistic for *m* permutations of the remaining tuples over the other locations. [*Anselin page 10, 2017*](https://geodacenter.github.io/docs/LA_multivariateGeary1.pdf)

I'm currently experimenting with randomization of zipped lists, although I'm not actually sure this will work! I hope to have some answers by the end of August. On a more positive note, the univariate local Geary is [fully operational with inference](https://github.com/jeffcsauer/GSOC2020/blob/master/functions/local_geary.py), and I am in the process of writing up a doc for the estimator.

## Experiment with the LOSH statistic as it relates to [Bivand and Wong's 2018 paper](https://link.springer.com/article/10.1007/s11749-018-0599-x)

This is the most ambitious of the additional tasks and is still underway. Most of the work is happening behind-the-scenes at the moment, but I will post general updates as they come! 