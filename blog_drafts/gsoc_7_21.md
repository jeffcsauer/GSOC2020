# Introduction

Two weeks are left in the Phase II of GSOC! It is a bit amazing - and terrifying - when one thinks of how fast this summer passed. I am incredibly grateful for the opportunity provided by GSOC and I hope to continue the summer productivity into the Fall semester. But before then, there remains plenty to do in GSOC! This blog post will act as a primary recap for what has happened over the course of Phase II. 

## Summarizing the work in a PR

I recently opened [PR#139](https://github.com/pysal/esda/pull/139) on the PySAL ESDA github repository. This PR adds the four new estimator functions as well as two demonstration notebooks. Each function is written in the form of a scikit-learn style estimator. PEP8 formatting has been applied to all of the functions, although a handful of lines are left long due to readability. Regarding the notebooks, each notebook briefly reviews the core math of the statistic and explains how a user might deploy the statistic in a geospatial workflow. The included examples match those found in [`R` `spdep::LOSH` and `spdep::LOSH.cs`](https://www.google.com/search?client=firefox-b-1-d&q=losh.cs) as well as the [GeoDa Local Join Counts tutorial](https://www.google.com/search?client=firefox-b-1-d&lei=ytAVX7L4EYKztAb23K9A&q=local%20indicators%20of%20spatial%20association&ved=2ahUKEwiyhYuJqtzqAhWCGc0KHXbuCwgQsKwBKAB6BAgPEAE&biw=1920&bih=938).  

I also wrote some straightforward .py test files that execute successfully [on my fork](https://github.com/jeffcsauer/esda). These new tests follow the lead of the existing tests in the PySAL testing framework, specifically `test_moran.py` for the LOSH function and `test_join_counts.py` for the local join count functions.

Now, there are some slight 'issues' (or considerations) that I highlight in the PR. I designate these as slight because the functionality and documentation of the estimators are working, yet there are some considerations that might delay the PR. Specifically:

- Local join count functions are currently using the default `n_jobs=1` to avoid issues fixed in potential [merges](#diff-1fd239570b4623d6ce23f0e268495633L510). Users may enter values other than 1. 
- One doctest dataset, `commpop`, appears to have a [dead link](https://geodacenter.github.io/data-and-lab/) at the moment. Having trouble loading this as a `libpysal.example`. Once this dataset is working the doctest can be updated. Currently hosting it on a personal Github link.

## Challenge: `numba` and multicore processes

During a recent GSOC-wide call, Dani pointed out that I had hard-coded the number of jobs (`n_jobs`) equal to 1. This `n_jobs` setting is similar to the `numba` [number of threads](https://numba.pydata.org/numba-doc/dev/reference/envvars.html#envvar-NUMBA_NUM_THREADS) argument. Essentially, it tells the `_crand()` engine how many computer cores to use.  

When using `n_jobs=1` the functions had increased in speed anywhere from 40-60% on average. One would presume that increasing the number of cores would directly improve this speed. Yet when `n_jobs` was set to `-1` (use all cores) or a specific number `>1`, I was experiencing quite perplexing errors regarding the size of input arrays. The specific error messages can be seen [here](https://github.com/jeffcsauer/GSOC2020/blob/master/validation/Validation_LJC_univariate.ipynb) under the 'Next steps' header. Levi explained that 'the issue (I think) is when unequally-sized jobs are sent around [to several cores], `hstack` can fail. Not sure why/how. I think using `row_stack` solved the similar issue I was seeing in the end.' What is quite interesting about the switch of these functions is `row_stack` has the same requirement that the input arrays must be of the same dimensions (see [Code #2 example here](https://www.geeksforgeeks.org/numpy-ma-row_stack-in-python/), yet `row_stack` does not create a `numba` error. This issue is fixed in a draft PR by [Levi](https://github.com/pysal/esda/pull/133/commits/68f623c268e29e0ad24236fcc10d88536f9df3bc#diff-1fd239570b4623d6ce23f0e268495633L510), but for the time being setting `n_jobs=1` will provide both increased speed and stable functionality to users. Once this issue is fixed users may set `n_jobs` equal to any number appropriate for their machine.

It is worth noting that there remains the possibility that using more cores may not improve speed. There are several reasons as to why this may not be so. **I am not a computer scientist so I cannot provide detailed reasoning - the following are a summary of what I have found via some quick online searches.**

1. A lay summary of reasons is available [here](https://www.forbes.com/2009/11/23/google-microsoft-programming-technology-cio-network-multicore-hardware.html#63abb07e6914). Not everything in the article applies to the current issue (i.e. `numba` parallelizes the operations), but considerations like lower clock speeds in individual cores, lack of (shareable) memory, and locked resources might. When I forced the operations to use multiple cores I actually saw a slowdown - this could be time wasted splitting up the job and bringing it back together.  
2. Another potential issue is that of thread switching, visualized [here](http://www.dabeaz.com/GIL/gilvis/fourthread.html). The linked visualizations are prepared by David M. Beazley and are meant to '...[O]bserve how thread switching gets more rapid as the number of CPUs increases. Zoom in to see all sorts of crazy GIL contention.' Given that each thread switch can lead to what is called 'computational overhead', we are wary of the dramatic increases in thread switching as the number of cores increases. Now this may not be the case with the `numba`-ized functions, but it could be! 
3. More discussion available on this [S/O post](https://stackoverflow.com/questions/9183476/how-to-do-the-same-calculations-faster-on-4-core-cpu-4-threads-or-50-threads?noredirect=1&lq=1), specifically on the topic of balancing cores versus threads.

## Going further 

Over the next couple of weeks I will have a meeting with Levi and Serge to discuss the PR and what changes might be made. There is certainly room for improvement on the demonstration `DOC` notebooks, and there may be some cleanup required in the functions. 

I've also started working on the [Local Gamma statistic](https://github.com/jeffcsauer/GSOC2020/blob/master/scratch/Issues_Scratch.ipynb) in order to practice open-source contributing outside the scope of GSOC. 

# GSOCTo Do

- Create some experimental notebooks testing the functions across small and large datasets
