# Introduction

This has been quite an interesting week in GSoC. I've made several improvements to the `losh()` function and have (fingers-crossed) figured out inference for the local join counts! Details below. 

## `losh()` improvements 

As mentioned on the previous blog post, Levi had raised several [issues-as-feedback](https://github.com/jeffcsauer/GSOC2020/issues) on my GSOC20202 workbook repository. These included:

- Using warnings rather than print statements
- Simplifying the $VarH_i$ calculation
- Switching to sparse weights to calculate row sums
- Breaking notebooks into .py function files

I have been able to address each of these issues through various pull requests ([identified in each now-closed Github issue](https://github.com/jeffcsauer/GSOC2020/issues?q=is%3Aissue+is%3Aclosed)). The most challenging of these fixes was the simplification of the $VarH_i$ calculation. The brief recap of this issue is that the $VarH_i$ calculation is not 'easy-to-read' code at first glance. Thus, we were looking for ways to break up the function such that it was a bit easier to digest (at least visually). One part of the calculation involves squaring the sum of the row weights. After some back-and-forth of trying to implement Levi's original solution, Levi realized that we were getting caught up between matrix multiplication versus exponentiation. We finally arrived at Levi's elegant solution...

```
w_squared = numpy.asarray(w.sparse.multiply(w.sparse).sum(axis=1)).flatten()
```

...which allows us to flexibly calculation the sum of squared row weights. The key advantage of this solution is that it avoids using list comprehension.  

I also appreciated the practice of [using the `warnings` module](https://github.com/jeffcsauer/GSOC2020/issues/4). This module is used throughout Python software development as is more regimented in the handling of invalid inputs. In my first solution, I was carrying out a basic string comparison to see if the user was inputting something like `chi-square` or `permutation`. If this user passed through something like `garbage`, the function would still run with a printed warning. However, this printed warning would not actually stop the function! Thus, with the new use of the `warnings` module I can now use `warnings.warn()` or `raise ValueError()` to either 1) offer a more recognizable warning or 2) block the function from running. For now, I am issuing a warning, although in the future this will likely need to change to `raise ValueError()` as it would eliminate the possibility of users accidentally accepting the p-values with an incorrect `a` input. 

## Inference for the local join count functions

Given my progress on setting up the `losh` and `local_join_count` functions, I have been feeling relatively positive about my progress in GSoC so far. The PySAL community has been incredibly supportive, and my mentors have pushed by my coding abilities and intellectual thought on the topics of spatial statistics. However, I the past week or so I was dreading one of the last major coding hurdles I would have to face: conditional permutation and inference.

I am not sure why this topic is so intimidating to me. Perhaps it is simply the sound of the words: *conditional permutation*. Or, more realistically, the fear likely comes from a similar place of other fears: misunderstanding. When I first learned about inference via permutation, I was in GEOG202: Introduction to Spatial Statistics at McGill University. We learned about the process in a simple point-and-click method in the GeoDa software. My youthful inexperience could not grasp what was going on as I repeatedly clicked 'Generate' and saw the distribution and accompanying p-values shift ever so slightly. 

Now, returning several years later to the topic and faced with its implementation I was of course hesitant. To begin my understanding of the method, I turned to Anselin and Li (2019) explanation: 

>A conditional permutation test as proposed in Anselin (1995) to compute a pseudo
p value for the LISA statistics can be constructed in the usual way. The general
principle, for those locations i where xi = 1 , is to carry out a series of random permutations of the remaining observations, while counting the times the number of
neighbors with xj = 1 equals or exceeds qi , the observed value of the join counts. In
practice, this is implemented by taking ki (the number of neighbors for i) draws without replacement from a set of N − 1 observations with K − 1 values of 1 for those
observations where xi = 1 . A pseudo p value can be computed as (v + 1) ∕ (r + 1) ,
where v is the number of times the neighbors have qi or more values equal 1, and r is
the number of permutations. The standard caveats apply (e.g., sensitivity to the number
of permutations, varying results depending on the random number sequence,
multiple comparisons, etc.).
It should be noted that, as a one-sided test, the conditional permutation approach
includes instances as rejecting the null of spatial randomness where there are more
than qi neighbors with xj = 1 in the computation of the pseudo p value. **Anselin and Li, 2019, p.7**

Breaking the above paragraph down, we recognize that we need to accomplish several steps:
1. Randomly permute values of 0 and 1 for the length of the input data (excluding $x_i$)
2. Count the nubmer of times the local join counts of the permuted values = 1 or exceeds the actual observed local join counts
3. Based on the number generated in Step 2, calculate a pseudo p-value as $(v+1) / (r+ 1)$

Now, to my benefit, conditional permutation has already been implemented in the `esda` `Moran_Local()` function. I was able to lift almost all of the code from a sub-function of `Moran_Local()` called `_crand()`. Because of the differences in how randomization is carried out across architectures and software, it difficult to check if these randomized p-values will be exactly the same. For now, I am doing a surface level comparison to the output from GeoDa to see if the values are within a certain range. I'm also getting some strange values of 0.001 where GeoDa reports `NaN`. I have scheduled an additional call with Levi tomorrow (6/22) where we will talk about these issues and carry out some code review. I hope to also walk through the origin of the `_crand()` sub-function as there is quite a bit going on in only a few line of code! Next week I will try to provide a line-by-line breakdown of what is going on. 

## GSOC To-Do list

- Implement inference for bivariate and multivariate local join counts
- Begin a summary notebook that demonstrates all of the functions for the next GSoC-wide call (July 3rd)
- Start thinking about how and what tests will be necessary for the functions
- PR for citations present in the new function docstrings to main pysal repository 
- Ask Pablo/panel team aobut these models and how they might be helpful for my own work!

