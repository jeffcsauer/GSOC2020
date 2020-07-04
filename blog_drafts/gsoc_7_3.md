# Introduction

What a whirlwind of a day to wrap up the week! This was a challenging week for GSoC. As I shared last week, I was feeling a bit down as I had encountered my first major hurdle in the project. This hurdle - inference for Case 1 of Bivariate Local Join Counts - was something that took my about a week of tinkering, reflecting, and trial-and-error to resolve. After an inspiring `PySAL` developer call I decided to spend the afternoon working on it and a solution finally emerged! The project now has completed all goals of Phase 1, namely:

- Reviewing literature on the spatial estimators
- Coding first attempts at each estimator (producing correct values, inference, and handling invalid inputs)
- Integrating mentor feedback

I'll share a bit more of the updates in the remainder of the blog post.

# Progress updates - a major push towards completion!

Using the table from last week I want to highlight some updates across all of the functions.

| Function              | Generating correct values | Generating correct inference | PEP8 styling | Docstrings | Tests      |
|-----------------------|----------------------------|-------------------------------|---------------|------------|------------|
| `LOSH`                | Yes                        | Yes                           | Yes (a few long lines)        | Yes     | Basic tests implemented |
| `Local_Join_Count`    | Yes                        | Yes                           | Yes (a few long lines)        | Yes     | Basic tests implemented |
| `Local_Join_Count_BV` | Yes*                       | Yes                           | Yes (a few long lines)        | Yes     | Basic tests implemented |
| `Local_Join_Count_MV` | Yes                        | Yes                           | Yes (a few long lines)        | Yes     | Basic tests implemented |

*As discussed in last week's post, there seems to be some interesting behavior with this function depending on input dataset. 

A summary of updates brought by this week include:
- Working inference for Bivariate Local Join Counts (Case 1) (!)
- PEP8 code formatting (with the exception of a few long supplementary lines)
- Formatted and edited docstrings
- Doctests for both a toy dataset and external reference. External reference uses the same datasets as the vignette for R `spdep::LOSH()` or the GeoDa online tutorials. The only exception is the doctest for `Local_Join_Count_MV()` as there is no multivariate example in the GeoDa tutorials, so I simply added a variable to the bivariate example. 

# Inference for `Local_Join_Count_BV` function

By a sheer stroke of luck and tinkering, we now have working inference for Bivariate Local Join Counts (Case 1). The trick was using the following operation in the _crand() function (function [here](https://github.com/jeffcsauer/GSOC2020/blob/master/functions/local_join_count_bv.py#L188)): 

```
...
tmp_z = z[idsi[rids[:, 0:wc[i]]]]
if case == "BJC":
    joins[i] = x[i] * (w[i] * tmp_z).sum(1)
...
```

This calculates the number of joins by multiplying a given `x[i]` value by `(w[i] * tmp_z).sum(1)` (sum of weights for `z_i`). When we compare the simulated p-values to those of GeoDa, we get a strong agreement and no strange values. This is a huge relief as this was a major obstacle that - for whatever reason - I could not wrap my head around. 

# Moving forward: getting to know `numba`!

Now that a basic inference method is working across all of the functions we can begin to take a look at Levi and Dani's `numba` endeavor. The TLDR of this initiative is that they are finding huge performance increases using the `numba` compiler, specifically increasing the efficiency of the conditional randomization procedure. I will need to take a look at some of the [`numba` documentation](http://numba.pydata.org/) and then [Levi's `esda` fork](https://github.com/ljwolf/esda/blob/moran-perf/esda/crand.py) to understand how they work together. 

# GSoC To Do

- Finalize docstrings/doctests/functions
- Talk to Levi about migrating each function over to `numba` 
- Start to prepare for a PR with esda by setting up a local testing environment
    - Need to be able to run `from esda import losh`, `from esda import Local_Join_Count`, etc. in local testing environment to fully check tests