# Introduction

The midpoint of summer is approaching. This coming week will mark a midway point in GSoC2020. As I move away from the hurdle of the last couple of weeks I am setting into a pace of fine-tuning the functions and creating robust notebooks to demonstrate their functionality. While bugs are nuanced errors are sure to arise, I am feeling good about the latter 2/3 of GSoC2020.

# Strategy: demonstration notebooks

I am working on two 'final' demonstration notebooks that showcase the LOSH and LJC functions. These two notebooks are going to focus less on the coding and more on the application of the functions. These types of notebooks are important as many users want to understand the substantive meaning of the functions (i.e. 'what does a high LOSH value mean?'). To assuage any user concerns about deploying the functions in Python I will use the same examples that appear in `R` `spdep::LOSH()` as well as GeoDa.

# Software development reflections

## Understanding tests

If you take a look at my [GSoC2020 proposal](https://docs.google.com/document/d/1WjHjy5Eyk4WG5QWfnsnhWg1r4-e09JXXCx0iaPphg6c/edit?usp=sharing), you will see that we are heading into Week 8. The docstrings and doctests are present in each function. I've also written a series of straightforward tests using `unittest`, available for view [here](https://github.com/jeffcsauer/GSOC2020/tree/master/tests). I implemented the functions and tests in a fork of `esda`. All the tests are passing at the moment and - if the functions remain stable over the next few weeks - I will likely submit a pull request. 

On the topic of tests, I was quite concerned with what makes a 'good' test. I read up on several sources from [here]() and [there](), but I wasn't really understanding the fundamental point of testing. In my mind, I understand the logic of a test. For example: 

```{python}
assert np.array_equal(ljc.LJC, [0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 3, 2, 2, 3, 3, 2])
```

The above test is checking the values of `ljc.LJC` against an array. Because the input dataset is fixed in its arrangement of 0s and 1s, we want to ensure thta that the functions returns the expected values. 

The question then arises: what if you have fundamentally misunderstood how to implement the function? Or perhaps you your reference values have been calculated incorrectly? If either of these were the case, your test becomes meaningless (or, even worse, permitting wrong values). When I brought up my concerns with Levi, he linked a [fitting picture](https://www.reddit.com/r/ProgrammerHumor/comments/3bg3ro/yay_all_unit_tests_passing/) which captured my concerns. Levi explained that we want to write tests that can capture large and obvious incorrect behavior. In this approach, tests are less about making sure that the function gets the correct values. They are more about checking that the right types of data are returned, or that inappropriate inputs raise errors. 

## Should the three `Local Join Count` functions be merged into a single file?

In the `esda` implementation of Moran, the various versions of Moran are implemented in a single file. This single-file approach reduces the total number of files and centralizes the location where future edits are to be made. On the other hand, keeping the functions separate might prevent editing errors. Given the partial overlap of the code in each function, I could imagine a scenario wherein a user carries out a Control+F operation and accidentally alters all three functions. 

There is a good discussion of the advantages on keeping the functions split in this [SO post](https://stackoverflow.com/questions/15580539/what-are-the-advantages-of-using-more-then-1-code-file-for-a-project-c/15580555#15580555). To briefly summarize, the core advantage of keeping the files independent are: 

- leveraging single-file compiler speed
- project organization and readability 
- facilitating 'code reuse'

The above post actually links to [this](https://www.gamedev.net/tutorials/_/technical/general-programming/organizing-code-files-in-c-and-c-r1798/) explanation posted on Gamedev.net.

# Challenge: migrating to `numba`

As a final push to improve the functions, Levi suggested that I try migrating the inference in my functions to use the new `numba`-based conditional randomization engine (detailed here on [Levi's `esda` fork](https://github.com/ljwolf/esda/blob/moran-perf/esda/crand.py)). The new engine is live on the primary PySAL `esda` submodule, so it is a worthwhile effort to migrate as many functions as possible.  

This new engine only applies to the local join count functions as of now. I spent some time on Friday migrating these functions. I encountered the most difficulty at the conceptual stage of implementation. For example, take a look at this `numba`-ized `_crand()` code for the univariate local join count: 

```
@_njit(fastmath=True)
def _ljc_uni(i, z, permuted_ids, weights_i, scaling):
    zi, zrand = _prepare_univariate(i, z, permuted_ids, weights_i)
    return zi * (zrand @ weights_i)
```

Although it's only a few lines, there is actually quite a bit going on. The first is the header:

```
@_njit(fastmath=True)
```

This 'reduces numerical rigour with view of gaining additional performance' ([source](https://numba.pydata.org/numba-doc/latest/user/performance-tips.html#fastmath)). Next, we feed in a series of...

```
def _ljc_uni(i, z, permuted_ids, weights_i, scaling):
```

is, zs, permuted_ids, weights_i, and scaling into a function...wait a second, what are we feeding these values into? Aha - now you see my conceptual hiccup. You see, function is actually a sub-routine of the larger [`_crand()` engine](https://github.com/ljwolf/esda/blob/moran-perf/esda/crand.py#L67). So these is, zs, permuted_ids, and weights are being passed to the `_crand()` function with an additional *statistical_function* (i.e. `stat_func`), which we defined above as `_ljc_uni`. I was stuck for quite some time on understanding which of these functions is actually called first. Understanding that `_ljc_uni` - or any of the other statistical functions in `esda` - are a sub-routine of `_crand()` is essential for writing your own. 

From there, we can generate a series of `zi` and `zrand` objects using yet another pre-written function called `_prepare_univariate`. This function handles the shuffling and creation of z-minus-i and randomly values for these i units.

```
zi, zrand = _prepare_univariate(i, z, permuted_ids, weights_i)
```

We can then do some quick math with these output products to get the local join counts:

```
return zi * (zrand @ weights_i)
``` 

This operation is taking a given value, `zi`, and multiplying it by the dot product of the randomized values times the weights matrix. Note that this final operation changes depending on what local statistic you are interested in. You can see the results at the end of each new [validation notebook](https://github.com/jeffcsauer/GSOC2020/tree/master/validation). 

# GSoC To Do

- Finalize demonstration notebooks comparing PySAL functionality to `R` `spdep::LOSH()` and GeoDa. 
- Ask Levi/Serge about a contributing a pull request.