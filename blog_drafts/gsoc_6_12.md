# Introduction

Slightly shorter blog post this week! More substantial update last week available [here](https://jeffcsauer.github.io/post/2020/06/07/gsoc-blog-2-coding-progress-revisions-and-more/). En route to North Carolina to help settle a relative's house. Listening to the audiobook of [*How To Be An Antiracist* by Ibram X. Kendi](https://www.ibramxkendi.com/how-to-be-an-antiracist-1) during the drive.  

## ESDA call

Had a ESDA-specific call with Serge this week on 6/12. We discussed how to fine-tune the LOSH function, which is now *mostly* operational. The function currently offers chi-square inference (demostration of code available [here](https://github.com/jeffcsauer/GSOC2020/blob/master/validation/Validation_LOSH.ipynb)). I also filled out the docstrings and reorganized the Github page so the functions are now clearly identifiable in a functions folder. 

Serge and I spent a bulk of the call talking about the advantages and disadvantages between implementing these estimators as functions versus classes. The primary disadvantage of functions is that they may store an excess number of objects used in the calculation of a given estimator. This can cause a range of problems beyond memory limitations, such as unwittingly overwriting existing objects or overwhelming the end-user with too many miscellaneous objects. Beyond keeping programs more organized, classes offer the advantage of inheritance wherein a later part of a program or function may use the output from a class. 

Serge and I also reviewed some of the [issues-as-feedback](https://github.com/jeffcsauer/GSOC2020/issues) Levi provided on my GSOC20202 workbook repository. I hope to address these in the next week or so. The issues include:
- Using warnings rather than print statements
- Simplifying the $VarH_i$ calculation
- Switching to sparse weights to calculate row sums
- Breaking notebooks into .py function files

I greatly appreciate this feedback as it directly addresses some of the weakest parts of my `losh` implementation. As a new user of PySAL, I was struggling to understand how to get some of the properties I wanted (such as rowsums) out of a weights object without using something like list comprehension. For example, this was my solution to get the row sums:

`Wrs = [round(np.sum(list(w[i].values()))) for i in range(len(y))]`

Levi pointed out that I can use sparse weights to achieve the same array in a much more intuitive fashion:

`row_sums = numpy.asarray(w.sparse.sum(axis=1)).flatten()`

Exactly what I was wanted to do but had no idea how!

## GSoC project progress - Wrapping up week 2

As previously stated, the LOSH function with chi-square inference is now *mostly* operational. I'm getting nearly all of the same results as the `R` `spdep::LOSH()` function, although there are a few slight differences due to a rounding of an integer. The LOSH function came together after a late-night session of coding to address concern [3] from last week's blog post - handling both standardized and unstandardized weights. Ended up implementing the fix by Bivand that expertly handles both situations in the following manner (found in the third line of the [`losh` $VarH_i$ calculation](https://github.com/jeffcsauer/GSOC2020/blob/master/validation/Validation_LOSH.ipynb)): 

`sum(y^2)`

For situations where the weights are row standardized, this simply squares a decimal and takes the sum. When the weights are not row standardized, it still squares the values, but because the values are now row standardized it simply does $1^2$ (i.e. 1) and then takes the sum. To me this is an amazingly clever solution as it eliminates the need to have additional code that checks what type of neighbors are passed to the function.

## GSOC To-Do list

- Reading and reviewing inference strategy for Local Join Counts. Will be working on implementing these in the coming week(s)! Following Anselin and Li (2019).
- Addressing issues raised by Levi to improve the `losh` function. 