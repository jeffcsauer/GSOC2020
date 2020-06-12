DRAFT

# Introduction

Brief blog post this week! More substantial update last week available [here](https://jeffcsauer.github.io/post/2020/06/07/gsoc-blog-2-coding-progress-revisions-and-more/). En route to North Carolina to help settle my grandparent's house and clear out the remaining junk. Listening to audiobook of *How To Be An Antiracist* by Ibram X. Kendi during the drive.  

## ESDA call

- call with Serge on 6/12
- discussed fine-tuning of LOSH function
- question was handling unstandardized weights, yet it seems that the LOSH function in R gives the same results (https://github.com/jeffcsauer/GSOC2020/blob/master/validation/data/Validation_LOSH_real_world.Rmd) regarldess of standardized vs unstandardized weights. 
- my modified function does not at the moment. figuring this out
- added more to docstrings, function running quite well imo

## GSoC project progress - Wrapping up week 2

LOSH function with chi-square inference now *mostly* operational. Available [here](https://github.com/jeffcsauer/GSOC2020/blob/master/validation/Validation_LOSH.ipynb). Validating the results against the `R` `spdep::LOSH()` function. 

Reading and reviewing inference strategy for Local Join Counts. Will be working on implementing these in the coming week! Following Anselin and Li (2019), we need the following elements to get calculate the p-values:
-
-
-


As previously mentioned, I posted a longer blog post last week. i would recommend reading that post to get an idea of the project progress!

