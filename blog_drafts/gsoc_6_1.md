# Coding, Progress, Revisions, and More

It has been a few weeks since the last GSoC post (5/22/2020). [A lot has happened in the United States during this time](https://www.theatlantic.com/newsletters/archive/2020/05/george-floyd-americas-racial-contract/612313/), and more is sure to come. In the following  post I will share my progress on the project and reflect on the challenges faced so far. 

## More community bonding

The plan of hosting calls on Friday has continued without a hitch. On 5/29 we held a PySAL GSoC-wide call to ask general questions that might pertain to all of the projects. Although each of the PySAL GSoC projects are quite different, there were quite a few questions whose answers were useful to all of us. For example, a student asked about how progress should be publicly tracked in each progress. A brief discussion was held, and it was ultimtely decided that students should raise issues in their PySAL sub-modules pertinent to their projects. Additionally, the PySAL GSoC mentors and mentees considered different ways of hosting student code (e.g. work journal repositories, GitHub projects, or several small pull requests). Ultimately, it was decided that each project could use whatever system they preferred - I am opting for a work journal repository). After the GSoC-wide call I hopped on a quick meeting with Levi to talk about steps for the coming weeks. Given my progress (described in more detail below), he suggested I circulate the current implementation of the `losh()` function to the GSoC mentors. 

On 6/5 the GSoC students attended the PySAL-wide developer meeting, a monthly meeting where all active contributors to the PySAL organization hop on a call to talk about pressing issues and development in the software. Around 9 developers attended the meeting, 12 including the GSoC student mentors. These meetings are useful in that they provide students with direct insight into the trajectory of a software suite. When looking at PySAL from the outside, an introductory user might be overwhelmed by the sheer number of submodules and moving parts (both of these translate into a large, sprawling codebase). However, when sitting on in the developer meetings you can start to associate developers with specific projects and understand how PySAL is compartmentalized. At this particular meeting the GSoC students introduced themselves and their projects. A history of conversations and notes from previous developer meetings are available [here](https://github.com/pysal/pysal/wiki/Developer-meetings). 

## GSoC project progress - Wrapping up week 1

An updated project timeline for my project is available [here](https://docs.google.com/document/d/1WjHjy5Eyk4WG5QWfnsnhWg1r4-e09JXXCx0iaPphg6c/edit#). The most notable changes are a modification of the dates associated with each phase. According to the project timeline, these are the Phase 1 objectives: 

> Phase I: Implementation of estimators (June 1 - July 3, 5 weeks)

> MILESTONE by end of week 5: Implement univariate and multivariate local join count statistics with inference. 

> **Week 1** - code preparation: I will review relevant literature on the spatial estimators proposed for implementation.(Anselin & Xi, 2019; Ord & Getis, 2012) I will pseudo-code initial attempts at the estimators with a structure that mirrors existing estimators within esda such as Local_Moran. I will discuss with my mentor(s) what type of inference should be prioritized for each estimator (i.e. both chi-square and bootstrap inference have been proposed for LOSH).(Ord & Getis 2012; Xu et al., 2014)
**Week 2-4** - coding: Code first attempts at each estimator. Submit draft to mentor(s) for feedback and identify priority areas for optimization. Attend PySAL developer meeting the first week of June.
**Week 5** - code review: Integrate mentor(s) feedback to optimize calculation of estimators and ensure code structure consistency with other esda estimators. I will establish expected efficiency using functions like timeit.

Due to some luck of starting early, I have already [pseudo-coded each of the estimators (LOSH, uni/bi/multivariate local join counts)](https://github.com/jeffcsauer/GSOC2020/tree/master/scratch). I am glad for the head start because there is now ample time for GSoC mentors to provide specific feedback on the behavior of the functions (e.g. LOSH, see below). 

In addition to the pseudo-code, I've [migrated each estimator](https://github.com/jeffcsauer/GSOC2020/blob/master/scratch/migration.ipynb) to a scikit-learn style function. This was a bit tricky as I still do not fully understand how the different parts of a scikit-learn style function feed into each other (i.e. why are some objects available in one level but not another). Given the progress, my mentors advised that we focus on fine-tuning one estimator before moving on to the others. We will continue to work on LOSH for the next week or so, and then we will move on to different forms of local join counts. This may slightly alter the GSoC project timeline. Specifically, I imagine that we may reallocate the time spent for docstrings and tests towards perfecting the functions. Indeed, many of the docstrings are already drafted!

Notes on the various papers relevant to my projects are available [here](https://github.com/jeffcsauer/GSOC2020/tree/master/review).

## A note on the LOSH overhaul 

It is worthwhile to highlight a brief example of a benefit to starting the project early. In my initial implementation of the LOSH statistic all seemed to be going well. Given that I had a working function, I decided to go ahead and compare my output to that of Dr. Roger Bivand, a lead maintainer of the `spdep` package in `R`. However, when comparing the output of my `losh()` function to `spdep::LOSH()`, I noticed that my $H_i$ values were equal to the raw residuals. After spending some time comparing my implementation to that of `spdep::LOSH()`, I two critical mistakes:

1. Extra division. I was carrying out an extra division of the spatial weights. I had forgot that my function already assumed the spatial weights were row-standardized, so I did not need to divide again by the number of neighbors for a given unit $i$.
2. Wrong version of the $H_i$ statistic. This relates to the many forms of the LOSH statistic that appears in the original Ord and Getis 2012 paper. I was implementing the unscaled version of the LOSH statistic when in fact the scaled version was of interest. This scaled version incorporates the neighborhood residuals into the denominator.

After changing to the scaled version, my `losh()` function was matching `spdep::LOSH()`. I shared my progress with the GSoC mentors and next steps were identified: [1] cleaning up some of the function to follow PEP-8 styling, [2] adding chi-square based inference to the statistic, and [3] handling unstandardized weights. I've accomplished most of [1] and [2], but I'm still working on [3].

Ultimately was able to match the output from Bivand's `spdep::LOSH.cs()` function! this is a strong indication that the function is behaving as expected (so far).

## The importance of inference

The remaining work for all of the functions are developing strategies of inference. Levi suggested I read Bivand and Wong's 2018 paper entitle 'Comparing implementations of global and local indicators of spatial association'. It is an interesting paper for a few different reasons, but what struck me the most is that it appears to be a meta-analysis of several implementations of common spatial statistics. This meta-analysis leads Bivand and Wong to new insights, namely that simulation-based p-values (rather than analytical) may present a systematic bias.

An a separate 'meta' note, this paper helps quantitative geography take a step towards the p-value crisis and the replication crisis. These topics have been covered extensively elsewhere, but there needs to be a greater understanding of these topics within the field of geography. There have been questions like this before - namely the modifiable areal unit problem (MAUP) and S-MAUP statistics - but there is room for investigation into how specific types of geographic structures might impact inference (e.g. lattice vs. political boundaries).

## More advances in my understanding of development

A few things have come up this week that are explicitly related to development choices. The first was what style the new estimators should take. [I raised this issue](https://github.com/pysal/esda/issues/118) on the ESDA submodule page and got some initial feedback. I've stuck with the sci-kit learn format and spent some time reviewing their [developer agreement](https://scikit-learn.org/stable/developers/develop.html).

Additionally, I spent some time carrying out PEP8 formatting in all of the functions. There is handy jupyter function that can check each cell for PEP8 formatting consistency. At the beginning of the notebook, run the following command: 

`%load_ext pycodestyle_magic`

This loads the code format checker, and we turn it on with the following command:

`%pycodestyle_on`

When you next run a cell, this will evaluate the cell and highlight lines that are inconsistent with PEP8 guidelines. Just make sure to run `%pycodestyle_off% when you want to start running the cells for their output!