# GSOC Updates

As Phase II comes to an end, I am happy to state that significant progress has been made on the original details of the GSOC ESDA Enhancements project. All proposed estimators have been implemented with docstrings, doctests, and stylized using the PEP8 style guide. A PR is open documenting these [proposed contributions](https://github.com/pysal/esda/pull/139). A GSOC call on 7/24 outlined several additional tasks to keep the contribution momentum going: 

- Updating LOSH inference to new the `__crand()` engine
 
    - The existing `losh()` estimator can (and maybe should!) be migrated to the new `__crand()` engine.  

- Draft univariate and multivariate local Geary estimator

    - An additional exploratory spatial statistic of interest is the Local Geary estimator, originally outlined in [Anselin 1995](https://www.google.com/search?client=firefox-b-1-d&q=anselin+1995).   

- Experiment with the LOSH statistic as it relates to [Bivand and Wong's 2018 paper](https://link.springer.com/article/10.1007/s11749-018-0599-x). 

# Strategy

I have begun implementing the local Geary statistics [here](https://github.com/jeffcsauer/GSOC2020/blob/master/review/Local_Geary_Workbook.ipynb). The local Geary values are matching the output from GeoDa, but there are a few additional steps that are still necessary. First and foremost is inference via conditional randomization. Second is the partitioning of the local Geary values and inference into cluster groups. 

After basic inference is up and running for the local Geary estimators I will likely switch to updating the new `__crand()` engine. Once this is working I will experiment with the statistic on some of the claims raised in Bivand and Wong 2018. Plenty to do in the next eight weeks or so!