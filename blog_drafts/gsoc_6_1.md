# DRAFT

# More community bonding!

GSOC-wide call, what did we talk about

hopped on a quick call with Levi to talk about steps for the coming weeks. will be circulating the losh function to the ESDA enhancemnets mentors 

# A note on the LOSH overhaul 

i was implementing the unscaled version of the statistic (eq X). however, when comparing to the output from Roger Bivand's `spdep::LOSH()` function my Hi values were equal to tthe raw residuals. After spending some times, I noticed I had made two mistakes:

1. I carrying out an 

# What have I been up to in addition to community bonding?

The past couple of weeks I have been working on pseudo-coding the implementation of LOSH, uni/bi/multi-variate LJC. I've now migrated all of them into functions with the structure of scikit-learn. 

- more details

# The importance of inference

working through Bivand and Wong 2018, 'Comparing implementations of global and local indicators of spatial association'. a beast of a paper! this could be a pathway into the p-value crisis for geography. needs to be a greater understanding about p-value sensitivity and develop more rigorous frameworks with perhaps strange quesitons. for example, how do different lattice structures impact p-value calculations? 

**More advances in my understanding of development**

scikit learn developer agreement: https://scikit-learn.org/stable/developers/develop.html

PEP8 formatting. mentioned in pysal/esda issues #118. I implemented this in my own `arcospy` package and I agree! working on this now.
