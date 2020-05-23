# Welcome to my first GSOC blog post!

Google Summer of Code (GSOC) 2020 has been going for a couple weeks now and I thought it would be fun to share what progress has been made so far. If you're just joining now you should know that the first month or so is the 'Community Bonding Period'. During this time the student developers get to know their community, establish expectations, and set recurring check-in meetings. 

# Initial community bonding

My main avenues of building rapport with the community are participating in the gitter channel and GSOC project calls. In the gitter channel I've linked a SO post about identifying [which CPU core is being used a multi-threaded task](https://stackoverflow.com/questions/20669881/identifying-a-processor-core-or-worker-id-parallel-python), as well as linked some large, (nearly!) ready-to-go spatial datasets for testing functions. In this latter case I learned that one of PySAL's contributors, Eli Knaap,  had already made the [https://open.quiltdata.com/b/spatial-ucr/packages/](https://open.quiltdata.com/b/spatial-ucr/packages/) and they are available to the public. I plan to use these in my own testing down the line! 

The first PySAL-wide GSOC call was 5/15/2020. It was a pleasure to chat with the other PySAL GSOC students [Mragank Shekhar](https://summerofcode.withgoogle.com/projects/?sp-search=pysal#5775104799145984) and [Pablo Estrada](https://summerofcode.withgoogle.com/projects/?sp-search=pysal#6472262816890880). Mragank will be focusing on raster integration with PySAL, and Pablo will be developing panel spatial econometric models. Both of these projects will offer some tremendous developments to PySAL and I look forward to getting to know them more over the summer. 

An ESDA enhancements call was held today on 5/22/2020. I updated the project mentors on my initial attempts to code up the estimators - I've already gotten a bit of a head start (see below!). This was a positive update and I spoke at length with Dr. Levi Wolf and Dr. Taylor Oshan as to the next steps for the project steps over the coming weeks. 

# What have I been up to in addition to community bonding?

I have created a [Github repository that will act as a work diary](https://github.com/jeffcsauer/GSOC2020). The repository is organized as follows:

```
GSOC2020
│   README.md
│   Sauer_GSOC2020_ESDA.pdf   
|
└───scratch                 # this is where my work on coding up the estimators
│   │   LOSH_workbook.ipynb
|   |   migration.ipynb
│   │   ...
└───review                  # this is where I review papers and key concepts
|   │   OLJC_workbook_paper_review.ipynb
|   │   ...
└───functions               # this is where the 'final' or 'polished' functions will reside
|   |   LOSH.py  
|   |   Local_Join_Counts.py
|   |   ...
└───notebooks               # notebooks demonstrating the functions on toy and real-world data
|   |   LOSH.ipynb
|   |   Local_Join_Counts.ipynb
|   |   ...
```

So far I have worked through the PySAL implementation of local join counts and Moran's I, and I carried out some local calculations of uni-/bi-/multi-variate local join counts and LOSH. In the coming weeks I will adapt each of these into a function with the following API format:

`function_name(vector_of_y_values, spatial_weights_object, parameter_1, parameter_2, ... )`

One of the key takeaways of the today's call was the structure of the function. Dr. Wolf explained that older ESDA functions in like `localmoran.py` tend to 'get heavy' due to duplicate caching of vectors and other internal attributes. More recent ESDA functions like `lee.py` use a scikit-learn structure, with others using scipy. Personally I am leaning towards a scipy approach! 

# What do I want to get out of GSOC?

**Conceptual takeaway(s) so far**

After already spending a couple of weeks really digging in 'under the hood' of what is going in the ESDA functions, I am starting to understand the intrinsic value of the concept that is spatial autocorrelation. We need to have a way to quantify the innate reasoning of what Tobler described as 'nearer things being more related than farther things' It has made me think of so many other metrics that might be useful (but might not have the catchiest name), like a 'cluster disruption test' that would assess the influence of removing a unit on local/global spatial autocorrelation. 

**Developer takeaway(s) so far**

First and foremost, with each passing day I gain a deeper appreciation for what PySAL represents. It is a federation of open source, user contributed packages that represents one of the most developed spatial analysis toolkits out there. There is an immense amount of work that is done each day to enhance existing PySAL code and identify relevant areas for expansion. By reviewing existing parts of the PySAL codebase I have been able to understand complex functions and implement tricks that would have never occurred to me. It really is an aspiring piece of work. It is an inspiring group of academics, industry professionals, 

Looking back as an undergrad eager to get involved in the research process I was so intimidated by tools like Github, notebooks, and, in general, python! Now I am quite fond of them! :grin: 

# GSOC, PySAL, Join Counts, LOSH - what are these terms?

- Google Summer of Code (GSOC): an annual program that grows each year that focuses on contributing funds towards open-source software development. More info [here](https://summerofcode.withgoogle.com/).
- The Python Spatial Analysis Library (PySAL): a federation of spatial analysis packages fueled by a passionate community of contributors. More info [here](https://pysal.org/). 
- Exploratory Spatial Data Analysis (ESDA): a technique in spatial analysis that focuses on revealing spatial autocorrelation in your data. A nice Towards Data Science post summarizing some concepts [here](https://towardsdatascience.com/what-is-exploratory-spatial-data-analysis-esda-335da79026ee). 
- OLJC: operational local join counts, a type of ESDA statistic from [Anselin and Li 2019](https://link.springer.com/article/10.1007%2Fs10109-019-00299-x). 
- LOSH: local indicators of spatial heteroskedasticity, a type of ESDA statistic from [Ord and Getis (2012)](https://link.springer.com/article/10.1007/s00168-011-0492-y).




