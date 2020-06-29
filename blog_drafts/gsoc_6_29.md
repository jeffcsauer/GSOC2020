# Introduction

We are rapidly heading towards the final week for the first phase of GSoC 2020! What a wild ride it has been - I can't believe it's only been around a month and a half! Brief blog post this week as there is lots to do. 

# Phase 1 Progress Recap

I've produced solid draft of all of the functions original identified in the proposal, with inference working for *nearly* all of them. Here is a table summarizing the progress made so far: 

| Function              | Generating correct values | Generating correct inference | PEP8 styling | Docstrings | Tests      |
|-----------------------|----------------------------|-------------------------------|---------------|------------|------------|
| `LOSH`                | Yes                        | Yes                           | Mostly        | Mostly     | Need to do |
| `Local_Join_Count`    | Yes                        | Yes                           | Mostly        | Mostly     | Need to do |
| `Local_Join_Count_BV` | Selectively                | No for Case 1, Yes for Case 2 | Mostly        | Mostly     | Need to do |
| `Local_Join_Count_MV` | Yes                        | Yes                           | Mostly        | Mostly     | Need to do |

# Update on inference for the `local_join_count` functions

Thanks to the input from Levi and Serge I was able to get inference up and running for the most of the `local_join_count` functions. Quite interestingly, I've run into a semi-large hiccup with `Local_Join_Count_BV`. I am able to produce the correct local join counts on some toy datasets, but I can only get the function to agree with GeoDa in select datasets. I'm investigating this issue in detail [here](https://github.com/jeffcsauer/GSOC2020/blob/master/validation/Understanding_BV_LJC.ipynb).

To add a bit of reflection to this post, I would share with readers that sometimes it is worthwhile to trust your gut instinct. As explained on the previous blog post, I had anticipated that inference would be the most challenging aspect of the implementation. This turned out to be exactly the case. Due to my lack of understanding and lack of confidence, finalizing the methods of inference has proven a bit taxing. I often get frustrated with myself and go for a jaunt in the woods of self-loathing. In pursuit of better health I am following Serge's advice to let things settle for a moment. Often these periods of brief respite will allow the solution to become apparent. 

# 'To-do' for the final week of Phase 1

For the final week, I would like to touch base with Levi and/or Serge about what is going on in the bivariate function. Relating to the above paragraph, I think I am at a point where I need their expertise in a close read of the code to figure out what is going on. 

I've also prepared a notebook demonstrating all of the functions [here](https://github.com/jeffcsauer/GSOC2020/blob/master/scratch/GSoC_Phase1_Demonstration.ipynb). Take a look if you're interested!

As we transition into the second phase I will need to start thinking about tests and finalizing the docstrings. 