# Welcome To Marcus Viscardi's Mess
A directory full of half-baked scripts and programing ideas

##Table of Contents:

1\. *assignReadsToGenesMultiprocessing*
> Trying various forms of [multiprocessing](https://docs.python.org/3.8/library/multiprocessing.html) in an effort to 
speed up read assignment. As well as, a major effort to utilize pandas dataframes in order to utilize C backend for
faster sorting/lookup/etc.

> Another option of using [Pandarallel](https://github.com/nalepae/pandarallel) could be cool, although current
**(April 9, 2020)** tests make vectorized processing (utilizing *pandas.apply()*) look faster that both
[Pandarallel](https://github.com/nalepae/pandarallel) and
[Muliprocessing.pool.Pool()](https://docs.python.org/3/library/multiprocessing.html#module-multiprocessing.pool)

>Large steps made in parsing in the \*.SAM and \*.allChrs.txt **(as of April 14, 2020)**, splitting them up based on
mapped chromosome, success with vectorize-ing the recoverMappedPortion function (may still be a bottleneck), and lots
of argParse library usage (mostly for debugging) **(as of April 15, 2020)**. Next step will be utilizing the
[pandas.Database.Merge()](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.merge.html)
functionality to utilize pandas back end for assigning the \*.SAM reads to the genes identified by chromosome index
in the \*.allChrs.txt annotations.

2\. *multiIndexingDataframes*
> Attempt at using [multiIndexing](https://pandas.pydata.org/docs/user_guide/advanced.html) with 
[pandas](https://pandas.pydata.org/docs/) dataframes to see if they'll be faster/lighter/easier-to-process

3\. *otherNoodling*
> Pretty much just random stuff that does not fit in elsewhere