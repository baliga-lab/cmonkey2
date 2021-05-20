Changes
=======

Version 1.2.12, 2021/05/20
--------------------------

  - support for MEME up to 5.3.3

Version 1.2.11, 2018/08/24
--------------------------

  - bug fix: resume feature

Version 1.2.10, 2018/07/16
--------------------------

  - bug fix: BSCM switch ratios reference

Version 1.2.9, 2018/06/20
-------------------------

  - bug fix: set enrichment some parameters not
    referenced correctly

Version 1.2.8, 2018/06/20
-------------------------

  - removed memory profiling debug (can cause problems)
  - more tolerant debug score writing

Version 1.2.7, 2018/06/19
-------------------------

  - support for MEME versions 4.11.4 and 4.12.0
  - run log data now stored in database for consistency

Version 1.2.6, 2017/02/03
-------------------------

  - database aspects now handled with sqlalchemy
  - support for MEME 4.11.x added

Version 1.2.5, 2017/02/03
-------------------------

  - cm2view integrates the updated sequence logo viewer

Version 1.2.4, 2017/01/17
-------------------------

  - dependency update, missing packages included

Version 1.2.3, 2016/11/17
-------------------------

  - tweaks to export and plotting
  - fix for clusters_per_row option

Version 1.2.2, 2016/09/06
-------------------------

  - added the cm2plot tool

Version 1.2.1, 2016/08/25
-------------------------

  - fix: numeric row names are converted to string
  - fix: DataMatrix.subtract_with_quantile() did not have any effect
  - fix: set enrichment, bonferroni cutoff denominator-numerator were
    in wrong positions, value of the minium reference matrix value

Version 1.2.0, 2016/08/08
-------------------------

Initial PyPI version
