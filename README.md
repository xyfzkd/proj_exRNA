# Introduction

# proj_exRNA
## mapping
`mapping.sh`
## mapping results visualization
### aim
* 统计不同RNA类型reads的比例并以饼图展示
* 统计不同RNA类型mapped reads的长度分布
### counts file scripts
my own `length_ratio.sh` creats files do not match format of other sample's counts results, but my plot scripts work on my counts data

to match format of other sample's counts result, we need another script `read.sh` to retrieve length counts.

* why must match?

No need, because they are used for plotting separately. However, if you plot a merged pic, using new script `read.sh` will be convenient.(length ingored, just do as `length_ratio.sh` )
### counts file plotting scripts
`ratio.ipynb`
`length.ipynb`
may not work well
## construct matix
`construct_mx.sh`

five files containing reads counts of different RNA have to be merged to a some-column-drop one provided in cnode.
## matix processing
...
## feature selection
...
