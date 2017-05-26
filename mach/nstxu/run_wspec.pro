.compile /u/mreinke/GENIE/mlr_functions.pro
.compile /u/mreinke/GENIE/load_cmod.pro
.compile /u/mreinke/GENIE/tsv_read.pro
.compile /u/mreinke/GENIE/mach/nstxu/w_spec_old.pro
loadct,12,/silent
xwplot
set_plot,'x'
w_spec
