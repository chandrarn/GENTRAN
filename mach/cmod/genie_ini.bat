; codes will use these environment variables for loading configs (ie vuv_lines, etc), so set it to default if unset
setenv,'MACHINE=cmod'
LOGNAME=getenv('LOGNAME')
MACHINE=getenv('MACHINE')
GENIE_PATH=getenv('GENIE_PATH') & if GENIE_PATH eq '' then begin & GENIE_PATH='/home/mlreinke/GENIE/' & setenv,'GENIE_PATH='+GENIE_PATH & endif

;!path='/usr/local/cmod/codes/efit/idl/:'+!path
!path=expand_path('+'+GENIE_PATH)+':'+!path

;!path='/usr/local/cmod/idl/GENIE/usr/mlreinke/:'+!path
!path=expand_path('+'+GENIE_PATH+'usr/'+LOGNAME)+':'+!path

; define genie routines to compile
gpros=[GENIE_PATH+'mlr_functions.pro', $
       GENIE_PATH+'load_'+MACHINE+'.pro', $
       GENIE_PATH+'mach/'+MACHINE+'/set_genie_env.pro', $
       GENIE_PATH+'csv_convert.pro', $
       GENIE_PATH+'csv_read.pro', $
       GENIE_PATH+'tsv_read.pro', $
       GENIE_PATH+'ssvdc.pro', $
       GENIE_PATH+'genie_help.pro', $
       GENIE_PATH+'general_plot.pro', $
       GENIE_PATH+'read_atomicmass_tables.pro', $
       GENIE_PATH+'shellinvert.pro', $
       GENIE_PATH+'zeff_neo.pro', $
       GENIE_PATH+'neutron_rate.pro', $
       '/usr/local/cmod/idl/quick_fit.pro',$
       GENIE_PATH+'qfit.pro', $
       GENIE_PATH+'read_gpfit.pro', $
       GENIE_PATH+'GENRAD/cooling_curves.pro', $
       GENIE_PATH+'GENRAD/read_ark_table.pro', $
       '/usr/local/cmod/idl/HIREXSR/hirexsr_load_data.pro', $
       GENIE_PATH+'GENRAD/johnson_hinnov_emissivity.pro', $
       GENIE_PATH+'GENRAD/sxr_vuv_molybdenum.pro', $
       GENIE_PATH+'GENRAD/read_hullac_data.pro', $
       GENIE_PATH+'GENTRAN/read_fracabund_tables.pro', $
       GENIE_PATH+'GENTRAN/ave_ionization.pro', $
       GENIE_PATH+'GENTRAN/read_adas_files.pro', $
       GENIE_PATH+'GENTRAN/read_loch_data.pro', $
       GENIE_PATH+'GENTRAN/read_mattioli_data.pro', $
       GENIE_PATH+'GENTRAN/read_putterich_data.pro', $
       GENIE_PATH+'GENTRAN/ionrec.pro', $
       GENIE_PATH+'GENTRAN/gentran.pro', $
       '/usr/local/cmod/codes/spectroscopy/strahl/idl/set_color_right.pro',$
       '/usr/local/cmod/codes/spectroscopy/strahl/idl/cw.pro',$
       '/usr/local/cmod/codes/spectroscopy/strahl/idl/psprint.pro',$
       '/usr/local/cmod/codes/spectroscopy/strahl/idl/bilin_interpol.pro',$
       GENIE_PATH+'GENTRAN/strahl/make_cmod_grid.pro', $
       GENIE_PATH+'GENTRAN/strahl/make_cmod_param.pro', $
       GENIE_PATH+'GENTRAN/strahl/make_cmod_pp.pro', $
       GENIE_PATH+'GENTRAN/strahl/run_cmod_strahl.pro', $
       '/usr/local/cmod/codes/spectroscopy/strahl/idl/read_strahl.pro',$
       GENIE_PATH+'genie_line.pro', $
       GENIE_PATH+'genie_spec.pro', $
       GENIE_PATH+'GENPOS/genpos.pro', $
       GENIE_PATH+'GENSPEC/genspec.pro', $
       GENIE_PATH+'GENPOS/genpos_utility.pro', $
       GENIE_PATH+'GENPOS/genpos_xraytomo.pro', $
       GENIE_PATH+'GENSPEC/genspec_utility.pro', $
       GENIE_PATH+'GENRAD/genrad_argon.pro', $
       GENIE_PATH+'GENRAD/genrad_xraytomo.pro', $
       GENIE_PATH+'GENRAD/genrad_hirexsr_thaco.pro', $
       GENIE_PATH+'IMPSPEC/impspec.pro' ]

; don't compile widgets for now since they can overwrite each other (ie some procedures may be named the same, etc)
;       GENIE_PATH+'widgets/w_axuv.pro', $
;       GENIE_PATH+'widgets/w_chromex.pro', $
;       GENIE_PATH+'widgets/w_impspec.pro', $
;       GENIE_PATH+'widgets/w_spec.pro', $
;       GENIE_PATH+'widgets/w_xtomo.pro' $
;       ]

; make the genie init script and run it...
openw,gini,'~/gini.bat',/get_lun
 npros=n_elements(gpros)
 for i=0,npros-1 do begin &$
  str='.compile '+gpros[i] &$
  printf,gini,str &$
 endfor
close,gini
@~/gini.bat
spawn,'rm ~/gini.bat'
makesym,10
print, 'GENIE loaded from: '+GENIE_PATH
set_genie_env
print, 'use helpme to view HTML-based documentation'

