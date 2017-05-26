;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Use to display syntax for my procedures
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


PRO helpme,compile=compile
	print, ''
	print, 'GENIE HELP'
	print, '(1) utilities'
	print, '(2) plotting/widgets'
	print, '(3) GENPOS'
	print, '(4) GENSPEC'
	print, '(5) GENRAD'
	print, '(6) GENTRAN'
	print, '(7) IMPSPEC'
	print, '(8) THACO'
	print, '(9) upper level routines'
	print, ''
	read, 'Choice?  ',var
	print, ''
	print, ''

	help_dir='/usr/local/cmod/idl/GENIE/'
	IF var EQ 1 THEN BEGIN
		IF logname() EQ 'mlreinke' AND keyword_set(compile) THEN BEGIN
			MK_HTML_HELP,help_dir+['mlr_functions.pro','oploterror.pro','makesym.pro','mpfit.pro','mpfitfun.pro'], $
			'/home/mlreinke/GENIE/doc/genie_utilities.html', title='GENIE - utilities'
			spawn,'firefox /home/mlreinke/GENIE/doc/genie_utilities.html',output,ereror
		ENDIF ELSE spawn, 'firefox /usr/local/cmod/idl/GENIE/doc/genie_utilities.html', output, error
	ENDIF

	IF var EQ 2 THEN BEGIN
		IF logname() EQ 'mlreinke' AND keyword_set(compile) THEN BEGIN
			MK_HTML_HELP,help_dir+['general_plot.pro','widgets/w_xtomo.pro','widgets/w_axuv.pro',$
				'widgets/w_spec.pro','widgets/w_impspec.pro'],'/home/mlreinke/GENIE/doc/genie_plotting.html', $
				title='GENIE - plotting/widgets'
			spawn, 'firefox /home/mlreinke/GENIE/doc/genie_plotting.html',output,error
		ENDIF ELSE spawn, 'firefox /usr/local/cmod/idl/GENIE/doc/genie_plotting.html', output, error
        ENDIF
	
	IF var EQ 3 THEN BEGIN
		IF logname() EQ 'mlreinke' AND keyword_set(compile) THEN BEGIN
			MK_HTML_HELP,help_dir+['genie_line.pro','shellinvert.pro','GENPOS/'+['genpos.pro','genpos_sphere.pro','genpos_utility.pro',$
				'genpos_xraytomo.pro']],'/home/mlreinke/GENIE/doc/genie_genpos.html', $
				title='GENIE - GENPOS'
			spawn, 'firefox /home/mlreinke/GENIE/doc/genie_genpos.html',output,error
		ENDIF ELSE spawn, 'firefox /usr/local/cmod/idl/GENIE/doc/genie_genpos.html', output, error
        ENDIF	

	IF var EQ 4 THEN BEGIN
		IF logname() EQ 'mlreinke' AND keyword_set(compile) THEN BEGIN
			MK_HTML_HELP,help_dir+['genie_spec.pro','GENSPEC/'+['genspec.pro','genspec_sphere.pro','genspec_utility.pro',$
				'genspec_chromex.pro']],'/home/mlreinke/GENIE/doc/genie_genspec.html', $
				title='GENIE - GENSPEC'
			spawn, 'firefox /home/mlreinke/GENIE/doc/genie_genspec.html',output,error
		ENDIF ELSE spawn, 'firefox /usr/local/cmod/idl/GENIE/doc/genie_genspec.html', output, error
        ENDIF	
	IF var EQ 5 THEN BEGIN
		IF logname() EQ 'mlreinke' AND keyword_set(compile) THEN BEGIN
			MK_HTML_HELP,help_dir+'GENRAD/'+['cooling_curves.pro','read_ark_table.pro','genrad_argon.pro',$
			'genrad_hirexsr_thaco.pro','genrad_mcp.pro','genrad.pro', 'genrad_xraytomo.pro',$
			'johnson_hinnov_emissivity.pro', 'read_hullac_data.pro','sxr_vuv_molybdenum.pro'],$
			'/home/mlreinke/GENIE/doc/genie_genrad.html', $
				title='GENIE - GENRAD'
			spawn, 'firefox /home/mlreinke/GENIE/doc/genie_genrad.html',output,error
		ENDIF ELSE spawn, 'firefox /usr/local/cmod/idl/GENIE/doc/genie_genrad.html', output, error
        ENDIF	
	IF var EQ 6 THEN BEGIN
		IF logname() EQ 'mlreinke' AND keyword_set(compile) THEN BEGIN
			MK_HTML_HELP,help_dir+'GENTRAN/'+['ave_ionization.pro','gentran.pro','ionrec.pro',$
			'read_adas_files.pro','read_fracabund_tables.pro','read_loch_data.pro',$
			'read_mattioli_data.pro','read_mazzotta_data.pro','read_putterich_data.pro','/strahl/run_cmod_strahl.pro',$
			'/strahl/make_cmod_grid.pro','/strahl/make_cmod_pp.pro','/strahl/make_cmod_param.pro'],$
			'/home/mlreinke/GENIE/doc/genie_gentran.html', $
				title='GENIE - GENTRAN'
			spawn, 'firefox /home/mlreinke/GENIE/doc/genie_gentran.html',output,error
		ENDIF ELSE spawn, 'firefox /usr/local/cmod/idl/GENIE/doc/genie_gentran.html', output, error
        ENDIF	
	IF var EQ 7 THEN BEGIN
		IF logname() EQ 'mlreinke' AND keyword_set(compile) THEN BEGIN
			MK_HTML_HELP,help_dir+['IMPSPEC/impspec.pro'], $
			'/home/mlreinke/GENIE/doc/genie_impspec.html', title='GENIE - IMPSPEC'
			spawn,'firefox /home/mlreinke/GENIE/doc/genie_impspec.html',output,ereror
		ENDIF ELSE spawn, 'firefox /usr/local/cmod/idl/GENIE/doc/genie_impspec.html', output, error
        ENDIF
	IF var EQ 8 THEN BEGIN
		IF logname() EQ 'mlreinke' AND keyword_set(compile) THEN BEGIN
			MK_HTML_HELP,'/usr/local/cmod/idl/HIREXSR/'+['hirexsr_bin_spectra.pro','hirexsr_calc_moments.pro',$
			'hirexsr_calc_profiles.pro','hirexsr_fit_ellipse.pro','hirexsr_fit_spectra.pro','hirexsr_tree_utilities.pro',$
			'hirexsr_shot_analysis_tools.pro','batch_hirexsr_run.pro','hirexsr_load_data.pro'], $
			'/home/mlreinke/GENIE/doc/thaco.html', title='THACO'
			spawn,'firefox /home/mlreinke/GENIE/doc/thaco.html',output,ereror
		ENDIF ELSE spawn, 'firefox /usr/local/cmod/idl/GENIE/doc/thaco.html', output, error
        ENDIF
	IF var EQ 9 THEN BEGIN
		IF logname() EQ 'mlreinke' AND keyword_set(compile) THEN BEGIN
			MK_HTML_HELP,help_dir+['neutron_rate.pro','zeff_neo.pro','qfit.pro','read_atomicmass_tables.pro','read_xray_data.pro'], $
			'/home/mlreinke/GENIE/doc/genie_general.html', title='GENIE - GENERAL'
			spawn,'firefox /home/mlreinke/GENIE/doc/genie_general.html',output,ereror
	ENDIF ELSE spawn, 'firefox /usr/local/cmod/idl/GENIE/doc/genie_general.html', output, error
        
     ENDIF
END
