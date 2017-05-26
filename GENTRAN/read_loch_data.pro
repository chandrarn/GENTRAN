FUNCTION read_loch_rec_data,debug=debug,load=load
	;file_path='/home/mlreinke/idl/impurities/data/adas/acd96_ar.dat'
	;save_path='/home/mlreinke/idl/impurities/data/adas/acd96_ar.sav'
	file_path='/usr/local/cmod/idl/atomic_physics/adas/acd96_ar.dat'
	save_path='/usr/local/cmod/idl/atomic_physics/adas/acd96_ar.sav'
	IF keyword_set(load) THEN BEGIN
		restore, save_path
		RETURN,rates
	ENDIF
	rates=read_acd_file(file_path)
	;save,rates,filename=save_path
	RETURN,rates
END

FUNCTION read_loch_ion_data,debug=debug,load=load
	;file_path='/home/mlreinke/idl/impurities/data/adas/scd96_ar.dat'
	;save_path='/home/mlreinke/idl/impurities/data/adas/scd96_ar.sav'
	file_path='/usr/local/cmod/idl/atomic_physics/adas/scd96_ar.dat'
	save_path='/usr/local/cmod/idl/atomic_physics/adas/scd96_ar.sav'
	IF keyword_set(load) THEN BEGIN
		restore, save_path
		RETURN,rates
	ENDIF
	rates=read_scd_file(file_path)
	;save,rates,filename=save_path
	RETURN,rates
END
