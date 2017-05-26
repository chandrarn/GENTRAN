PRO read_gpfit,path,data,debug=debug
	data=-1
	id=ncdf_open(path)
	ncdf_attget,id,/global,'shot',shot
	ncdf_attget,id,/global,'t_min',t1
	ncdf_attget,id,/global,'t_max',t2
	data={shot:shot,t1:t1,t2:t2}
	info=ncdf_inquire(id)
	vars=strarr(info.nvars)
	FOR i=0,info.nvars-1 DO BEGIN	;read variable names to determine type of file
		inq=ncdf_varinq(id,i)
		vars[i]=inq.name
	ENDFOR	
	inq=ncdf_varinq(id,0)	;assume zeroth index is the radial since there's so many versions
	ncdf_varget,id,inq.name,rho
	rvar=inq.name
	pname='err'
	tmp=where(vars EQ 'T_e')
	IF tmp[0] NE -1 THEN BEGIN
		pname='temp'
		lname='alTe'
		ncdf_varget,id,'T_e',prof
		ncdf_varget,id,'err_T_e',err
		ncdf_varget,id,'a_LT_e',lprof
		ncdf_varget,id,'err_a_LT_e',lprof_err
		units='keV'
	ENDIF
	tmp=where(vars EQ 'n_e')
	IF tmp[0] NE -1 THEN BEGIN
		pname='dens'
		lname='alne'
		ncdf_varget,id,'n_e',prof
		ncdf_varget,id,'err_n_e',err
		ncdf_varget,id,'a_Ln_e',lprof
		ncdf_varget,id,'err_a_Ln_e',lprof_err
		units='10!u20!n [m!u-3!n]'
	ENDIF
	IF pname EQ 'err' THEN BEGIN
		stop
		print, 'ERROR - unknown profile stored in file'
		RETURN
	ENDIF
	data=create_struct(pname,prof,'rho',rho,'err',err,lname,lprof,'alerr',lprof_err,data,'radius',rvar,'units',units)
	IF keyword_set(debug) THEN stop
	ncdf_close,id
END

PRO load_gpfit_time,path,shot,npro,data
	
	;load initial profile
	read_gpfit,path+'gpfit_te_'+num2str(shot,1)+'_0',tedata
	nrho=n(tedata.rho)+1
	ntime=npro
	etemp=fltarr(nrho,ntime)
	edens=fltarr(nrho,ntime)
	eterr=fltarr(nrho,ntime)
	ederr=fltarr(nrho,ntime)
	time=fltarr(ntime)
	dt=fltarr(ntime)
	rho=tedata.rho
	
	FOR i=0,npro-1 DO BEGIN
		read_gpfit,path+'gpfit_te_'+num2str(shot,1)+'_'+num2str(i,1),tedata
		read_gpfit,path+'gpfit_ne_'+num2str(shot,1)+'_'+num2str(i,1),nedata
		time[i]=(tedata.t1+tedata.t2)/2
		dt[i]=(tedata.t2-tedata.t1)/2
		etemp[*,i]=tedata.temp
		eterr[*,i]=tedata.err
		edens[*,i]=nedata.dens
		ederr[*,i]=nedata.err

	ENDFOR
	data={etemp:etemp,edens:edens,rho:rho,time:time,dt:dt,eterr:eterr,ederr:ederr,shot:tedata.shot,radius:tedata.radius,etemp_units:tedata.units,edens_units:nedata.units}

END
