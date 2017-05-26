PRO run_gpfit,shot,t1,t2,etemp=etemp,edens=edens,output=output,config=config,debug=debug,remove=remove,bounds=bounds,inter=inter,starts=starts,all=all,xplot=xplot
	IF NOT keyword_set(etemp) AND NOT keyword_set(edens) THEN RETURN
	IF NOT keyword_set(starts) THEN starts=8

	IF NOT keyword_set(config) THEN BEGIN
		is_cts=1
		is_ets=1
		is_gpc=1
		is_gpc2=1
		xr=[0,1.05]
		xvar='r/a'
		tree='EFIT20'
  		config={is_cts:is_cts,is_ets:is_ets,is_gpc:is_gpc,is_gpc2:is_gpc2,xr:xr,xvar:xvar,tree:tree}
      ENDIF ELSE BEGIN
		is_cts=config.is_cts
		is_ets=config.is_ets
		is_gpc=config.is_gpc
		is_gpc2=config.is_gpc2
		xr=config.xr
		xvar=config.xvar
		tree=config.tree
	ENDELSE
	IF NOT keyword_set(xplot) THEN xplot=xr

	IF keyword_set(etemp) THEN BEGIN
		signal='Te'
		system=''
		IF is_ets THEN system+=' ETS'
		IF is_cts THEN system+=' CTS'
		IF is_gpc THEN system+=' GPC'
		IF is_gpc2 THEN system+=' GPC2'
	ENDIF

	IF keyword_set(edens) THEN BEGIN
		signal='ne'
		system=''
		IF is_ets THEN system+=' ETS'
		IF is_cts THEN system+=' CTS'
	ENDIF

	bndstr=''
	IF NOT keyword_set(bounds) THEN bounds=[0,14.5,0.1,0.6,0.05,0.3,0,0.2,0.85,1.05]	
	FOR i=0,n(bounds) DO bndstr+=' '+num2str(bounds[i],dp=2)

	IF NOT keyword_set(inter) THEN intstr=' --no-interaction' ELSE intstr=''
	
	rmvstr=''
	IF keyword_set(remove) THEN FOR i=0,n(remove) DO rmvstr+=' '+num2str(remove[i],1)

	IF NOT keyword_set(output) THEN output='/home/'+logname()+'/gpfit/gpfit_'+signal+'_'+num2str(shot,1)

	IF t1 EQ t2 THEN tstr=' -t '+num2str(t1,dp=4) ELSE tstr=' --t-min '+num2str(t1,dp=4)+' --t-max '+num2str(t2,dp=4)
	
	spawnstr='fit_profile --signal '+signal+' --shot '+num2str(shot,1)+tstr+' --system '+system+' --coordinate '+xvar+$
		' --EFIT-tree '+tree+' --bounds '+bndstr+' --x-min '+num2str(xr[0],dp=2)+' --x-max '+num2str(xr[1],dp=2)+' -o '+output+$
		' --plot-idxs --random-starts '+num2str(starts,1)+' --x-lim '+num2str(xplot[0],dp=2)+' '+num2str(xplot[1],dp=2)
	
	IF rmvstr NE '' THEN spawnstr+=' --remove-points '+rmvstr
	IF keyword_set(all) THEN spawnstr+=' --all-points'
	IF intstr NE '' THEN spawnstr+=intstr
	IF keyword_set(debug) THEN stop
	spawn,spawnstr,err

END
