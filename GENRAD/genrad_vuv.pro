PRO load_impspec_data,shot,elem,line,br,time,status=status,sig=sig
	path='\SPECTROSCOPY::TOP.IMPSPEC.'+strupcase(elem)+'.LINE'+num2str(line,1)+':BR'
	mdsopen,'spectroscopy',shot
	br=mdsvalue('_data='+path,/quiet,status=status)
	time=mdsvalue('dim_of(_data,0)',/quiet)
	sig=mdsvalue('dim_of(_data,1)',/quiet)
	mdsclose,'spectroscopy',shot
END

PRO genrad_vuv_molybdenum,shot,data,br,xeus=xeus,loweus=loweus,mcp=mcp,emiss=emiss,wave=wave,chg=chg
	ntime=n(data.time)+1
	nrho=n(data.rho)+1
	IF keyword_set(mcp) THEN RETURN
	IF NOT keyword_set(xeus) THEN loweus=1.0
	IF keyword_set(xeus) THEN pos=[2.561, .2158, 0.196,0.1136]
	IF keyword_set(loweus) THEN pos=[2.561,-.2158, 0.196,-0.1136]
	genpos_pos_reform,pos,[0.44,1.0,-0.6,0.6] 
	
	FOR i=0,ntime-1 DO BEGIN
		sxr_vuv_molybdenum,rotate(data.csden[*,*,i],4),data.temp[*,i],data.dens[*,i],wave,chg,iemiss
		IF i EQ 0 THEN BEGIN
			nlines=n(wave)+1
			emiss=fltarr(nrho,ntime,nlines)
                ENDIF
		FOR j=0,nlines-1 DO emiss[*,i,j]=iemiss[j,*]
        ENDFOR
	glines=[3,9,11]	
	ngood=n(glines)+1
	br=fltarr(ntime,nlines)
	FOR j=0,ngood-1 DO BEGIN
		i=glines[j]
		br[*,i]=genpos_line_brtau(pos,emiss[*,*,i],data.rmaj,data.time,shot,data.time,rmap=rmap)
        ENDFOR
END

PRO genrad_vuv_tungsten,shot,data,br,xeus=xeus,loweus=loweus,mcp=mcp,emiss=emiss,wave=wave,chg=chg

END
