FUNCTION mcp_pos,shot,jack=jack
	IF NOT keyword_set(jack) THEN BEGIN
		mdsopen, 'spectroscopy',(shot)
		jack=mdsvalue('\SPECTROSCOPY::TOP.VUV:JACK_POS')
	        mdsclose, 'spectroscopy',(shot)
	ENDIF

	ang=ACOS((72.375*72.375+87.5*87.5-(.9536*jack+27.457)*(.9536*jack+27.457))/(2.0*72.375*87.5))-.280392
	pivot=[3.738,.3315]
	;adjust z0 so that r0 = 1.1.  This will make it so that the default n_s will work
	zo=pivot[1]+tan(-1*ang)*(pivot[0]-1.1)

	output=[1.1,zo,0.0,ang]
	RETURN, output
END

FUNCTION mcp_genrad_fb_emiss,shot,t,rmaj,temp,dens,lam,Ip,z,debug=debug,gpv=gpv,nz=nz
	IF NOT keyword_set(nz) THEN nz=dens
	h=6.626068e-34	;m^2 kg/s
	e=1.60218e-19	;J/eV
	c=2.9979e8	;m/s
	;E=hc/lam

	pos=mcp_pos(shot)

	IF NOT keyword_set(zeff) THEN zeff=1.0
	emiss=3.3e-32*dens/1.0e6*dens/1.0e6*z^4*1.15/temp^(1.5)*exp(-1.0*(ang2ev(lam)-Ip)/temp)	;W/cm^3/eV/str from PoP v7 pg 4052
	emiss*=h*c/e/lam^2*1.0e6*4.0*!pi	;change to W/m^3/Ang
	emiss*=lam/(h*c)			;change to ph/s/m^3/Ang
	stop
	spec=lineint_spectrum(pos,emiss,rmaj,shot,t)
	
	output=spec
	IF keyword_set(debug) THEN stop
	RETURN,output
END


PRO kr_comp_spec,lam,spec,load=load,ylog=ylog,plot=plot
	save_path='/home/mlreinke/idl/mcpher/kr_comp_spec.dat'

	shot=1080116000+[9,10,11,14]
	t1=[1.0,0.90,0.88,1.0]
	t2=[1.2,0.95,0.91,1.2]

	IF NOT keyword_set(load) THEN BEGIN
		spec=fltarr(1024,n(shot)+1)
		lam=fltarr(1024,n(shot)+1)
		FOR i=0,n(shot) DO BEGIN
			mdsopen, 'spectroscopy',shot[i]
			cnts=mdsvalue('\SPECTROSCOPY::TOP.VUV.ANALYSIS:MCP_BRI_MOD',/quiet, status=status)
			t=mdsvalue('dim_of(\SPECTROSCOPY::TOP.VUV.ANALYSIS:MCP_BRI_MOD,0)',/quiet, status=status)
			lam_i = mdsvalue('dim_of(\SPECTROSCOPY::TOP.VUV.ANALYSIS:MCP_BRI_MOD,1)',/quiet, status=status)
			mdsclose, 'spectroscopy', shot[i]

			ilow=ipt(t,t1[i])
			ihigh=ipt(t,t2[i])
			spec[*,i]=sum_array(cnts[*,ilow:ihigh],/i)/(ihigh-ilow+1.0)
			lam[*,i]=lam_i
			save,spec,lam,filename=save_path
		ENDFOR
	ENDIF ELSE restore, save_path 	
	spec[*,3]*=1.015
	spec[*,2]*=0.92
	spec[*,1]*=0.8
	spec[*,0]*=0.62
	lam[*,0]-=0.25
	lam_low=[215.5,182,162,147.25]
	lam_high=[259,223,203,188.5]
	IF keyword_set(plot) THEN BEGIN
		IF keyword_set(ylog) THEN ylow=0.3 ELSE ylow=0
		plot, [0],[0],xr=[145,260],yr=[ylow,13],/xsty,xtit='Wavelength [Ang]',ytit='Brightness [AU]',ylog=ylog,/ysty
		color=[0,30,100,200]
		FOR i=0,n(shot) DO BEGIN
			tmp=where(lam[*,i] GE lam_low[i] AND lam[*,i] LE lam_high[i])
			oplot, lam[tmp,i],spec[tmp,i],color=color[i]
			IF i EQ 0 THEN BEGIN
				lamtot=lam[tmp,n(shot)-i]
				spectot=spec[tmp,n(shot)-i]
			ENDIF ELSE BEGIN
				tmp=where(lam[*,n(shot)-i] GT max(lamtot) AND lam[*,n(shot)-i] LE lam_high[n(shot)-i])
				lamtot=[lamtot,lam[tmp,n(shot)-i]]
				spectot=[spectot,spec[tmp,n(shot)-i]]
			ENDELSE
		ENDFOR
		;plot,lamtot,spectot,thick=3
		restore,save_path
		save,spec,lam,spectot,lamtot,filename=save_path
	ENDIF

	lam=lamtot
	spec=spectot
	
END

PRO mcp_genrad_profiles,shot,t1,t2,csden,cserr,temp,temperr,dens,denserr,rhovec,qmin=qmin,qmax=qmax,w=w,nz=nz,plotwin=plotwin,zeff=zeff,adas=adas,out=out
	IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0
	eta=2.89e20
	z=n(csden[*,0])
	IF NOT keyword_set(zeff) THEN zeff=1.0
	IF NOT keyword_set(nz) THEN nz=1.0e17
	IF NOT keyword_set(qmin) THEN BEGIN
		CASE z OF
			18 : qmin=10
			36 : qmin=0
		ENDCASE
	ENDIF
	IF NOT keyword_set(qmax) THEN BEGIN
		CASE z OF
			18 : qmax=z
			36 : qmax=z
		ENDCASE
	ENDIF
		
	plotwin+=1
	IF keyword_set(ps) THEN BEGIN
		xsize=7.0
		ysize=7.0*800/1400.0
		ls=1.25
	ENDIF ELSE BEGIN
		xsize=1400.0
		ysize=800.0
		ls=2.0
	ENDELSE
	IF NOT keyword_set(ps) THEN BEGIN
		device, window_state=var
		IF var[plotwin] EQ 0 THEN window,plotwin,xsize=xsize,ysize=ysize,xpos=1610,ypos=670,title='output profiles,'+num2str(plotwin) $
			ELSE wset,plotwin
	ENDIF ELSE BEGIN
		device, xsize=xsize, ysize=ysize, /inches
	ENDELSE

	;get efit_data
	mdsopen,'analysis',shot
	rmid=mdsvalue('\efit_rmid')
	efit_t=mdsvalue('dim_of(\efit_rmid)')
	mdsclose,'analysis',shot
	i1=ipt(efit_t,t1)
	i2=ipt(efit_t,t2)
	ro=min(sum_array(rmid[i1:i2,*],/j)/(i2-i1+1.0))
	a=max(sum_array(rmid[i1:i2,*],/j)/(i2-i1+1.0))-ro

	mdsopen, 'spectroscopy',shot
	cnts=mdsvalue('\SPECTROSCOPY::TOP.VUV.ANALYSIS:MCP_BRI_MOD',/quiet, status=status)
	t=mdsvalue('dim_of(\SPECTROSCOPY::TOP.VUV.ANALYSIS:MCP_BRI_MOD,0)',/quiet, status=status)
	lam = mdsvalue('dim_of(\SPECTROSCOPY::TOP.VUV.ANALYSIS:MCP_BRI_MOD,1)',/quiet, status=status)
	mdsclose, 'spectroscopy', shot
	cnts*=eta
	pos=mcp_pos(shot)

	i1=ipt(t,t1)
	i2=ipt(t,t2)
	mcpspec=sum_array(cnts[*,i1:i2],/i)/(i2-i1+1.0)
	maxloc=maxloc(mcpspec)
	wmin=2.0
	tmp=where(mcpspec GE 0.5*mcpspec[maxloc] AND lam GE lam[maxloc]-wmin/2.0 AND lam LE lam[maxloc]+wmin/2.0)
	w=(max(lam[tmp])-min(lam[tmp]))/(2.0*sqrt(2.0*alog(2)))
	;w=0.001
	specem=calc_specemiss(lam,csden*nz,temp,dens,qmin=qmin,qmax=qmax,w=w,adas=adas,thr=thr,cserr=cserr,temperr=temperr,denserr=denserr,emerr=emerr)
	thspec=lineint_spectrum(pos,specem,rhovec*a+ro,shot,0.5*(t1+t2))
	err_thspec=lineint_spectrum(pos,emerr,rhovec*a+ro,shot,0.5*(t1+t2))

	bremem=fltarr(n(lam)+1,n(rhovec)+1)
	FOR i=0,n(rhovec) DO BEGIN
		bremem[*,i]=1.89e-28*zeff*temp[i]^0.182*(dens[i]/1.0e6)*(dens[i]/1.0e6)/(sqrt(temp[i])*lam^2)*exp(-12.4000/(temp[i]*lam)) ;watts/cm^3/Ang
		bremem[*,i]*=lam/(1.989e-17)	;watts to ph/s at given wavelength
	ENDFOR
	bremem*=1.0e6	;ph/s/m^3/Ang
	brspec=lineint_spectrum(pos,bremem,rhovec*a+ro,shot,0.5*(t1+t2))
	;thspec+=brspec

	plotwin+=1
	openwin,plotwin
	tmp=where(mcpspec GT 0)
	order=sort(mcpspec[tmp])
	back=mean(mcpspec[tmp[order[0:4]]])
	;mcpspec-=back
	plot,lam,mcpspec/max(mcpspec),yr=[0,1.1],/ysty,xr=[min(lam[tmp]),max(lam[tmp])],/xsty,xtit='Wavelength [Ang]',ytit='Normalized Brigthness'
	;oplot,lam,thspec,color=200
	tmp=indgen(170)*6
	oplot,lam,thspec/max(thspec),color=200
	;oploterror,lam[tmp],thspec[tmp]/max(thspec),err_thspec[tmp]/max(thspec),color=200,errcolor=200,psym=3
	;oplot,lam,(thspec-err_thspec)/max(thspec),color=200,linestyle=2.0
	;oplot,lam,(thspec+err_thspec)/max(thspec),color=200,linestyle=2.0

	stop
	out={lam:lam,mcp:mcpspec,th:thspec}
END
