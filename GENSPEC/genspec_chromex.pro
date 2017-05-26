FUNCTION chromex_kbot_pos,ch=ch
	vertex=[83.30,-48.72]/100.0
	ch=indgen(18)+1
	rpts=[56.80,54.91,50.20,fltarr(10)+46.83,46.65,46.25,45.0,44.1,44.05]/100.0
	zpts=-1.0*[50.80,49.68,48.36,46.72,45.08,43.43,41.77,40.08,38.35,36.58,34.75,32.86,30.9,28.75,26.34,23.20,20.01,17.2]/100.0
	n_ch=n(ch)+1
	pos=fltarr(4,n_ch)
	FOR i=0,n_ch-1 DO BEGIN
		rz2=[rpts[i],zpts[i]]
		rz1=vertex
		pos[[0,1],i]=rz1
		pos[[2,3],i]=genpos_cyl2tok(rz1,rz2,1.0)
	ENDFOR
	RETURN,pos
END

FUNCTION chromex_abot_pos,ch=ch
	vertex=[33.06,-75.79]/100.0
	ch=[7,8,9,10,11,12,13,14,15,16]
	rpts=[62.89,62.72,62.56,62.40,62.20,62.02,61.85,61.66,61.44,61.25]/100.0
	zpts=-1.0*[57.04,56.01,54.95,53.86,52.76,51.61,50.42,49.18,47.92,46.61]/100.0
	n_ch=n(ch)+1
	pos=fltarr(4,n_ch)
	FOR i=0,n_ch-1 DO BEGIN
		rz1=[rpts[i],zpts[i]]
		rz2=vertex
		pos[[0,1],i]=rz1
		pos[[2,3],i]=genpos_cyl2tok(rz1,rz2,1.0)
	ENDFOR
	RETURN,pos
		
END
FUNCTION chromex_ktop_pos,ch=ch
	vertex=[74.68,42.56]/100.0
	ch=indgen(12)+1
	rpts=[0.4611,0.4680,0.4683,0.4683,0.4838,0.5113,0.5374,0.5645,0.5792,0.6343,0.6589,0.6835,0.7081]
	zpts=-1.0*[0.2555,0.2906,0.3575,0.4403,0.4789,0.4686,0.4942,0.5027,0.5827,0.4137,0.3961,0.3768,0.3661]
	n_ch=n(ch)+1
	pos=fltarr(4,n_ch)
	FOR i=0,n_ch-1 DO BEGIN
		rz2=[rpts[i],zpts[i]]
		rz1=vertex
		pos[[0,1],i]=rz1
		pos[[2,3],i]=genpos_cyl2tok(rz1,rz2,1.0)
        ENDFOR
	coefs=poly_fit(ch,pos[3,*],1)
	pos[3,*]=coefs[1]*ch+coefs[0]
	RETURN,pos
		
END


FUNCTION chromex_aside_pos,ch=ch
	vertex=[126.49,23.81]/100.0
	ch=indgen(13)+1
	rpts=fltarr(13)+0.441
	zpts=[0.428,0.388,0.344,0.308,0.266,0.228,0.183,0.147,0.107,0.067,0.027,-0.0131,-0.0532]
	n_ch=n(ch)+1
	pos=fltarr(4,n_ch)
	FOR i=0,n_ch-1 DO BEGIN
		rz2=[rpts[i],zpts[i]]
		rz1=vertex
		pos[[0,1],i]=rz1
		pos[[2,3],i]=genpos_cyl2tok(rz1,rz2,1.0)
	ENDFOR
	RETURN,pos	
END

FUNCTION chromex_pos,shot=shot,kbot=kbot,ktop=ktop,abot=abot,aside=aside
	IF keyword_set(shot) THEN BEGIN
		mdsopen,'spectroscopy',shot
		peri=mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.SETTINGS:CFG_PERISCOP')
		fiber=mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.SETTINGS:CFG_PER_FIBR')
		mdsclose,'spectroscopy',shot
		tmp=where(peri EQ 'K_BOTTOM')
		IF tmp[0] NE -1 THEN kbot=fiber[tmp]
		tmp=where(peri EQ 'K_TOP')
		IF tmp[0] NE -1 THEN ktop=fiber[tmp]
		tmp=where(peri EQ 'A_BOTTOM')
		IF tmp[0] NE -1 THEN abot=fiber[tmp]
		tmp=where(peri EQ 'A_SIDE')
		IF tmp[0] NE -1 THEN aside=fiber[tmp]

        ENDIF 
	kb=chromex_kbot_pos(ch=kbch)
	kt=chromex_ktop_pos(ch=ktch)
	ab=chromex_abot_pos(ch=abch)
	as=chromex_aside_pos(ch=asch)
	pos=-1
	IF keyword_set(kbot) THEN BEGIN
		FOR i=0,n(kbot) DO BEGIN
			index=where(kbch EQ kbot[i])
			IF index[0] NE -1 THEN BEGIN
				IF pos[0] EQ -1 THEN pos=kb[*,index] ELSE pos=[[pos],[kb[*,index]]]
			ENDIF
                ENDFOR
	ENDIF
	IF keyword_set(ktop) THEN BEGIN
		FOR i=0,n(ktop) DO BEGIN
			index=where(ktch EQ ktop[i])
			IF index[0] NE -1 THEN BEGIN
				IF pos[0] EQ -1 THEN pos=kt[*,index] ELSE pos=[[pos],[kt[*,index]]]
			ENDIF
		ENDFOR
	ENDIF
	IF keyword_set(abot) THEN BEGIN
		FOR i=0,n(abot) DO BEGIN
			index=where(abch EQ abot[i])
			IF index[0] NE -1 THEN BEGIN
				IF pos[0] EQ -1 THEN pos=ab[*,index] ELSE pos=[[pos],[ab[*,index]]]
			ENDIF
		ENDFOR
        ENDIF
	IF keyword_set(aside) THEN BEGIN
		FOR i=0,n(aside) DO BEGIN
			index=where(asch EQ aside[i])
			IF index[0] NE -1 THEN BEGIN
				IF pos[0] EQ -1 THEN pos=as[*,index] ELSE pos=[[pos],[as[*,index]]]
			ENDIF
		ENDFOR
        ENDIF


	RETURN,pos
END

;R=74.68cm, z=42.56cm
;views terminate at roughly R=0.441m, z=[0.428 0.388 0.344 0.308 0.266 0.228 0.183 0.147 0.107 0.067 0.027 -0.0131 -0.0532];
;1-4(inner divertor),5-9(dome),10-13(upper outer divertor) 

;8/4/2014 - modified the /raw load to use the logic of W_SPEC_NG
PRO chromex_load_spec,shot,int,lam,ch,t,lab=lab,raw=raw,per=per,fib=fib
	mdsopen,'spectroscopy',shot
	IF NOT keyword_set(raw) THEN BEGIN
		int=mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.ANALYSIS:ANAL_SPECTRA')
		lam=mdsvalue('dim_of(\SPECTROSCOPY::TOP.CHROMEX.ANALYSIS:ANAL_SPECTRA,0)')
		ch=mdsvalue('dim_of(\SPECTROSCOPY::TOP.CHROMEX.ANALYSIS:ANAL_SPECTRA,1)')
		t=mdsvalue('dim_of(\SPECTROSCOPY::TOP.CHROMEX.ANALYSIS:ANAL_SPECTRA,2)')
	ENDIF ELSE BEGIN
		int=mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.ANALYSIS:RAW_SPECTRA')
		x=size(int)
		lam = mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.CALIBRATION:LAMBDA_TRUE')
		IF total(finite(lam)) EQ 0 OR max(lam EQ 1023) THEN BEGIN	;use W_SPEC_NG approach from GENIE/load_cmod.pro
			lam=fltarr(x[2],x[3])
	 		lambda_ctr=mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.SETTINGS:LAMBDA_CTR',/quiet) ; [nm]
			grating=mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.SETTINGS:GRATING',/quiet) ; [grooves/mm]
			CASE grating OF
				600  : disp = 360.0/grating ; [angstrom/pixel]
				1200 : disp = 390.0/grating ; [angstrom/pixel]
				1800 : disp = 375.0/grating ; [angstrom/pixel]
				ELSE : disp = 360.0/grating
	  		ENDCASE
			FOR i=0,x[3]-1 do lam[*,i] = (findgen(x[2]) - floor(x[2]/2.))*disp + lambda_ctr*10
			lam/=10.0	;convert to nm
		ENDIF
		dt=mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.CAMAC:DELTA_TIME')
		tstart=mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.CAMAC:TIME_START')
		t=tstart+dt*findgen(x[1])+dt/2
		;lam=mdsvalue('dim_of(\SPECTROSCOPY::TOP.CHROMEX.ANALYSIS:RAW_SPECTRA,0)')
		;ch=mdsvalue('dim_of(\SPECTROSCOPY::TOP.CHROMEX.ANALYSIS:RAW_SPECTRA,1)')
		;t=mdsvalue('dim_of(\SPECTROSCOPY::TOP.CHROMEX.ANALYSIS:RAW_SPECTRA,2)')
	ENDELSE
	per=mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.SETTINGS:CFG_PERISCOP')
	fib=mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.SETTINGS:CFG_PER_FIBR')
	mdsclose,'spectroscopy',shot
	lab=strarr(n(per)+1)
	FOR i=0,n(per) DO lab[i]=per[i]+'_'+num2str(fib[i],1)
END


PRO chromex_fitbalmer_back,xfit,yfit,xback,yback,i_pk=i_pk,lam_pk=lam_pk,i_low=i_low
	low_n=5
	high_n=13
	lam_n=make(low_n,high_n,high_n-low_n+1)
	lam_pk=ev2ang(13.6056*(1/4.0-1.0/lam_n^2))/10.0
	n_pk=n(lam_pk)+1
	i_pk=ivec(xfit,lam_pk)
	xback=fltarr(n_pk)
	yback=fltarr(n_pk)
	i_low=fltarr(n_pk)
	FOR i=1,n(lam_n) DO BEGIN
		i_low[i-1]=minloc(yfit[i_pk[i]:i_pk[i-1]])+i_pk[i]
		xback[i-1]=xfit[i_low[i-1]]
		yback[i-1]=yfit[i_low[i-1]]
	ENDFOR
	xback[n_pk-1]=xfit[2]
	yback[n_pk-1]=median(yfit[0:4])
END

FUNCTION chromex_despike,y,f=f
	IF NOT keyword_set(f) THEN f=0.5e-3
	n=n(y)+1
	y_filt=y
	FOR i=0,n-3 DO BEGIN
		sub=y_filt[i:i+2]
		IF sub[2]-median(sub) GT f THEN y_filt[i+2]=mean(sub[0:1])			
	ENDFOR
	RETURN,y_filt
END
			

FUNCTION chromex_findbr,int,lam,pk,dx=dx,plot=plot
	IF NOT keyword_set(dx) THEN dx=1.0
	cntr_max=10
	n_time=n(int[0,*])+1
	n_pk=n(pk)+1
	br=fltarr(n_time,n_pk)
	x=lam
	FOR i=0,n_time-1 DO BEGIN
		y=int[*,i]
		FOR j=0,n_pk-1 DO BEGIN
			ilow=ipt(x,pk[j]-dx/2.0)
			ihigh=ipt(x,pk[j]+dx/2.0)
			local_max=maxloc(y[ilow:ihigh])+ilow
			cntr=0
			br_tmp=fltarr(cntr_max)
			br_tmp[0]=y[local_max]
			WHILE cntr LT cntr_max-1 DO BEGIN
				br_tmp[cntr+1]=br_tmp[cntr]+y[local_max+cntr]+y[local_max-cntr]
				cntr+=1
			ENDWHILE
			stop
		ENDFOR
	ENDFOR
			
END

PRO chromex_aside_imp,shot
	chromex_load_spec,shot,int,lam,ch,t,lab=lab
	tmp=where(strlowcase(lab) EQ 'a_side_0')
	IF tmp[0] EQ -1 THEN RETURN
	int=rotate(int[*,*,tmp],4)
	lam=lam[*,tmp]
	pk_loc=[384.68,402.47]
	x=lam
	y=int[*,20]

	stop
END

PRO chromex_fitbalmer,shot,load=load

	i=19
	j=9
	save_path='/home/mlreinke/idl/chromex/tmp.dat'


	IF NOT keyword_set(load) THEN BEGIN
		chromex_load_spec,shot,int,lam,ch,t	
		save,int,lam,ch,t,filename=save_path
	ENDIF ELSE restore, save_path
	yfit=reform(int[i,*,j])
	xfit=lam[*,j]
	
	chromex_fitbalmer_back,xfit,yfit,xback,yback,lam_pk=lam_pk,i_pk=i_pk,i_low=i_low
	spline_p,xback,yback,xb_spl,yb_spl

	openwin,0
	plot,xfit,yfit,/xsty,/ylog
	oplot,xback,yback,color=200,psym=8
	oplot,xb_spl,yb_spl,color=200
	yfit_sub=yfit-interpol(yb_spl,xb_spl,xfit)

	
	fit_param=fltarr(3,n(i_low)-1)
	FOR i=1,n(i_low)-1 DO BEGIN
		x=xfit[i_low[i]:i_low[i-1]]
		y=yfit_sub[i_low[i]:i_low[i-1]]
		x_off=i_pk[i]
		weights=fltarr(n(x)+1)+1.0
		ini=[1.0,x[maxloc(y)],max(y)]
		out=curvefit(x,y,weights,ini,function_name='lorentz_fit',/noderiv)
		fit_param[*,i-1]=ini
		openwin,1
		plot,x,y
		oplot,x,out,color=200
		stop
	ENDFOR

	openwin,1
	plot, xfit,yfit,/xsty,xr=[370,390]
	oplot, xfit,yfit_sub,color=200	
	
END


PRO chromex_imp_profile,shot,br,ch,t,line=line,per=per,shift=shift
	IF NOT keyword_set(per) THEN per='abot'
	IF NOT keyword_set(line) THEN line='mo1'
	IF NOT keyword_set(shift) THEN shift=0.0

	CASE per OF 
		'abot' :BEGIN
			prefix='A'
			suffix='BOTTOM'
		END
		
		'kbot' :BEGIN
			prefix='K'
			suffix='BOTTOM'
		END
		
		'ktop' :BEGIN
			prefix='K'
			suffix='TOP'
		END
	ENDCASE

	CASE line OF

		'mo1' 	: 	dlam=[385.6,387.0]
		'w1'	: 	dlam=[400.4,401.3]
		'n2'	:	dlam=[399.0,400.1]
		'n3'	:	dlam=[437.4,438.4]
		'd7'	:	dlam=[394.8,398.65]
		'd5' 	:	dlam=[433.4,434.9]
		ELSE : dlam=-1

	ENDCASE
	IF dlam[0] EQ -1 THEN BEGIN
		print, 'no spectral line found'
		RETURN
	END

	CHROMEX_LOAD_SPEC,shot,int,lam,ch,t,lab=lab
	perch=where(strmatch(lab,prefix+'_'+suffix+'*') EQ 1)
	int=int[*,*,perch]
	lam=lam[*,perch]
	lab=lab[perch]
	nch=n(perch)+1
	ch=intarr(nch)
	FOR i=0,nch-1 DO BEGIN
		str=lab[i]
		str=strsplit(str,'_',/extract)
		ch[i]=int(str[2])
	ENDFOR
	
	ntime=n(t)+1
	br=fltarr(nch,ntime)
	order=sort(ch)
	ch=ch[order]
	
	FOR i=0,nch-1 DO BEGIN
		FOR j=0,ntime-1 DO BEGIN
			int_tmp=reform(int[j,*,i])
			lam_tmp=lam[*,i]+shift
			klow=ipt(lam_tmp,dlam[0])
			kup=ipt(lam_tmp,dlam[1])	
			br[order[i],j]=total(int_tmp[klow:kup])-0.5*(int_tmp[klow]+int_tmp[kup])
		ENDFOR
	ENDFOR
	

END


PRO chromex_mp564_plots,shot,tau,line=line,k=k,a=a,shift=shift
	load_xpt_data,shot,rxpt,zxpt,txpt
	chromex_imp_profile,shot,kbr,kch,t,line=line,per='kbot',shift=shift
	chromex_imp_profile,shot,abr,ach,t,line=line,per='abot',shift=shift
	kpos=chromex_kbot_pos(ch=kpch)
	apos=chromex_abot_pos(ch=apch)
	kpos=kpos[*,kch-1]
	apos=apos[*,ach-1-apch[0]]
	pos=[[kpos],[apos]]
	br=[kbr,abr]
	IF keyword_set(k) THEN BEGIN
		pos=kpos
		br=kbr
	ENDIF
	IF keyword_set(a) THEN BEGIN
		pos=apos
		br=abr
	ENDIF
	npos=n(pos[0,*])+1

	;line_path_plots,pos,/bspline,shot=shot,tpt=tau,/div,thick=10
	
	rx=rxpt[ipt(txpt,tau)]/100.0
	zx=zxpt[ipt(txpt,tau)]/100.0
	
	zpos=fltarr(npos)
	FOR i=0,npos-1 DO BEGIN
		s=line_s(pos[*,i],r=rx)
		zpos[i]=line_z(s[1],pos[*,i])
	ENDFOR

	stop


	;chromex_mp564_plots,1100903021,1.1,line='n2',/k
	;plot,(zpos-zx)*100.0,br[*,ipt(t,1.1)]*100.0,psym=-8,xtit='Z-Z!lX!n [cm]',ytit='Line Brightness [AU]',chars=1.2
	;chromex_mp564_plots,1100903021,1.1,line='n3',/k
	;oplot,(zpos-zx)*100.0,br[*,ipt(t,1.1)]*100.0,psym=-8,color=200
END
