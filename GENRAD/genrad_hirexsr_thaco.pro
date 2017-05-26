;load emissivity, brightness and POS vectors

PRO hirexsr_load_thaco_bright,shot,group,br,pos,tau

	CASE group OF
		1 : BEGIN	;w+n>=3
			pospath='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.MOMENTS.W:POS'
			cpath='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.FITS.WN3:COEFS'
			labels='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.FITS.WN3:LABELS'
			lo=labels
		END
		2 : ;x+y
		3 : BEGIN ;qra
			pospath='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.MOMENTS.Z:POS'
			cpath='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.FITS.ZJK:COEFS'
			labels='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.FITS.ZJK:LABELS'
			lo=['q','r','a']
		END

		4 : BEGIN ;zjk
			pospath='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.MOMENTS.Z:POS'
			cpath='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.FITS.ZJK:COEFS'
			labels='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.FITS.ZJK:LABELS'
			lo=['z','j','k']
		END
		5 : lya1,2
		ELSE : RETURN
	ENDCASE

	mdsopen,'spectroscopy',shot
	pos=mdsvalue(pospath)
	coefs=mdsvalue(cpath)
	time=mdsvalue('dim_of('+cpath+',1)')
	labels=mdsvalue(labels)
	mdsclose,'spectroscopy',shot
	
	;select time range of interest
	tmp=where(time NE -1)
	tau=time[tmp]
	coefs=coefs[*,tmp,*]
	
	;select spatial coords of interest (assume only one CHMAP)
	tmp=where(pos[0,*,0] NE -1)
	pos=pos[*,tmp]
	coefs=coefs[tmp,*,*]

	ntime=n(tau)+1
	nch=n(pos[0,*])+1
	
	br=fltarr(nch,ntime)
	FOR i=0,n(lo) DO BEGIN
		index=where(labels EQ lo[i])
		FOR j=0,nch-1 DO BEGIN
			FOR k=0,ntime-1 DO br[j,k]+=coefs[j,k,index*3]*coefs[j,k,index*3+2]*sqrt(2.0*!pi)
		ENDFOR
	ENDFOR

END


;update of genrad to use THACO-based data 
PRO hirexsr_genrad_thaco_bright,shot,t1,t2,csden,cserr,temp,temperr,dens,denserr,rhovec,plotwin=plotwin,del_i=del_i,back=back,phase=phase,out=out
	IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0
	IF NOT keyword_set(plotwin) THEN plotwin=10

	CASE group OF 
		1 : ;wn3
		2 : ;xy
		3 : BEGIN ;qra
			wl_roi=[3.9780,3.9875]
		END
		4 : BEGIN ;zjk
			wl_roi=[3.9875,4.0000]
		END
		5 : BEGIN ;lya
			wl_roi=[3.7265,3.743]
		END
		ELSE : RETURN
	ENDCASE

	;calculate emissivity profiles from radiation modeling
	emiss=ar_xray_emiss(csden,dens,temp,wl_roi=wl_roi,cserr=cserr,sigte=temperr,signe=denserr,emerr=emerr,csemiss=csem)
	time=0.5*(t1+t2)

	;get efit_data
	mdsopen,'analysis',shot
	rmid=mdsvalue('\efit_rmid')
	efit_t=mdsvalue('dim_of(\efit_rmid)')
	mdsclose,'analysis',shot
	i1=ipt(efit_t,t1)
	i2=ipt(efit_t,t2)
	ro=min(sum_array(rmid[i1:i2,*],/j)/(i2-i1+1.0))
	a=max(sum_array(rmid[i1:i2,*],/j)/(i2-i1+1.0))-ro


	tit=[' ',' ','XCIS: q+r+a Line Group','XCIS: z+j+k Line Group','XCIS: Ly-'+n2g('alpha')+' Line Group']
	plotwin+=1
	IF keyword_set(ps) THEN BEGIN
		xsize=7.0
		ysize=7.0*800/1100.0
		ls=1.5
	ENDIF ELSE BEGIN
		xsize=1100.0
		ysize=800.0
		ls=2.0
	ENDELSE
	IF NOT keyword_set(ps) THEN BEGIN
		device, window_state=var
		IF var[plotwin] EQ 0 THEN window,plotwin,xsize=xsize,ysize=ysize,xpos=1610,ypos=670,title='output profiles,'+num2str(plotwin) $
			ELSE wset,plotwin
	ENDIF ELSE BEGIN
		d_old=!d
		device, xsize=xsize, ysize=ysize, /inches
	ENDELSE

	hirexsr_load_thaco_bright,shot,group[j],br,pos,tau
	tmp=where(tau GE t1 AND tau LE t2)
	brm=sum_array(br[*,tmp],/i)/(n(tmp)+1.0)
	brx=line_br(pos,emiss,rhovec*a+ro,[time],shot,time)
		
	stop
			
	IF keyword_set(ps) THEN device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm
END	

;+
;NAME:
;	HIREXSR_THACO_DELTALAM_EMISS
;
;PURPOSE:
;	This function will invert the brightness specificed by a wavelength interval rather than an specific line which is
;	already done by THACO
;
;CALLING SEQUENCE:
;	result=HIREXSR_THACO_DELTALAM_EMISS(shot,lr,line,rho)
;
;INPUTS:
;	shot	LONG	shot number
;	lr	FLTARR	[lam0,lam1] of the wavelength region to sum brightness [Ang]
;	line	INT	THACO line number (0=wn3, 1=xy, 2=zjk, 3=lya)
;	rho	FLTARR	[nrho] of the radial grid over which to do inversion [r/a]
;
;OPTIONAL INPUTS:
;	tr	FLTARR	[t1,t2] of the time interval to truncate, set t1=t2 to force and IPT call
;	eps	FLT	of the 2nd deriv smoothing [DEFAULT=1.0]
;	eta	FLT	of the edge zero weighting [DEFAULT=100.0]
;	good	INTARR	[nch] of 1's and 0's indicating the use of channel in the inversion
;	tht	INT	of the THACO tree number to use
;	etree	STRING	of the EFIT tree to use for equilibrium [DEFAULT: ANALYSIS]
;	jch	INT	of the bin to use for characteristic spectra [DEFAULT ~ on-axis]
;	
;KEYWORD PARAMETERS:
;	sine	/sine will add an m=1 sine component to the inversion
;
;OUTPUTS:
;	result	STRUC	containing a wide amount of pre/post-inversion data
;		*.em	 FLTARR [nrho,ntime] emissivity profile [X/m^3]
;		*.ems	 FLTARR [nrho,ntime] m=1 sine emissivity profile [X/m^3]
;		*.rho	 FLTARR [nrho] r/a radial scaling
;		*.time	 FLTARR [ntime]	time points [sec]		
;		*.emerr	 FLTARR [nrho,ntime] uncertainty in *.em [X/m^3]
;		*.emserr FLTARR [nrho,ntime] uncertainty in *.ems [X/m^3]
;		*.br	 FLTARR [nch,ntime] brightness profile [X/m^2]
;		*.brchk	 FLTARR [nch,ntime] brightness check from derived emissivity [X/m^2]	
;		*.brerr	 FLTARR [nch,ntime] uncertainty in *.br [X/m^2]
;		*.pos	 FLTARR [4,nch] of the POS vectors
;		*.good	 INTARR [nch] of the channels selected for inversion
;		*.lam	 FLTARR	[nlam] of the wavelenths for charac.spectra [Ang]
;		*.specbr FLTARR	[nlam] of the spectral brightness [X/m^2/Ang]
;		*.sig 	 FLTARR	[nlam] of the uncertainty in *.specbr [X/m^2/Ang]
;
;PROCEDURE:
;
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke (4/18/12) adapated from a variety of THACO tools
;
;-

FUNCTION hirexsr_thaco_deltalam_emiss,shot,lr,line,rho,tr=tr,eps=eps,eta=eta,good=good,tht=tht,etree=etree,sine=sine,jch=jch

	IF min(lr) LT 3.93 THEN h=1 ELSE h=0	;set h=1 (h-like spectra) 

	;load HIREXSR data
	hirexsr_load_avespec,shot,specbr,lam,sig,resid,tau,nave,h=h,tht=tht
	hirexsr_load_line_pos,shot,line,pos,tpos=tpos,tht=tht
	hirexsr_load_fits,shot,line,coefs,ave,double,labels,tht=tht

	;determine what baseline was used in the fit
	ncoefs=n(coefs[0,0,*])+1
	CASE ncoefs MOD 3 OF
		0 : BEGIN
			IF ncoefs EQ 3 THEN BEGIN
				n_line=1
				basen=3
			ENDIF ELSE BEGIN		;quadradic baseline
				n_line=ncoefs/3-1
				basen=2
			ENDELSE
		END
		1 : BEGIN		;constant baseline
			n_line=ncoefs/3
			basen=0
		END

		2 : BEGIN		;linear baseline
			n_line=ncoefs/3
			basen=1	
		END
	ENDCASE

	IF NOT keyword_set(tr) THEN tr=[tau[0],max(tau)]	;calculate emissivity for all tau's				
	IF tr[0] EQ tr[1] THEN tpts=ipt(tau,tr[0]) ELSE tpts=where(tau GE tr[0] AND tau LE tr[1])

	IF NOT keyword_set(eps) THEN eps=1.0
	IF NOT keyword_set(eta) THEN eta=1.0e2
	time=tau[tpts]
	ntime=n(tpts)+1
	nrho=n(rho)+1

	dx=0.01		;hard coded for zjk which is the only one with non-constannt baseline
	x0=3.94		;hard coded for zjk
	FOR k=0,ntime-1 DO BEGIN
		i=tpts[k]
		pindex=last(where(tpos LE tau[i]))	;find the proper POS index
		ipos=pos[*,*,pindex]
		tmp=where(ipos[2,*] GE 0)
		nch=n(tmp)+1
		IF k EQ 0 THEN BEGIN			;initialize variables (probably will break if tpos > 0)
			bright=fltarr(nch,ntime)
			brerr=fltarr(nch,ntime)
			brchk=fltarr(nch,ntime)
			emiss=fltarr(nrho,ntime)
			emerr=fltarr(nrho,ntime)
			emsine=fltarr(nrho,ntime)
			emserr=fltarr(nrho,ntime)
			IF NOT keyword_set(jch) AND h EQ 0 THEN jch=nch/2
			IF NOT keyword_set(jch) AND h EQ 1 THEN jch=1
		ENDIF

		FOR j=0,nch-1 DO BEGIN
			order=sort(lam[j,i,*])
			x=lam[j,i,order]
			y=specbr[j,i,order]
			ysig=sig[j,i,order]
			
			;remove baseline from COEFS
			p=reform(coefs[j,i,*])
			CASE basen OF 
				3 : base=0
				2 : base=p[n(p)]+p[n(p)-1]*(x-x0)/dx+p[n(p)-2]*(x-x0)^2/dx^2	;
				1 : base=p[n(p)]+p[n(p)-1]*(x-x0)/dx				; 
				0 : base=p[n(p)]
				ELSE :	
                        ENDCASE
			y-=base
			
			;truncate dataset to delta lambda
			lpts=where(x GE lr[0] AND x LE lr[1])
			x=x[lpts]
			y=y[lpts]
			ysig=sqrt(ysig[lpts]^2+min(ysig[lpts])^2)	;adds baseline subtraction error

			IF h EQ 1 THEN BEGIN
 				m=where(labels EQ '4d')
				pmo=p[3*m:3*m+2]
				ymo=pmo[0]*exp(-(x-pmo[1])^2/(2.0*pmo[2]^2))
				y-=ymo
			ENDIF

			;calculated moments
			bright[j,k]=int_tabulated(x,y,/double)
			FOR m=0,n(lpts)-1 DO brerr[j,k]+=(0.5*(x[m+1]-x[m]))^2*ysig[m]^2+(0.5*(x[m+1]-x[m]))^2*ysig[m+1]^2
			IF k EQ 0 AND j EQ jch THEN BEGIN
				jlam=lam[j,i,order]
				jspecbr=specbr[j,i,order]
				jerr=sig[j,i,order]
			ENDIF
                ENDFOR
		brerr=sqrt(brerr)
		
		iu=fltarr(n(tmp)+1)+4.0*!pi			;this puts voxel in units of length so brightness can be inverted rather than power
		ivox=genpos_pos2voxel_matrix(ipos[*,tmp],iu,shot,tpts=tau[i],rho_vec=rho,r_ap=r_ap,n_s=n_s,tree=etree)
		IF keyword_set(sine) THEN im1svox=genpos_pos2voxel_matrix(ipos[*,tmp],iu,shot,tpts=tau[i],rho_vec=rho,r_ap=r_ap,n_s=n_s,m=1,/sine,tree=etree)

		;combine voxel matrices
		IF NOT keyword_set(good) THEN BEGIN
			good=intarr(nch)+1
                ENDIF
		chk=where(good EQ 1)
		a=max(ivox)
		IF keyword_set(sine) THEN BEGIN
			ivoxtot=[[ivox[chk,*]],[im1svox[chk,*]]] 
			nprof=2
                ENDIF ELSE BEGIN 
			ivoxtot=ivox[chk,*]
			nprof=1
                ENDELSE

		;invert data
		out=genpos_profile_invert(bright[chk,k],ivoxtot,double(eps[0]*a^2),brchk=br0,eta=double(eta[0]*a^2),err=brerr[chk,k],nprof=nprof)
		emiss[*,k]=out[0:nrho-1]
		emerr[*,k]=br0.inverr[0:nrho-1]
		brchk[chk,k]=br0.mom
		IF keyword_set(sine) THEN BEGIN
			emsine[*,k]=out[nrho:*]
			emserr[*,k]=br0.inverr[nrho:*]
		ENDIF
        ENDFOR
	
	output={em:emiss,ems:emsine,rho:rho,time:time,emerr:emerr,emserr:emserr,br:bright,brchk:brchk,brerr:brerr,pos:ipos[*,tmp],good:good,lam:jlam,specbr:jspecbr,sig:jerr}
	IF keyword_set(debug) THEN stop
	RETURN,output
END	

PRO hirexsr_genrad_thaco_profiles,shot,t1,t2,csden,cserr,temp,temperr,dens,denserr,rhovec,plotwin=plotwin,tht=tht,hgood=hgood,hegood=hegood
	IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0

	IF NOT keyword_set(plotwin) THEN plotwin=10
	;setup ROI's for emissivity profiles
	wl_roi=[[3.9875,3.9985],$	;k+j+z
		[3.9780,3.9875],$	;q,r,a
		[3.7265,3.7405]]	;lya1+lya2+mo	
	h_roi=[0,0,1]
	line=[2,2,3]
	normi=[1,1,1]
	normval=fltarr(3)
	tit=['IXCS: z+j+k Line Group','IXCS: q+r+a Line Group','IXCS: Ly-'+n2g('alpha')+' Line Group']
	label=['z','z','lya1']
	nspec=n(wl_roi[0,*])+1

	;get efit_data
	mdsopen,'analysis',shot
	rmid=mdsvalue('\efit_rmid')
	efit_t=mdsvalue('dim_of(\efit_rmid)')
	mdsclose,'analysis',shot
	i1=ipt(efit_t,t1)
	i2=ipt(efit_t,t2)
	ro=min(sum_array(rmid[i1:i2,*],/j)/(i2-i1+1.0))
	a=max(sum_array(rmid[i1:i2,*],/j)/(i2-i1+1.0))-ro
	
	;select time points
	IF NOT keyword_set(phase) THEN BEGIN
		t=(t1+t2)/2.0
		tr=[t,t]
	ENDIF ELSE BEGIN
		;hirexsr_sawtooth_phase,shot,t1,t2,ph,ph_time,dph=0.65
		;ipt=ipt(ph,phase)
		;IF ipt[0] EQ -1 THEN BEGIN
		;	IF phase LT min(ph) THEN ipt=0
		;	IF phase GT max(ph) THEN ipt=n(ph)
		;ENDIF
		;t=ph_time[ipt]
		;tr=[t,t]
		;del_i=1
	ENDELSE

	;calc predicted and measured emissivity/brightness for each wl_roi and plot
	FOR j=0,nspec-1 DO BEGIN
		hirexsr_load_line_pos,shot,line[j],pos,tpos=tpos,tht=tht
		nch=n(where(pos[2,*] GT 0))+1
		jch=0
		IF h_roi[j] THEN BEGIN 		;setup for H-like inversion
			rho=make(0,1.0,25)
			sine=0
			IF NOT keyword_set(hgood) THEN BEGIN
				good=intarr(nch)+1
			ENDIF ELSE good=hgood
			chnorm=[3,4,5]
			q=17
                ENDIF ELSE BEGIN		;setup for He-like inversions
			rho=make(0,1.0,40)
			sine=1
			IF NOT keyword_set(hegood) THEN BEGIN
				good=intarr(nch)+1
				npm=nch/24
				good[[1*npm-[1,2,3],2*npm-[1,2,3],8*npm-[1,2,3],14*npm-[1,2,3]]]=0
                        ENDIF ELSE good=hegood
			chnorm=12*npm+[-1,-2,0,1,2]
			q=16
		ENDELSE

		;calculate emissivity and brightness profiles from experiment
		exp=hirexsr_thaco_deltalam_emiss(shot,wl_roi[*,j],line[j],rho,tr=tr,eps=eps,eta=eta,good=good,tht=tht,etree=etree,sine=sine,jch=jch)		

		;calculate emissivity and brightness profiles from simulations
		emiss=ar_xray_emiss(csden,dens,temp,wl_roi=wl_roi[*,j],cserr=cserr,sigte=temperr,signe=denserr,emerr=emerr,csemiss=csem)
		xpos=exp.pos
		genpos_pos_reform,xpos,[0.44,1.0,-0.6,0.6]
		br=genpos_line_br(xpos,emiss,rhovec,[0.5*(t1+t2)],shot,0.5*(t1+t2),emerr=emerr,/rho,brerr=brerr)
		th={em:emiss,emerr:emerr,rho:rhovec,csem:csem,time:[0.5*(t1+t2)],br:br,brerr:brerr}

		plotwin+=1
		IF keyword_set(ps) THEN BEGIN
			xsize=7.0
			ysize=7.0*800/1100.0
			ls=1.5
		ENDIF ELSE BEGIN
			xsize=1100.0
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
		pos1=[0.075,0.075,0.6,0.95]
		pos2=[0.075,0.1,0.6,0.3]
		pos3=[0.7,0.5,0.975,0.975]
		pos4=[0.7,0.35,0.975,0.5]
		pos5=[0.65,0.1,0.975,0.275]

		IF normi[j] EQ 1 THEN BEGIN
			norm=mean(exp.br[chnorm]/th.br[chnorm])
			normval[j]=norm
			;expr=string(mean(exp.em))
			;expr=strsplit(expr,'+',/extract)
			;expr=expr[1]
			;pow=10.0^float(expr)
			pow=1.0
		ENDIF	
		th.br*=norm
		th.brerr*=norm
	
		;emissivity comparison
		maxpt=max(th.em*norm/pow) > max(exp.em/pow)
		lab=num2str(shot,1)+' t='+num2str(t,dp=2)
		IF normi[j] THEN normstr='NORMALIZATION: '+num2str(norm,dp=2) ELSE normstr=' '
		plot, [0],[0],pos=pos1,yr=[0,maxpt*1.05],ytit='Emissivity [AU]',/ysty,xr=[0,max(rhovec)],/xsty,chars=0.5*ls,$
			tit=tit[j],xtit='r/a'
		oploterror,rhovec,th.em*norm/pow,th.emerr*norm/pow,color=200,errcolor=200
		oplot, rhovec,csem[q,*]*norm/pow,color=30
		oplot, rhovec,csem[q-1,*]*norm/pow,color=30,linestyle=2.0
		oplot, rhovec,csem[q+1,*]*norm/pow,color=30,linestyle=4.0
		oplot, [0.9,1.0],0.6*maxpt*[1.0,1.0],color=30,linestyle=0
		oplot, [0.9,1.0],0.4*maxpt*[1.0,1.0],color=30,linestyle=2
		oplot, [0.9,1.0],0.2*maxpt*[1.0,1.0],color=30,linestyle=4
		xyouts,0.91,0.62*maxpt,'EXC',color=30
		xyouts,0.91,0.42*maxpt,'ION',color=30
		xyouts,0.91,0.22*maxpt,'REC',color=30

		oploterror,exp.rho,exp.em/pow,exp.emerr/pow
		xyouts,1.0,maxpt*0.7,normstr,chars=0.4*ls,orient=90

		;emissivity residual
;		residual=genspec.em[*,i]/pow-interpol(emiss[*,i]*norm/pow,rhovec,genspec.rho)
;		residual_err=sqrt((genspec.emerr[*,i]/pow)^2+(interpol(emerr[*,i]*norm/pow,rhovec,genspec.rho))^2)
;		max=max(residual) > 0
;		min=min(residual) < 0
;		plot, [0],[0],pos=pos2,/noerase,yr=[min,max]*1.05,xtit='r/a',ytit=n2g('Delta'),xr=[0,max(rhovec)],/xsty,/ysty,chars=0.5*ls
;		oploterror,genspec.rho,residual,residual_err,psym=8,symsize=0.5*ls
;		oplot,[0,1],[0,0],linestyle=2.0

		;brightness
		tmp=where(exp.good EQ 1)
		ch=indgen(n(exp.good)+1)
		maxpt=max(exp.br[tmp]+exp.brerr[tmp]) > max(th.br+th.brerr)
		plot, [0],[0],pos=pos3,/noerase,yr=[0,maxpt*1.05],/ysty,xr=[0,max(ch)],/xsty,ytit='Counts',chars=0.5*ls,xtickname=replicate(' ',10)
		makesym,9
		oploterror,ch,exp.br,exp.brerr,psym=8,symsize=0.5*ls
		makesym,10
		oplot,ch[tmp],exp.br[tmp],psym=8,symsize=0.5*ls
		oplot,ch[tmp],exp.brchk[tmp],color=100
		oploterror,ch,th.br,th.brerr,color=200,errcolor=200

		;brightness residual
		fit_residual=exp.br-exp.brchk
		chk_residual=exp.br-th.br
		err_residual=sqrt(exp.brerr^2+th.brerr^2)
		max=max(fit_residual[tmp]+err_residual[tmp]) > 0 > max(chk_residual[tmp]+err_residual[tmp])
		min=min(fit_residual[tmp]-err_residual[tmp]) < 0 < min(chk_residual[tmp]-err_residual[tmp])
		plot, [0],[0],pos=pos4,/noerase,xr=[0,max(ch)],/xsty,yr=[min,max]*1.05,/ysty,xtit='CH #',ytit=n2g('Delta'),chars=0.5*ls
		oploterror,ch[tmp],fit_residual[tmp],exp.brerr[tmp],psym=8,color=100,errcolor=100,symsize=0.5*ls
		oploterror,ch[tmp],chk_residual[tmp],err_residual[tmp],psym=8,color=200,errcolor=200,symsize=0.5*ls
		oplot,[0,50],[0,0],linestyle=2.0
		xyouts,max(ch)*1.05,min*0.5,lab,chars=0.3*ls,orient=90

		;spectra
		xr=[min(exp.lam),max(exp.lam)]
		plot, [0],[0],pos=pos5,/noerase,xr=xr,/xsty,yr=[1,max(exp.specbr)]*1.1,/ylog,chars=0.5*ls,xtit=n2g('lambda')+' [mAng]'
		ct=[1,1,12,12]
		col=[200,200,120,120]
		tmp=where(exp.lam GE wl_roi[0,j] AND exp.lam LE wl_roi[1,j])
		fillx=[wl_roi[0,j],exp.lam[tmp],wl_roi[1,j],wl_roi[0,j]]
		filly=[1,exp.specbr[tmp],1,1]
		tmp=where(filly LT 1)
		IF tmp[0] NE -1 THEN filly[tmp]=1
		loadct,ct[0],/s
		polyfill,fillx,filly,color=col[0]
		loadct,12,/s
		oplot,exp.lam,exp.specbr
	ENDFOR
END

;+
;NAME:
;	HIREXSR_GENRAD_MOLY
;
;PURPOSE:
;	This procedure loads experimental Mo data from HIREXSR and computes
;	the equivilent using data from an impurity transport simulation
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 1/18/13
;	1/28/15		M.L. Reinke - modifed to allow time-evolving simulated data to be computed if data.time > 1 element
;
;-
PRO hirexsr_genrad_moly,shot,t1,t2,csden,cserr,temp,temperr,dens,denserr,rhovec,tht=tht,line=line,eptr=eptr,efit=efit,sptr=sptr,data=data,nobr=nobr,index=index
	print, 'starting'
	IF NOT keyword_set(index) THEN index=0
	IF keyword_set(data) THEN BEGIN
		csden=data.csden
		cserr=data.cserr
		temp=data.temp
		temperr=data.terr
		dens=data.dens
		denserr=data.derr
		rhovec=data.rho
		time=data.time
		x=size(rhovec)
		y=size(csden)
		IF y[1] NE x[1] THEN BEGIN
			csden=rotate(csden,4)
			cserr=rotate(cserr,4)
		ENDIF
	ENDIF
	IF keyword_set(index) THEN BEGIN
		csden=data.csden[*,*,index]
		cserr=data.cserr[*,*,index]
		temp=data.temp[*,index]
		temperr=data.terr[*,index]
		dens=data.dens[*,index]
		denserr=data.derr[*,index]
		time=time[index]
	ENDIF
	IF NOT keyword_set(tht) THEN tht=0
	IF NOT keyword_set(line) THEN line=[4]
	IF shot GT 1120821000 AND shot LT 1121003000 THEN BEGIN		;allow for high-Te configuration based on shot range
		line=8
	ENDIF
	nrho=n(rhovec)+1
	nlines=n(line)+1
	ntime=n(time)+1

	tit='XICS: Mo!u32+!n'
	label='4d'

	;get efit_data if not included as optional input
	IF NOT keyword_set(efit) THEN BEGIN
		mdsopen,'analysis',shot
		rmid=mdsvalue('\efit_rmid')
		efit_t=mdsvalue('dim_of(\efit_rmid)')
		ipsin=mdsvalue('dim_of(\efit_rmid,1)')
		mdsclose,'analysis',shot
		i1=ipt(efit_t,t1)
		i2=ipt(efit_t,t2)
		ro=min(sum_array(rmid[i1:i2,*],/j)/(i2-i1+1.0))
		a=max(sum_array(rmid[i1:i2,*],/j)/(i2-i1+1.0))-ro
		irmid=reform(sum_array(rmid[i1:i2,*],/j))/(i2-i1+1.0)
		irho=(irmid-ro)/a
		efit={rmid:irmid,rho:irho,psin:ipsin,ro:ro,a:a}
		;print, 'EFIT data loaded'
        ENDIF
	
	
        ;load emissivity and brightness profiles from experiment if not included as optional inputs
	IF NOT keyword_set(eptr) THEN BEGIN
		eptr=ptrarr(nlines,/allocate_heap)
		hirexsr_load_line_pos,shot,line[0],pos,tpos=tpos,tht=tht		;load POS vector
		tmp=where(pos[1,*] NE -1)
		pos=pos[*,tmp]

		hirexsr_load_mlintptr,shot,line,lint,tau,tht=tht			;load brightness profile
		i1=ipt(tau,t1)
		i2=ipt(tau,t2)
		ilint=*lint[i1]
		bright=ilint[*,8]
		brerr=ilint[*,9]^2
		FOR i=i1+1,i2 DO BEGIN
			ilint=*lint[i]
			bright+=ilint[*,8]
			brerr+=ilint[*,9]^2
                ENDFOR
		heap_free,lint
		bright/=(i2-i1+1.0)
		brerr=sqrt(brerr)/(i2-i1+1.0)

		hirexsr_load_profile,shot,line[0],prof,err,psin,tau,tht=tht	;load emissivity profile
		i1=ipt(tau,t1)
		i2=ipt(tau,t2)
	
		FOR i=i1,i2 DO BEGIN		;check to make sure radial grid is constant over time
			chk=total(psin[*,i]-psin[*,i1])
			IF chk NE 0 THEN BEGIN
				print, 'time-evolving radial grid for THT='+num2str(tht,1)+', - incompatible'
				RETURN
	                ENDIF
		ENDFOR
	
		;check for m=1 term, isolate m=0 only
		npsi=n(psin[*,i1])+1
		IF psin[0] EQ psin[npsi/2] THEN ism=1 ELSE ism=0
		CASE ism OF
			0 : BEGIN
				emiss=prof[*,i1:i2,0]
				emerr=err[*,i1:i2,0]
				emrho=interpol(irho,ipsin,psin[*,i1])		;convert to r/a grid
                	END
			1 : BEGIN
				emiss=prof[0:npsi/2-1,i1:i2,0]
				emerr=err[0:npsi/2-1,i1:i2,0]
				emrho=interpol(irho,ipsin,psin[0:npsi/2-1,i1])	;convert to r/a grid	
	                END
		ENDCASE
		emiss=sum_array(emiss,/i)/(i2-i1+1.0)
		emerr=sqrt(sum_array(emerr^2,/i))/(i2-i1+1.0)
		ied={emiss:emiss,emerr:emerr,rho:emrho,bright:bright,brerr:brerr,pos:pos,shot:shot,t1:t1,t2:t2,line:line,tht:tht,tit:tit,label:label,enorm:0.0,bnorm:0.0}
		*eptr[0]=ied
		;print, 'Measured Emissivity and Brightness Profiles Loaded'
        ENDIF

	;calculate emissivity and brightness profiles from simulations
	;sxr_vuv_molybdenum,csden,temp,dens,wave,chrg,emsim
	;index=where(wave EQ 3.74 AND chrg EQ 32)		;isolate the relevant line
	;emsim=reform(emsim[index,*])
	IF NOT keyword_set(sptr) THEN sptr=ptrarr(nlines,/allocate_heap)
	IF size(*sptr[0],/type) NE 0 THEN BEGIN						;restore from sdata
		sdata=*sptr[0]
		pos=sdata.pos
        ENDIF ELSE BEGIN
		hirexsr_load_line_pos,shot,line[0],pos,tpos=tpos,tht=tht		;load POS vector
		tmp=where(pos[1,*] NE -1)
		pos=pos[*,tmp]
	ENDELSE
	FOR k=0,ntime-1 DO BEGIN
		IF k EQ 0 THEN BEGIN
			emsim=fltarr(nrho,ntime)
			eserr=fltarr(nrho,ntime)
		ENDIF
		emsim[*,k]=nelike_mo_4d(rotate(csden[*,*,k],4),temp[*,k],dens[*,k])
		eserr[*,k]=emsim[*,k]*0.0						;no error propigation in SXR_VUV_MOLYBDENUM right now
	ENDFOR
	IF keyword_set(nobr) THEN BEGIN
		nch=n(pos[0,*])+1
		brsim=fltarr(nch)+1.0
		bserr=fltarr(nch)
		print, 'Simulated Brightness Profile Skipped'
        ENDIF ELSE BEGIN
		xpos=pos
		genpos_pos_reform,xpos,[0.44,1.0,-0.6,0.6]
		IF ntime NE 1 THEN brsim=genpos_line_brtau(xpos,emsim,rhovec,time,shot,time,brerr=bserr,/rho) ELSE brsim=genpos_line_br(xpos,emsim,rhovec,time,shot,time,/rho,brerr=bserr)
		print, 'Simulated Brightness Profile Computed'
	ENDELSE
	isd={emiss:emsim,emerr:eserr,rho:rhovec,bright:brsim,brerr:bserr,pos:pos,shot:shot,time:time,t1:t1,t2:t2,enorm:0.0,bnorm:0.0}
	*sptr[0]=isd

END

;+
;NAME:
;	HIREXSR_GENRAD_ARGON
;
;PURPOSE:
;	This procedure loads experimental Mo data from HIREXSR and computes
;	the equivilent using data from an impurity transport simulation
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 1/18/13
;	1/26/15		M.L. Reinke - modifed to allow time-evolving simulated data to be computed if data.time > 1 element
;
;-

PRO hirexsr_genrad_argon,shot,t1,t2,csden,cserr,temp,temperr,dens,denserr,rhovec,tht=tht,line=line,eptr=eptr,efit=efit,sptr=sptr,data=data,nobr=nobr
	IF keyword_set(data) THEN BEGIN
		csden=data.csden
		cserr=data.cserr
		temp=data.temp
		temperr=data.terr
		dens=data.dens
		denserr=data.derr
		rhovec=data.rho
		time=data.time
		x=size(rhovec)
		y=size(csden)
		IF y[1] NE x[1] THEN BEGIN		;assume data.csden is [space, cs, time] (STRAHL output)
			csden=rotate(csden,4)		;if GENTRAN, rotate from [cs,space] to [space,cs]
			cserr=rotate(cserr,4)
		ENDIF
	ENDIF
	IF NOT keyword_set(tht) THEN tht=0
	IF NOT keyword_set(line) THEN line=[2,3]
	nrho=n(rhovec)+1
	nlines=n(line)+1
	ntime=n(time)+1


	;get efit_data if not included as optional input
	IF NOT keyword_set(efit) THEN BEGIN
		mdsopen,'analysis',shot
		rmid=mdsvalue('\efit_rmid')
		efit_t=mdsvalue('dim_of(\efit_rmid)')
		ipsin=mdsvalue('dim_of(\efit_rmid,1)')
		mdsclose,'analysis',shot
		i1=ipt(efit_t,t1)
		i2=ipt(efit_t,t2)
		ro=min(sum_array(rmid[i1:i2,*],/j)/(i2-i1+1.0))
		a=max(sum_array(rmid[i1:i2,*],/j)/(i2-i1+1.0))-ro
		irmid=reform(sum_array(rmid[i1:i2,*],/j))/(i2-i1+1.0)
		irho=(irmid-ro)/a
		efit={rmid:irmid,rho:irho,psin:ipsin,ro:ro,a:a}
        ENDIF 
	
	
        ;load emissivity and brightness profiles from experiment if not included as optional inputs
	IF NOT keyword_set(eptr) THEN BEGIN
		eptr=ptrarr(nlines,/allocate_heap)
		FOR j=0,nlines-1 DO BEGIN
			hirexsr_load_line_pos,shot,line[j],pos,tpos=tpos,tht=tht		;load POS vector
			tmp=where(pos[1,*] NE -1)
			pos=pos[*,tmp]
			CASE line[j] OF 
				2 : BEGIN
					tit='XICS: Ar!u16+!n'
					label='z'
				END
				3 : BEGIN
					tit='XICS: Ar!u17+!n'
					label='lya'
				END
			ENDCASE


			hirexsr_load_mlintptr,shot,line[j],lint,tau,tht=tht		;load brightness profile
			i1=ipt(tau,t1)
			i2=ipt(tau,t2)
			ilint=*lint[i1]
			bright=ilint[*,8]
			brerr=ilint[*,9]^2
			FOR i=i1+1,i2 DO BEGIN
				ilint=*lint[i]
				bright+=ilint[*,8]
				brerr+=ilint[*,9]^2
                	ENDFOR
			heap_free,lint
			bright/=(i2-i1+1.0)
			brerr=sqrt(brerr)/(i2-i1+1.0)

			hirexsr_load_profile,shot,line[j],prof,err,psin,tau,tht=tht	;load emissivity profile
			i1=ipt(tau,t1)
			i2=ipt(tau,t2)
		
			FOR i=i1,i2 DO BEGIN		;check to make sure radial grid is constant over time
				chk=total(psin[*,i]-psin[*,i1])
				IF chk NE 0 THEN BEGIN
					print, 'time-evolving radial grid for THT='+num2str(tht,1)+', - incompatible'
					RETURN
		                	ENDIF
			ENDFOR
	
			;check for m=1 term, isolate m=0 only
			npsi=n(psin[*,i1])+1
			IF psin[0] EQ psin[npsi/2] THEN ism=1 ELSE ism=0
			CASE ism OF
				0 : BEGIN
					emiss=prof[*,i1:i2,0]
					emerr=err[*,i1:i2,0]
					emrho=interpol(irho,ipsin,psin[*,i1])		;convert to r/a grid
        	        	END
				1 : BEGIN
					emiss=prof[0:npsi/2-1,i1:i2,0]
					emerr=err[0:npsi/2-1,i1:i2,0]
					emrho=interpol(irho,ipsin,psin[0:npsi/2-1,i1])	;convert to r/a grid	
		                END
			ENDCASE
			emiss=sum_array(emiss,/i)/(i2-i1+1.0)
			emerr=sqrt(sum_array(emerr^2,/i))/(i2-i1+1.0)
			ied={emiss:emiss,emerr:emerr,rho:emrho,bright:bright,brerr:brerr,pos:pos,shot:shot,t1:t1,t2:t2,line:line,tht:tht,tit:tit,label:label,enorm:0.0,bnorm:0.0}
			*eptr[j]=ied
		ENDFOR
        ENDIF

	;calculate emissivity and brightness profiles from simulations
	IF NOT keyword_set(sptr) THEN sptr=ptrarr(nlines,/allocate_heap)
	FOR j=0,nlines-1 DO BEGIN
		IF size(*sptr[j],/type) NE 0 THEN BEGIN						;restore from sdata
			sdata=*sptr[j]
			pos=sdata.pos
                ENDIF ELSE BEGIN
			hirexsr_load_line_pos,shot,line[j],pos,tpos=tpos,tht=tht		;load POS vector
			tmp=where(pos[1,*] NE -1)
			pos=pos[*,tmp]
		ENDELSE
		FOR k=0,ntime-1 DO BEGIN
			;print, j,k,nlines,ntime
			IF k EQ 0 THEN BEGIN
				emsim=fltarr(nrho,ntime)
				eserr=fltarr(nrho,ntime)
			ENDIF
			CASE line[j] OF 
				2 : BEGIN
					calc_ar_line_rates,'z',temp[*,k]/1.0e3,q,rates
					emsim[*,k]+=dens[*,k]*csden[*,q[0],k]*rates[0,*]
					emsim[*,k]+=dens[*,k]*csden[*,q[1],k]*rates[1,*]
					emsim[*,k]+=dens[*,k]*csden[*,q[2],k]*rates[2,*]
					eserr[*,k]=emsim[*,k]*0.0					
				END
				3 : BEGIN	
					rates=reform_ark_data([3.73105,3.73105],/load)			;hard coded for the lya1 transition in the ArK database 
					emsim[*,k]+=dens[*,k]*csden[*,16,k]*interpol(rates.ion,alog10(rates.temp),alog10(temp[*,k]/1.0e3))
					emsim[*,k]+=dens[*,k]*csden[*,17,k]*interpol(rates.exc,alog10(rates.temp),alog10(temp[*,k]/1.0e3))
					emsim[*,k]+=dens[*,k]*csden[*,18,k]*interpol(rates.rec,alog10(rates.temp),alog10(temp[*,k]/1.0e3))
					eserr[*,k]=emsim[*,k]*0.0						;no error propigation right now
                	        END
			ENDCASE
		ENDFOR
		IF keyword_set(nobr) THEN BEGIN
			nch=n(pos[0,*])+1
			brsim=fltarr(nch)+1.0
			bserr=fltarr(nch)
			print, 'Simulated Brightness Profile Skipped'
		 ENDIF ELSE BEGIN
			xpos=pos
			genpos_pos_reform,xpos,[0.44,1.0,-0.6,0.6]	
			IF ntime NE 1 THEN brsim=genpos_line_brtau(xpos,emsim,rhovec,time,shot,time,brerr=bserr,/rho) ELSE brsim=genpos_line_br(xpos,emsim,rhovec,time,shot,time,/rho,brerr=bserr)
			print, 'Simulated Brightness Profile Computed'
                ENDELSE
		isd={emiss:emsim,emerr:eserr,rho:rhovec,bright:brsim,brerr:bserr,pos:pos,shot:shot,time:time,t1:t1,t2:t2,enorm:0.0,bnorm:0.0}
		*sptr[j]=isd
	ENDFOR
    END


PRO hirexsr_genrad_resid,eptr,sptr,rhomean=rhomean,eres=eres,bres=bres
	IF NOT keyword_set(rhomean) THEN rhomean=0.05
	x=size(eptr)
	nlines=x[1]	;see how many lines are stored
	IF NOT keyword_set(eres) THEN eres=ptrarr(nlines,/allocate_heap)
	IF NOT keyword_set(bres) THEN bres=ptrarr(nlines,/allocate_heap)
	FOR i=0,nlines-1 DO BEGIN
		edata=*eptr[i]
		sdata=*sptr[i]
		edata.enorm=mean(edata.emiss[where(edata.rho LE rhomean)])
		sdata.enorm=mean(sdata.emiss[where(sdata.rho LE rhomean)])
		sdata.bnorm=1.0
		edata.bnorm=edata.enorm/sdata.enorm
		*eptr[i]=edata
		*sptr[i]=sdata
		res=edata.emiss/edata.enorm-interpol(sdata.emiss/sdata.enorm,sdata.rho,edata.rho)
		err=sqrt((edata.emerr/edata.enorm)^2+(interpol(sdata.emerr/sdata.enorm,sdata.rho,edata.rho))^2)
		ier={res:res,rho:edata.rho,err:err}
		*eres[i]=ier
		res=edata.bright/edata.bnorm-sdata.bright/sdata.bnorm
		err=sqrt((edata.brerr/edata.bnorm)^2+(sdata.brerr/sdata.bnorm)^2)
		ibr={res:res,ch:indgen(n(edata.bright)+1)+1,err:err}
		*bres[i]=ibr
	ENDFOR
END

;+
;NAME:
;	HIREXSR_GENRAD_PLOT
;
;PURPOSE:
;	Used to plot comparisons of measured and simulated brightness
;	and emissivity profiles from HIREXSR
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 1/18/13
;
;-

PRO hirexsr_genrad_plot,eptr,sptr,plotwin=plotwin,ncase=ncase,eres=eres,bres=bres
	x=size(eptr)
	nlines=x[1]	;see how many lines are stored
	IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0
	IF NOT keyword_set(plotwin) THEN plotwin=10
	IF keyword_set(ps) THEN BEGIN
		xsize=7.0
		ysize=7.0*800/1100.0
		ls=1.5
	ENDIF ELSE BEGIN
		xsize=1100.0
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
	pos1=[0.09,0.3,0.55,0.95]
	pos2=[0.09,0.1,0.55,0.3]
	pos3=[0.64,0.3,0.98,0.95]
	pos4=[0.64,0.1,0.98,0.3]

	hirexsr_genrad_resid,eptr,sptr,rhomean=rhomean,eres=eres,bres=bres	;compute normalization and residuals
	nstr='MEAN < '+num2str(rhomean,dp=1)
	emax=1.0
	bmax=0.0
	nch=0
	ermax=0.0
	brmax=0.0
	FOR i=0,nlines-1 DO BEGIN						;compute axis scaling for plots
		edata=*eptr[i]
		sdata=*sptr[i]
		emax=max(edata.emiss/edata.enorm) > emax
		emax=max(sdata.emiss/sdata.enorm) > emax
		bmax=max(sdata.bright/sdata.bnorm) > bmax
		nch=max(n(edata.bright)+1) > nch
		resid=*eres[i]
		ermax=max(abs(resid.res+resid.err)) > max(abs(resid.res-resid.err)) > ermax
		resid=*bres[i]
		brmax=max(abs(resid.res+resid.err)) > max(abs(resid.res-resid.err)) > brmax
	ENDFOR
	
	;plot emissivity comparison
	normstr='NORMALIZATION: '+nstr
	tit='SHOT: '+num2str(edata.shot,1)+' '+num2str(edata.t1,dp=2)+' < t < '+num2str(edata.t2,dp=2)
	plot, [0],[0],pos=pos1,yr=[0,emax*1.05],ytit='Norm. Emissivity',/ysty,xr=[0,1.0],/xsty,chars=0.5*ls,$
		tit=tit,xtickname=replicate(' ',6),xticks=5,xminor=4
	FOR i=0,nlines-1 DO BEGIN
		edata=*eptr[i]
		sdata=*sptr[i]
		oploterror,sdata.rho,sdata.emiss/sdata.enorm,sdata.emerr/sdata.enorm,color=200,errcolor=200
		oploterror,edata.rho,edata.emiss/edata.enorm,edata.emerr/edata.enorm
	ENDFOR
	xyouts,0.95,emax*0.5,normstr,chars=0.4*ls,orient=90
	xyouts,0.1,emax*0.15,'MEASURED',chars=0.8*ls
	xyouts,0.1,emax*0.25,'SIMULATED',chars=0.8*ls,color=200

	;plot emissivity residual
	plot, [0],[0],pos=pos2,/noerase,yr=[-ermax,ermax]*1.05,xtit='r/a',ytit=n2g('Delta'),xr=[0,1.0],/xsty,/ysty,chars=0.5*ls,xticks=5,xminor=4
	FOR i=0,nlines-1 DO BEGIN
		resid=*eres[i]
		oploterror,resid.rho,resid.res,resid.err,psym=4,symsize=0.5*ls
	ENDFOR
	oplot,[0,1],[0,0],linestyle=2.0

	;plot brightness comparison
	plot, [0],[0],pos=pos3,/noerase,yr=[0,bmax*1.5],/ysty,xr=[0,nch+1],/xsty,ytit='Norm. Bright',chars=0.5*ls,xtickname=replicate(' ',10),tit=edata.tit
	FOR i=0,nlines-1 DO BEGIN
		edata=*eptr[i]
		sdata=*sptr[i]
		ch=indgen(n(edata.bright)+1)+1
		oploterror,ch,edata.bright/edata.bnorm,edata.brerr/edata.bnorm,psym=3,symsize=0.8*ls
		oploterror,ch,sdata.bright/sdata.bnorm,sdata.brerr/sdata.bnorm,color=200,errcolor=200
	ENDFOR

	;plot brightness residual
	plot, [0],[0],pos=pos4,/noerase,yr=[-brmax,brmax]*1.03,xr=[0,max(ch)+1],/xsty,/ysty,xtit='CH #',ytit=n2g('Delta'),chars=0.5*ls
	FOR i=0,nlines-1 DO BEGIN
		resid=*bres[i]
		oploterror,resid.ch,resid.res,resid.err,psym=4,symsize=0.5*ls
	ENDFOR
	oplot,[0,50],[0,0],linestyle=2.0
	
END
