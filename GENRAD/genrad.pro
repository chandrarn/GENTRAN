
;+
;NAME: 
;	CALC_SPECEMISS
;
;-

FUNCTION calc_specemiss,lam,csden,temp,dens,qmin=qmin,qmax=qmax,w=w,thr=thr,adas=adas,cserr=cserr,temperr=temperr,denserr=denserr,emerr=emerr
	IF NOT keyword_set(thr) THEN thr=0.0001
	
	dlam=[min(lam),max(lam)]
	IF NOT keyword_set(w) THEN w=0.01*(dlam[1]-dlam[0])
	nrho=n(temp)+1
	nlam=n(lam)+1
	z=n(csden(*,0))
	IF NOT keyword_set(cserr) THEN cserr=csden*0.0
	IF NOT keyword_set(temperr) THEN temperr=temp*0.0
	IF NOT keyword_set(denserr) THEN denserr=dens*0.0

	em=dblarr(nlam,nrho)
	emerr=dblarr(nlam,nrho)
	IF NOT keyword_set(adas) THEN out=calc_spec(z,temp,dlam,qmin=qmin,qmax=qmax,terr=temperr) ELSE out=adas_pec_spec(z,temp,dens,dlam)

	;reduct to meaningful number of lines	
	int=sum_array(out.sigv,/i)
	int/=max(int)
	tmp=where(int GT thr)
	lam_o=out.lam[tmp]
	q=out.q[tmp]
	sigv=out.sigv[tmp,*]
	sigv_err=out.err[tmp,*]
	nline=n(lam_o)+1

	FOR i=0,nrho-1 DO BEGIN
		print, i
		FOR j=0,nline-1 DO BEGIN
			em[*,i]+=sigv[j,i]/(w*sqrt(2.0*!pi))*exp(-(lam-lam_o[j])^2/(2.0*w^2))*csden[q[j],i]*dens[i]
			emerr[*,i]+=(em[*,i]/sigv[j,i])^2*sigv_err[j,i]^2+(em[*,i]/csden[q[j],i])^2*(cserr[q[j],i])^2+(em[*,i]/dens[i])^2*denserr[i]^2
		ENDFOR
		emerr[*,i]=sqrt(emerr[*,i])
	ENDFOR
	output=em
	IF keyword_set(debug) THEN stop
	RETURN,output
END


FUNCTION calc_plc,csplc,fq,te=te,nel=nel,debug=debug,contplc=contplc

	z=fq.z	
	temp=fq.temp
	IF NOT keyword_set(te) THEN te=temp
	IF NOT keyword_set(nel) THEN nel=fq.dens[0]

	ntemp=n(te)+1
	IF n(csplc.dens) EQ 0 THEN csdpt=0 ELSE csdpt=ipt(csplc.dens,nel)
	IF n(contplc.dens) EQ 0 THEN contdpt=0 ELSE contdpt=ipt(contplc.dens,nel)

	plc=fltarr(ntemp)
	FOR i=0,ntemp-1 DO BEGIN
		FOR j=0,z DO BEGIN
			plc[i]+=interpol(csplc.plc[j,*,csdpt],csplc.temp,te[i])*interpol(fq.fq[j,*],fq.temp,te[i])
			IF keyword_set(contplc) THEN plc[i]+=interpol(contplc.plc[j,*,contdpt],contplc.temp,te[i])*$
				interpol(fq.fq[j,*],fq.temp,te[i])
		ENDFOR
	ENDFOR
	IF keyword_set(contplc) THEN cont_path=contplc.path ELSE cont_path='NONE'	

	output={plc:plc,temp:te,dens:[nel],z:z,csplc:csplc.path,contplc:cont_path,ion:fq.ion,rec:fq.rec}
	IF keyword_set(debug) THEN stop
	RETURN,output
END


;+
;NAME:
;	TCRASH
;	
;PURPOSE:
;	This function determines the times of sawteeth crash by looking at a core ECE channel (GPC_T0)
;
;CALLING SEQUENCE:
;	result=TCRASH(shot,t1,t2)
;	
;INPUTS:
;	shot	LONG	shot number
;	t1	FLOAT	lower time bound [sec]
;	t2	FLOAT 	upper time bound [sec]
;
;OPTIONAL INPUTS:
;	dt	FLOAT	lower bound on the length of a sawtooth [sec] (DEFAULT: 0.005)
;	thr	FLOAT	bound on the Te time derivative [keV/s] (DEFAULT: -1.0)
;
;KEYWORD PARAMETERS:
;	plot	/plot displays an IDL plot showing the Te data and the identified crashes
;
;OUTPUTS:
;	result	FLTARR	[ncrash] of the times [sec] of the arb. number of sawtooth crashes found
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke (2009?) 
;
;-

;finds the time point of the sawtooth crash
FUNCTION tcrash,shot,t1,t2,dt=dt,thr=thr,plot=plot
	IF NOT keyword_set(dt) THEN dt=0.005	;10 ms minimum sawtooth period 
	IF NOT keyword_set(thr) THEN thr=-1.0

	mdsopen,'electrons',shot
	te=mdsvalue('\ELECTRONS::TOP.ECE.GPC_RESULTS.TE:GPC_T0')
	t=mdsvalue('dim_of(\ELECTRONS::TOP.ECE.GPC_RESULTS.TE:GPC_T0)')
	mdsclose,'electrons',shot

	tmp=where(t GE t1 AND t LE t2)
	t=t[tmp]
	te=te[tmp]
	ncrash=floor((t2-t1)/dt)+1
	tcr=fltarr(ncrash*10.0)
	dte=deriv(t,te)*1.0e-3
	tmp=1
	cntr=0
	IF keyword_set(plot) THEN plot, t,te,xtit='Time [sec]',ytit='T!le!n [keV]',chars=1.3,$
		tit=n2g('Delta')+'t='+num2str(dt*1.0e3)+' [ms]   thr='+num2str(thr,dp=2)+' [keV/ms]'
	WHILE min(dte) LE thr AND tmp[0] NE -1 DO BEGIN
		check=intarr(n(t)+1)+1
		mpt=minloc(dte)
		tcr[cntr]=t[mpt]
		cntr+=1
		tmp=where(t GE t[mpt]-dt/2.0 AND t LE t[mpt]+dt/2.0)
		IF tmp[0] NE -1 THEN BEGIN
			check[tmp]=0
			tmp=where(check EQ 1)
			dte=dte[tmp]
			te=te[tmp]
			t=t[tmp]
		ENDIF
	ENDWHILE
	tcr=tcr[where(tcr NE 0)]
	tcr=tcr[sort(tcr)]
	IF keyword_set(plot) THEN FOR i=0,n(tcr) DO oplot,tcr[i]*[1,1],[0,10.0],linestyle=2.0,color=200

	output=tcr
	RETURN,output
END		

