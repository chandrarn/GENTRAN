;translates the output from quick_fits.pro
;written by A. Dominguez (9/10)
;generalized to use Te - M.L. Reinke (9/21/10)

pro qfit_ng,shot,dens,temp,rmaj,time,data=data,plot=plot,verb=verb
	if not(keyword_set(ne_deg)) then ne_deg=2
	if not(keyword_set(te_deg)) then te_deg=2

	quick_fit_ng,shot,data,plot=plot,verb=verb	;run quick fit

	;reorganize output data [space,time] arrays
	time        = data.tdat			;time in [sec]
	rmaj        = transpose(data.RR)
	dens = transpose(data.nefit)/1.e20 ; density in 10^20 [m-3]
	temp = transpose(data.tefit)/1000. ; temperature in [keV]
end

FUNCTION tanhfn,x,r,w
	; a hyperbolic tangent function for splicing together the piecewise functions
	return,( 1+tanh( (x-r)/w*!pi ) )/2.
END

FUNCTION mtanh,x,p
	; modified tanh for boundary fit, c.f. R. J. Groebner, Nucl. Fusion. 44 (2004) 204-213
	;p=[offset, amplitude, center, width]
	return, p[0] + p[1]*( 1+tanh( (x-p[2])/(-p[3])*!pi) )/2.
END

;+
;NAME:
;	MPOLYPROF
;
;PURPOSE:
;	This function calculates a polynomial based piecewise profile function 
;	for plasma parameters such as temperature and/or density. The format 
;	is made to be used with MPFIT.
;
;CALLING_SEQUENCE:
;	result=MPOLYPROF(x,p)
;
;INPUTS:
;	x	FLTARR	[n_x] of the points for which to calculate the profile
;	p	FLTARR 	[2+deg] where nb will determine the order of the baseline fit
;			 p[0,*] = [r0,r1, c0, c1, c2, ... cn] ; boundary fit
;			           minimum and maximum radii, polynomial coefficients
;			 p[1,*] = [r0,r1, c0, c1, c2, ... cn] ; outer core (q>1)
;			 p[2,*] = [r0,r1, c0, c1, c2, ... cn] ; internal core (q<=1)
;	deg	INT	(optional) defaults to 2. degree of polynomials to fit to
;			each section.
;	tr_wid	FLT	(optional) [m] defaults to 1cm=0.01m. Transition width between 
;			piecewise polynomials
;
;OUTPUTS:
;	result:	FLTARR  [n_x] of the discretized profile
;
;PROCEDURE:
;	
;	Profiles are fitted using a set of three polynomial fits to data:
;	 1. Edge/boundary from the separatrix to the limiter boundary
;	 2. Outer core (q> 1)
;	 3. Inner core (q<=1)
;	These fits are assembled piecewise using tanh functions with a transition
;	width defined by tr_wid
;
;MOFICATION HISTORY:
;	Written by: 	A. N. James, influenced by Yunxing Ma's quick_fit routine
;-
FUNCTION mpolyprof,x,p
	nseg=3
	nd=size(p,/n_el)/nseg-2 ; number of polynomial coefficients
	nx=size(x,/n_el)

	p=reshape(p,[nseg,2+nd])
	y=fltarr(nx)
	w=0.01 ;[m]

	; construct polynomial pieces for each segment
	fpoly=fltarr(nseg,nx)
	for i=0,nseg-1 do $
	 for j=0,nd-1 do $
	  fpoly[i,*]=fpoly[i,*] + p[i,2+j]*x^j

	; suppress segments outside their region of validity
	fseg=fpoly
	for i=0,nseg-1 do $
	 case i of
		0 : $ ; suppress only on left for first segment (boundary)
			fseg[i,*]=fseg[i,*] *tanhfn(x,p[nseg,1], w)
		nseg-1 : $ ; suppress only on right for last segment (inner core)
			fseg[i,*]=fseg[i,*] *tanhfn(x,p[nseg,0],-w)
		else : $
			fseg[i,*]=fseg[i,*] *tanhfn(x,p[nseg,0], w) $
					    *tanhfn(x,p[nseg,1],-w)
	endcase

	f=total(fseg)

	return,y

; i.e.	x=make(0,1.2,1000)
;	p[[0.,0.6,3000.,0,0] $
; 	  [.6,0.9,2000.,-200/.6,0], $
;	  [.9,1.2,1800,
;	y=mpolyprof(x,p)
;	window,/free & plot,x,y
END


Pro quick_fit_ng,shot,d,ets=ets,efit_tree=efit_tree,plot=plot,verb=verb

str_shot = strtrim(string(shot),2)
if keyword_set(efit_tree) eq 0 then efit_tree='analysis'
if keyword_set(verb) eq 0 then verb=0

; retrieve efit data
efit_good=efit_check(times,ntimes,shot=shot,tree=efit_tree)
if ntimes eq 0 then message, 'No good EFIT points'
mdsopen, efit_tree,shot
	t_efit_bdry=mdsvalue('dim_of(\efit_geqdsk:nbbbs)',quiet=(verb eq 0),status=s0) ;Efit Boundary time
	t_ip = mdsvalue('dim_of(\'+efit_tree+'::efit_aeqdsk:aout)',quiet=(verb eq 0))  ;Plasma current time
	t_efit=mdsvalue('dim_of(\EFIT_AEQDSK:rout)',quiet=(verb eq 0),status=s1)       ;Efit time
	r_mag = mdsvalue('\'+efit_tree+'::efit_aeqdsk:rmagx',quiet=(verb eq 0),status=status)/100. ;Mag axis in m
	r_out = mdsvalue('\'+efit_tree+'::efit_aeqdsk:rout',quiet=(verb eq 0),status=status) /100. ;R0 in m
	a_out = mdsvalue('\'+efit_tree+'::efit_aeqdsk:aout',quiet=(verb eq 0)) /100. ; a0 in m
	q_95  = mdsvalue('\'+efit_tree+'::efit_aeqdsk:q95',quiet=(verb eq 0)) ;q95
	chord_04 = .01*mdsvalue('\'+efit_tree+'::efit_aeqdsk:rco2v[3,*]',quiet=(verb eq 0)) ;[m] length of interferometer chord through plasma volume
	r_sep = r_out+a_out      ;Rsep in meter 

	if efit_tree eq 'analysis' $
	 then bt_0   = mdsvalue('\ANALYSIS::TOP:EFIT.RESULTS.A_EQDSK:BT0', quiet=(verb eq 0)) $ ;[T]
	 else bt_0   = mdsvalue('\'+efit_tree+'::TOP.RESULTS.A_EQDSK:BT0', quiet=(verb eq 0))  ;[T]

	nt=n_elements(t_efit) 	;number of Efit time
	nbdry=n_elements(t_efit_bdry) ;number of Efit boundary time
	if (s0 and s1) then begin
		if nt ne nbdry then begin	        	
			print,'Number of boundary time points '+num2str(nbdry,1) $
			      +' is not equal to number of "rout" times '+num2str(nt,1)
			tmp=where(t_efit_bdry eq t_efit)
			t_efit=t_efit[tmp]
		endif		
		t_efit = t_efit[efit_good]
		r_sep = r_sep[efit_good]  ;[m]
	endif else message,'Could not get EFIT data'
mdsclose, efit_tree,shot


;build data arrays

mdsopen,'electrons',shot
	tdat  =  mdsvalue('dim_of(\ELECTRONS::TOP.YAG_NEW.RESULTS.PROFILES:R_MID_T)',quiet=(verb eq 0),status=stat_ts)
	nt  =  size(tdat,/n_el)
	if stat_ts then begin

		; start with one assumed boundary point
		zeros=reform(fltarr(nt),[nt,1])
		rlim=0.90 ; limiter location
		nlim=0.1 ; limiter density
		tlim=10. ; limiter temperature
		rdat=rlim+zeros
		nedat=0.1+zeros ; [x1e20/m^3]
		nedat_err=0.1+zeros ; x1e20/m^3]
		tedat=10.+zeros ;[eV]
		tedat_err=10.+zeros ;[eV]

		; core TS data
		rmid_cts = mdsvalue('\ELECTRONS::TOP.YAG_NEW.RESULTS.PROFILES:R_MID_T', quiet=(verb eq 0)) ;[m]
		ne_cts = mdsvalue('\ELECTRONS::TOP.YAG_NEW.RESULTS.PROFILES:NE_RZ', quiet=(verb eq 0)) ;[/m^3]
		ne_cts_err = mdsvalue('\ELECTRONS::TOP.YAG_NEW.RESULTS.PROFILES:NE_ERR', quiet=(verb eq 0)) ;[m^3]
		te_cts = mdsvalue('\ELECTRONS::TOP.YAG_NEW.RESULTS.PROFILES:TE_RZ', quiet=(verb eq 0))*1000. ;[eV]
		te_cts_err = mdsvalue('\ELECTRONS::TOP.YAG_NEW.RESULTS.PROFILES:TE_ERR', quiet=(verb eq 0))*1000. ;[eV]
		; Get Core TS channel Mask for fitting. If masked the channel will not be used in fitting
		mask = mdsvalue('\ELECTRONS::TOP.YAG_NEW.RESULTS.PROFILES:FIT_MASK',quiet=(verb eq 0),status=st_mask)
		if st_mask then begin 
			zvert = mdsvalue('\calib:z_nom',quiet=(verb eq 0))
			zsort = mdsvalue('\ELECTRONS::TOP.YAG_NEW.RESULTS.PROFILES:Z_SORTED',quiet=(verb eq 0))
			zmask = zvert[mask-1]
			ptr = intarr(n_elements(zsort))
			for j=0,n_elements(zmask)-1 do begin 
				i_mask = where(zsort eq zmask[j], nmask)
				ptr[i_mask] = 1
			endfor
			i_nomask = where(ptr eq 0, n_nomask)
			if n_nomask ne 0 then begin 
				rmid_cts = rmid_cts[*,i_nomask]
				ne_cts = ne_cts[*,i_nomask]
				ne_cts_err = ne_cts_err[*,i_nomask]
				te_cts = te_cts[*,i_nomask]
				te_cts_err = te_cts_err[*,i_nomask]
			endif
		endif

		rdat=[[rdat],[rmid_cts]]
		nedat=[[nedat],[ne_cts]]
		nedat_err=[[nedat_err],[ne_cts_err]]
		tedat=[[tedat],[te_cts]]
		tedat_err=[[tedat_err],[te_cts_err]]

		; edge TS data
		if keyword_set(ets) then begin
			rmid_ets = mdsvalue('\ELECTRONS::TOP.YAG_EDGETS.RESULTS:RMID', quiet=(verb eq 0)) ;[m] 
			ne_ets = mdsvalue('\ELECTRONS::TOP.YAG_EDGETS.RESULTS:NE', quiet=(verb eq 0)) ;[/m^3]
			ne_ets_err = mdsvalue('\ELECTRONS::TOP.YAG_EDGETS.RESULTS:NE:ERROR', quiet=(verb eq 0)) ;[/m^3]
			te_ets = mdsvalue('\ELECTRONS::TOP.YAG_EDGETS.RESULTS:TE', quiet=(verb eq 0)) ;[eV]
			te_ets_err = mdsvalue('\ELECTRONS::TOP.YAG_EDGETS.RESULTS:TE:ERROR', quiet=(verb eq 0)) ;[eV]

			rdat=[[rdat],[rmid_ets]]
			nedat=[[nedat],[ne_ets]]
			nedat_err=[[nedat_err],[ne_ets_err]]
			tedat=[[tedat],[te_ets]]
			tedat_err=[[tedat_err],[te_ets_err]]
		endif

		for i=0,nt-1 do begin
			tmp=sort(rdat[i,*]) ; sort by radius
			rdat[i,*]=rdat[i,tmp]
			nedat[i,*]=nedat[i,tmp]
			nedat_err[i,*]=nedat_err[i,tmp]
			tedat[i,*]=tedat[i,tmp]
			tedat_err[i,*]=tedat_err[i,tmp]
		endfor

		; add a bogus point at R=0
		;rdat=[[zeros],[rdat]]
		;nedat    =[[total(nedat[*,0:2]    ,2)/3],[nedat]]
		;nedat_err=[[total(nedat_err[*,0:2],2)/3],[nedat_err]]
		;tedat    =[[total(tedat[*,0:2]    ,2)/3],[tedat]]
		;tedat_err=[[total(tedat_err[*,0:2],2)/3],[tedat_err]]

	endif else begin 
		print,'No TS data for this shot -Exit fitting'
	endelse

	;TCI data for density cross calibration
	ttci  = mdsvalue('dim_of(\ELECTRONS::TOP.TCI.RESULTS:NL_04)', quiet=(verb eq 0), status=stat_tci)
	if stat_tci then nl04  = mdsvalue('\ELECTRONS::TOP.TCI.RESULTS:NL_04', quiet=(verb eq 0))/1.e20 ;[x1e20/m^3]

mdsclose,'electrons', shot


;sync timebases
tmp = where(tdat ge min(t_efit) and tdat le max(t_efit), nfit)
if nfit eq 0 then begin 
	print,'No TS data for good efit times -Exit fitting'
	return
endif else begin 
	tdat = tdat[tmp]
	 nt=n_elements(tdat)
	rdat=rdat[tmp,*]
	nedat=nedat[tmp,*]
	nedat_err=nedat_err[tmp,*]
	tedat=tedat[tmp,*]
	tedat_err=tedat_err[tmp,*]

	;interpolate EFIT outputs onto TS time axis 
	r_mag = interpol(r_mag, t_ip, tdat)
	r_sep = interpol(r_sep,t_efit, tdat)
	r_out = interpol(r_out, t_ip, tdat)
	a_out = interpol(a_out, t_ip, tdat)
	q_95  = interpol(q_95,t_ip, tdat)
	bt_0  = interpol(bt_0,t_ip, tdat)
	t_efit= tdat

	;Interpolate TCI data onto TS time axis
	nl04  = interpol(nl04, ttci, tdat)
	chord_04 = interpol(chord_04, t_ip, tdat)
endelse
 
; output variables
nx   = 100
RR   = dblarr(nt,nx)  ;Std radius axis that the output profiles use
Rrho = dblarr(nt,nx)
nefit= dblarr(nt,nx)
tefit= dblarr(nt,nx)
rped = dblarr(nt)
c0   = dblarr(nt,4)
c1   = dblarr(nt,4)

; Profile Fitting Loop
if keyword_set(plot) && plot ge 2 then window,0
for i=0, nt-1 do begin 

	R0  = r_out[i]  ;midplane geometric center of LCFS 
	a0  = a_out[i]  ;midplane minor radius
	rsep= r_sep[i]  ;magnetic separatrix
	Rx  = r_mag[i]  ;magnetic axis
	q95 = q_95[i]   ;q95
	nl  = nl04[i]   ;nl04
	cd04 = chord_04[i]   ;chord_04
	if keyword_set(verb) then print,  'i='+num2str(i,1) $
					+' t='+num2str(tdat[i])+'s' $
					+' R0='+num2str(R0,dp=2)+'m' $
					+' a0='+num2str(a0,dp=2)+'m' $
					+' rsep='+num2str(rsep,dp=2)+'m' $
					+' Rx='+num2str(Rx,dp=2)+'m' $
					+' q95='+num2str(q95,dp=2) $
					+' nl='+num2str(nl,dp=2) $
					+' cd04='+num2str(cd04,dp=2)


	for ifit=0,1 do begin
		; Select data, and exclude non-positive points and points beyond limiter from fitting
		case ifit of
		  0 : begin
			; Ne DATA PREPARATION
			r = reform(rdat[i,*])
			n = reform(nedat[i,*])
			nerr = reform(nedat_err[i,*])
			tmp = where(n ge nlim and r le rlim)
			if keyword_set(verb) then print,' ne : max='+num2str(max(n))+' min='+num2str(min(n))
		      end
		  1 : begin
			; Ne DATA PREPARATION
			r = reform(rdat[i,*])
			n = reform(tedat[i,*])
			nerr = reform(tedat_err[i,*])
			tmp = where(n ge tlim and r le rlim)
			if keyword_set(verb) then print,' te : max='+num2str(max(n))+' min='+num2str(min(n))
		      end
		endcase
		if tmp[0] ne -1 then begin
			nel_tmp=n_elements(tmp)
			nel_r=n_elements(r)
			if keyword_set(verb) && nel_tmp lt nel_r then print,'  removing '+num2str(nel_r-nel_tmp,1)+' negative or external data points'
			r=r[tmp]
			n=n[tmp]
			nerr=nerr[tmp]
		endif else begin
			if keyword_set(verb) then print,'  no non-negative data inside limiter, skipping fits'
			continue
		endelse
		
		nr=n_elements(r)
		RR0=min(r) < rx
		RR1=rlim
		RRs=rsep < rlim
		RR[i,*]=make(RR0,RR1,nx)
		Rrho[i,*] = (RR[i,*]-RR0)/(RRs-RR0)

		; Sometimes some channels errors=0.0, set to 15% then
		tmp = where(nerr eq 0.0)
		if tmp[0] ne -1 then begin
			nerr[tmp] = 0.15*n[tmp]
			if keyword_set(verb) then print,'  changed error for following points from zero to 15% : '+strjoin(string(tmp,format='(i4)'),' ')
		endif


		; PROFILE FITTING SECTION

		; fit outside plasma boundary
		if rsep lt rlim $
		 then rped[i]=rsep-0.02 $
		 else rped[i]=rlim-0.02
		tmp = where(r ge rped[i],ntmp)
		nmin=6
		if ntmp lt nmin then begin
			tmp=indgen(nmin)-nmin+nr
			if keyword_set(verb) then print,'  insufficient data outside of boundary, using outermost '+num2str(nmin,1)+' points for fit'
		endif
		dum=min(r[tmp],imn,subscript_max=imx) & imn=imn+tmp[0] & imx=imx+tmp[0]
		tmpb=tmp

		c0[i,*]=[n[imx], n[imn]-n[imx], rped[i], max(r[tmp])-min(r[tmp])] ; initial guess [offset, amplitude, center, width]
		parinfo = replicate({value:0., fixed:0, limited:[1,1],limits:[0.,0], tied:''}, n_elements(c0[i,*]))
		 parinfo.value=reform(c0[i,*])
		 parinfo.limits=[[n[imx],		n[imn]], $
		 		 [n[imn]-n[imx],	5.0*(n[imn]-n[imx])], $
				 [rped[i],		rlim], $
				 [0.1*(r[imx]-r[imn]),	r[imx]-r[imn]] ]

		c0[i,*]=mpfitfun('mtanh',r[tmp],n[tmp],nerr[tmp],c0[i,*], $
				 perror=error,parinfo=parinfo,status=status,niter=niter)
		if keyword_set(verb) then print,'  mpfit status='+num2str(status,1)+' coeffs='+strjoin(string(transpose(c0[i,*]),format='(e8.2)'),' ')
		if status eq 0 then begin & print,[transpose(parinfo.value),[parinfo.limits]] & stop & endif
		nfit0=mtanh(RR[i,*],c0[i,*])

		if keyword_set(plot) && plot ge 2 then begin
			yr=[min(n),max(n)]
			plot,[0,0],[0,0],xr=[min(rdat),max(rdat)],yr=[min(n),max(n)],title='shot '+num2str(shot,1)+' t='+num2str(tdat[i],dp=3)+'s'
			oploterr,r[tmp],n[tmp],nerr[tmp]
			oplot,RR[i,*],nfit0
			oplot,[rx,rx],yr & oplot,rped[[i,i]],yr
		endif

		; fit inside plasma boundary
		tmp = where(r le rped[i],ntmp)
		dum=min(r[tmp],imn,subscript_max=imx) & imn=imn+tmp[0] & imx=imx+tmp[0]
		tmpc=tmp

		;npoly=3
		;c1=poly_fit(r[tmp],n[tmp],npoly,measure_errors=nerr[tmp],/double,yfit=yfit1)
		;nfit1=poly(RR[i,*],c1)
		c1[i,*]=[n[imx], n[imn]-n[imx], mean(r[tmp]), max(r[tmp])-min(r[tmp])] ; initial guess [offset, amplitude, center, width]
		parinfo = replicate({value:0., fixed:0, limited:[1,1],limits:[0.,0], tied:''}, n_elements(c1[i,*]))
		 parinfo.value=reform(c1[i,*])
		 parinfo.limits=[[n[imx],		n[imn]], $
		 		 [n[imn]-n[imx],	5.0*(n[imn]-n[imx])], $
				 [r[imn],		r[imx]], $
				 [0.1*(r[imx]-r[imn]),	r[imx]-r[imn]] ]

		c1[i,*]=mpfitfun('mtanh',r[tmp],n[tmp],nerr[tmp],c1[i,*], $
				 perror=error,parinfo=parinfo,status=status,niter=niter)
		if keyword_set(verb) then print,'  mpfit status='+num2str(status,1)+' coeffs='+strjoin(string(transpose(c1[i,*]),format='(e8.2)'),' ')
		if status eq 0 then begin & print,[transpose(parinfo.value),[parinfo.limits]] & stop & endif
		nfit1=mtanh(RR[i,*],c1[i,*])

		; splice together the fits
		wspl=0.01 ; splice width
		nfit = nfit1 ;* tanhfn(RR[i,*],rped[i],-wspl) ; core contribution
			;+ nfit0 * tanhfn(RR[i,*],rped[i], wspl) ; boundary contribution

		if keyword_set(plot) && plot ge 2 then begin
			oploterr,r[tmp],n[tmp],nerr[tmp]
			oplot,RR[i,*],nfit1
			oplot,RR[i,*],nfit
			;stop
		endif

		case ifit of
		 0 : nefit[i,*]=nfit
		 1 : tefit[i,*]=nfit
		endcase
	endfor
endfor


if keyword_set(plot) then begin
	if tmp[0] ne -1 then begin
		window,/free

		pos=[.1,.1,.9,.5]
		plot,[0,0],[0,0],xr=[min(RR),max(RR)],yr=[0,max(nedat)*1.15],pos=pos,ytitle='n!Ie!N [/m!E3!N]',xtitle='R [m]'
		 if size(dene,/type) ne 7 $
		  then for i=0,nt-1 do begin
			colr=round(1+(i+1.)/nt*253)
			oplot,RR[i,*],nefit[i,*],color=colr
			errplot,rdat[i,*],nedat[i,*]-nedat_err[i,*],nedat[i,*]+nedat_err[i,*],color=colr
		       endfor

		pos=[.1,.5,.9,.9]
		plot,[0,0],[0,0],xr=[min(RR),max(RR)],yr=[0,max(tedat)*1.15],pos=pos,ytitle='T!Ie!N [eV]',/noerase,title='shot '+num2str(shot,1)+' thomson data and fits'
		 if size(te,/type) ne 7 $
		  then for i=0,nt-1 do begin
			colr=round(1+(i+1.)/nt*253)
			oplot,RR[i,*],tefit[i,*],color=colr
			errplot,rdat[i,*],tedat[i,*]-tedat_err[i,*],tedat[i,*]+tedat_err[i,*],color=colr
		       endfor

		for i=0,nt-1 do $
		 if i mod (nt/10) eq 0 then begin
			 col=round(1+(i+1.)/nt*254)
			 xyouts,.8,.8-round(10.*i/nt)*0.02,'t='+num2str(tdat[i],dp=3)+'s',/normal,color=col
		 endif
	endif else print,'no data before t=0 to plot'
endif


d={tdat:tdat,rdat:rdat,nedat:nedat,nedat_err:nedat_err,tedat:tedat,tedat_err:tedat_err, $ ; ts data
   r_mag:r_mag,r_sep:r_sep,r_out:r_out,a_out:a_out,q_95:q_95,bt_0:bt_0,nl04:nl04,chord_04:chord_04, $ ; efit data
   RR:RR,Rrho:Rrho,nefit:nefit,tefit:tefit,rped:rped,c0:c0,c1:c1} ; fit data

	;c = poly_fit(x, y, deg, /double, measure_errors=yerr, yfit=ypoly, chisq=chi,status=stat, covar=cov)

end





