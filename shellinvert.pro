function foldmatrix,npts,rm,rcent
; Return a matrix that represents the difference between corresponding 
; values when folded about rcent. Such that foldmatrix#emiss gives
; a vector of difference values. Subsequently a minimization of the
; norm of the matrix constitutes a symmetrization about the presumed
; center position.
fm=fltarr(npts,npts)
rdiff=rm-rcent
wneg=where(rdiff lt 0.)
if(wneg(0) ne -1) then begin
    for j=0,n_elements(wneg)-1 do begin
        rcorr= -rdiff(wneg(j))
        nrcorr=where((rdiff ge rcorr ) and (shift(rdiff,1) le rcorr))
        nrcorr=nrcorr(0)
        fm(wneg(j),wneg(j))=fm(wneg(j),wneg(j))-1.
        denom=rdiff(nrcorr)-rdiff(nrcorr-1)
        fm(wneg(j),nrcorr)=fm(wneg(j),nrcorr)+(rdiff(nrcorr)-rcorr)/denom
        fm(wneg(j),nrcorr-1)=fm(wneg(j),nrcorr-1)-(rdiff(nrcorr-1)-rcorr)/denom
;        print,wneg(j),rdiff(wneg(j)),nrcorr,rdiff(nrcorr),rdiff(nrcorr-1)
    endfor
endif
;print,fm,form='(17f4.1)'
return,fm
end

;+
;NAME:
;	SHELLINVERT
;
;PURPOSE:
;	This procedure inverts brightness profiles into emissivity profiles using the
;	typicall shell inversion procedure.  Optional smoothing using 2nd derivative
;	minimization is also included.
;
;CALLING SEQUENCE:
;	SHELLINVERT,bright,rtang,emiss,radii
;
;INPUTS:
;	bright:		FLTARR [nrtang,ntime] array of brightness data
;	rtang:		FLTARR [nrtang] vector of tangency radii for brightness data
;
;OPTIONAL INPUTS:
;	good:		INTARR [nrtang] vetor of 1's for good channels and 0 for bad channels [DEFAULT all good]
;	rout:		FLT value, in same units as rtang, where there is no plasma [DEFAULT 0.93] 
;	eps:		FLT weighting constant for smoothing [DEFAULT = 0 - no smoothing]
;	npts:		INT number of points to create emissivity profile when eps > 0 [DEFAULT = 34]
;	br_err:		FLTARR [nrtang] of the absolute error in units of bright
;	ff:		INTARR [nrtang] of "fudge factors" to be sent to CALC_GEOM to account for partical occlusion of a channel by
;				structures like the inner wall
;	drt:		FLTARR [nrtang] of the calibration uncertainties in RTANG
;
;KEYWORD PARAMETERS:
;	verb:		/verb to print out non-critical information to terminal
;	debug:		/debug to stop code before END to access local variables
;
;OUTPUTS:
;	emiss:		FLTARR [npts,ntime] array of emissivity data
;	radii:		FLTARR [npts] vector of radii to index emissivity data
;
;OPTIONAL OUTPUTS:
;	svdout:		STRUC holds the u,v,w data returned from SSVDC.PRO
;	brchk:		STRUC holds the br and rad data using L#emiss as a check on the inversion
;	em_err:		FLTARR [npts] of the error in the emissivity values (if br_err is non-zero)
;	rp:		FLTARR [npts] of the radii at rtang+drt
;	rm:		FLTARR [npts] of the radii at rtang-drt
;	
;RESTRICTIONS:
;	 MLR_FUNCTIONS must be compiled before SHELLINVERT.
;	This file can usually be found in /home/mlreinke/idl/general/mlr_fucntions
;
;PROCEDUE:
;	The workings of this procedure were generalized from the generic geometry 
;	of the shell inversion.  The smoothing algorithm was taken from Ian Hutchinson's
;	BOLO_SVD.PRO code as was code to calculate the L matrix.  
;	
;	See IHH 2/20/2004 note for the description of the math.
;	***as of 5/10/2011 this has moved to an LA_INVERT(/double) scheme***
;
;EXAMPLE:
;	SHELLINVERT,bright,rtang,emiss,radii,npts=40,rout=0.95,eps=0.1
;	
;	This will take the (m x n) brightness array and creat a (40 x n)		
;	emissivity array with a 40 point radial scale that is smoothed.
;
;MODIFICATION HISTORY
;	Adapted from:	IHH BOLO_SVD.PRO
;	Written By:	ML Reinke 7-30-05
;	10-08-05:	ML Reinke - updated some /verb outputs, no
;                                   functional changes
;	2-16-09:	ML Reinke - added functionality to calculate  error direction w/o multiple
;                                   inversions using propigation of error
;	2-12-10:	ML Reinke - added the ff optional input to account for inner wall occlusion.
;	5-10-11:	ML Reinke - moved inversion from SVD to LA_INVERT, made error calculation
;                                   done the same as GENPOS_PROFILE_INVERT. Added drt tools.
;  
;-


PRO shellinvert,bright,rtang,emiss,radii,beta=beta,good=good,npts=npts,eps=eps,redge=redge,svdout=svdout,brchk=brchk,verb=verb,debug=debug,$
                br_err=br_err,em_err=em_err,drt=drt,rp=rp,rm=rm,ff=ff,fmweight=fmweight,axis=axis


	;bright must be indexed as (r,t)
	ntime=n(bright[0,*])+1
        n_rtang=n(rtang)+1
        brin=bright
        rtin=rtang
        IF keyword_set(br_err) THEN errin=br_err
        IF NOT keyword_set(br_err) THEN br_err=-1.0				;set error equal to -1.0 if not set
        x=size(br_err)
        IF x[0] EQ 1 AND br_err[0] NE -1 THEN BEGIN				;if error is set 1d then write constant over time
               	err=fltarr(x[1],ntime)
                FOR i=0,ntime-1 DO err[*,i]=br_err 
                br_err=err
        ENDIF
        IF NOT keyword_set(ff) THEN ff=fltarr(n_rtang)+1.0			;goes through to calc_geom for adjusting length matrix due to obstruction
        IF NOT keyword_set(redge) THEN redge=0.93		 		;last radial point
        IF NOT keyword_set(npts) THEN npts=34 					;number of points to define emissivity profile
	IF NOT keyword_set(eps) THEN eps=0					;default smoothing
        
        
	;truncate dataset based on good, if called
	IF keyword_set(good) THEN BEGIN
		IF n(good) NE n(bright[*,0]) THEN BEGIN
			print, 'good vector does not match sizing of brightness'
			RETURN
		ENDIF
		IF keyword_set(verb) THEN print, ' Removing Bad Channels'
		tmp=where(good EQ 1)
		bright=bright[tmp,*]
                IF br_err[0] NE -1 THEN br_err=br_err[tmp,*]
                rtang=rtang[tmp]
                ff=ff[tmp]
                n_rtang=total(good)
	ENDIF

	;add outermost data point on rtang and fix brightness=0
	rtang=[rtang,redge]
        ff=[ff,1.0]
	bright=[bright,transpose(fltarr(ntime))]
        IF br_err[0] NE -1 THEN br_err=[br_err,transpose(fltarr(ntime))]			;set error in brightness zero to be equal to zero
        IF keyword_set(verb) THEN print, ' Edge point at '+num2str(redge,dp=3)+' [m]'
        IF keyword_set(debug) THEN stop

	;verify order for shell inversion technique
	order=sort(rtang)
	rtang=rtang[order]
        ff=ff[order]
	bright=bright[order,*]
        IF br_err[0] NE -1 THEN br_err=br_err[order,*]
	IF keyword_set(verb) THEN print, ' Brightness array formatted, radpts = '$
		+num2str(n(bright[*,0])+1,1)+' tpts= '+num2str(n(bright[0,*])+1,1)

	;calculate length matrix for shell inversion
	IF eps EQ 0 THEN BEGIN
            	IF keyword_set(verb) THEN print, ' Performing exact inversion, no smoothing'
		npts=n(rtang)+1			;override npts so that L matrix is square
		calc_geom,rtang,npts,radii,L,ff=ff
		ssvdc,L,w,u,v
		w_recip=fltarr(npts,npts)
		FOR i=0,npts-1 DO IF w[i] NE 0 THEN w_recip[i,i]=(1.0/w[i])
		IF abs(min(w)/max(w)) LT 1.0e-5 THEN print, ' Inversion is ill-conditioned, |w_min/w_max| = '$
			+num2str(abs(min(w)/max(w)))
		invL=v#w_recip#transpose(u)
		emiss=invL#bright
	ENDIF ELSE BEGIN
            	IF keyword_set(verb) THEN print, 'Performing inversion with smoothing'
		IF n_rtang GT npts THEN npts=n_rtang		;check to make sure resolution isn't downgraded
		calc_geom,rtang,npts,radii,L,beta=beta,ff=ff
                IF keyword_set(drt) THEN BEGIN
                    calc_geom,rtang+drt,npts,rp,Lp
                    calc_geom,rtang-drt,npts,rm,Lm
                ENDIF

		;create the identity matrix
		ident = fltarr(npts,npts) & for i = 0,npts-1 do ident(i,i)=1.		
	
		;create the second derivative matrix
		d = -2.*ident
		for i = 1,npts-2 do begin
			d(i-1,i) = 1.
			d(i+1,i) = 1.
		endfor
		
		;handle the corners of the second derivative matrix
		;This boundary condition is just double weighting the second deriv.
		d(0,0) = 1. 
		d(0,1) = -2. 
		d(0,2) = 1. 
		d(1,0) = 1.
		
		;This bc is equivalent to taking bright(npts)=0.
		d(npts-2,npts-1) = 1.

 		; The foldmatrix for symmetrization.
          	IF keyword_set(fmweight) THEN BEGIN
                	IF NOT keyword_set(axis) THEN rcent=.69 ELSE rcent=axis
                        fm=foldmatrix(npts,radii,rcent)
                ENDIF
                        
		;create the least-squares + second derivative matrix
		IF keyword_set(fmweight) THEN mat = (transpose(L)#L + eps*transpose(d)#d+fmweight*transpose(fm)#fm) ELSE $
                   	mat = (transpose(L)#L + eps*transpose(d)#d)

		;do the singular value decomposition
		inv_matrix=la_invert(mat,/double)
		emiss=inv_matrix#(transpose(L)#bright)
                IF br_err[0] NE -1 THEN BEGIN
                    	IF keyword_set(verb) THEN print, 'calculating inversion error'
                    	errmatrix=inv_matrix#transpose(L)
                	em_err=fltarr(npts,ntime)
                        FOR i=0,npts-1 DO BEGIN
                        	FOR j=0L,ntime-1 DO em_err[i,j]=sqrt(total((errmatrix[i,*]*br_err[*,j])^2))	;calculate error via propigation
                        ENDFOR
                ENDIF ELSE em_err=fltarr(npts,ntime)
        ENDELSE

        brchk={br:L#emiss,rad:rtang,br_err:br_err}
        rtang=rtin
        bright=brin
        IF br_err[0] NE -1 THEN br_err=errin
	IF keyword_set(debug) THEN stop
        IF keyword_set(verb) THEN print, ' Inversion successful, returning'
END

PRO shellinvert_error,bright_in,err_br,rtang_in,emiss,radii,num=num,good=good,npts=npts,eps=eps,redge=redge,verb=verb,debug=debug
	
	IF NOT keyword_set(num) THEN num=1000
        shellinvert,bright_in[*,0],rtang_in,em_i,radii,good=good,npts=npts,eps=eps,redge=redge,verb=verb
        emiss=fltarr(n(em_i)+1,num)
        n_br=n(bright_in[*,0])+1
        IF n(err_br) EQ 0 THEN err_br=fltarr(n_br)+err_br
        FOR i=0,num-1 DO BEGIN
        	bright_err=bright_in+randomn(seed,n_br)*err_br
                shellinvert,bright_err,rtang_in,em_i,radii,good=good,npts=npts,eps=eps,redge=redge,verb=verb
                emiss[*,i]=em_i
        ENDFOR
END

;CALC_GEOM was taken from /home/hutch/work/bolo/foil/bolo_svd.pro
;and is included here to make shellinvert more portable
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro calc_geom,rt,npts,rm,s,beta=beta,ff=ff
;
;	this is a program to generate the chord lengths
;	needed for inverted the chord brightness data to get
;	local emissivities
;
;	rt(n_rt) = array of radii for the chordal views of the array
;	npts = number of points used in the inversion
;	s(npts,n_rt) = array of chord lengths between rings
;	re(npts+1) = array of ring radii, i.e. region ring edges.
;	rm(npts) = array of midpoint radii between rings
;       
;
;	create the radius array - from rt(0) - dr/2  to rt(max) + dr/2
;	where dr is half the 'average' distance between data chords
;
n_rt = n_elements(rt)
dr = (rt(n_rt-1) - rt(0))/(npts-1)
re = findgen(npts+1)*dr + rt(0) - dr/2
rm = (re + shift(re,1))/2. & rm = rm(1:*)
s = fltarr(npts,n_rt) & s(*,*) = 0.
IF NOT keyword_set(beta) THEN beta=fltarr(npts)
IF NOT keyword_set(ff) THEN ff=fltarr(npts)+1.0
;
;	determine chord lengths through each slice
;
for j = 0,npts-1 do begin				
	for i = 0,n_rt-1 do begin
		IF beta[j] EQ 0 THEN s(j,i) = (sqrt((re(j+1)^2-rt(i)^2)>0.)-sqrt((re(j)^2-rt(i)^2)>0.))*ff[i] ELSE BEGIN $
                  s[j,i]=exp(beta[j]*sqrt(((max(rm)+dr/2.0)^2-rt[i]^2) > 0.0))*1.0/beta[j]*(sinh(beta[j]*sqrt(((rm[j]+dr/2.0)^2-rt[i]^2) > 0.0))-$
                                                                                            sinh(beta[j]*sqrt(((rm[j]-dr/2.0)^2-rt[i]^2) > 0.0)))
                                                                                
                ENDELSE
	endfor
endfor
;		reflect path length about tangency point
s = 2.*s
s = transpose(s)
; print, s
return
end


;SSVDC was taken from /home/hutch/work/bolo/foil/ssvdc.pro
;and is included here to make shellinvert more portable
pro ssvdc,a,s,u,v
; Single precision SVDecomposition. Uses a LINPACK  routine.
; No tolerance, please. See USER10:[HORNE]DSVDC.TXT for documentation. 
; GT 11/19/93. IHH 30 Aug 2001.
;c     dsvdc is a subroutine to reduce a double precision nxp matrix x
;c     by orthogonal transformations u and v to diagonal form.  the
;c     diagonal elements s(i) are the singular values of x.  the
;c     columns of u are the corresponding left singular vectors,
;c     and the columns of v the right singular vectors.
; That is, u^T a v = diag(s), u diag(s) v^T = a.

ss=size(a)
if (ss(0) ne 2) then begin
        print,'Array must be 2 dimensional, not',ss(0),'dimensional.'
        return
endif
; Correct dimensions, let's go.

x=float(a)    ; array must be float precision
ldx=long(ss(1))
n=long(ss(1))
p=long(ss(2))
ldu=n
ldv=p
work=fltarr(n)
job=long(11)
s=fltarr(min([n+1,p]))
e=fltarr(p)
u=fltarr(ldu,n)
v=fltarr(ldv,p)
info=long(111)
;err=call_external('/home/hutch/idld/linpack/liblinpackidl.so','ssvdc',x,ldx,n,p,s,e,u,ldu,v,ldv,work,job,info)
err=call_external('/usr/local/cmod/lib/liblinpackidl.so','ssvdc',x,ldx,n,p,s,e,u,ldu,v,ldv,work,job,info)
;print,'INFO = ',info
return
end


;BOLO_POWER_PSI1_SVD was taken from /home/hutch/work/bolo/foil/bolo_svd.pro
;and is included here to make shellinvert more portable
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro bolo_power_psi1_svd,shot,rm,time,emiss,ptime,tpow,ok,ipow,rmidp,all=all
;
;	program to calculate the total radiated power
;	
;	inputs/outputs to this program are:
;		shot = shot number
;		rm = array of radii associated with the emissivity
;		time = time array associated with the emissivity
;		emiss = array of radiation emissivities
;
;		ptime = time array associated with the radiated power
;		tpow = total radiated power from psi < 1
;               ipow = power radiated in volume element shell
;               rmidp= volume element radius
;
;	EFIT data is required from analysis tree.
;
; Document by comments IHH Mar 05

mdsopen,'analysis',shot
ok=1
efit_time = mdsvalue('dim_of(\analysis::efit_fitout:volp)',/quiet,status=ok)
if not ok then print,'problem with EFIT times, returning' else begin 
    good_indices=efit_check(good_times,shot=shot)
    if (good_indices(0) eq -1) then begin
        print,'bolo_power_psi1_svd: No good efit times'
        ok=0
    endif
    all=1.0
    IF NOT keyword_set(all) THEN efit_time = efit_time(where(efit_time lt 1.8)) ELSE efit_time=efit_time[good_indices]
        
    rmid = mdsvalue('\analysis::efit_rmid',/quiet,status=ok)
    if not ok then print,'problem with EFIT radii, returning' else begin
        evol = mdsvalue('\analysis::efit_fitout:volp',/quiet,status=ok)
        if not ok then print,'problem with EFIT volume, returning'
    endelse
    

endelse
if(not ok) then return
mdsclose,'analysis',shot

nefit = n_elements(efit_time)
print,'Efit times:',nefit
;ptime = findgen((nefit*4)+2)*.005 + efit_time(0)
; Instead of assuming efit is on 20ms timebase, assume just that it is uniform.
ptime = findgen((nefit*4)+2)*(efit_time(1)-efit_time(0))/4. + efit_time(0)
npt = n_elements(ptime)
nshot=shot
; This number is built in to the choice of 33x33 efits. nefitradii=33
nefitradii=n_elements(rmid(0,*))

;print,'nefitradii, rmid',nefitradii,rmid(0,*)
; locate is an mdsplus-specific idl routine.
emissv = emiss(*,locate(time,ptime))
; Here emissv is (nrm,npt) on the emissivity mesh.

nrm = n_elements(rm)
z = fltarr(nrm) & z(*) = -.01
if shot gt 1000101001 then z(*) = -.0265

;print,'efit_time=',efit_time
rmidrz = efit_rz2rmid(rm,z,efit_time,rho=1,shot=shot,PPP=1,check=1)
; Presumably rmidrz is (nefit,nrm), the rho=1 tells to use normalized rho.
; Whereas rmid is (nefit,nefitradii)
volum = fltarr(npt,nefitradii)
tpow = fltarr(npt)
emissp = fltarr(npt,nefitradii)
rmidp = fltarr(npt,nefitradii)
rmidrzp = fltarr(nrm,npt)
rho = fltarr(nefit,nefitradii)
rhop = fltarr(npt,nefitradii)

;Sloppy fix in case last time slice is screwed up
;BY 04
if rmid(nefit-1,0) eq 0. then begin
print,'Last time step is bad - set equal to penultimate time step'
rmid(nefit-1,*)=rmid(nefit-2,*)
endif
;end of sloppy fix
; rmid was got back from tree in real space. Normalize it.
for i = 0,nefit-1 do rho(i,*) = (rmid(i,*)-rmid(i,0))/(rmid(i,nefitradii-1)-rmid(i,0))
for i = 0,nefitradii-1 do rhop(*,i) = spline(efit_time,rho(*,i),ptime)
rhop(*,0) = .05
; Typical maximum evol value is 0.88. This is m^3.
for i = 0,nefitradii-1 do volum(*,i) = spline(efit_time,evol(*,i),ptime)
for i = 0,nefitradii-1 do rmidp(*,i) = spline(efit_time,rmid(*,i),ptime)
; Here we get the rmidrz on the ptime so rmidrz is (npt,nrm)
for i = 0,nrm-1 do rmidrzp(i,*) = spline(efit_time,rmidrz(i,*),ptime)
; Here emissp is (npt,nrm) but the interpol function gets the emission
; value at values of rmidrzp equal to rhop.
; For interpol to work, the values of rmidrz must be monotonic.
; But they are not if we go past the origin. 
; Experiments seem to show that it just gives the outer values.
; This might work correctly therefore, but ought to be fixed.
; Also, unless rmid is already normalized to run from 0 to 1 as rhop
; does, then this interpolation is incorrect.
; Actually, according to documentation, they are rho-space if rho=1.
for i = 0,npt-1 do emissp(i,*) = interpol(emissv(*,i),rmidrzp(*,i),rhop(i,*))

emissp(*,0) = emissp(*,1) ;take care of interpolation problems at last point


for i = 0,npt-1 do begin
	for j = 0,nefitradii-2 do begin
		tpow(i) = tpow(i) + (volum(i,j+1)-volum(i,j))*.5*$
			(emissp(i,j+1)>0+emissp(i,j)>0)
	endfor
endfor
;
ipow = fltarr(npt,nefitradii)
for i = 0,npt-1 do begin
;	ipow(i,0) = (volum(i,1)-volum(i,0))*.5*(emissp(i,1)>0+emissp(i,0)>0)
	ipow(i,0) = 0
	for j = 1,nefitradii-2 do begin
		ipow(i,j+1) = ipow(i,j) + (volum(i,j+1)-volum(i,j))*.5*$
			(emissp(i,j+1)>0+emissp(i,j)>0)
	endfor
endfor

end
