function get_spred_2006,shot,grating
; 20120709 - ANJ modified to turn off 'extra background' subtraction which causes negative data in disruptions
; retrieves spred data for shots taken after the multi-integration time upgrade to the data acquisition
;
;inputs
;   shot = shot # (long integer)
;   grating = grating id (either 1 or 2) integer
;        1 is the short wavelength grating from 100 to 295 Angstroms
;        2 is the long wavelength grating from ~250 to 1050 Angstroms
;
; keywords  err is returned as 0 if no errs are encountered
;
; returns structure a containing 
;   a.raw(1024,nf)  the raw spred data for 1024 pixels and nf frames
;   a.data(1024,nf) the normalized spred data for 1024 pixels and nf frames
;   a.ierr  error information,0 is ok
;   a.t(nf) the end of each integration period for nf frames
;    a.dt(nf)  the integration period for each frame
;   a.w(1024)  the wavelength for each pixel
;   a.mcp_hv  high voltage on micrchannel plate (volts)
;   a.mcp_hv  high voltage on mcp to fluorescent screen gap (volts)

COMMON PTDATA_COM,IARRAY,RARRAY,ASCII,INT16,INT32,REAL32,REAL64


if shot lt 123500 then begin
 a=get_spred_2000_raw(shot,grating)
 return,a
endif
n_pixels=1024L

if grating eq 1 then point='frostgrat1'
if grating eq 2 then point='frostgrat2'
  if shot ge 124439 and shot le 124533 then begin  ;screw up in digitizer inputs 15-23 May 2006
     if grating eq 1 then point='frostgrat2'
     if grating eq 2 then point='frostgrat1'
  endif
if grating ne 1 and grating ne 2 then begin
   print,'error in the entry for grating--terminating'
   a={ierr:err}
   return,a
endif

gadat2,t,d,point,shot,err=err,ical=0,/header ; header keyword returns real32 by common block
nn=n_elements(d)
d=[0,0,d(0:nn-3)] ; shift in start trigger with new acquisiition software starting with for shots > 124000
if err ne 0 then begin
   print,'error in reading spred data for grating ',grating
   a={ierr:err}
   return,a
endif
tstart=real32(4)/1000.  ;background data taken before tstart  (tstart = 0.0 prior to shot 127000)

;assign gap and mcp voltages
mcp_hv=real32(6)
if grating eq 2 then mcp_hv=real32(7)
gap_hv=real32(10)
if grating eq 2 then gap_hv=real32(11)

nf=n_elements(d)/n_pixels

raw=reform(d,n_pixels,nf)
if shot ge 124439 and shot le 124533 then begin   ;screw up in digitizer inputs 15-23 May 2006
  raw1=raw
  for i=0,511 do begin
   raw1(2*i,*)=raw(2*i+1,*)
   raw1(2*i+1,*)=raw(2*i,*)
  endfor
  raw=raw1
endif
data=float(raw)   ;this array will be normalized later

gadat2,t,times,'frosttimes',shot,err=err
if err ne 0 then begin
   print,'error in reading frame times'
   a={ierr:err}
   return,a
endif
dt=fltarr(nf)
;for i=1,nf-1 do dt(i-1)=times(i)-times(i-1)
for i=0,nf-1 do dt(i)=times(i+1)-times(i)

;time=times
time=times(1:nf)

;set coefficients for 3rd order polynomial fit for pixel vs wavelength
;coefficients from R. Wood 5/19/94
cwps=[291.205,-0.21918,4.73785e-6,8.62097e-9];for high resolution short wavelength grating
cwpl=[172.7,0.938225,1.46646e-4,-3.78789e-8];for low resolution long wavelength grating
if grating eq 2 then cwp=cwpl else cwp=cwps
;setup wavelength array
w=fltarr(n_pixels)
for j=0,n_pixels-1 do w(j)=cwp(0) + cwp(1)*j + cwp(2)*j*j + cwp(3)*j*j*j


;Now begin calculation of normalized data
;first determine integration periods used after t=0

;In this section ignore differences of a few clock
;periods (a few microseconds) in the integration period

dta=round(1.e5*dt)/1.e5 ;round integration period to nearest 0.01 ms
k0=where(time gt 2.e-6)   ;again we're igoring differences of a few ticks

dt1=dta(k0(0))
tint_per=dt1
k2=where(dta ne dt1 and time gt 2.e-6)
 if k2(0) lt 0 then goto,jumpper
 dt2=dta(k2(0))
 tint_per=[tint_per,dt2]
k3=where(dta ne dt1 and dta ne dt2 and time gt 2.e-6)
 if k3(0) lt 0 then goto,jumpper
 dt3=dta(k3(0))
 tint_per=[tint_per,dt3]
k4=where(dta ne dt1 and dta ne dt2 and dta ne dt3 and time gt 2.e-6)
 if k4(0) lt 0 then goto,jumpper
 dt4=dta(k4(0))
 tint_per=[tint_per,dt4]

jumpper:
n_per=n_elements(tint_per)

;now find background frames and subtract from data frames
; background frames are taken before tstart or at end of shot

;data=float(raw)
for i=0,n_per-1 do begin

 kb=where(dta eq tint_per(i) and time le tstart+2.e-6)
   if kb(0) lt 0 then begin
    if dta(nf-1) eq tint_per(i) then begin
     print,'******************************'
     print,'using frame at end of shot for backgrond'
     print,'for integration period = tint_per(i)'
     print,'******************************'
     kb=[nf-1]
     goto,jumpob1
    endif else begin
    if i gt 0 then begin
     print,'**************************'
     print,'using previous background for new background'
     print,'**************************'
     goto,jumpob
    endif
     print,'********************************'
     print,'using maximum pre t=0 integration period for'
     print,'background for integration period = tint_per(i)'
     print,'********************************'
     kpz=where(time le tstart+2.e-6)
     npz=n_elements(kpz)
     kpz=kpz(1:npz-1)
     tmax=max(dta(kpz),mpz)
     kb=[kpz(mpz)]
    endelse
   endif 
 jumpob1:
 bkg=total(raw(*,kb),2)/n_elements(kb)
 jumpob:
 kd=where(dta eq tint_per(i) and time gt tstart+2.e-6)
 for j=0,n_elements(kd)-1 do data(*,kd(j))=data(*,kd(j))-bkg
 for j=0,n_elements(kb)-1 do data(*,kb(j))=data(*,kb(j))-bkg
endfor
;print,'done bkg subtraction'

;eliminate pre t=0 filler frames
case n_per of
 1 : kk1=where(dta eq tint_per(0))
 2 : kk1=where(dta eq tint_per(0) or dta eq tint_per(1))
 3 : kk1=where(dta eq tint_per(0) or dta eq tint_per(1) or dta eq tint_per(2))
 4 : kk1=where(dta eq tint_per(0) or dta eq tint_per(1) or dta eq tint_per(2) or tint_per(3))
 else : print,'too many cases!'
endcase

raw=raw(*,kk1)
data=data(*,kk1)
time=time(kk1)
dt=dt(kk1)
n_frames=n_elements(kk1)

if grating eq 1 then begin
  ; subtract "scattered light" and neutron/gamma background using channels which are beyond the edge of the microchannel plate
  ;for l=0,n_frames-1 do data(*,l) = data(*,l)-total(data(949:1018,l))/70. ; ANJ commented this line 20120706

  ; correct for odd/even pixel gain difference
  z=indgen(512)*2
  data(z,*)=data(z,*)*1.0945
endif

;digitizer was changed after shot 81049
if shot ge 81050 then data=data/4.

a={shot:shot,$
   grating:grating,$
   raw:raw,$
   data:data,$
   w:w,t:time,dt:dt,$
   ierr:err,$
   mcp_hv:mcp_hv,$
   gap_hv:gap_hv }
return,a
end



PRO vuv_load_spec,shot,specbr,lam,time,short=short,long=long,status=status,sigbr=sigbr
	IF keyword_set(short) THEN BEGIN
	 print,'loading data from short grating for shot '+num2str(shot,1)
	 spr=get_spred_2006(shot,1)
	 arrsiz=size(spr.data)
	  n_frames=arrsiz(2)
	  n_pixels=arrsiz(1)

	;retrieve inverse sensitivity curve
         cal_file='/u/west/idl/spred/isler97.dat'
	 openr,cal_lun,cal_file,/get_lun
	  is_short = fltarr(n_pixels)
	   is_long = fltarr(n_pixels)
	  readf,cal_lun,is_short
	  readf,cal_lun,is_long
         free_lun,cal_lun
				
	 integ_time = spr(0).t(1) - spr(0).t(0)
         for i=0,n_frames-1 do spr(0).data(*,i) = spr(0).data(*,i) $
                *is_short*2.0*!pi/integ_time/1.0e13/ (spr(0).mcp_hv/1000.0)^(11.51)  ;hv gain 

         specbr=spr(0).data
         lam=spr.w
	 time=spr.t
	ENDIF

	IF keyword_set(long) THEN BEGIN
	 print,'loading data from long grating for shot '+num2str(shot,1)
	 spr=get_spred_2006(shot,2)
	  arrsiz=size(spr.data)
	   n_frames=arrsiz(2)
	   n_pixels=arrsiz(1)

	;retrieve inverse sensitivity curve
         cal_file='/u/west/idl/spred/isler97.dat'
	 openr,cal_lun,cal_file,/get_lun
	  is_short = fltarr(n_pixels)
	   is_long = fltarr(n_pixels)
	  readf,cal_lun,is_short
	  readf,cal_lun,is_long
         free_lun,cal_lun
				
	 integ_time = spr(0).t(1) - spr(0).t(0)
         for i=0,n_frames-1 do spr(0).data(*,i) = spr(0).data(*,i) $
               *is_long*1.56*!pi/integ_time/1.0e13/ (spr(0).mcp_hv/1000.0)^(11.51)  ;hv gain 

         specbr=spr(0).data
         lam=spr.w
	 time=spr.t
	ENDIF
	
	IF shot eq 136180 THEN BEGIN
	 print,' DANGER !!!! - ANJ horrible (sort of) hack ahead for calculating spectral sigma'
	 tmp=where(time gt 3.1)
	 sigbr=dblarr(n(lam)+1,n(time)+1)
	 for i=0,n(lam) do begin
	  sigbr[i,*]=stddev(specbr[i,tmp],/double)
	 endfor
	ENDIF
END

;+
;NAME:
;	READ_ISPEC_PATH
;
;PURPOSE:
;	This function returns the path to the ISPEC file for 
;	a given impurity, Z.
;
;CALLING SEQUENCE:
;	result=READ_ISPEC_PATH(z)
;	
;INPUTS:
;	z	INT	atomic number of impurity 
;
;OUTPUTS:
;	result	STRING	path
;
;PROCEDURE:
;	Currently, all ISPEC files are kept in /home/mlreinke/idl/vuv/
;	during the initial testing and development of IMPSPEC
;
;	If the Z is not on file, then path='none' will be returned
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 3/2012
;	6/8/12		M.L. Reinke - modified the ispec paths to usr /usr/local/cmod/idl/GENIE/IMPSPEC/ispec/
;
;-

FUNCTION read_ispec_path,z,ver=ver
	GENIE_PATH=getenv('GENIE_PATH') 
	if GENIE_PATH eq '' then GENIE_PATH='/usr/local/cmod/idl/GENIE/'

	if not keyword_set(ver) then ver=1
	 if ver eq 2 then ver='2' else ver=''

	CASE z OF 
		5  : path=GENIE_PATH+'IMPSPEC/ispec_d3d/B.ispec'+ver
		7  : path=GENIE_PATH+'IMPSPEC/ispec_d3d/N.ispec'+ver
		8  : path=GENIE_PATH+'IMPSPEC/ispec_d3d/O.ispec'+ver
		9  : path=GENIE_PATH+'IMPSPEC/ispec_d3d/F.ispec'+ver
		10 : path=GENIE_PATH+'IMPSPEC/ispec_d3d/Ne.ispec'+ver
		18 : path=GENIE_PATH+'IMPSPEC/ispec_d3d/Ar.ispec'+ver
		42 : path=GENIE_PATH+'IMPSPEC/ispec_d3d/Mo.ispec'+ver
		ELSE : begin
			print,'z='+num2str(z,1)+' not known to ispec. cannot find path!'
			path='none'
		end
	ENDCASE
	RETURN,path
END

;+
;NAME:
;	READ_ISPEC_FILE
;
;PURPOSE:
;	This function reads an ASCII file in the "ISPEC" format that
;	can be used in the IMPSPEC line-fitting code.
;
;CALLING SEQUENCE:
;	result=READ_ISPEC_FILE(path)
;
;INPUTS:
;	path	STRING	of the path to the ISPEC file
;
;OPTIONAL INPUTS:
;	z	INT	of the atomic number of interest (calls READ_ISPEC_PATH)
;
;OUTPUTS:
;	result	STRUC	"ISPEC" file formatted for use in other IMPSPEC codse
;		*.LINE#		STRUC	information describing how to fit LINE number "#"
;			*.DLAM		FLTARR	[lam0,lam1] to truncate wavelength range
;			*.SPEC		STRING	spectrometer identifier ('short','long','chromex')
;			*.INST		FLOAT	maximum line width for (< 0 specifies a lorentzian, > 0 gaussian)
;			*.ILINES	INT	number of spectral lines in the group (index=0 is line of interest)
;			*.LAM		FLTARR	[ilines] center wavelengths
;			*.Z		INTARR	[ilines] atomic number
;			*.LABEL		STRARR	[ilines] line label in spectroscopy notation
;			*.ISO		STRARR	[ilines] isoelectronic element symbol
;		*.ELEM		STRING	elemental symbol
;		*.Z		INT	atomic number
;		*.NLINES	INT	number of lines specified in the file
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 3/12
;
;-
PRO readfn,lun,line,str=str ; reads with reads if str is set, readf if not set
 if keyword_set(str) then begin
  line=lun[0]
  if n(lun) ne 0 then lun=lun[1:n(lun)] else lun=''
 endif else readf,lun,line
END
FUNCTION eoffn,lun,str=str
 if keyword_set(str) then return,(lun[0] eq '') and (n(lun) eq 0) else return,eof(lun)
END

FUNCTION read_ispec_file,path,z=z,str=str
	IF keyword_set(z) THEN $
         if z ne 0 $
          then path=read_ispec_path(z) $
          else str=1
	IF path[0] EQ 'none' THEN RETURN,-1
	IF keyword_set(str) THEN lun=strsplit(path,10B,/extract) ELSE openr,lun,path,/get_lun
	line=strarr(1)

	;HEADER
	;-------------------------
	readfn,lun,line,str=str & tmp=strsplit(line,',',/extract)
	 elem=tmp[0]
	 z=int(tmp[1])
	readfn,lun,line,str=str & tmp=strsplit(line,'=',/extract)
	 nlines=int(tmp[1])
	readfn,lun,line,str=str
	WHILE eoffn(lun,str=str) NE 1 $
              AND strmatch(line,'*#DATA#*',/fold_case) EQ 0 $
              DO readfn,lun,line,str=str
	if keyword_set(str) then data_start=lun else point_lun,-lun,data_start ; -ve lun gets pointer for data_start

	;DATA
	;-------------------------
	FOR i=0,nlines-1 DO BEGIN
		if keyword_set(str) then lun=data_start else point_lun,lun,data_start ; +ve lun sets pointer to data_start return to beginning of data
		WHILE eoffn(lun,str=str) NE 1 $
                      AND strmatch(line,'line='+num2str(i,1),/fold_case) EQ 0 $
                      DO readfn,lun,line,str=str ;search for line=i
		readfn,lun,line,str=str & tmp=strsplit(line,'=',/extract) & spec=tmp[1] ; first line after 'line=...' is spec
		readfn,lun,line,str=str & tmp=strsplit(line,',',/extract) & dlam=float(tmp) ; second is wl span
		readfn,lun,line,str=str & tmp=strsplit(line,'=',/extract) & inst=float(tmp[1]) ; third is instrument function
		readfn,lun,line,str=str & tmp=strsplit(line,'=',/extract) & ilines=int(tmp[1]) ; fourth is number of lines
		 lam=fltarr(ilines)
		 iz=intarr(ilines)
		 label=strarr(ilines)
		 iso=strarr(ilines)
		FOR j=0,ilines-1 DO BEGIN
			readfn,lun,line,str=str & tmp=strsplit(line,',',/extract)
			 lam[j]=float(tmp[0])
			 label[j]=tmp[1]
			 iso[j]=tmp[2]
			 tmp=strsplit(label[j],'  ',/extract)
			 iz[j]=read_atomic_charge(tmp[0])
                ENDFOR
	 	;shift=0.5			;allowable shift in angstroms (put in ISPEC file?)
		shift=1.0			;allowable shift in angstroms (put in ISPEC file?)
		result=execute('line'+num2str(i,1)+'={dlam:dlam,spec:spec,inst:inst,shift:shift,ilines:ilines,lam:lam,z:iz,label:label,iso:iso}')
	ENDFOR
	ispec={elem:elem,z:z,nlines:nlines}
	FOR i=0,nlines-1 DO BEGIN
		j=nlines-1-i
		result=execute("ispec=create_struct('line"+num2str(j,1)+"',line"+num2str(j,1)+",ispec)")
        ENDFOR
	if not keyword_set(str) then begin
	 close,lun
	 free_lun,lun
	endif

	return,ispec
END

FUNCTION read_ispec2_file,path,z=z,str=str,verb=verb
	if not keyword_set(verb) then verb=0
	IF keyword_set(z) THEN $
         if z ne 0 $
          then path=read_ispec_path(z,ver=2) $
          else str=1
	IF path[0] EQ 'none' THEN RETURN,-1
	IF keyword_set(str) THEN lun=strsplit(path,10B,/extract) ELSE openr,lun,path,/get_lun
	line=strarr(1)

	;HEADER
	readfn,lun,line,str=str & tmp=strsplit(line,',',/extract) ; 1st line is element and z
	 elem=tmp[0]
	 z=int(tmp[1])
	 nlines=-1
	
	;DATA
	WHILE eoffn(lun,str=str) NE 1 DO BEGIN
	 readfn,lun,line,str=str 
	 if strmatch(strmid(line,0,5),'spec=') then begin
	  if verb gt 0 then print,line
	  
	  nlines=nlines+1
	  ilines=0
	  lam=0
	  iz=0
	  label=''
	  iso=''
	  inst=1.0 ;[A] default instrument (ie minimum) width
	  shift=1.0 ;[A] default wavelength shift for group of lines
	  
	  tmp=strsplit(line,  ',',/extract) & dlam=[float(tmp[1]),float(tmp[2])]
	  tmp2=strsplit(tmp[0],'=',/extract) & spec=tmp2[1]

	  if n(tmp) gt 2 then tmp=tmp[2:*]
	  while n(tmp) gt 0 do begin
	  	tmp=tmp[1:*]
	  	tmp2=strsplit(tmp[0],'=',/extract)
	  	case tmp2[0] of
	  	 'inst' : inst=float(tmp2[1])
	  	 'shift' : shift=float(tmp2[1])
	  	 'else' : print,' !!! unable to parse ispec line: '+tmp[0]
	  	endcase
	  endwhile
	  if verb gt 0 then begin
	   print,' inst='+num2str(inst,dp=3)
	   print,' shift='+num2str(shift,dp=3)
	  endif
	  
	 endif else if nlines ne -1 and line ne '' then begin
	  if verb gt 0 then print,' '+line
	  tmp=strsplit(line,',',/extract)
	  if ilines eq 0 then begin
	   lam=float(tmp[0])
	   label=tmp[1]
	   if n(tmp) gt 1 then iso=tmp[2]
	   tmp=strsplit(label,' ',/extract)
	    iz=read_atomic_charge(tmp[0])
	  endif else begin
	   lam=[lam,float(tmp[0])]
	   label=[label,tmp[1]]
	   if n(tmp) gt 1 then iso=[iso,tmp[2]]
	   tmp=strsplit(label[ilines],' ',/extract)
	    iz=[iz,read_atomic_charge(tmp[0])]
	  endelse
	  ilines=ilines+1
	 endif else if nlines ne -1 and line eq '' then begin
	  result=execute('line'+num2str(nlines,1)+'={dlam:dlam,spec:spec,inst:inst,shift:shift,ilines:ilines,lam:lam,z:iz,label:label,iso:iso}')
	 endif
	ENDWHILE
	ispec={elem:elem,z:z,nlines:nlines}
	FOR i=0,nlines-1 DO BEGIN
		j=nlines-1-i
		addstr="'line"+num2str(j,1)+"',line"+num2str(j,1)
		result=execute("ispec=create_struct("+addstr+",ispec)")
        ENDFOR
	if not keyword_set(str) then begin
	 close,lun
	 free_lun,lun
	endif

	return,ispec
END

;+
;NAME:
;	ISPEC2STRING
;
;PURPOSE:
;	This function converts the ISPEC ASCII file into a single
;	string for storage into the tree.
;
;CALLING SEQUENCE:
;	result=ISPEC2STRING(path)
;
;INPUTS:
;	path	STRING	of the path to the ISPEC file
;
;OPTIONAL INPUTS:
;	z	INT	of the atomic number of interest (calls READ_ISPEC_PATH)
;
;OUTPUTS:
;	result	STRING	of the ISPEC files with string(10B) inserted
;			for carriage returns.
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 3/12
;
;-

FUNCTION ispec2string,path,z=z
	IF keyword_set(z) THEN path=read_ispec_path(z)
	openr,lun,path,/get_lun
	line=strarr(1)
	readf,lun,line
	string=line+string(10B)
	WHILE eof(lun) NE 1 DO BEGIN
		readf,lun,line
		string=string+line+string(10B)
	ENDWHILE
	close,lun
	free_lun,lun
	RETURN,string[0]
END

;+
;NAME:
;	GAUSSIAN_FITS
;
;PURPOSE:
;	This function calculates the gaussian line profile for the sum of 
;	an arbitrary number of gaussians plus an optional DC offset.  The format
;	is made to be used with MPFIT
;
;CALLING_SEQUENCE:
;	result=GAUSSIAN_FITS(x,p)
;
;INPUTS:
;	x	FLTARR	[n_x] of the points for which to calculate the spectra
;	p	FLTARR 	[n_gauss*3+nb] where nb will determine the order of the baseline fit
;
;OUTPUTS:
;	result:	FLTARR  [n_x] of the sum of all the gaussians at each point
;
;PROCEDURE:
;	Gaussians are specified as a*exp(-(x-b)^2/(2*c^2)).  If the parameter values of p
;	are a double precision array, then the result is returned as a double.
;	
;	IF n_elements(p) MOD 3 = 0 then the last three are assume to be a quadraditc baseline
;	IF n_elements(p) MOD 3 = 1 then the last three are assume to be a constant baseline
;	IF n_elements(p) MOD 3 = 2 then the last three are assume to be a linear baseline
;
;MOFICATION HISTORY:
;	Written by: 	copied from hirexsr_fit_spectra.pro in THACO
;	3/12		modified the baseline reference to use the
;			first element of the x scaling vector
;
;-

FUNCTION gaussian_fits, x, p,base=base
	CASE n_elements(p) MOD 3 OF
		0 : BEGIN
			IF n_elements(p) EQ 3 THEN BEGIN ;no baseline
				n_line=1
				basen=3
			ENDIF ELSE BEGIN		;quadradic baseline
				n_line=n_elements(p)/3-1
				basen=2
			ENDELSE
		END
		1 : BEGIN		;constant baseline
			n_line=n_elements(p)/3
			basen=0
		END

		2 : BEGIN		;linear baseline
			n_line=n_elements(p)/3
			basen=1	
		END
	ENDCASE
			
	type=size(p,/type)
	IF type EQ 5 THEN L = dblarr(3,n_line) ELSE L=fltarr(3,n_line)
	FOR i = 0, n_line-1 DO FOR j = 0,2 DO L[j,i] = p[3*i+j]
	
	nx=n(x)+1
	y=fltarr(nx)
	FOR i = 0, n_line-1 DO y += L[0,i]*exp(-(x-L[1,i])^2/(2.*L[2,i]^2))
	
	x0=x[0]	
	CASE basen OF 
		3 : base = 0
		2 : base = p[n(p)] + p[n(p)-1]*(x-x0) + p[n(p)-2]*(x-x0)^2	;
		1 : base = p[n(p)] + p[n(p)-1]*(x-x0)				; 
		0 : base = p[n(p)]
		ELSE :	
	ENDCASE
	y+=base

	RETURN, y
END

;+
;NAME:
;	LORENTZIAN_FITS
;
;PURPOSE:
;	This function calculates the lorentzian line profile for the sum of 
;	an arbitrary number of lines plus an optional DC offset.  The format
;	is made to be used with MPFIT
;
;CALLING_SEQUENCE:
;	result=LORENTZIAN_FITS(x,p)
;
;INPUTS:
;	x	FLTARR	[n_x] of the points for which to calculate the spectra
;	p	FLTARR 	[n_gauss*3+nb] where nb will determine the order of the baseline fit
;
;OUTPUTS:
;	result:	FLTARR  [n_x] of the sum of all the lines at each point
;
;PROCEDURE:
;	Lorentzians are specified as a/pi*c/2/((x-b)^2+(c/2)^2), for the amplitude a (height of the line, not area), center b, and FWHM c.  
;       If the parameter values of p are a double precision array, then the result is returned as a double.
;	
;	IF n_elements(p) MOD 3 = 0 then the last three are assume to be a quadraditc baseline
;	IF n_elements(p) MOD 3 = 1 then the last three are assume to be a constant baseline
;	IF n_elements(p) MOD 3 = 2 then the last three are assume to be a linear baseline
;
;MOFICATION HISTORY:
;	Written by: 	A. N. James 20120712
;                       copied from hirexsr_fit_spectra.pro in THACO
;	3/12		modified the baseline reference to use the
;			first element of the x scaling vector
;
;-

FUNCTION lorentzianfits, x, p,base=base
	CASE n_elements(p) MOD 3 OF
		0 : BEGIN
			IF n_elements(p) EQ 3 THEN BEGIN ;no baseline
				n_line=1
				basen=3
			ENDIF ELSE BEGIN		;quadradic baseline
				n_line=n_elements(p)/3-1
				basen=2
			ENDELSE
		END
		1 : BEGIN		;constant baseline
			n_line=n_elements(p)/3
			basen=0
		END

		2 : BEGIN		;linear baseline
			n_line=n_elements(p)/3
			basen=1	
		END
	ENDCASE
			
	type=size(p,/type)
	IF type EQ 5 THEN L = dblarr(3,n_line) ELSE L=fltarr(3,n_line)
	FOR i = 0, n_line-1 DO FOR j = 0,2 DO L[j,i] = p[3*i+j]
	
	nx=n(x)+1
	y=fltarr(nx)
	;FOR i = 0, n_line-1 DO y += L[0,i]*exp(-(x-L[1,i])^2/(2.*L[2,i]^2))
	FOR i = 0, n_line-1 DO y += L[0,i]/!PI*L[2,i]/2/( (x-L[1,i])^2 + (L[2,i]/2)^2 )
	
	x0=x[0]	
	CASE basen OF 
		3 : base = 0
		2 : base = p[n(p)] + p[n(p)-1]*(x-x0) + p[n(p)-2]*(x-x0)^2	;
		1 : base = p[n(p)] + p[n(p)-1]*(x-x0)				; 
		0 : base = p[n(p)]
		ELSE :	
	ENDCASE
	y+=base

	RETURN, y
END

;+
;NAME:
;	LORENTZIAN_INTEG FITS
;
;PURPOSE:
;	This function calculates the lorentzian line profile for the sum of 
;	an arbitrary number of lines plus an optional DC offset.  The format
;	is made to be used with MPFIT. Each point is calculated as the integral 
;       of a lorentzian from midway to neighboring points.
;
;CALLING_SEQUENCE:
;	result=LORENTZIAN_INTEG_FITS(x,p)
;
;INPUTS:
;	x	FLTARR	[n_x] of the points for which to calculate the spectra
;	p	FLTARR 	[n_gauss*3+nb] where nb will determine the order of the baseline fit
;
;OUTPUTS:
;	result:	FLTARR  [n_x] of the sum of all the lines at each point
;
;PROCEDURE:
;	Lorentzians are specified as a/pi*c/2/((x-b)^2+(c/2)^2), for the amplitude a (height of the line, not area), center b, and FWHM c.
;	The integral is then a difference of arctangents:
;         int_x1^x2 lor dx = a/pi*( atan((x1-b)*2/c) - atan((x2-b)*2/c) )
;       If the parameter values of p are a double precision array, then the result is returned as a double.
;	
;	IF n_elements(p) MOD 3 = 0 then the last three are assume to be a quadraditc baseline
;	IF n_elements(p) MOD 3 = 1 then the last three are assume to be a constant baseline
;	IF n_elements(p) MOD 3 = 2 then the last three are assume to be a linear baseline
;
;MOFICATION HISTORY:
;	Written by: 	A. N. James 20120712
;                       copied from hirexsr_fit_spectra.pro in THACO
;	3/12		modified the baseline reference to use the
;			first element of the x scaling vector
;
;-

FUNCTION lorentzian_integ_fits, x, p,base=base
	CASE n_elements(p) MOD 3 OF
		0 : BEGIN
			IF n_elements(p) EQ 3 THEN BEGIN ;no baseline
				n_line=1
				basen=3
			ENDIF ELSE BEGIN		;quadradic baseline
				n_line=n_elements(p)/3-1
				basen=2
			ENDELSE
		END
		1 : BEGIN		;constant baseline
			n_line=n_elements(p)/3
			basen=0
		END

		2 : BEGIN		;linear baseline
			n_line=n_elements(p)/3
			basen=1	
		END
	ENDCASE
			
	type=size(p,/type)
	IF type EQ 5 THEN L = dblarr(3,n_line) ELSE L=fltarr(3,n_line)
	FOR i = 0, n_line-1 DO FOR j = 0,2 DO L[j,i] = p[3*i+j]
	
	nx=n(x)+1
	y=fltarr(nx)
	dx=x(1:n(x)-1)-x(0:n(x)-1)
	dx=[dx[1],dx]

	;FOR i = 0, n_line-1 DO y += L[0,i]*exp(-(x-L[1,i])^2/(2.*L[2,i]^2))
	;FOR i = 0, n_line-1 DO y += L[0,i]/!PI*L[2,i]/2/( (x-L[1,i])^2 + (L[2,i]/2)^2 )
	FOR i = 0, n_line-1 DO y += L[0,i]*L[2,i]/dx*( atan((x+dx/2-L[1,i])/L[2,i]) - atan((x-dx/2-L[1,i])/L[2,i]) )
;         int_x1^x2 lor dx = a/pi*( atan((x1-b)*2/c) - atan((x2-b)*2/c) )
	
	x0=x[0]	
	CASE basen OF 
		3 : base = 0
		2 : base = p[n(p)] + p[n(p)-1]*(x-x0) + p[n(p)-2]*(x-x0)^2	;
		1 : base = p[n(p)] + p[n(p)-1]*(x-x0)				; 
		0 : base = p[n(p)]
		ELSE :	
	ENDCASE
	y+=base

	RETURN, y
END

;+
;NAME:
;	ISPEC_FIT_LINE
;
;PURPOSE:
;	This function is a lower level code that actually does the
;	fitting of a spectrum based on the LINE information from the ISPEC file
;
;CALLING SEQUENCE:
;	result=IMPSPEC_FIT_LINE(spec,lam,sig,line)
;
;INPUTS:
;	spec	FLTARR	[nlam] of the spectral brightness vs. wavelength
;	lam	FLTARR	[nlam] of the wavelength values
;	sig	FLTARR	[nlam] of the photon statistic uncertainties in spec
;	line	STRUC	substructure of the ISPEC file (see READ_ISPEC_FILE)
;			(*.ilines specifies the # of lines in the fit)
;
;OPTIONAL INPUTS:
;	fitz	INTARR	[nz] of the atomic numbers to include in fit DEFAULT: all listed
;
;KEYWORD PARAMETERS:
;	plot	/plot will display a graphical output of the data, seed and fit
;
;OUTPUTS:
;	result	FLTARR	[ilines*3+2) of the fit coefficients along
;			linear background
;
;OPTIONAL OUTPUTS:
;	error	FLTARR	[ilines*3+2] of the uncertainty in fit coefficients
;
;PROCEDURE:
;	The line.dlam array defines the subset of spectrum that is
;	used in the fit.  Non-linear least-squares fitting is done using
;	MPFIT using the following constraints.
;		line intensitity and DC offset > 0
;		line width > 0 and < line.inst
;		all line widths are the same
;		all lines are tied to the main peak, shifting together
;		constrains shift to < 0.5 angstroms
;		
;	Those impurities which are not specified with the
;	fitz optional input have their peaks fixed at 0.0, thus removing them
;	from the fit.
;	
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 3/12
;
;-

FUNCTION impspec_fit_line,spec,lam,sig,line,error=error,plot=plot,status=status,fitz=fitz

	tmp=where(lam GE line.dlam[0] AND lam LE line.dlam[1])
	x=lam[tmp]
	y=spec[tmp]
	ysig=sig[tmp]
	
	chk=where(ysig LE 0)
	IF chk[0] NE -1 THEN ysig[chk]=abs(mean(y))
	
	isat=where(y+ysig le 0) ; throw out negative values
	if isat[0] ne -1 then begin
	 ;x=x[tmp]
	 ;y=y[tmp]
	 ;ysig=ysig[tmp]
	 ysig[isat]=2*max(y)
	endif

	n_lines=line.ilines
	L=fltarr(3,n_lines)
	 L[2,*]=line.inst ; width
	 L[1,*]=line.lam ; line center
	 L[0,*]=max(y) ; amplitude of all lines start at the max measured level
	 ;L[0]=max([ 0, y(ipt(x,L[1])) ]) ; amplitude of only first line
	 ;L[0,*]=max([ 0, y(ipt(x,L[1])) ]) ; amplitude of all lines start at that of first line

	; prepare initial estimate
	estimate=fltarr(n_lines*3+2)
	FOR i = 0,n_lines-1 DO $
	 FOR j = 0,2 DO $
	  estimate[3*i+j] = L[j,i]
	
	base_line=min(y)
	estimate[n_lines*3+1]=max([ 0, base_line ]) ; set the base estimate
	
	parinfo = replicate({ $
	 value:estimate, $
	 fixed:0, $
	 limited:[0,0],$
	 limits:[0.,0],$
	 tied:''$
	 }, n_elements(estimate) )

	; set line bounds
	for i=0,n_lines-1 do begin
	 parinfo[3*i].limited=[1,  1         ] ; intensity bounds
	 parinfo[3*i].limits =[0.0,max(y+ysig)]
	 IF keyword_set(fitz) THEN BEGIN
	  chk=where(fitz EQ line.z[i])
	  IF chk[0] EQ -1 THEN $
	   parinfo[3*i].fixed=1	;if not in the fitted-z list then fix background line to 0.0
	 ENDIF

	 ;limit the wavelength within a shift to prevent gross wander
	 if i eq 0 then begin
	  parinfo[3*i+1].limited=[1,1]
	  parinfo[3*i+1].limits =line.lam[i]+line.shift*[-1,1]
	 endif else begin
	  IF line.lam[0]-line.lam[i] LT 0 $ ;tie all wavelength shifts together
	   THEN parinfo[3*i+1].tied = 'P(1)+'+strtrim(abs(line.lam[0]-line.lam[i]),1) $
	   ELSE parinfo[3*i+1].tied = 'P(1)-'+strtrim(abs(line.lam[0]-line.lam[i]),1)
	 endelse

	 parinfo[3*i+2].limited=[1,        1          ] ; width bounds
	 parinfo[3*i+2].limits =[line.inst,5*line.inst]
	 parinfo[3*i+2].tied = 'P(2)' ; tie all widths together
	endfor
	
	;fixes DC offset to positive definite
	parinfo[n_lines*3+1].limited=[ 1,   0                ]
	parinfo[n_lines*3+1].limits =[ 0.0, max([0,mean(y)]) ]

	linetype=1 
	case linetype of
	 0 : coefs = mpfitfun('gaussian_fits', x,y,ysig, estimate,perror=error,parinfo=parinfo,status=status,niter=niter,maxiter=2000,quiet=(not keyword_set(verb)) )
	 1 : coefs = mpfitfun('lorentzian_integ_fits', x,y,ysig, estimate,perror=error,parinfo=parinfo,status=status,niter=niter,maxiter=2000,quiet=(not keyword_set(verb)) )
	endcase
	
	IF keyword_set(plot) AND status NE 0 THEN BEGIN
		plot,x,y,/xst,/ysty,xtit='Wavelength [Ang]',ytit='Spectral Brightness [AU]',yr=[0,max(y)*1.1],tit=line.label[0],/nodata
		oploterror,x,y,ysig,psym=3
		xplot=make(x[0],last(x),1000)
		case linetype of
		 0 : begin
		  oplot,xplot,gaussian_fits(xplot,estimate),color=100,linestyle=1
		  oplot,xplot,gaussian_fits(xplot,coefs),color=200
		  FOR i=0,n_lines-1 DO oplot,xplot,gaussian_fits(xplot,coefs[i*3:(i+1)*3-1]),color=30,linestyle=2
		 end
		 1 : begin
		  oplot,xplot,lorentzian_integ_fits(xplot,estimate),color=100,linestyle=1
		  oplot,xplot,lorentzian_integ_fits(xplot,coefs),color=200
		  FOR i=0,n_lines-1 DO oplot,xplot,lorentzian_integ_fits(xplot,coefs[i*3:(i+1)*3-1]),color=30,linestyle=2
		 end
		endcase
		oplot,xplot,coefs[n_lines*3+1]+(xplot-xplot[0])*coefs[n_lines*3],color=30,linestyle=2
        ENDIF
	if status eq 0 then stop, 'mpfitfun error 0: improper input parameters'
	RETURN,coefs
END

;+
;NAME:
;	IMPSPEC
;
;PURPOSE:
;	This is a high level procedure which takes in a shot number
;	and an ISPEC file and computes fit and calculates line brightness
;
;CALLING SEQUENCE:
;	IMPSPEC,shot,ispec,br,coefs
;
;INPUTS:
;	shot		LONG	shot number
;	ispec		STRUC	formatted as per READ_ISPEC_FILE
;	
;OPTIONAL INPUTS:
;	kline		INT	select one line from ISPEC to run (use -1 to select 0th line)
;	fitz		INTARR	of the Z's to include in fit (sent to IMPSPEC_FIT_LINE)
;
;KEYWORD PARAMETERS:
;	plot	/plot sent to IMPSOEC_FIT_LINE to plot results of fit
;	verb	/verb not used yet
;	debug	/debug stops the code at the end of each line fitting
;
;OUTPUTS:
;	br	STRUC	containing the line brightness
;		*.LINE#		FLTARR	[ntime,3] of the brightnes [*,0],
;					time [*,1] and uncertainty [*,2] 
;		*.ELEM		STRING	elemental symbol
;		*.Z		INT	atomic number
;		*.NLINES	INT	number of lines specified in the file		
;	coefs	STRUC	containing the fit coefficients
;		*.LINE#		FLTARR	[line.iline*3+2,3] of the
;					spectral fit coefficients
;		*.ELEM		STRING	elemental symbol
;		*.Z		INT	atomic number
;		*.NLINES	INT	number of lines specified in the file
;
;OPTIONAL OUTPUTS:
;	short	STRUC	containing the short data, used for repeated calls
;	long	STRUC	containing the long data, used for repeated calls
;	
;PROCEDURE:
;	
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 3/12
;
;-

PRO impspec,shot,ispec,br,coefs,tlim=tlim,short=short,long=long,plot=plot,debug=debug,kline=kline,fitz=fitz,verb=verb
	
	nlines=ispec.nlines
	instr=strarr(nlines)
	FOR i=0,nlines-1 DO BEGIN
		j=nlines-1-i
		instr[i]=ispec.(j).spec
	ENDFOR
	IF total(where(instr EQ 'short')) NE -1 THEN need_short=1 ELSE need_short=0
	IF total(where(instr EQ 'long')) NE -1 THEN need_long=1 ELSE need_long=0

	;load spectroscopy data if needed
	IF need_short AND NOT keyword_set(short) THEN BEGIN
		vuv_load_spec,shot,specbr,lam,time,/short,sigbr=sigbr
		 if keyword_set(tlim) then tsel=where( (time ge tlim[0]) and (time le tlim[1]) ) $
		 		      else tsel=indgen( n(time)+1 )
		 short={specbr:specbr(*,tsel),sig:sigbr(*,tsel),lam:lam,time:time(tsel),nt:n(tsel)+1}
        ENDIF
	IF need_long AND NOT keyword_set(long) THEN BEGIN
		vuv_load_spec,shot,specbr,lam,time,/long,sigbr=sigbr
		 if keyword_set(tlim) then tsel=where( (time ge tlim[0]) and (time le tlim[1]) ) $
		 		      else tsel=indgen( n(time)+1 )
		 long={specbr:specbr(*,tsel),sig:sigbr(*,tsel),lam:lam,time:time(tsel),nt:n(tsel)+1}
        ENDIF
	br   ={elem:ispec.elem,z:ispec.z,nlines:ispec.nlines}
	coefs={elem:ispec.elem,z:ispec.z,nlines:ispec.nlines}
	IF NOT keyword_set(kline) THEN BEGIN
		start=0
		stop=nlines-1
        ENDIF ELSE BEGIN
		IF iline EQ -1 THEN iline=0
		start=nlines-1-kline
		stop=nlines-1-kline
	ENDELSE
	FOR i=start,stop DO BEGIN
		j=nlines-1-i
		if keyword_set(verb) then print,'running impspec for j='+num2str(j,1)+' and '+instr[i]+' grating...'
		CASE instr[i] OF 	
			'short' : instrument=short
			'long'  : instrument=long
			else    : print,'unknown instrument: '+instr[i]
		ENDCASE
		jcoefs=fltarr(3*ispec.(j).ilines+2,instrument.nt,2)
		jbr=fltarr(instrument.nt,3)
		jstatus=intarr(instrument.nt)
		FOR k=0,instrument.nt-1 DO BEGIN
			kcoefs=impspec_fit_line(instrument.specbr[*,k],instrument.lam,instrument.sig[*,k],ispec.(j),plot=plot,status=status,error=kerror,fitz=fitz)
			if keyword_set(verb) and (verb gt 1) then begin
				print,'k=',k
				print,' status=',status
				print,' kerror=',kerror
				print,' kcoefs=',kcoefs
			endif
			IF status GT 0 THEN BEGIN
			;IF status EQ 6 OR status EQ 7 THEN BEGIN
				jcoefs[*,k,0]=kcoefs
				jcoefs[*,k,1]=kerror
				jbr[k,0]=sqrt(2*!pi)*kcoefs[0]*kcoefs[2]
				jbr[k,2]=jbr[k,0]*sqrt((kerror[0]/kcoefs[0])^2+(kerror[2]/kcoefs[2])^2)
			ENDIF
			jstatus[k]=status
		ENDFOR
		jbr[*,1]=instrument.time
		
		result=execute("coefs=create_struct('line"+num2str(j,1)+"',jcoefs,coefs)")	
		result=execute(   "br=create_struct('line"+num2str(j,1)+"',jbr,br)")	
		IF keyword_set(debug) THEN stop
	ENDFOR
END

;+
;NAME:
;	RUN_IMPSPEC
;	
;PURPOSE:
;	This procedure runs IMPSPEC and stores the data in the tree
;
;CALLING SEQUENCE:
;	RUN_IMPSPEC,shot,z
;
;INPUTS:
;	shot	LONG 	shot number
;	z	INT	atomic number of impurity to look at
;
;OPTIONAL INPUTS:
;	fitz	INTARR	of the other impurities to include in the fit (sent to IMPSPEC)
;
;OUTPUTS:
;	all data output is stored into the \SPECTROSCOPY::TOP.IMPSPEC tree
;
;PROCEDURE:
;	The ISPEC file is loaded from the tree using READ_ISPEC_TREE
;
;MODIFICAION HISTORY:
;	Written by:	M.L. Reinke - 3/12
;
;-

PRO run_impspec,shot,z,tlim=tlim,fitz=fitz,plot=plot,verb=verb
; eg run_impspec,136180,10,tlim=[3.0,3.02],/plot,verb=2
	IF NOT keyword_set(fitz) THEN fitz=0
	IF NOT keyword_set(verb) THEN verb=2

	;load ISPEC configuration file
	;ispec=read_ispec_tree(shot,z) ; only use trees for saving data, not config!
	;ispec=read_ispec_file(path,z=z) & help,/str,ispec ; parses a file from disk
	ispec=read_ispec2_file(path,z=z,verb=verb) & help,/str,ispec
	ispec_str=ispec2string(path) & print,ispec_str
	print,'loaded ispec file from '+path
	
	;run IMPSPEC
	impspec,shot,ispec,br,coefs,verb=verb,tlim=tlim,fitz=fitz,plot=plot

	; eventually might want to save to a tree instead of disk?
;	mdsopen,'spectroscopy',shot
;	mdsput,'\SPECTROSCOPY::TOP.IMPSPEC.'+ispec.elem+':NLINES','build_with_units($,"")',ispec.nlines
;	FOR i=0,ispec.nlines-1 DO BEGIN
;		path='\SPECTROSCOPY::TOP.IMPSPEC.'+ispec.elem+'.LINE'+num2str(i,1)
;		mdsput,path+':BR','build_signal(build_with_units($1,"ph/s/m^2"),*,build_with_units($2,"seconds"),build_with_units($3,"ph/s/m^2"))',$
;			br.(i)[*,0],br.(i)[*,1],br.(i)[*,2]
;		mdsput,path+':COEFS','build_signal(build_with_units($1," "),*,build_with_units($2,"seconds"),build_with_units($3," "))',$
;			coefs.(i)[*,*,0],br.(i)[*,1],coefs.(i)[*,*,1]
;		mdsput,path+':DLAM','build_with_units($,"Angstroms")',ispec.(i).dlam
;		mdsput,path+':LAM','build_with_units($,"Angstroms")',ispec.(i).lam
;		mdsput,path+':LABEL','build_with_units($," ")',ispec.(i).label
;        ENDFOR
;	mdsclose,'spectroscopy',shot

	fnam='impspec_'+num2str(shot,1)+'.dat'
	save,ispec,ispec_str,br,coefs,filename=fnam & print,'saved output to '+fnam
END

