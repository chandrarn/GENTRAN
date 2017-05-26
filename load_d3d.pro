; this file is used to configure w_spec and impspec for an experiment
; data retrieval routines should also be specified here

; machine specific configs should be implemented through this file, not 
; hard-coded into w_spec and impspec, so that the codes are more easily 
; portable to other experiments.

PRO load_d3d,machine
line= {nam:['C III',	'C VI',	'B IV',	'Ni XXVI',	'Fe XXIV'], $ ; default lines to track
      lam0:[964.5,	179.4,	116.5,	163.3,		190.2], $
      lam1:[989.5,	185.0,	119.6,	167.2,		193.6], $
      spec:['long',	'short','short','short',	'short'],$
      mult:[0.1,	1.0,	1.0,	1.0,		1.0], $
      plot:[1,		1,	1,	1,		1]}

line_Ne= {nam:['Ne I',	'Ne III','Ne VIII',	'Ne IX','Ne X'], $ ; default lines to track
          lam0:[620.8,	827.2,	 281.3,		264.6,	126.7], $
          lam1:[632.4,	838.5,	 285.7,		270.8,	126.7], $
          spec:['long',	'long',	 'long',	'long',	'short'],$
          mult:[1,	0.1,	 1.0,		2.0,	1.0], $
          plot:[1,	1,	 1,		1,	1]}

machine={name:'d3d', $
	 inst:['short','long','divspec'], $ ; spectroscopic instruments
	 load:[1,      1,     1     ], $
	 plot:[1,      1,     1     ], $
	 nch: [1,      1,     16    ], $
	 timetr:['ip','dens','te','beams','ech'], $ ; time traces
	 line:line, $
	 elem:'D,B,C,Fe,Ni,Cr' $ ; default elements to overplot
	}

END

function load_d3d_shot
	; retrieve the current shot
	return,mdscur_shot('d3d')
end

function get_spred,shot,grating
; adapted from Phil West's /u/west/idl/spred/public_version/get_spred_2006.pro
; 20120709 - ANJ modified to turn off 'extra background' subtraction which causes negative data in disruptions
;
;retrieves spred data for shots taken after the multi-integration time upgrade to the data acquisition
;
;inputs
;   shot = shot # (long integer)
;   grating = grating id (either 1 or 2) integer
;        1 is the short wavelength grating from 100 to 295 Angstroms
;        2 is the long wavelength grating from ~250 to 1050 Angstroms

; keywords  err is returned as 0 if no errs are encounter
;returns structure a containing 
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
if err ne 0 then begin
   print,'error in reading spred data for grating ',grating
   a={ierr:err}
   return,a
endif

;assign gap and mcp voltages
if grating eq 1 then begin
 mcp_hv=real32(6)
 gap_hv=real32(10)
endif else begin
 mcp_hv=real32(7)
 gap_hv=real32(11)
endelse

tstart=real32(4)/1000. ;[s] background data taken before tstart  (tstart = 0.0 prior to shot 127000)
nn=n_elements(d)
nf=n_elements(d)/n_pixels
d=[0,0,d(0:nn-3)]    ;shift in start trigger with new acquisiition software starting with for shots > 124000
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
if grating eq 2 $
 then cwp=[172.7,   0.938225,1.46646e-4,-3.78789e-8] $ ; for low resolution, long wavelength grating
 else cwp=[291.205,-0.21918, 4.73785e-6, 8.62097e-9]   ; for high resolution, short wavelength grating
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
 if n_elements(kb) gt 1 then bkg=total(raw(*,kb),2)/n_elements(kb) $
                        else bkg=raw(*,kb)
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
 else : print,'too many cases! ',n_per
endcase

raw=raw(*,kk1)
data=data(*,kk1)
time=time(kk1)
dt=dt(kk1)
n_frames=n_elements(kk1)

if grating eq 1 then begin
;  subtract "scattered light" and neutron/gamma background using channels
;  which are beyond the edge of the microchannel plate
  ;for l=0,n_frames-1 do data(*,l) = data(*,l)-total(data(949:1018,l))/70. ; ANJ commented this line 20120706

;  correct for odd/even pixel gain difference
  z=indgen(512)*2
  data(z,*)=data(z,*)*1.0945
endif

;digitizer was changed from 12-bit to 14-bit after shot 81049
if shot ge 81050 then data=data/4.
; daq routines were only saving 12-bits until shot 149898 (~15 years later...)
; but digitizer was correctly scaling, so data would 'wrap' around 4096

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

PRO load_d3d_short,shot,d
	 spr=get_spred(shot,1)
	 if spr.ierr eq 0 $
	  then status=1 $ ; success
	  else status=0 ; error
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
				
	 integ_time = spr.t(1) - spr.t(0)
         for i=0,n_frames-1 do spr.data(*,i) = spr.data(*,i) $
                *is_short*2.0*!pi/integ_time/1.0e13/ (spr.mcp_hv/1000.0)^(11.51)  ;hv gain 

	msg=' VMCP='+num2str(spr.mcp_hv,1)+',VGAP='+num2str(spr.gap_hv,1)

	d={t:spr.t, $
	   l:spr.w,llabel:n2g('lambda')+' ['+string(197B)+']', $
	   d:spr.data,dlabel:'[1e13 photons/cm!E2!Nsr/'+string(197B)+'/s]', $
	   title:'SPRED short grating',s:status,p:status,index:0,msg:msg}
END

PRO load_d3d_long,shot,d
	 spr=get_spred(shot,2)
	 if spr.ierr eq 0 $
	  then status=1 $ ; success
	  else status=0 ; error
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
				
	 integ_time = spr.t(1) - spr.t(0)
         for i=0,n_frames-1 do spr.data(*,i) = spr.data(*,i) $
               *is_long*1.56*!pi/integ_time/1.0e13/ (spr.mcp_hv/1000.0)^(11.51)  ;hv gain 

	msg=' VMCP='+num2str(spr.mcp_hv,1)+',VGAP='+num2str(spr.gap_hv,1)

	d={t:spr.t, $
	   l:spr.w,llabel:n2g('lambda')+' ['+string(197B)+']', $
	   d:spr.data,dlabel:'[1e13 photons/cm!E2!Nsr/'+string(197B)+'/s]', $
	   title:'SPRED long grating',s:status,p:status,index:0,msg:msg}
END

PRO load_d3d_ip,shot,d
	pt='Ip'
	gadat, time, ip, pt, shot, err=err, /alldata ; err=0 on success
	if err eq 0 then status=1 else status=0 ; stutes=1 on success
	d={t:time/1000.0, $
	   d:ip/1.0e6,dlabel:'[MA]', $
	   pt:pt,title:pt,oplot:0,s:status,yr:[0.0,0.0]}
END

PRO load_d3d_dens,shot,d
	pt='denr0'
	gadat, time, dens, pt, shot, err=err, /alldata
	if err eq 0 then status=1 else status=0 ; stutes=1 on success
	d={t:time/1000.0, $
	   d:dens/1.e13,dlabel:'[1e13/cc?]', $
	   pt:pt,title:pt,oplot:0,s:status,yr:[0.0,0.0]}
END

PRO load_d3d_te,shot,d
	pt='ece40'
	gadat, time, te, pt, shot, err=err, /alldata
	if err eq 0 then status=1 else status=0 ; stutes=1 on success
	d={t:time/1000.0, $
	   d:te,dlabel:'[keV]', $
	   pt:pt,title:pt,oplot:0,s:status,yr:[0.0,0.0]}
END

PRO load_d3d_beams,shot,d
	pt='beams'
	gadat, time, beams, pt, shot, err=err, /alldata
	if err eq 0 then status=1 else status=0 ; stutes=1 on success
	d={t:time/1000.0, $
	   d:beams,dlabel:'[MW]', $
	   pt:pt,title:pt,oplot:0,s:status,yr:[0.0,0.0]}
END

PRO load_d3d_ech,shot,d
	pt='echpwrc'
	gadat, time, ech, pt, shot, err=err, /alldata
	if err eq 0 then status=1 else status=0 ; stutes=1 on success
	d={t:time/1000.0, $
	   d:ech,dlabel:'[MW]', $
	   pt:pt,title:pt,oplot:1,s:status,yr:[0.0,0.0]}
END


pro load_d3d_divspec,shot,d ; multichannel divertor spectrometer
	a=get_divspec_anj(shot)
	t=a.t

	dch=a.track.d/2.0^16 ; [nl nt nch]
	lch=a.track.w ; [nl nch]

	dim=size(dch,/dim)
	 nl=dim[0]
	 nt=dim[1]
	 nch=dim[2]

	tmp=where(t lt 0)
	if tmp[0] ne -1 then tmp=tmp[where(tmp gt 0)]
	if tmp[0] ne -1 then begin
	 bg=total(reform(dch[*,tmp,*]),2)/n_elements(tmp)
	 for i=0,nt-1 do dch[*,i,*]=reform(dch[*,i,*])-bg
	endif

	ch_nam=strarr(nch)

	d=total(dch,3)/nch
	l=total(lch,2)/nch

	c={slit:a.slit,dt:a.integ_t,cwl:a.wave_set,grmm:a.grmm}
	cr=string(10B)
	msg=cr+'  lambda_ctr='+num2str(c.cwl,dp=1) $
	   +cr+'  grating='+num2str(c.grmm,1) $
	   +cr+'  int_time='+num2str(c.dt,1) $
	   +cr+'  slit='+num2str(c.slit,dp=1)

	status = a.ierr eq 0

	d={t:t, $
	   l:l,llabel:n2g('lambda')+' ['+string(197B)+']', $
	   d:d,dlabel:'[AU]', $
	   dch:dch,lch:lch,ch_nam:ch_nam, $
	   c:c, $
	   title:'divspec',s:status,p:status,index:0,msg:msg }
end






;*******************************
FUNCTION Defined, var
;*******************************
;~ Returns 1 if VAR is defined, and 0 otherwise. [Programming]
	On_Error, 2
	if n_params() NE 1 $
	 then message, /noname, "Usage: Defined(var)"

	s = size(var)
	return, (s[s[0]+1] NE 0)
end  ;; defined


;+
;**************************************************
FUNCTION int2bin, $
data, nibble=nibble, nb = nb, bitarray=bitarray, formatStr=formatStr
;**************************************************
;~ Converts integer-type DATA into its binary representation [Conversion, String]
;
;  Data must be an integer-type scalar or 1D array.
;
;  By default, a string is returned consisting of the binary
;  representation of DATA, in groups of eight bits, most significant
;  bit leftmost.  If NIBBLE is set, the bits are instead displayed in
;  groups of 4.
;
;  You can restrict the result to NB bytes (or, if NIBBLE is set, NB nibbles) by
;  setting NB to a positive integer value. (This option is notavailable currently
;  if BITARRAY is set.)
;
;  If BITARRAY is set, a byte array is returned, in which each
;  element is the corresponding binary bit of DATA, BUT IN REVERSE ORDER.
;  The array has 8, 16, 32 or 64 elements, when DATA is a byte, integer,
;  integer long or long-long, respectively.
;
;  You can create a string representation of fixed-point integers, e.g., "1100.01"
;  by supplying FORMATSTR in the form 'bN.M', where N is the total number of binary
;  bits returned and M is the number of bits displayed to the right of 
;  the binary point.
;  The result is padded with leading 0's to exactly N binary digits,
;  if necessary.
;
;  Examples:
;
;  help, Int2Bin('8f'x)
;   <Expression>    STRING    = '00000000 10001111'
;
;  print, int2bin('8f'x, nb=1)
;   10001111
;
;  help, Int2Bin(fix('8f'x), /bitarray)
;    <Expression>    BYTE      = Array[16]
;  print, Int2Bin(fix('8f'x), /bitarray)
;    0   0   0   0   0   0   0   0   1   0   0   0   1   1   1   1
;
;  print, int2bin('8f'x, format='b5.1')
;   0111.1
;-

;  Modification history:
;  Dec 30 1993  W. Rigby.  Original by C. D. Pike, 7-Oct-93, with my minor mods.
;  ...
;  Jan 14 2005  stripped down before posting.
 
On_Error, 2
MsgPrefix = "[Int2Bin] "
maxBits = 64

if n_params() NE 1 $
 then Message, /noname, "Usage: Int2Bin, data {[/nibble || /bitarray] [nb]} || {formatStr=}"

doBitArray = Keyword_Set(BITARRAY)
doNibble = Keyword_Set(NIBBLE)
doFormatStr= Defined(formatStr)
doNB = Defined(nb)

;;if (CheckArgument(data, "integerType", errMsg)) then $
;;  Message, /noname, MsgPrefix + "DATA" + errMsg
;;
;;if (CheckArgument(n, /optional, "integerType", /positive, errMsg)) then $
;;  Message, /noname, MsgPrefix + "N" + errMsg
;;
;;if (CheckArgument(formatStr, /optional, "string", /scalar, errMsg)) then $
;;  Message, /noname, MsgPrefix + "FORMATSTR" + errMsg
;;
;;if (CheckArgument(nb, /optional, "integertype", /scalar, minvalue=1, 
;;maxValue=8, errMsg)) then $
;;  Message, /noname, MsgPrefix + "NB" + errMsg

if doNibble AND doBitArray $
 then Message, /noname, MsgPrefix + "BITARRAY can't be set when NIBBLE is."

if doNB and doBitArray $
 then Message, /noname, MsgPrefix + "The keyword NB not supported when BITARRAY is set"

if doFormatStr then begin
	if doBitArray $
	 then Message, /noname, MsgPrefix + "BITARRAY can't be set when FORMATSTR is supplied."
	if doNibble $
	 then Message, /noname, MsgPrefix + "NIBBLE can't be set when FORMATSTR is supplied."
	if Defined(nB) $
	 then Message, /noname, MsgPrefix + "NB can't be set when FORMATSTR is supplied."

	msg = "FORMATSTR must be in the form 'bn.m' with (n GT 0) and (0 LE m LE n)"

	if StrLowCase(StrMid(formatStr, 0, 1)) NE 'b' $
	 then Message, /noname, MsgPrefix + msg

	p = StrPos(formatStr, ".")  ;;
	if p LT 2 then Message, /noname, MsgPrefix + msg
	n = StrLen(formatStr)
	if p EQ n-1 then Message, /noname, MsgPrefix + msg

	width = StrMid(formatStr, 1, p-1)
	;;   if (MyNot(StrIsNumber(width, /int))) then $
	;;     Message, /noname, MsgPrefix + msg

	nBinary = StrMid(formatStr, p+1, n-1-(p+1)+1)
	;;   if (MyNot(StrIsNumber(nBinary, /int))) then $
	;;     Message, /noname, MsgPrefix + msg

	width = fix(width)
	nBinary = fix(nBinary)

	 if ((width LE 0) OR (nbinary LT 0)) then Message, /noname, MsgPrefix + msg

	 if width LT nBinary then Message, /noname, MsgPrefix + msg
	 doBitArray = 0   ;; force a string result
endif

;;nd = ndim(data)
nd = (size(data))[0]
case nd of
 0 : ;; scalars OK
 1 : ;; 1D arrays OK
 else : Message, /noname, MsgPrefix + "DATA must be a scalar or a one-dimensional array."
endcase

; make a mask for each of the 64 possible bit positions
; mask[0] has the MSB set and mask[maxBits-1] has bit 0 set.
mask = Rotate(2ULL^indgen(maxBits), 2)
case nd of
 0 : out = ((ulong64(data) AND mask) NE 0)
 1 : begin
	ndata = n_elements(data)
	nmask = n_elements(mask)

	result = strarr(ndata)
	out = ulon64arr(ndata, nmask)       ;; out[i] is 1 if bit MSB-i is set
	for i = 0, nmask-1 do out[*,i] = ((ulong64(data) AND mask[i]) NE 0) ; thus out[0] <--> MSB, out[ 
     end
endcase

;  trim output depending on input type
switch size(data, /tname) of 
 'UNDEFINED' : Message, /noname, MsgPrefix + "DATA is undefined!"
 
 'BYTE': begin
		ii1 = maxbits-1    ;; last element in OUT --> bit 0 in DATA, so we want
		ii0 = ii1 - 8 + 1  ;; the last 8 bits of OUT
		if doBitArray then begin
		 case nd of
		  0 : result = byte(out[ii0:ii1])
		  1 : result = byte(out[*,ii0:ii1])
		 endcase
		endif else begin
		 if doNibble $
		  then fmtstr = '(2(4I1, 1x))' $
		  else fmtstr = '(8I1)'

		 case nd of
		  0 : result = String(format=fmtstr, out[ii0:ii1])
		  1 : for i = 0, ndata-1 do $
		       result[i] = String(format=fmtstr, out[i, ii0:ii1])
		 endcase
		endelse
		break
	end

 'UINT' :
 'INT': begin
		ii1 = maxbits-1
		ii0 = ii1 - 16 + 1   ;; we want the last 16 bits of OUT
		if doBitArray then begin
		 case nd of
		  0 :  result = byte(out[ii0:ii1])
		  1 :  result = byte(out[*, ii0:ii1])
		 endcase
		endif else begin
		 if doNibble $
		  then fmtstr = '(4(4I1, 1x))' $
		  else fmtstr = '(2(8I1, 1x))'

		 case nd of
		  0 : result = String( format=fmtstr, out[ii0:ii1])
		  1 : for i = 0, ndata-1 do $
		       result[i] = String( format=fmtstr, out[i,ii0:ii1] )
		 endcase
		endelse
		break
	end

 'ULONG' :
 'LONG': begin
		ii1 = maxbits-1
		ii0 = ii1 - 32 + 1   ;; we want the last 32 bits of OUT

		if doBitArray then begin
		 case nd of
		  0 : result = byte(out[ii0:ii1])
		  1 : result = byte(out[*,ii0:ii1]) ; just to be explicit about what's going on
		 endcase
		endif else begin
		 if doNibble $
		  then fmtstr = '(8(4I1, 1x))' $
		  else fmtstr = '(4(8I1, 1x))'
	
		 case (nd) of
		  0 : result = String(format=fmtstr, out[ii0:ii1] )
		  1 : for i = 0, ndata-1 do $
		       result[i] = String(format=fmtstr, out[i,ii0:ii1])
		 endcase
		endelse
		break
	end

 'LONG64' :
 'ULONG64' : begin
		ii1 = maxbits-1
		ii0 = ii1 - 64 + 1 ;; i.e., 0
		if doBitArray then begin
		 case nd of
		  0 : result = byte(out[ii0:ii1])
		  1 : result = byte(out[*,ii0:ii1]) ; just to be explicit about what's going on
		 endcase
		endif else begin
		 if doNibble $
		  then fmtstr = '(16(4I1, 1x))' $
		  else fmtstr = '(8(8I1, 1x))'

		 case (nd) of
		  0 : result = String( format=fmtstr, out)
		  1 : for i = 0, ndata-1 do $
		       result[i] = String( format=fmtstr, out)
		 endcase
		endelse
		break
	end

 else : Message, /noname, MsgPrefix + "This can't happen."  ;; new integer data type?
endswitch

if doNB then begin
 if doBitArray then begin
  ;; easy to add, but no need currently
  Message, /noname, MsgPrefix + "This can't happen."
 endif else begin
  ;; result is a string, with bytes (or nibbles) separated by spaces
  nr = n_elements(result)        ;; = n_elements(data)

  for j = 0, nr-1 do begin
   ;; strArr = SplitString(result[j])   ;; split at whitespace
   strArr = StrSplit(result[j], " ", /extract)
   n = n_elements(strArr)            ;; number of bytes or nibbles in the string array

   if (nb LT n)  then begin
    r = strArr[n-1]         ;; least significant byte or nibbles
    for i = 1, nb-1 do r = strArr[n-1-i] + " " + r  ;; add bytes or nibbles from right to left
    result[j] = r   ;; OK to index scalars (a=4, a[0] = 3) -- since when?
   endif ;; otherwise just return the whole string
  endfor

 endelse ; doBitArray
endif ; doNB

if doFormatStr then begin
 ;; put the bitarray into n.m format. remember: lsb is in the last element of RESULT
 ;; remove all the (byte or nibble) spaces first -- too messy to try to keep them.
 result = StrCompress(result, /remove_all)

 ;; trim or pad to WIDTH bits
 n = StrLen(result)
if width LE n $
 then temp = StrMid(result, n-width, width) $
 else temp = string(make_array(n, value=(byte('0'))[0]))
 ;    temp = MakeString(width-n, char='0') + result
 n = width

 result = StrMid(temp, 0, width-nbinary) + "." + StrMid(temp, n-nbinary, nbinary)
endif

return, result

end   ;  Int2Bin


PRO split4ascii,binary32,ascii2
; from: /d2/divertor/mds/split4ascii.pro
;split the 32-bit output of Int2Bin into four separate bytes.
;convert 8-bit binary pattern of each byte into an ascii character
;each pair of ascii characters forms one ROI label.
byte4  = strarr(4)
ascii4 = strarr(4)
ascii2 = strarr(2)
;binary = '00110001 00110110 01010000 01010011'
byte4 = strsplit( binary32,/extract)
for i=0,3 do begin 
  bin2dec,byte4(i),out,/quiet
  ascii4(i)=string(byte(out))
endfor
ascii2(0) = ascii4(0) + ascii4(1)
ascii2(1) = ascii4(2) + ascii4(3)

return
end


pro bin2dec,inp,out,quiet=quiet
;from: /u/mclean/idl/utils/bin2dec.pro 
;+
; Project     : SOHO - CDS     
;                   
; Name        : BIN2DEC
;               
; Purpose     : Convert binary representation to decimal integer.
;               
; Explanation : The binary representation of a decimal number is converted
;               to a decimal integer and can be displayed or returned or 
;               both or neither.
;               
; Use         : IDL> bin2dec, binary [, decimal, /quiet]
;    
; Inputs      : binary - the binary representation to convert. It can either
;                        be a string of zeros and ones or an array with each
;                        element a zero or one.
;                        eg bin2dec,'01010101'    or
;                           bin2dec,['1','0','1','0','1','0','1','0']    or
;                           bin2dec,[1,0,1,0,1,0,1,0]
;                        The MSB is assumed to come first
;
; Opt. Inputs : None
;               
; Outputs     : See below
;               
; Opt. Outputs: decimal - the decimal integer equivalent of the input.
;               
; Keywords    : quiet - unless given the decimal number is printed to the
;                       terminal
;
; Calls       : None
;               
; Restrictions: Input must be a string or an array of integers or strings.
;               
; Side effects: None
;               
; Category    : Utils, Numerical
;               
; Prev. Hist. : None
;
; Written     : C D Pike, RAL, 7-Oct-93
;               
; Modified    : 
;
; Version     : Version 1, 7-Oct-93
;-            

;  initialise output and find size of input
out = 0L
n = n_elements(inp)

;  is input a single string? If so break it up.
; replaced calls to /u/mclean/idl/utils/datatype.pro with idl size intrinsic
if n eq 1 and size(inp,/type) eq 7 then begin ; 7 is string
   n = strlen(inp)
   inp1 = intarr(n)
   for i=0,n-1 do inp1(i) = strmid(inp,i,1)
   inp = inp1
endif

if n gt 1 and size(inp,/type) eq 7 then inp = fix(inp);  is input a multiple string?
x = reverse(byte(inp));  switch array around for convenience
for i=0,n-1 do out = out + long(x(i))*2L^i;  calculate integer
if not keyword_set(quiet) then print,out;  if not silenced then report

end



function grating_root,x
; find root of grating equation to determine theta_i (x) corresponding to 
; angle of rotation for selected wavelength

; pass grating variables from calc_mds_wave in COMMON block
	common wave_block,wave_set,a,B
	return,a*(sin(x)+sin(x+B))-wave_set
end

pro calc_mds_wave,wave_calc,npixels,shot
; from: /d2/divertor/mds/calc_mds_wave.pro
; wavelength calibration routine for MDS, Multichannel Divertor Spectrometer
; 
; G. McKee, 8/13/96, General Atomics
; N. Brooks, 1/19/98, General Atomics - generalized to track with npixels<770
; N. Brooks, 1/22/07  modify for smaller pixel size of new CCD
; N. Brooks  6/18/07  reversal of pixel order after repair of CCD

common wave_block,wave_set,a,B  ;reinstall today

; spectrometer parameters, all double precision.
; Consider chip dimensions as possible source of uncertainty (or others)
a=1.0d/1200.0d*1d7	; grating groove separation in Angstroms (1200 g/mm)
B=9.96d*!pi/180.0	; separation angle of optical axes of focussing/collimating mirrors, in radians (9.96 deg)
f=1330.0d		; focal length, in mm
N=double(npixels)	; number of pixels, floating, double prec.

if shot gt 127425 $ ; d=width of CCD chip, in m
 then d=npixels*0.013d $;PI VersArray 1024BFT camera CCD chip	
 else d=npixels*0.0225d ;Wright camera CCD chip	

;wave_cal = [4471.48,5015.68,5875.62,6678.15,7065.19] 
;deuterium  [4339, 4860,6560,9545]
;argon      [8115, 9123 ]

	; Wavelengths of observed HeI lines in Angstroms

;disp_cal = [0.13196,0.13019,0.12703,0.12366,0.12189]
	; dispersion, Angstroms/pix
	; Calculated disp. at each of the above HeI lines. These lines are used
	; since they can be easily observed with plasma or even-better with 
	; helium glow discharge spectra to determine wavelength offset.
	; Dispersions are calculated from grating equation:
	;		m(lambda)=a(sin(theta_i)+sin(theta_i+B))
	; with B the angular separation of the focussing/collimating mirrors
	; (optical axes) as seen from the grating, theta_o=theta_i+B

wave_cal = [ 3451.3, 4647.4, 6578.05, 9123.0 ]  ; tokamak data from 01/25/07
;offset_cal =[-0.33,-0.96,-2.56,-3.84,-4.32]	; offset of HeI lines
offset_cal =[  -2.0,   -4.0,   -6.6,   -11.8 ]	; offset of BII, CIII, CII, ArI
line_coeff=poly_fit(wave_cal,offset_cal,1)	; linear fit to offset
offset=wave_set*line_coeff(1)+line_coeff(0)	; calculated wavelength offset

; shot-dependent (data acquired in 1997) correction to wavelength array:
; Note: this is a constant correction, though a linear correction may be
; necessary once a sufficient data base is acquired (GM, 5/97).

x=[0.00,0.25,0.5] ;initial guess vector for solution of 'grating_root' function

; calculate rotation angle of grating, theta_i, for given wavelength setting,
; from grating equation
theta_i=fx_root(x,'grating_root',/double,tol=1e-6)

; calculate wavelength array to be loaded into MDS data structure
; during repair of CCD at PI, chip was flipped. Correct wavelength display
; requires reversal of pixel order 
if shot lt 129000 $
 then wave_calc=a*(sin(theta_i)+sin(theta_i+B+d/f*(dindgen(npixels)/N-0.5d))) + offset $ ; before repair
 else wave_calc=a*(sin(theta_i)+sin(theta_i+B+d/f*( ( N - 1 - dindgen(npixels) )/N - 0.5d))) + offset ; after repair

return
end




function get_divspec_anj,shot		
; copied from ~mclean/idl/get_divspec_adam.pro

common PTDATA_COM,IARRAY,RARRAY,ASCII,INT16,INT32,REAL32,REAL64
common wave_block,wave_set,a,B
common inst_prof, ifun, offset, he_prof, he_offset, da_prof, da_offset, $
	ne_prof,ne_offset,ampratio

;retrieve data from the ADPCCD pointname for the divertor spectrometer 
print,'Retrieving CCD data for shot: '+strtrim(shot,2)

pointname='DIVSPEC' ;pointname='ADPMDS'
ical=0  ; retrieve digitizer counts
;get_data,pointname,shot,ical,time,data,rc
gadat2,time,data,pointname,shot,err=err,ical=ical,/header ; header keyword returns real32 by common block

if err gt 0 then begin
 print,'Error from PTDATA: '+strtrim(err,2)
 print,'No CCD data retrieved for: '+strtrim(shot,2)
 return,{ierr:1l,err:err}
endif
ntracks = fix(real32[21])	; number of tracks

label = strarr(ntracks)
for i=0,ntracks/2 - 1 do begin 
	binary32 = int2bin( ascii[6 + i] )
	split4ascii,binary32,ascii2

	; Check for little endian. If so, reverse retrieved label and comment strings.
	; Required at least since Hydra was decommissioned.
	machine_little_endian = (BYTE(1,0,1))[0]
	if machine_little_endian then ascii2 = string(reverse(byte(reverse(ascii2))))

	j = 2*i
	label(j) = ascii2[0]
	k = 2*i + 1
	label(k) = ascii2[1]
endfor

comment=''
for i = 44,58 do begin  ;15 words or 60 characters will fit in title space
	binary32 = int2bin( ascii[i] )
	split4ascii,binary32,ascii2

	; Check for little endian. If so, reverse retrieved label and comment strings.
	; Required at least since Hydra was decommissioned.
	machine_little_endian = (BYTE(1,0,1))[0]
	if machine_little_endian then ascii2 = string(reverse(byte(reverse(ascii2))))
  
	comment = comment + ascii2[0] + ascii2[1]
endfor

ntimes = int16[4]  ;number of frames
npixels_total = int16[5] ; number of total pixels per frame
npixels = npixels_total/ntracks
slit = real32[8] ; input slit width in microns
wave_set = real32[7]	; wavelength setting on MDS (used for wave cal.)

d = reform(data,npixels,ntracks,ntimes)  ; data values in digitizer counts
t = real32[23:23+ntimes-1]  ; times of frame readouts in seconds
                            ; t(0) = .100s, where fr(0) contains big dark noise
			    ; from last ten minutes
t_sync_interval = t[1] - t[0]
IF shot GE 141810 $
 THEN t_trigger =  -10.000 $
 ELSE t_trigger =  -0.050 
t_delay = real32[23+ntimes]  ;offset between 6C (-50ms) or 70 (-10s) trigger 
                ;and first sync pulse 
		;t_delay does not respond to changes in its value in GUI   
exp_t = real32[6]	; frame time interval (exposure time)
t = t + t_trigger  + exp_t/2.  ;t_center of each frame 

dum = {d:fltarr(npixels,ntimes),w:fltarr(npixels),dispersion:0.0,label:''}
track = replicate(dum,ntracks)

;dispersion = 0.1238	; approx. dispersion. Will need to be calculated
			; as a function of wavelength
p =findgen(npixels)

; get wavelength calibration for MDS
calc_mds_wave,wave_calc,npixels,shot
	a=1.0d/1200.0d*1d7	; grating groove separation in Angstroms (1200 g/mm)
	grmm=1d7/a ; [/mm] grooves per mm 

; calculate average dispersion over spectrum to pass on to main routine
dispersion = (wave_calc(n_elements(wave_calc)-1)-wave_calc(0)) $
		/ (float(n_elements(wave_calc))-1.0)

; assuming 770 pixel arrangement on MDS_CCD4, but leave general
for i=0,ntracks-1 do begin ;fill data arrays
  track[i].d = d[*,i,*]
  track[i].w = wave_calc	; all wave arrays the same for the moment...
  track[i].label = label[i]
  track[i].dispersion = dispersion
;                               perform background subtraction if triggered
                                ;with 6A (t = -10 sec) and first readout
                                ;  delayed until -0.300 s. 
   
;Replace last two
;statements with nested do loops below. It subtracts the average of 
;frames 2-5 from all frames.
;
;  bkgd = fltarr(npixels,ntracks)
;  if ( t_delay ne 0.0 ) then begin
;    for i=0,ntracks-1 do begin
;      bkgd(*,i) = total( d(*,i,1:4),3 )/4.
;      track(i).label = label(i)
;      for fr = 0,ntimes-1 do $
;        track(i).d = d(*,i,fr) - bkgd(*,i) 
;    endfor
;  endif else begin track(i).d = d(*,i,*)
;  endelse
endfor

; construct track structure using assembled read-in data
; fibers are physically aligned on MDS input slit with 7 radial fibers,
; then the DIMES fiber, then the 6 toroidal fibers. Arrange so that DIMES
; fiber is located between data indices 2 & 3, where it is actually mapped 
; on the floor of the divertor.

pos=fltarr(ntracks)	; define major radial position array
if_read_file='/d2/divertor/mds/calib/ifun_5015.dat'
openr,iflun,if_read_file,/get_lun
	readf,iflun,if_shot
	date=''
	readf,iflun,date
	readf,iflun,if_track
	readf,iflun,if_frame
	readf,iflun,if_pixels
	readf,iflun,centroid

	ifun = fltarr(if_pixels)
	readf,iflun,ifun
close,iflun	

offset = float(if_pixels-1.0)/2.0 - centroid

; import helium profile data calculated by R. Isler considering fine structure
; and Zemman splitting (unpolarized) at 2.0 Tesla (like 8/2/96 XP)

restore,'/d2/divertor/mds/calib/HeII_profs.dat'
	; calculate the helium profile offset from the centroid pixel index
	he_centroid=fltarr(8)
	for i=0,7 do he_centroid(i)=total(findgen(15)*he_prof(*,i))/total(he_prof(*,i))
	he_offset=(float(n_elements(he_prof(*,0)))-1.0)/2.0-he_centroid

; same for D_alpha, B=2.1 T., variable: da_prof
restore,'/d2/divertor/mds/calib/D_alpha_profs.dat'
	; calculate the offset from the centroid pixel index
	da_centroid=fltarr(8)
	for i=0,7 do da_centroid(i)=total(findgen(25)*da_prof(*,i))/total(da_prof(*,i))
	da_offset=(float(n_elements(da_prof(*,0)))-1.0)/2.0-da_centroid

; ...and for NeI (6402)
restore,'/d2/divertor/mds/calib/NeI_profs.dat'
	ne_centroid=fltarr(8)
	for i=0,7 do ne_centroid(i)=total(findgen(21)*ne_prof(*,i))/total(ne_prof(*,i))
	ne_offset=(float(n_elements(ne_prof(*,0)))-1.0)/2.0-ne_centroid

a={	ierr:0,	$
	shot:shot,	$
	ntimes:ntimes,	$
	t:t,		$
	integ_t:exp_t,$
	slit:slit,	$
	comment:comment,$
	npixels:npixels,$
	track:track,$
	pos: pos,$
	maxdata:0.0,$
	wave_set:wave_set,$
	grmm:grmm,$
	ibacksub:0,$
	ixtalk:0,$
	iflatfield:0,$
	icalib:0,$
	ineutron:0}

a.maxdata=max(a.track.d)
print,'DivSpec data retrieved for shot: '+strtrim(shot,2)
return,a
end
;                   each 9 or 10 digit number word in base 10 must be
;                   converted to binary. The binary counterpart
;                   consists of four bytes in reverse order, where
;                   each byte represents an ASCII character.

;ascii array
; 0     128       number of words in array following first two elements
; 1     128       number of words in array following first two elements
; 2  1145656915   DIVS  name of experiment (2 words)                           
; 3  1346716448   PEC_
; 4   538980658   __12  shot number (2 words)  
; 5   926102582   7346
; 6   859067442   34T2      labels of sixteen ROIs (64 words)
; 7  1412650036   T3T4 
; 8  1412781110   T5T6
; 9  1429296434   U1U2
;10  1429427508   U3U4
;11  1429558582   U5U6
;12  1278299186   L1L2
;13  1278430260   L3L4
;14   538976288
;15   538976288
;16   538976288
;17   538976288
;18   538976288
;19   538976288
;20   538976288
;21   538976288
;22   538976288
;23   538976288
;24   538976288
;25   538976288
;26   538976288
;27   538976288
;28   538976288
;29   538976288
;30   825699658   18-J
;31  1634610480   an-0
;32   924856630   7 13
;33   976367674   :16:
;34   859390778   58C:
;35  1547912004   \CCD
;36   541286765    Cam
;37  1701994844   era\
;38  1131376230   Conf
;39  1768381507   ig\C
;40  1128538180   CD D
;41  1701208437   efau 
;42  1819553379   lt.c
;43  1718034464   fg



;int16 array

; 0     128       number of words in array following first two elements
; 1     128       number of words in array following first two elements
; 2       2       CCD mode
; 3      16       number of ROIs
; 4       5       number of frames
; 5   16384       total number of pixels in frame
; 6       0
; 7       0


;real32 array

; 0     512.      number of words in array following first two elements
; 1     512.      number of words in array following first two elements
; 2       2.      CCD mode
; 3       0.
; 4       5.      number of frames
; 5       0.
; 6       0.050   exposure time
; 7    9100.      wavelength
; 8     300.      slit width
; 9-13    0. 
;14       0.      begin X (row index) 
;15       0.      begin Y (column index)
;16    1023.      end   X
;17    1023.      end   Y
;18       1.      bin   X (sum pixels in row)
;19      64.      bin   Y (sum columns)
;20   16384.      total number of pixels in single readout
;21      16.      number of ROIs
;22       0.      number of timeouts
;23       0.050   time of first readout
;24       0.100   time of second readout
;23:nframes-1   
;         0.150   times of subsequent readouts which can be stored, but only
;        12.750   n_frames have data
;23+nframes      
;         0.      delay time relative to synchronous DIII-D trigger at
;                   -10 sec
       

;old ccd setup for 14 tracks
; [56,45,44,45,44,44,42,38,36,36,37,36,36,37] of 22.5 microns/row
