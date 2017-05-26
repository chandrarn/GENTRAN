; this file is used to configure w_spec and impspec for an experiment
; data retrieval routines should also be specified here

; machine specific configs should be implemented through this file, not 
; hard-coded into w_spec and impspec, so that the codes are more easily 
; portable to other experiments.

PRO load_nstx,machine 
line_std={nam:['Mo XXXII','Ar XV','Ca XVII','Fe XXIII'], $ ; default lines to track
	 lam0:[126.1,     220.4,  192.3,    132.4     ], $
	 lam1:[129.2,     222.7,  193.6,    133.5     ], $
	 spec:['loweus',  'loweus','loweus','loweus'  ], $
	 mult:[1.0,       1.0,    1.0,     1.0        ], $
	 plot:[1,         1,      1,       1          ]}

line_Mo ={nam:['Mo XXXII','Mo I','C III','Li I'], $ ; shot 139375, xp1041
	 lam0:[126.1,     3785.,  4062., 3908. ], $
	 lam1:[129.2,     3815.,  4078., 3921. ], $
	 spec:['loweus',  'dims','vips2','vips2'], $
	 mult:[1.0,       10.0,   1.0,   1.0        ], $
	 plot:[1,         1,      1,     1          ]}

machine={name:'nstx', $
	 inst:['short', 'loweus','xeus','dims','vips2'], $ ; spectroscopic instruments
	 load:[1,       1,       1,     1,     1      ], $
	 plot:[1,       1,       1,     1,     1      ], $
	 nch: [1,       1,       1,     19,    10     ], $
	 timetr:['totpwr','dalpha','dens','beams','hhfw'], $ ; time traces
	 line:line_Mo,$
	 elem:'Li,C,Cu,Mo,W' $ ; default elements to overplot
	}

END

function load_nstx_shot
	; retrieve the current shot
	return,mdscur_shot('nstx')
end

PRO load_nstx_spred,shot,d,short=short,verb=verb
        mdsopen,'passivespec',shot
         t=mdsvalue('\PASSIVESPEC::TOP.SPRED.TIMES',/quiet,status=tstatus)

         d=mdsvalue('_sig=\PASSIVESPEC::TOP.SPRED.RAWDATA.SPECTRA',/quiet,status=dstatus) ; [nwl nch nt]
          if keyword_set(short) $
	   then d=reform(d[*,0,*]) $
	   else d=reform(d[*,1,*])
	  if size(d,/n_dim) lt 2 then begin
		nl=size(d,/n_el)
		d=reform(d,nl,1)
		nt=1
		dim=size(d,/dim)
	  endif else begin
		dim=size(d,/dim)
	 	nl=dim[0]
	 	nt=dim[1]
	  endelse

	; background subtraction
	tmp=where(t lt 0)
	if tmp[0] ne -1 then tmp=tmp[where(tmp gt 0)] ; first frame tends to be garbage?
	if tmp[0] ne -1 then begin
	 bg=total(d[*,tmp],2)/n_elements(tmp)
	 for i=0,nt-1 do d[*,i]=reform(d[*,i])-bg
	endif

         ;l=mdsvalue('\PASSIVESPEC::TOP.SPRED.WAVELENGTH',/quiet,status=lstatus) ; doesn't work...
         C=mdsvalue('\PASSIVESPEC::TOP.SPRED.POLY_COEFFS',/quiet,status=lstatus) ; doesn't work...
          if keyword_set(short) $
	   then l=poly(indgen(nl),C[0:3]) $
	   else l=poly(indgen(nl),C[4:7])

         vph=mdsvalue('\PASSIVESPEC::TOP.SPRED.PH_VOLTAGE',/quiet)
         vmcp=mdsvalue('\PASSIVESPEC::TOP.SPRED.MCP_VOLTAGE',/quiet)
	mdsclose

	c={vph:vph,vmcp:vmcp}
	cr=string(10B)
	msg=cr+'  VPH='+num2str(vph,1) $
           +cr+'  VMCP='+num2str(vmcp,1)

	if tstatus gt 0 then begin
	 if n_elements(t) gt nt then t=t[0:nt-1]
	endif else t=indgen(nt-1)

	status= dstatus and tstatus and lstatus

	if keyword_set(short) $
	 then title='SHORT' $
	 else title='LONG'

	d={t:t, $
	   l:l,llabel:n2g('lambda')+' ['+string(197B)+']', $
	   d:d,dlabel:'[AU]', $
	   title:title,c:c,s:status,p:status,index:0,msg:msg}
if keyword_set(verb) then help,/str,d
END

PRO load_nstx_short,shot,d,verb=verb
	load_nstx_spred,shot,d,/short,verb=verb
END

PRO load_nstx_long,shot,d,verb=verb
	load_nstx_spred,shot,d,verb=verb
END

PRO load_nstx_loweus,shot,d,verb=verb
        mdsopen,'pspec_pc',shot
         d=mdsvalue('\PSPEC_PC::TOP.LOWEUS.SPECTRA',/quiet,status=dstatus) ; [nwl nch nt]
         l=mdsvalue('\PASSIVESPEC::TOP.LOWEUS.WAVELENGTHS',/quiet,status=lstatus)
         t=0.0
	mdsclose

	status = dstatus mod 2 && lstatus mod 2

	d={t:t, $
	   l:l,llabel:n2g('lambda')+' ['+string(197B)+']', $
	   d:d,dlabel:'[AU]', $
	   dstatus:dstatus,lstatus:lstatus, $
	   title:'LOWEUS',s:status,p:status,index:0,msg:''}
if keyword_set(verb) then help,/str,d
END

PRO load_nstx_xeus,shot,d,verb=verb
        mdsopen,'pspec_pc',shot
         d=mdsvalue('\PSPEC_PC::TOP.XEUS.SPECTRA',/quiet,status=dstatus) ; [nwl nch nt]
         l=mdsvalue('\PASSIVESPEC::TOP.XEUS.WAVELENGTHS',/quiet,status=lstatus)
         t=0.0
	mdsclose

	status = dstatus mod 2 && lstatus mod 2

	d={t:t, $
	   l:l,llabel:n2g('lambda')+' ['+string(197B)+']', $
	   d:d,dlabel:'[AU]', $
	   dstatus:dstatus,lstatus:lstatus, $
	   title:'XEUS',s:status,p:status,index:0,msg:''}
if keyword_set(verb) then help,/str,d
END

PRO load_nstx_dims,shot,d,verb=verb
        mdsopen,'passivespec',shot
         dtmp=mdsvalue('\PASSIVESPEC::TOP.DIMS.RAWDATA.SPECTRA',/quiet,status=dstatus) ; [nwl nch nt]
	 if dstatus mod 2 ne 1 then dtmp=reform(fltarr(1,1,1),1,1,1)
	 dim=size(dtmp,/dim)
	  nl=dim[0]
	  nch=dim[1]
	  nt=dim[2]

         t=mdsvalue('\PASSIVESPEC::TOP.DIMS.DETECTOR.TIMES',/quiet,status=tstatus)

	 ch_nam=strarr(nch)
	 for i=0,nch-1 do $
	  ch_nam[i]=mdsvalue('\PASSIVESPEC::TOP.DIMS.DETECTOR.FIBER_MAP.BIN_'+string(i+1,format='(i02)'),/quiet)

         grmodel=mdsvalue('\PASSIVESPEC::TOP.DIMS.SPECTROGRAPH.GRATING.MODEL',/quiet)
         grmm=mdsvalue('\PASSIVESPEC::TOP.DIMS.SPECTROGRAPH.GRATING.GRMM',/quiet,status=lstatus1)
         cwl=mdsvalue('\PASSIVESPEC::TOP.DIMS.SPECTROGRAPH.CWL',/quiet,status=lstatus2)
         cwlnote=mdsvalue('\PASSIVESPEC::TOP.DIMS.SPECTROGRAPH.CWLNOTE',/quiet)
         entr_slit=mdsvalue('\PASSIVESPEC::TOP.DIMS.SPECTROGRAPH.ENTR_SLIT',/quiet)
         exp_time=mdsvalue('\PASSIVESPEC::TOP.DIMS.DETECTOR.EXP_TIME',/quiet)
	 filt_model=mdsvalue('\PASSIVESPEC::TOP.DIMS.SPECTROGRAPH.FILTER.MODEL',/quiet)
        mdsclose

	c={grmodel:grmodel,grmm:grmm,cwl:cwl,cwlnote:cwlnote,entr_slit:entr_slit,exp_time:exp_time,filt_model:filt_model}
	cr=string(10B)
	msg=cr+'  GRMM='+num2str(grmm,1) $
           +cr+'  CWL='+num2str(cwl,1) $
           +cr+'  SLIT='+num2str(entr_slit,1)

	dch=reform(fltarr(nl,nt,nch),nl,nt,nch)
	for i=0,nch-1 do dch[*,*,i]=dtmp[*,i,*]

	; background subtraction
	if tstatus then begin
		tmp=where(t lt 0)
		if tmp[0] ne -1 then tmp=tmp[where(tmp gt 0)]
		if tmp[0] ne -1 then begin
		 bg=total(reform(dch[*,tmp,*]),2)/n_elements(tmp)
		 for i=0,nt-1 do dch[*,i,*]=reform(dch[*,i,*])-bg
		endif
	endif

	d=total(dch,3)/nch

	if lstatus1 and lstatus2 then begin
	 gratings=[[3600, 0.012],$
	           [2400, 0.018],$
	       	   [1800, 0.02],$
	       	   [1200, 0.03],$
	     	   [600,  0.06],$
	           [300,  0.12],$
	           [150,  0.24],$
	           [75,   0.48]]
	 wlpp=gratings[where(grmm eq gratings[*,1]),1] ; [nm] wavelength per pixel
	 wl1 = cwl - wlpp*nl/2
	 l = 10.*(wl1 + findgen(nl)*wlpp) ; [A]
	endif else l=findgen(nl)
	lch=fltarr(nl,nch)
	for i=0,nch-1 do lch[*,i]=l

	status= dstatus mod 2 && tstatus mod 2 && lstatus1 mod 2 && lstatus2 mod 2

	d={t:t, $
	   l:l,llabel:n2g('lambda')+' ['+string(197B)+']', $
	   d:d,dlabel:'[AU]', $
	   dch:dch,lch:lch,ch_nam:ch_nam, $
	   c:c, dstatus:dstatus,tstatus:tstatus,lstatus1:lstatus1,lstatus2:lstatus2,$
	   title:'dims',s:status,p:status,index:0,msg:msg }
if keyword_set(verb) then help,/str,d
END

PRO load_nstx_vips2,shot,d,verb=verb
        mdsopen,'passivespec',shot
         dtmp=mdsvalue('\PASSIVESPEC::TOP.VIPS2.RAWDATA.SPECTRA',/quiet,status=dstatus) ; [nwl nch nt]
	 dim=size(dtmp,/dim)
	  nl=dim[0]
	  nch=dim[1]
	  nt=dim[2]

         t=mdsvalue('\PASSIVESPEC::TOP.VIPS2.TIMES',/quiet,status=tstatus)

	 ch_nam=strarr(nch)
	 for i=0,nch-1 do $
	  ch_nam[i]=strtrim( mdsvalue('\PASSIVESPEC::TOP.VIPS2.BIN_'+strtrim(string(i+1,format='(i)'),2),/quiet), 2)

         grating=mdsvalue('\PASSIVESPEC::TOP.VIPS2.GRATING',/quiet,status=lstatus1)
         cwl=mdsvalue('\PASSIVESPEC::TOP.VIPS2.CENTERWAVLEN',/quiet,status=lstatus2)
         slit=mdsvalue('\PASSIVESPEC::TOP.VIPS2.ENT_WID_SLIT',/quiet)
         exp_time=mdsvalue('\PASSIVESPEC::TOP.VIPS2.EXPOSURE_TIM',/quiet)
        mdsclose

	c={grating:grating,cwl:cwl,slit:slit,exp_time:exp_time}
	cr=string(10B)
	msg=cr+'  GRATING='+num2str(grating,1) $
           +cr+'  CWL='+num2str(cwl,1) $
           +cr+'  SLIT='+num2str(slit,1)

	dch=reform(fltarr(nl,nt,nch),nl,nt,nch)
	for i=0,nch-1 do dch[*,*,i]=dtmp[*,i,*]

	; background subtraction
	if tstatus then begin
		tmp=where(t lt 0)
		if tmp[0] ne -1 then tmp=tmp[where(tmp gt 0)]
		if tmp[0] ne -1 then begin
		 bg=total(reform(dch[*,tmp,*]),2)/n_elements(tmp)
		 for i=0,nt-1 do dch[*,i,*]=reform(dch[*,i,*])-bg
		endif
	endif else bg=reform(fltarr(nl,nch),nl,nch)

	d=total(dch,3)/nch

	if lstatus1 and lstatus2 then begin
	 gratings=[0.012,0.030,0.064]
	 ;wlpp=gratings[grating] ; [nm] wavelength per pixel
	 wlpp=gratings[grating-1] ; [nm] wavelength per pixel
	 l = 10.*(cwl + (findgen(nl)-nl/2.)*wlpp)
	endif else l=findgen(nl)
	lch=fltarr(nl,nch)
	for i=0,nch-1 do lch[*,i]=l

	status= dstatus mod 2 && tstatus mod 2 && lstatus1 mod 2 && lstatus2 mod 2

	d={t:t, $
	   l:l,llabel:n2g('lambda')+' ['+string(197B)+']', $
	   d:d,dlabel:'[AU]', $
	   dch:dch,lch:lch,ch_nam:ch_nam, $
	   c:c, dstatus:dstatus,tstatus:tstatus,lstatus1:lstatus1,lstatus2:lstatus2,$
	   title:'vips2',s:status,p:status,index:0,msg:msg }
if keyword_set(verb) then help,/str,d
END

PRO load_nstx_ornloo,shot,d

; load from: \NSTX::TOP.PASSIVESPEC.MINI_SPECS.ORNL_OO_1

END


PRO load_nstx_totpwr,shot,d
	pt='TOTPWR'
        mdsopen,'passivespec',shot
         d=mdsvalue('_sig=\PASSIVESPEC::TOP.BOLOM_ARRAY.ANALYSIS.TOTPWR',/quiet,status=status)
         t=mdsvalue('dim_of(_sig,0)',/quiet)
        mdsclose

	status = (status mod 2) eq 1

	d={t:t, $
	   d:d,dlabel:'[?]', $
	   pt:pt,title:pt,oplot:0,s:status,yr:[0.0,0.0]}
END

PRO load_nstx_divbolom,shot,d
	pt='DIV_BOLOM'
        mdsopen,'passivespec',shot
         d=mdsvalue('_sig=\PASSIVESPEC::TOP.DIV_BOLOM.ANALYSIS:IRRADIANCE',/quiet,status=status)
         t=mdsvalue('dim_of(_sig,0)',/quiet)
        mdsclose

	status = (status mod 2) eq 1

	d={t:t, $
	   d:d,dlabel:'[?]', $
	   pt:pt,title:pt,oplot:0,s:status,yr:[0.0,0.0]}
END

PRO load_nstx_dalpha,shot,d
	pt='DALPHA'
        mdsopen,'passivespec',shot
         d=mdsvalue('_sig=\PASSIVESPEC::TOP.FILTERED_VIS.BAYC_OPIPE:DALPHA_EIES',/quiet,status=status)
         t=mdsvalue('dim_of(_sig,0)',/quiet)
        mdsclose

	status = (status mod 2) eq 1

	d={t:t, $
	   d:d,dlabel:'[?]', $
	   pt:pt,title:pt,oplot:0,s:status,yr:[0.0,0.0]}
END

PRO load_nstx_coreusxr,shot,d
	pt='GTOP.CHORD_07'
        mdsopen,'usxr',shot
         d=mdsvalue('_sig=\USXR::TOP.'+pt,/quiet,status=status)
         t=mdsvalue('dim_of(_sig,0)',/quiet)
        mdsclose

	status = (status mod 2) eq 1

        IF status THEN BEGIN
                bl=mean(d[0:ipt(t,-0.01)])
                d-=bl
        ENDIF

	d={t:t, $
	   d:d,dlabel:'[?]', $
	   pt:pt,title:pt,oplot:0,s:status,yr:[0.0,0.0]}
END

PRO load_nstx_dens,shot,d
	pt='dens'
        ;mdsopen,'rf',shot
        ; d=mdsvalue('_sig=\RF::TOP.ORNLREFL:REFL_AVGDENS',/quiet,status=status)
        mdsopen,'activespec',shot
         d=mdsvalue('_sig=\ACTIVESPEC::TOP.MPTS.OUTPUT_DATA.BEST:FIT_NE',/quiet,status=status)
         t=mdsvalue('dim_of(_sig,0)',/quiet)
        mdsclose

	status = (status mod 2) eq 1

	if status then begin
		dims=size(d,/dimensions)
		 nt=dims[0]
		 nr=dims[1]
		d=total(d,2)/nr
	endif

	d={t:t, $
	   d:d,dlabel:'[1e20/m!E3!N]', $
	   pt:pt,title:pt,oplot:0,s:status,yr:[0.0,0.0]}
END

PRO load_nstx_beams,shot,d
	pt='BEAMS'
        mdsopen, 'nbi', shot
         d=mdsvalue('_sig=\NBI::TOP.ANALYSIS:P_INJ',/quiet,status=status) ;[MW]
         t=mdsvalue('dim_of(_sig,0)',/quiet)
        mdsclose

	status = (status mod 2) eq 1

	d={t:t, $
	   d:d,dlabel:'[MW]', $
	   pt:pt,title:pt,oplot:0,s:status,yr:[0.0,0.0]}
END

PRO load_nstx_hhfw,shot,d
	pt='HHFW'
        mdsopen, 'rf', shot
         d=mdsvalue('_sig=\RF::TOP.HHFW.POWERSIGS:HHFW_POWER',/quiet,status=status)/1.0e6 ; [MW]
         t=mdsvalue('dim_of(_sig,0)',/quiet)
        mdsclose

	status = (status mod 2) eq 1

	d={t:t, $
	   d:d,dlabel:'[?]', $
	   pt:pt,title:pt,oplot:1,s:status,yr:[0.0,0.0]}
END

PRO load_nstx_ech,shot,d
	pt='ECH'
        mdsopen, 'rf',shot
         d=mdsvalue('_sig=\RF::TOP.ECH:ECHPOWER',/quiet,status=status)
         t=mdsvalue('dim_of(_sig,0)',/quiet)
        mdsclose

	status = (status mod 2) eq 1

	d={t:t, $
	   d:d,dlabel:'[?]', $
	   pt:pt,title:pt,oplot:1,s:status,yr:[0.0,0.0]}
END

