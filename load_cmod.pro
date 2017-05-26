; this file is used to configure w_spec and impspec for an experiment
; data retrieval routines should also be specified here

; machine specific configs should be implemented through this file, not 
; hard-coded into w_spec and impspec, so that the codes are more easily 
; portable to other experiments.

; a note regarding MDSplus statuses: 
;  if the lowest significant bit is 1 (ie the status is 'odd' ), success  is indicated
;  if the lowest significant bit is 0 (ie the status is 'even'), an error is indicated
; in logical tests, IDL behaves the same way (ie 21 is true, 20 is false, etc)

PRO load_cmod,machine
line_std={nam:['D I (9)','F II',   'Mo I', 'D I (8)' , 'D I (7)'  ,'N II', 'D I (6)', 'B II'],$ ; default lines to track
	 lam0:[3825.1,     3841.4,   3856.4,  3871.3,   3946.9,    3989.2,  4084.9, 4116.3], $
	 lam1:[3842.5,     3856.9,   3870.8,  3905.6,   3987.7,    3997.2,  4112.9,  4121.3], $
	 spec:['chromex',  'chromex','chromex', 'chromex','chromex', 'chromex','chromex','chromex'  ], $
	 mult:[1.0,       1.0,    1.0,    1.0,      1.0,      1.0,	1.0,	1.0], $
	 plot:[1,         1,        0,      1,        1,        0,	1,	1         ],$
	 ct:[ 12,	12,	12,	12,	12,	12,	39,	39],$
	 col:[30,	90,	50,	100,	120,	150,	250,	210]}
machine={name:'cmod', $
	 inst:['hada','mcp','xeus','loweus','helike','hlike','chromex'], $ ; spectroscopic instruments
	 load:[0,     0,    0,     0,       0,       0,      1        ], $
	 plot:[1,     0,    1,     1,       0,       0,      1        ], $
	 nch: [1,     1,    1,     1,       16,      16,     16       ], $
	 timetr:['2pid','corextomo','nl04','rfpow','lhpow'], $ ; time traces
	 line:line_std,$
	 elem:'D,B,N,F,Mo' $ ; default elements to overplot
	}

END

function load_cmod_shot
	; retrieve the current shot
	shot=mdscur_shot('cmod')
	IF strpos(string(shot),'000') NE -1 THEN shot=-1 ;if shot is a 000 default to a functional shot
	return,shot
end

PRO load_cmod_mcp,shot,d
        mdsopen,'spectroscopy',shot
         d=mdsvalue('_sig=\SPECTROSCOPY::TOP.VUV.ANALYSIS:MCP_BRI_MOD',/quiet,status=status)
         l=mdsvalue('dim_of(_sig,1)',/quiet)
         t=mdsvalue('dim_of(_sig,0)',/quiet)
        mdsclose,'spectroscopy',shot
	pos=[0,0,0,0]
	status = status mod 2

	d={t:t, $
	   l:l,llabel:n2g('lambda')+' ['+string(197B)+']', $
	   d:d/10.0,dlabel:'[AU]', pos:pos,$
	   title:'mcp',s:status,p:status,index:0,msg:''}
END

PRO load_cmod_hada,shot,d
; from /usr/local/cmod/codes/spectroscopy/hada/hada_128_read_quiet.pro
        mdsopen,'spectroscopy',shot
	 ; spectrum
	 hasiga=mdsvalue('data(\SPECTROSCOPY::TOP.HALPH_DALPH.HARDWARE:HA_DA_A14:INPUT_2)',status=stat_a)
	 hasig_no_a=mdsvalue('data(\SPECTROSCOPY::TOP.HALPH_DALPH.HARDWARE:HA_DA_A14:INPUT_1)',status=stat_noa)
	 asize=n_elements(hasiga)
	  nl=272
	  nt=asize/nl

	 ; timing
	 digrate=mdsvalue('\SPECTROSCOPY::TOP:EARLS_MPB.CHANNEL_0:P1',status=stat_dig)
	 toffset=mdsvalue('\SPECTROSCOPY::TOP:EARLS_MPB.CHANNEL_1:EDGES',status=stat_toff) & toffset=toffset[0]

	 ; wavelength
	 slope=mdsvalue('\SPECTROSCOPY::TOP.HALPH_DALPH.ANALYSIS:SLOPE',status=stat_slop)
	 offset=mdsvalue('\SPECTROSCOPY::TOP.HALPH_DALPH.ANALYSIS:OFFSET',status=stat_loff)
        mdsclose,'spectroscopy',shot

	; timing
	stat_t=stat_dig mod 2 && stat_toff mod 2 && nt gt 0
	if stat_t then begin
		digperiod=1./digrate
		hatime=(indgen(nt)+1)*digperiod+toffset
		hatimetree=hatime-.5/digrate
		t=hatimetree
	endif else t=-1

	; wavelength
	stat_l=stat_slop mod 2 && stat_loff mod 2 
	if stat_l then begin
		index=reverse(indgen(nl/2)*2)
		indexx=indgen(nl/2)
		lamda=slope*indexx+offset  ;use the calibration & keep track of shifts during shot
		l=lamda
	endif else l=-1

	; spectrum
	stat_d=stat_a mod 2 && stat_noa mod 2 
	if stat_d && stat_t then begin
		hasiga=hasiga[0:nt*nl-1]
		hasig_no_a=hasig_no_a[0:nt*nl-1]
		hasigra=reform(hasiga,nl,nt)
		hasigr_no_a=reform(hasig_no_a,nl,nt)
		aara=hasigra[index,*]
		aar_no_a=hasigr_no_a[index,*]
		aari=aara
		satval=4.9
		tmp=where(aara le satval)
		if tmp[0] ne -1 then aar_no_a[tmp]=0.
		tmp=where(aara gt satval)
		if tmp[0] ne -1 then aara[tmp]=0.
		aari=aara+10.*(aar_no_a)
		aario=aari
		woffset=[indgen(5)+14,indgen(5)+120]
		spectra=fltarr(nl/2,nt)
		for j=0,nt-1 do begin	
		 aveoffset=mean(aario[woffset,j])
		 aaofff=aario[*,j]-aveoffset
		 spectra[*,j]=aaofff
		end
		d=spectra

		tmp=where(t lt 0)
		if tmp[0] ne -1 then begin
		 bg=total(d[*,tmp],2)/n_elements(tmp)
		 for i=0,nt-1 do d[*,i]=d[*,i]-bg
		endif

		tmp=where(t le 2.0)
		if tmp[0] ne -1 then begin
		 t=t[tmp]
		 d=d[*,tmp]
		endif
	endif else d=-1
	pos=[0,0,0,0]

	lamd=656.1032
	lamh=656.2793

	status=(stat_a mod 2) && (stat_noa mod 2) && (stat_dig mod 2) && (stat_toff mod 2) && (stat_slop mod 2) && (stat_loff mod 2)

	d={t:t, $
	   l:l*10,llabel:n2g('lambda')+' ['+string(197B)+']', $
	   d:d,dlabel:'[AU]', pos:pos,$
	   title:'hada',s:status,p:status,index:0,msg:''}
END

PRO load_cmod_xeus,shot,d
        mdsopen,'spectroscopy',shot
         d=mdsvalue('_sig=\SPECTROSCOPY::TOP.XEUS:SPEC',/quiet,status=status)
         l=mdsvalue('dim_of(_sig,0)',/quiet)
         t=mdsvalue('dim_of(_sig,1)',/quiet)
	 sig=mdsvalue('dim_of(_sig,2)',/quiet)
        mdsclose,'spectroscopy',shot
	pos=[0,0,0,0]

	; hack to fix NaN's... need to talk with Matt about why these are in here...
	tmp=where( finite(sig) eq 0 )
	if tmp[0] ne -1 then sig[tmp]=0

	status = status mod 2

	d={t:t, $
	   l:l,llabel:n2g('lambda')+' ['+string(197B)+']', $
	   d:d,dlabel:'[AU]', sig:sig, pos:pos,$
	   title:'xeus',s:status,p:status,index:0,msg:''}
END

PRO load_cmod_loweus,shot,d
        mdsopen,'spectroscopy',shot
         d=mdsvalue('_sig=\SPECTROSCOPY::TOP.LOWEUS:SPEC',/quiet,status=status) ; dim[nl nt]
         l=mdsvalue('dim_of(_sig,0)',/quiet)
         t=mdsvalue('dim_of(_sig,1)',/quiet)
	 sig=mdsvalue('dim_of(_sig,2)',/quiet) ; dim[nl nt], sigma errorbars
        mdsclose,'spectroscopy',shot
	pos=[0,0,0,0]

	; hack to fix NaN's... need to talk with Matt about why these are in here...
	tmp=where( finite(sig) eq 0 )
	if tmp[0] ne -1 then sig[tmp]=0

	status = status mod 2

	d={t:t, $
	   l:l,llabel:n2g('lambda')+' ['+string(197B)+']', $
	   d:d,dlabel:'[AU]', sig:sig, pos:pos,$
	   title:'loweus',s:status,p:status,index:0,msg:''}
END

PRO load_cmod_helike,shot,d
        mdsopen,'spectroscopy',shot
         dtmp=mdsvalue('_sig=\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.SPEC:SPECBR',/quiet,status=status)
         ttmp=mdsvalue('dim_of(_sig,1)',/quiet)
         ltmp=mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.SPEC:LAM',/quiet)
        mdsclose,'spectroscopy',shot

	status = status mod 2

        IF status THEN BEGIN
		; hirex channels can change dynamically in time, but we'll assume that they don't
		; and take dimensions from the initial time slice
		ich=where(dtmp[*,0,0] ne -1) & dtmp=dtmp[ich,*,*] & ltmp=ltmp[ich,*,*]
		it= where(ttmp        ne -1) & dtmp=dtmp[*,it,*]  & ltmp=ltmp[*,it,*]  & ttmp=ttmp[it]
		il= where(dtmp[0,0,*] ne -1) & dtmp=dtmp[*,*,il]  & ltmp=reform(ltmp[*,0,il])

		dim=size(dtmp,/dim)
		 nch=dim[0]
		 nt=dim[1]
		 nl=dim[2]

		dch=fltarr(nl,nt,nch)
		ch_nam=strarr(nch)
		for i=0,nch-1 do $
		 dch(*,*,i)=rotate(reform(dtmp(i,*,*)),4)/1.0e3
		 lch=rotate(ltmp,4)
		 t=ttmp

                nch0=fix(nch/2) ; assume CHMAP is not time evolving and take the middle (core) ch
                l=reform(lch[*,nch0])
                d=reform(dch[*,*,nch0])
		pos=fltarr(4,nch0)
	ENDIF ELSE BEGIN
		nan=!Values.F_NAN
		d=nan
		l=nan
		t=nan
		dch=nan
		lch=nan
		ch_nam=nan
		pos=[0,0,0,0]
        ENDELSE

	d={t:t, $
	   l:l,llabel:n2g('lambda')+' ['+string(197B)+']', $
	   d:d,dlabel:'[AU]', pos:pos,$ ; dims: [l,t]
	   dch:dch,lch:lch,ch_nam:ch_nam, $
	   title:'helike',s:status,p:status,index:0,msg:''}
END

PRO load_cmod_hlike,shot,d
        mdsopen,'spectroscopy',shot
         dtmp=mdsvalue('_sig=\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HLIKE.SPEC:SPECBR',/quiet,status=status) ; dims: [ch,t,l]
         ttmp=mdsvalue('dim_of(_sig,1)',/quiet)
         ltmp=mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HLIKE.SPEC:LAM',/quiet)
        mdsclose,'spectroscopy',shot
	
	status = status mod 2

        IF status THEN BEGIN
		; hirex channels can change dynamically in time, but we'll assume that they don't
		; and take dimensions from the initial time slice
		ich=where(dtmp[*,0,0] ne -1) & dtmp=dtmp[ich,*,*] & ltmp=ltmp[ich,*,*]
		it= where(ttmp        ne -1) & dtmp=dtmp[*,it,*]  & ltmp=ltmp[*,it,*]  & ttmp=ttmp[it]
		il= where(dtmp[0,0,*] ne -1) & dtmp=dtmp[*,*,il]  & ltmp=reform(ltmp[*,0,il])

		dim=size(dtmp,/dim)
		 nch=dim[0]
		 nt=dim[1]
		 nl=dim[2]

		dch=fltarr(nl,nt,nch)
		ch_nam=strarr(nch)
		for i=0,nch-1 do $
		 dch(*,*,i)=rotate(reform(dtmp(i,*,*)),4)/1.0e3
		 lch=rotate(ltmp,4)
		 t=ttmp

                nch0=fix(nch/2) ; assume CHMAP is not time evolving and take the middle (core) ch
                l=reform(lch[*,nch0])
                d=reform(dch[*,*,nch0])
		pos=fltarr(4,nch0)
	ENDIF ELSE BEGIN
		nan=!Values.F_NAN
		d=nan
		l=nan
		t=nan
		dch=nan
		lch=nan
		ch_nam=nan
		pos=[0,0,0,0]
        ENDELSE


	d={t:t, $
	   l:l,llabel:n2g('lambda')+' ['+string(197B)+']', $
	   d:d,dlabel:'[AU]', pos:pos,$ ; dims: [l,t]
	   dch:dch,lch:lch,ch_nam:ch_nam, $
	   title:'hlike',s:status,p:status,index:0,msg:''}
END


FUNCTION chromex_kbot_pos,ch=ch,vertex=vertex,rpts=rpts,zpts=zpts
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

FUNCTION chromex_abot_pos,ch=ch,vertex=vertex,rpts=rpts,zpts=zpts
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
FUNCTION chromex_ktop_pos,ch=ch,vertex=vertex,rpts=rpts,zpts=zpts
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

PRO load_cmod_chromex,shot,d
        mdsopen,'spectroscopy',shot
 	 dtmp=mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.ANALYSIS:RAW_SPECTRA',/quiet,status=dstatus)

	 t0=mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.CAMAC:TIME_START',/quiet,status=tstatus1)
	 dt=mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.CAMAC:DELTA_TIME',/quiet,status=tstatus2) ; actual measured time
	 ;dt=mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.SETTINGS:NOM_EXP_TIME',/quiet,status=tstatus2) ; actual time ends up being a bit longer...

         lch=mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.CALIBRATION:LAMBDA_TRUE', /quiet, status=lstatus)*10.0 ; [angstroms]
         px_offset=mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.CALIBRATION:PX_OFFSET', /quiet, status=pxstatus)

	 periscop=	mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.SETTINGS:CFG_PERISCOP',/quiet)
	 per_fibr=	mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.SETTINGS:CFG_PER_FIBR',/quiet)
	 bin_no=	mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.SETTINGS:BIN_NO', /quiet)
	 lambda_ctr=	mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.SETTINGS:LAMBDA_CTR',/quiet) ; [nm]
	 grating=	mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.SETTINGS:GRATING',/quiet) ; [grooves/mm]
	 slit=		mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.SETTINGS:SLIT_WIDTH',/quiet) ; [um]
        mdsclose,'spectroscopy',shot
	
	dstatus 	= dstatus  mod 2
	tstatus1	= tstatus1 mod 2
	tstatus2	= tstatus2 mod 2
	lstatus 	= lstatus  mod 2
	pxstatus	= pxstatus mod 2

	ch_nam=strarr(n(periscop)+1)
	for i=0,n(periscop) do ch_nam[i]=periscop[i]+'_'+num2str(per_fibr[i],1)

	c={periscop:periscop,per_fibr:per_fibr,bin_no:bin_no,lambda_ctr:lambda_ctr,grating:grating,slit:slit} ; config structure
	cr=string(10B)
	msg=cr+' .lambda_ctr='+num2str(lambda_ctr,dp=1) $
	   +cr+' .grating='+num2str(grating,1) $
	   +cr+' .slit='+num2str(slit,dp=1)

	dim=size(dtmp,/dim)
	 nt=dim[0]
	 nl=dim[1]
	 nch=dim[2]

	dch=fltarr(nl,nt,nch)
	rzview=fltarr(4,nch)
	kbot=chromex_kbot_pos(ch=kbch,vertex=kbv,rpts=kbr,zpts=kbz)
	ktop=chromex_ktop_pos(ch=ktch,vertex=ktv,rpts=ktr,zpts=ktz)
	abot=chromex_abot_pos(ch=abch,vertex=abv,rpts=abr,zpts=abz)
	FOR i=0,nch-1 DO BEGIN
		CASE c.periscop[i] OF
			'A_BOTTOM' : BEGIN
				index=where(abch EQ c.per_fibr[i])
				IF index[0] NE -1 THEN rzview[*,i]=[abv,abr[index],abz[index]]
                        END
			'K_BOTTOM' : BEGIN
				index=where(kbch EQ c.per_fibr[i])
				IF index[0] NE -1 THEN rzview[*,i]=[kbv,kbr[index],kbz[index]]
                        END
			'K_TOP' : BEGIN
				index=where(ktch EQ c.per_fibr[i])
				IF index[0] NE -1 THEN rzview[*,i]=[ktv,ktr[index],ktz[index]]
                        END
			ELSE : rzview[*,i]=-1.0
		ENDCASE
	ENDFOR
	for i=0,nch-1 do dch[*,*,i]=rotate(dtmp(*,*,i),4)

	if tstatus1 && tstatus2 $
	 then t = t0 + dt*findgen(nt) $
	 else begin
	  print,'chromex: timing data not available'
	  t = -1
	 endelse

	; calculate noise level
	;  noise figures for Princeton Instruments 1024B excelon3 camera:
	;  bit_depth=16 ;bits
	;  well_depth=80000 ;e-
	;  readout_noise=10 ;e-rms/pixel (Nr=8 counts rms)
	;  dark_current=0.002 ;e-/pixel/s (D<1 count/s)
	;  lQe=[300 350 400 500 600 700 800 900 1000] ;nm
	;   Qe=[4   60  76  88  94  92  80  48  10  ] ;%
	; noise = sqrt( P*Qe*t + D*t + Nr^2 ) = sqrt( digitizer counts BEFORE background subtraction )
	sigch = sqrt(dch)

	; normalize signals to bit depth
	dch=dch/2.^16
	sigch=sigch/2.^16

	; select frames before 0 for background and noise levels
	tmp=where( t lt 0-(t[1]-t[0]) ) 
	if tmp[0] ne -1  AND N_ELEMENTS(tmp) GT 2 then begin
                tmp=where(tmp NE 0 AND tmp NE 1) ; Ignore the first 2 frames, they are no good due to charge collecting.
		bg=total(dch[*,tmp,*],2)/n_elements(tmp)
		for i=0,nt-1 do dch[*,i,*]=dch[*,i,*]-bg
	endif

	if lstatus ne 1 or sum(finite(lch) eq 0) ne 0 then begin
	  print,'chromex: exact wavelength calibration not available, using estimate'
	  case c.grating of
		600  : disp = 360.0/c.grating ; [angstrom/pixel]
		1200 : disp = 390.0/c.grating ; [angstrom/pixel]
		1800 : disp = 375.0/c.grating ; [angstrom/pixel]
		else : disp = 360.0/c.grating
	  endcase
	  print, '  dispersion for '+num2str(c.grating,1)+'/mm grating is approximately '+num2str(disp,dp=4)+'nm/pixel'
	  lch = fltarr(nl,nch)
	  for i=0,nch-1 do lch[*,i] = (findgen(nl) - floor(nl/2.))*disp + c.lambda_ctr*10 ;[angstrom]
	endif

	l=total(lch,2)/nch
	d=total(dch,3)/nch
	sig=total(sigch,3)/nch

	if pxstatus then begin
	  ; the pixel offset is invariant with wavelength or grating, so calculate wavelength offset
	  l_offset = ( l[fix(nl/2)] - l[fix(nl/2)-1] )  * px_offset ;[angstrom]
	  print, ' adding wavelength offset correction of '+num2str(l_offset)+'A to lambda'
	  l   = l   + l_offset
	  lch = lch + l_offset
	endif

	status = dstatus && tstatus1 && tstatus2 && lstatus

	d={t:t, $
	   l:l,llabel:n2g('lambda')+' ['+string(197B)+']', $
	   d:d,dlabel:'[AU]', rzview:rzview,sig:sig,$
	   dch:dch,lch:lch,sigch:sigch,ch_nam:ch_nam, $
	   c:c, dstatus:dstatus,tstatus1:tstatus1,tstatus2:tstatus2,lstatus:lstatus,$
	   title:'chromex',s:status,p:status,index:0,msg:msg }
END

PRO load_cmod_cxrs,shot,d
; from /home/rmchurch/IDLWorkspace/cxrs_analysis/cxrs_separate.pro
	mdsopen, 'dnb', shot
	    data1 = mdsvalue('\dnb::top.mit_cxrs.photon1.analysis:raw_data',status=stat_d1) ; dim[nl nch nt]
	    data2 = mdsvalue('\dnb::top.mit_cxrs.photon2.analysis:raw_data',status=stat_d2)
	    ;this allows working when only one PHOTONMAX camera is working
	    if size(data1,/tname) eq 'STRING' then data1=fltarr(size(data2,/dim))
	    if size(data2,/tname) eq 'STRING' then data2=fltarr(size(data1,/dim))
	    ;setup
	    periscope1 = mdsvalue('\dnb::top.mit_cxrs.photon1.setup:periscope')
	    periscope2 = mdsvalue('\dnb::top.mit_cxrs.photon2.setup:periscope')
	    peri_fibr1 = mdsvalue('\dnb::top.mit_cxrs.photon1.setup:peri_fiber')
	    peri_fibr2 = mdsvalue('\dnb::top.mit_cxrs.photon2.setup:peri_fiber')
	    radius1=mdsvalue('\dnb::top.mit_cxrs.photon1.setup:radius')
	    radius2=mdsvalue('\dnb::top.mit_cxrs.photon2.setup:radius')
	    ;calibration
	    pixmin1 = mdsvalue('\dnb::top.mit_cxrs.photon1.calibration:pix_min')
	    pixmax1 = mdsvalue('\dnb::top.mit_cxrs.photon1.calibration:pix_max')
	    pixmin2 = mdsvalue('\dnb::top.mit_cxrs.photon2.calibration:pix_min')
	    pixmax2 = mdsvalue('\dnb::top.mit_cxrs.photon2.calibration:pix_max')
	    wave_c1 = mdsvalue('\dnb::top.mit_cxrs.photon1.calibration:wavelength')
	    wave_c2 = mdsvalue('\dnb::top.mit_cxrs.photon2.calibration:wavelength')
	    wave_a1 = mdsvalue('\dnb::top.mit_cxrs.photon1.analysis:wavelength',/quiet)
	    wave_a2 = mdsvalue('\dnb::top.mit_cxrs.photon2.analysis:wavelength',/quiet)
	    lin1 = mdsvalue('\dnb::top.mit_cxrs.photon1.calibration:lin_disp')
	    lin2 = mdsvalue('\dnb::top.mit_cxrs.photon2.calibration:lin_disp')
	    inst1 = mdsvalue('\dnb::top.mit_cxrs.photon1.calibration:inst_funct')
	    inst2 = mdsvalue('\dnb::top.mit_cxrs.photon2.calibration:inst_funct')
	    fltr1 = mdsvalue('\dnb::top.mit_cxrs.photon1.calibration:fltr_funct',/quiet)
	    fltr2 = mdsvalue('\dnb::top.mit_cxrs.photon2.calibration:fltr_funct',/quiet)
	    fltrw1 = mdsvalue('dim_of(\dnb::top.mit_cxrs.photon1.calibration:fltr_funct)',/quiet)
	    fltrw2 = mdsvalue('dim_of(\dnb::top.mit_cxrs.photon2.calibration:fltr_funct)',/quiet)
	    ;position. created FY09
	    R_pol_sig=mdsvalue('\dnb::top.mit_cxrs.position:R_pol',/quiet)
	    R_pol_bck=mdsvalue('\dnb::top.mit_cxrs.position:R_pol_bck',/quiet)
	    R_tor_sig=mdsvalue('\dnb::top.mit_cxrs.position:R_tor_sig',/quiet)
	    R_IW_sig=mdsvalue('\dnb::top.mit_cxrs.position:R_IW_SIG',/quiet)
	    R_IW_bck=mdsvalue('\dnb::top.mit_cxrs.position:R_IW_BCK',/quiet)
	    R_IW_pol_sig=mdsvalue('\dnb::top.mit_cxrs.position:R_IW_POL_SIG',/quiet)
	    R_IW_pol_bck=mdsvalue('\dnb::top.mit_cxrs.position:R_IW_POL_BCK',/quiet)
	    R_out_puff=mdsvalue('\dnb::top.mit_cxrs.position:R_out_puff',/quiet)
	    Z_pol_sig=mdsvalue('\dnb::top.mit_cxrs.position:Z_pol',/quiet)
	    Z_pol_bck=mdsvalue('\dnb::top.mit_cxrs.position:Z_pol_bck',/quiet)
	    Z_tor_sig=mdsvalue('\dnb::top.mit_cxrs.position:Z_tor_sig',/quiet)
	    Z_IW_sig=mdsvalue('\dnb::top.mit_cxrs.position:Z_IW_SIG',/quiet)
	    Z_IW_bck=mdsvalue('\dnb::top.mit_cxrs.position:Z_IW_BCK',/quiet)
	    Z_IW_pol_sig=mdsvalue('\dnb::top.mit_cxrs.position:Z_IW_POL_SIG',/quiet)
	    Z_IW_pol_bck=mdsvalue('\dnb::top.mit_cxrs.position:Z_IW_POL_BCK',/quiet)
	    Z_out_puff=mdsvalue('\dnb::top.mit_cxrs.position:Z_out_puff',/quiet)
	   
	    times = mdsvalue('\dnb::top.mit_cxrs.datx:times')
	    kin1 = mdsvalue('\DNB::TOP.MIT_CXRS.PHOTON1.KINETICS:SUBTIME', /quiet, status = kin1_ok)
	    kin2 = mdsvalue('\DNB::TOP.MIT_CXRS.PHOTON2.KINETICS:SUBTIME', /quiet, status = kin2_ok)
	    if not kin1_ok then kin1 = 0
	    if not kin2_ok then kin2 = 0
	    if kin1 gt 0 then data1 = mdsvalue('\dnb::top.mit_cxrs.photon1.kinetics:kin_data')
	    if kin2 gt 0 then data2 = mdsvalue('\dnb::top.mit_cxrs.photon2.kinetics:kin_data')

	    ;ADC Setup. If none (only created for FY2009+), use current default.    
	    gain1=mdsvalue('\dnb::top.mit_cxrs.photon1.setup:gain')
	    gain2=mdsvalue('\dnb::top.mit_cxrs.photon2.setup:gain')
	    gain_counts1=mdsvalue('\dnb::top.mit_cxrs.photon1.setup:gain_counts',status=statadc,/quiet)
	    if statadc mod 2 then begin
		gain_counts1=mdsvalue('\dnb::top.mit_cxrs.photon1.setup:gain_counts')
		gain_counts2=mdsvalue('\dnb::top.mit_cxrs.photon2.setup:gain_counts')
		read_noise1=mdsvalue('\dnb::top.mit_cxrs.photon1.setup:read_noise')
		read_noise2=mdsvalue('\dnb::top.mit_cxrs.photon2.setup:read_noise')
	    endif else begin
		mdsopen,'dnb',-1
		gain_counts1=mdsvalue('\dnb::top.mit_cxrs.photon1.setup:gain_counts')
		gain_counts2=mdsvalue('\dnb::top.mit_cxrs.photon2.setup:gain_counts')
		read_noise1=mdsvalue('\dnb::top.mit_cxrs.photon1.setup:read_noise')
		read_noise2=mdsvalue('\dnb::top.mit_cxrs.photon2.setup:read_noise')
	    endelse
	mdsclose

        if size(times,/tname) eq 'STRING' then begin
	 if times[0] eq '*' then print,'No cxrs timing available for this shot'
         times=float(times)
        endif

        ;code to remove pesky misplaced readouts
        t1=times[*,0]
        t2=times[*,1]
        diff1=t1[1:*]-t1[0:n_elements(t1)-2]
        diff2=t2[1:*]-t2[0:n_elements(t2)-2]

        ;get frames that are off. add a +1 since diff offshifted. this assumes most frames will be correct
        inds1=where(diff1 lt median(diff1)*0.9,count1)+1
        inds2=where(diff2 lt median(diff2)*0.9,count2)+1
        if count1 gt 0 then t1[inds1]=!values.f_nan
        if count2 gt 0 then t2[inds2]=!values.f_nan
        len=n_elements(t1)<n_elements(t2)
        tind1=fltarr(len)
        tind2=fltarr(len)
        good=0
        for i=0,len-1 do begin
            ind=where(abs(t1[i]-t2) lt 0.0005,count)
            if count gt 0 then begin
                tind1[good]=where(t1[i] eq times[*,0])
                tind2[good]=where(t2[ind[0]] eq times[*,1])
                good=good+1
            endif
        endfor
        tind1=tind1[0:good-1]
        tind2=tind2[0:good-1]
        
        if n_elements(mintime) eq 0 then mintime=0.0
        if n_elements(maxtime) eq 0 then maxtime=max(times)
        if mintime ge maxtime then print,'Mintime is greater than maxtime, program not set to handle, stopping'
	
        ;get times within min and max time. Assume
        start1=value_locate(times[*,0],mintime)+1
        start2=value_locate(times[*,1],mintime)+1
        stop1=value_locate(times[*,0],maxtime)
        stop2=value_locate(times[*,1],maxtime)
        inds1=where(tind1 ge start1 and tind1 le stop1,ct1)
        inds2=where(tind2 ge start2 and tind2 le stop2,ct2)
        if ct1 eq -1 or ct2 eq -1 then begin
            message,'No frames in between min and max time, program not set to handle, stopping'
        endif
        tind1=tind1[inds1]
        tind2=tind2[inds2]

        ;check both times arrays, take the max of the mins. This avoids getting the background
        ;of misplaced trigger frames, which have high counts
        backstart1=min(times[2:*,0])
        backstart2=min(times[2:*,1])
        backstart=max([backstart1,backstart2])
        backstop=-0.01
        backind1=where(times[*,0] gt backstart and times[*,0] lt backstop,ctback1)
        backind2=where(times[*,1] gt backstart and times[*,1] lt backstop,ctback2)
        if total(tind1[0] lt backind1) gt 0 or total(tind2[0] lt backind2) gt 0 then begin
            message,'Background frames in between min and max time, program not set to handle, stopping'
        endif

        ;code to take out saturation
        ind1=where(data1 eq 65535,count1)
        ind2=where(data2 eq 65535,count2)

        ;background subtraction
        back1=total(data1[*,*,backind1],3)/ctback1
        back2=total(data2[*,*,backind2],3)/ctback2
        
        for i=0,n_elements(data1[0,0,*])-1 do begin
            data1[*,*,i]=data1[*,*,i]-back1
            data2[*,*,i]=data2[*,*,i]-back2
        endfor

        times=times[tind1,0]
        data1=data1[*,*,tind1]
        data2=data2[*,*,tind2]

    frames=n_elements(times)

    wave1 = n_elements(wave_a1) gt 512? wave_a1 : wave_c1
    wave2 = n_elements(wave_a2) gt 512? wave_a2 : wave_c2
    if n_elements(wave_a1) lt 512 then print, "Using non-calibrated wavelength on Photon 1"
    if n_elements(wave_a2) lt 512 then print, "Using non-calibrated wavelength on Photon 2"
    
    range = indgen(120)
    
    if n_elements(fltrw1) eq 512 then fltrw1 = wave1
    if n_elements(fltrw2) eq 512 then fltrw2 = wave2
    
    ;bright_factor,shot,brightFact1,brightFact2
    
    peri=['POL','POL_BCK','TOR_SIG','IW_SIG','IW_BCK','IW_POL_SIG','IW_POL_BCK','OUT_PUFF']

    for k=0,n_elements(peri)-1 do begin
    
        whichbins1 = where(periscope1 eq peri[k], ct1)
        whichbins2 = where(periscope2 eq peri[k], ct2)
        
        N = ct1 + ct2
        if N eq 0 then begin
            print,'No '+peri[k]+' periscope on this shot'
            continue
        endif
        peri_data = fltarr(120,N,frames)
        R = fltarr(N)
        Z = fltarr(N)
        wave = fltarr(120,N)
        pixels = fltarr(120,N)
        gain_counts = fltarr(N)
        read_noise = fltarr(N)
        inst = fltarr(9,N)
        lin_disp = fltarr(N)
        bright = fltarr(N)
        peri_bright_data = fltarr(120,N,frames)
        view = intarr(N)

        if ct1 ne 0 then begin
            for j = 0, ct1-1 do begin
                bin = peri_fibr1[whichbins1[j]]-1
                view[j] = bin+1
                row = whichbins1[j] mod 18
                col = whichbins1[j]/18
                trash=min(abs(wave1[*,row,col]-4944.67),pixm)
                inds=pixm-60+range
                R[j]=peri_R[bin]
                Z[j]=peri_Z[bin]
                pixels[*,j]=inds
                wave[*,j] = wave1[inds,row,col]
                lin_disp[j]=lin1[row,col]
                if string(fltr1[0]) eq '*' then begin
                    fltr1 = fltarr(512,18,3)+1.0
                endif else fltrShift = interpol(fltr1[inds,row,col],fltrw1[inds,row,col],wave[*,j])
                bright[j]=max(brightFact1[inds,row,col])
                brightFact=brightFact1[inds,row,col]/bright[j]
                for i = 0, frames-1 do begin
                    peri_data[*,j,i] = data1[inds,row,i]/fltrShift*brightFact
                    peri_bright_data[*,j,i] = data1[inds,row,i]/fltrShift*brightFact1[inds,row,col]
                endfor
                ;this assumes readout rate always 5 MHz (change in future if not true)
                gain_counts[j]=gain_counts1[1,gain1-1]
                read_noise[j]=read_noise1[1,gain1-1]
                inst[*,j] = reform(inst1[row,col,*])
                ;inst[[1,2,4,5,7,8],j] = inst[[1,2,4,5,7,8],j]*lin1[r,c]
                inst[[2,5,8],j] = abs(inst[[2,5,8],j])
                inst[[0,3,6],j] = inst[[0,3,6],j]/total(inst[[0,3,6],j]*sqrt(2*!pi)*inst[[2,5,8],j]);renormalize just in case
            endfor
        endif
        
        if ct2 ne 0 then begin
            for n = 0, ct2-1 do begin
                bin = peri_fibr2[whichbins2[n]]-1
                j = n + ct1
                view[j] = bin+1
                row = whichbins2[n] mod 18
                col = whichbins2[n]/18
                trash=min(abs(wave2[*,row,col]-4944.67),pixm)
                inds=pixm-60+range
                R[j]=peri_R[bin]
                Z[j]=peri_Z[bin]
                wave[*,j] = wave2[inds,row,col]
                pixels[*,j]=inds
                lin_disp[j]=lin2[row,col]
                if string(fltr2[0]) eq '*' then begin
                    fltr2 = fltarr(512,18,3)+1.0
                endif else fltrShift = interpol( fltr2[inds,row,col],fltrw2[inds,row,col],wave[*,j])
                bright[j]=max(brightFact2[inds,row,col])
                brightFact=brightFact2[inds,row,col]/bright[j]
                for i = 0, frames-1 do begin
                    peri_data[*,j,i] = data2[inds,row,i]/fltrShift*brightFact
                    peri_bright_data[*,j,i] = data2[inds,row,i]/fltrShift*brightFact2[inds,row,col]
                endfor
                ;this assumes readout rate always 5 MHz (change in future if not true)
                gain_counts[j]=gain_counts2[1,gain2-1]
                read_noise[j]=read_noise2[1,gain2-1]
                inst[*,j] = reform(inst2[row,col,*])
                ;TODO decide pixel or wavelength
                ;inst[[1,2,4,5,7,8],j] = inst[[1,2,4,5,7,8],j]*lin2[r,c]
                inst[[2,5,8],j] = abs(inst[[2,5,8],j])
                inst[[0,3,6],j] = inst[[0,3,6],j]/total(inst[[0,3,6],j]*sqrt(2*!pi)*inst[[2,5,8],j])
            endfor
	endif

        order=sort(view)
        R = R[order]
        Z = Z[order]
        wave = wave[*,order]
        pixels = pixels[*,order]
        inst_funct = inst[*,order]
        lin_disp = lin_disp[order]
        bright = bright[order]
        gain_counts = gain_counts[order]
        read_noise = read_noise[order]
        view = view[order]
        peri_params={R:R,Z:Z,theta:theta,wave:wave,pixels:pixels,inst_funct:inst_funct,lin_disp:lin_disp,bright:bright,gain_counts:gain_counts,read_noise:read_noise,view:view}
        peri_data = peri_data[*,order,*]
        peri_bright_data = peri_bright_data[*,order,*]
        
        if k gt 0 then params=create_struct(params,peri[k],peri_params) else params=create_struct(peri[k],peri_params)
        if k gt 0 then data=create_struct(data,peri[k],peri_data) else data=create_struct(peri[k],peri_data)
        if k gt 0 then bright_data=create_struct(bright_data,peri[k],peri_bright_data) else bright_data=create_struct(peri[k],peri_bright_data)
        
    endfor

END



PRO load_cmod_llnlcid,shot,d,plot=plot
	mdsopen,'video_daq',shot
	 frames=mdsvalue('\VIDEO_DAQ::TOP.LLNL1.FRAMES') 
	 dim=size(frames,/dim) ;[nx ny nfr]
	  nx=dim[0]
	  ny=dim[1]
	  nframes=dim[2]
	 npre=mdsvalue('\VIDEO_DAQ::TOP.LLNL1.NPRE')
	 frames_t=1/30.*(findgen(nframes) - npre)

	 filt_info=mdsvalue('\VIDEO_DAQ::TOP.LLNL1.FILTER_INFO')
	 filt_cwl=mdsvalue('\VIDEO_DAQ::TOP.LLNL1.FILTER_CWL')
	 filt_wid=mdsvalue('\VIDEO_DAQ::TOP.LLNL1.FILTER_WIDTH')

	 gain=mdsvalue('\VIDEO_DAQ::TOP.LLNL1.GAIN')
	 gates_node='\'+mdsvalue('\VIDEO_DAQ::TOP.LLNL1.GATES_NODE')
	  gates_tree=strmid(gates_node,1,strpos(gates_node,'::')-1)
	  mdsopen,gates_tree,shot
	   gates=mdsvalue(gates_node)
	   gates_t=mdsvalue('dim_of('+gates_node+')')
	   if gates_t[0] ne '*' then ngt=n_elements(gates_t) & gates_t=gates_t[0:ngt-2] ; hack to fix different lengths
	  mdsclose,gates_tree,shot
	mdsclose,'video_daq'

	fields=bytarr(nx,ny,2*nframes)
	fields_t=fltarr(2*nframes)
	fields_dt=(frames_t[1]-frames_t[0])/2
	ih=indgen(nx)
	iv0=indgen(ny/2)*2
	ivfr=floor(indgen(480)/2)*2
	for ii=0,nframes-1 do begin
	 for ifield=1,2 do begin
	  it=2*ii + ifield-1
	  ;iv=iv0 + ifield-1
	  iv=ivfr + ifield-1

	  fields[*,*,it]=rotate(frames[ih,iv,ii],7)
	  fields_t[it]=frames_t[ii] + (ifield-2)*fields_dt
	 endfor
	endfor

	if keyword_set(plot) then begin
	 window & colorbar
	 for it=0,2*nframes-1 do begin
	  tv,fields[*,*,it]
	  xyouts,0,0,num2str(shot,1)+' '+num2str(fields_t[it])+'s'
	  wait,fields_dt
	 endfor
	end

	d={frames:frames,frames_t:frames_t, $
	   fields:fields,fields_t:fields_t, $
	   npre:npre, gain:gain, gates_node:gates_node, gates:gates,gates_t:gates_t $
	   }
	help,/str,d & stop
END


PRO load_cmod_2pid,shot,d
	pt='TWOPI_DIODE'
        mdsopen,'spectroscopy',shot
         d=mdsvalue('_sig=\SPECTROSCOPY::TOP.BOLOMETER:TWOPI_DIODE',/quiet,status=status)*2.75/1.0e3          ;define a 2.75 prad-2pi "fudge-factor" [MW]
         t=mdsvalue('dim_of(_sig,0)',/quiet)
        mdsclose,'spectroscopy',shot

	status = status mod 2

	d={t:t, $
	   d:d,dlabel:'[MW]', $
	   pt:pt,title:pt,oplot:0,s:status,yr:[0.0,0.0]}
END

PRO load_cmod_corextomo,shot,d
	pt='XTOMO'
        IF NOT keyword_set(ch) THEN ch=19
        IF ch LT 10 THEN ch_str='0' ELSE ch_str=''
        IF NOT keyword_set(array) THEN array=3
        mdsopen,'xtomo',shot
         d=mdsvalue('\XTOMO::TOP.BRIGHTNESSES.ARRAY_'+num2str(array,1)+':CHORD_'+ch_str+num2str(ch,1),/quiet,status=status)
         t=mdsvalue('dim_of(\XTOMO::TOP.BRIGHTNESSES.ARRAY_'+num2str(array,1)+':CHORD_'+ch_str+num2str(ch,1)+')',/quiet)
        mdsclose,'xtomo',shot

	status = status mod 2

        IF status THEN BEGIN
                bl=mean(d[0:ipt(t,-0.01)])
                d-=bl
        ENDIF

	d={t:t, $
	   d:d,dlabel:'[?]', $
	   pt:pt,title:pt,oplot:0,s:status,yr:[0.0,0.0]}
END

PRO load_cmod_nl04,shot,d
	pt='nl04'
        mdsopen,'electrons',shot
         d=mdsvalue('\ELECTRONS::TOP.TCI.RESULTS:NL_04',/quiet,status=status)*1.0e-20
         t=mdsvalue('dim_of(\ELECTRONS::TOP.TCI.RESULTS:NL_04)',/quiet)
        mdsclose,'electrons',shot

	status = status mod 2

	d={t:t, $
	   d:d,dlabel:'[1e20/m!E3!N]', $
	   pt:pt,title:pt,oplot:0,s:status,yr:[0.0,0.0]}
END

PRO load_cmod_rfpow,shot,d
	pt='ICRH'
        mdsopen, 'rf', shot
         d=mdsvalue('_sig=\rf::RF_POWER_NET',/quiet,status=status)
         t=mdsvalue('dim_of(_sig,0)',/quiet)
        mdsclose,'rf',shot

	status = status mod 2

	d={t:t, $
	   d:d,dlabel:'[MW]', $
	   pt:pt,title:pt,oplot:0,s:status,yr:[0.0,0.0]}
END

PRO load_cmod_lhpow,shot,d
	pt='LH'
        mdsopen, 'lh',shot,/quiet
         d=mdsvalue('_sig=\LH::TOP.RESULTS:NETPOW',/quiet,status=status)
         t=mdsvalue('dim_of(_sig,0)',/quiet)
        mdsclose,/quiet
        if n(d) ne n(t) then begin
         nlh=min([ n(d), n(t) ])
         t=t[0:nlh]
         d=d[0:nlh]
         print,' fixed unequal data and time vector lengths for lh power'
        endif
	if n(d) lt 1 or n(t) lt 1 then status=0

	status = status mod 2

	if status eq 0 then begin
	 d=-1
	 t=-1
	endif

	d={t:t, $
	   d:d/1.e3,dlabel:'[MW]', $
	   pt:pt,title:pt,oplot:1,s:status,yr:[0.0,0.0]}
END

