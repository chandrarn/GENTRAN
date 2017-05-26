PRO load_dalpha_profile,shot,br,rt,time,peri=peri,fiber=fiber
	nd1=6
	nd2=6
	np1=6
	nch=nd1+nd2+np1
	mdsopen,'dnb',shot
	time=mdsvalue('dim_of(\DNB::TOP.MIT_CXRS.DALPHA:PD1_1:BRIGHT,0)')
	ntime=n(time)+1
	xbr=fltarr(nch,ntime)
	peri=strarr(nch)
	fiber=intarr(nch)-1
	k=0
	FOR j=1,nd1 DO BEGIN
		ibr=mdsvalue('\DNB::TOP.MIT_CXRS.DALPHA:PD1_'+num2str(j,1)+':BRIGHT',status=status,/quiet)
		IF status THEN BEGIN
			xbr[k,*]=ibr
			peri[k]=mdsvalue('\DNB::TOP.MIT_CXRS.DALPHA:PD1_'+num2str(j,1)+':PERI')
			fiber[k]=mdsvalue('\DNB::TOP.MIT_CXRS.DALPHA:PD1_'+num2str(j,i)+':PERI_FIBER')
                ENDIF ELSE peri[k]='none'
		k+=1
        ENDFOR

	FOR j=1,nd2 DO BEGIN
		ibr=mdsvalue('\DNB::TOP.MIT_CXRS.DALPHA:PD2_'+num2str(j,1)+':BRIGHT',status=status,/quiet)
		IF status THEN BEGIN
			xbr[k,*]=ibr
			peri[k]=mdsvalue('\DNB::TOP.MIT_CXRS.DALPHA:PD2_'+num2str(j,1)+':PERI')
			fiber[k]=mdsvalue('\DNB::TOP.MIT_CXRS.DALPHA:PD2_'+num2str(j,i)+':PERI_FIBER')
                ENDIF ELSE peri[k]='none'
		k+=1
        ENDFOR
	FOR j=1,np1 DO BEGIN
		ibr=mdsvalue('\DNB::TOP.MIT_CXRS.DALPHA:PMT1_'+num2str(j,1)+':BRIGHT',status=status,/quiet)
		IF status THEN BEGIN
			xbr[k,*]=ibr
			peri[k]=mdsvalue('\DNB::TOP.MIT_CXRS.DALPHA:PMT1_'+num2str(j,1)+':PERI')
			fiber[k]=mdsvalue('\DNB::TOP.MIT_CXRS.DALPHA:PMT1_'+num2str(j,i)+':PERI_FIBER')
                ENDIF ELSE peri[k]='none'
		k+=1
        ENDFOR	
	mdsclose,'dnb',shot
	tmp=where(fiber NE -1 AND peri EQ 'POL_BCK')
	ngood=n(tmp)+1
	xbr=xbr[tmp,*]
	peri=peri[tmp]
	fiber=fiber[tmp]
	order=(sort(fiber))
	br=fltarr(ngood,ntime)
	FOR i=0,ngood-1 DO br[i,*]=xbr[order[i],*]
	fiber=fiber[order]

END
