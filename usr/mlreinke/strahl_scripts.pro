PRO run_gridding_scan
	k=[0.5,1.0,1.5,2.0,2.5,3.0]
	nrho=[50,100,200,300]
	shot=1101209004
	z=18
	t1=1.0
	t2=1.2
	nk=n(k)+1
	nnrho=n(nrho)+1
	FOR i=0,nk-1 DO BEGIN
		FOR j=nnrho-1,nnrho-1 DO BEGIN
			path=0
			print, k[i],nrho[j]
			RUN_CMOD_STRAHL,shot,z,t1,t2,path=path,k=k[i],nrho=nrho[j]
			savepath='/home/mlreinke/idl/impurities/strahl/gridding_scan_'+num2str(shot,1)+'_k'+num2str(i,1)+'_nrho'+num2str(j,1)+'.dat'
			spawn,'cp '+path+' '+savepath			
			print, savepath
                ENDFOR
	ENDFOR
END

PRO plot_gridding_scan,load=load
	IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0
	k=[0.5,1.0,1.5,2.0,2.5,3.0]
	nrho=[50,100,200,300]
	nk=6
	nnrho=4
	nmax=400
	shot=1101209004
	z=18
	t1=1.0
	t2=1.2
	savepath='/home/mlreinke/idl/impurities/strahl/gridding_scan_'+num2str(shot,1)+'.dat'
	IF NOT keyword_set(load) THEN BEGIN
		FOR i=0,nk-1 DO BEGIN
			FOR j=0,nnrho-1 DO BEGIN
				path='/home/mlreinke/idl/impurities/strahl/gridding_scan_'+num2str(shot,1)+'_k'+num2str(i,1)+'_nrho'+num2str(j,1)+'.dat'
				RUN_CMOD_STRAHL,shot,z,t1,t2,path=path,/nostrahl,data=idata
				IF i EQ 0 AND j EQ 0 THEN BEGIN
					psin=fltarr(nmax,nk,nnrho)
					ntime=n(idata.time)
					nz=fltarr(nmax,ntime,nk,nnrho)
					n0=fltarr(nmax,ntime,nk,nnrho)
					etemp=fltarr(nmax,nk,nnrho)
					edens=fltarr(nmax,nk,nnrho)
					time=idata.time
                                ENDIF
				irho=n(idata.rho)+1.0
				psin[0:irho-1,i,j]=idata.psin
				etemp[0:irho-1,i,j]=idata.temp[*,0]
				edens[0:irho-1,i,j]=idata.temp[*,0]
				FOR q=0,z DO nz[0:irho-1,*,i,j]+=idata.csden[*,q,*]
				FOR q=0,0 DO n0[0:irho-1,*,i,j]+=idata.csden[*,q,*]

			ENDFOR
		ENDFOR			
		save,psin,nz,n0,etemp,edens,time,filename=savepath
        ENDIF ELSE restore, savepath
	
	IF keyword_set(ps) THEN BEGIN
		device, xsize=5.0,ysize=5.5,/inches
		ls=1.0
		!p.thick=6
		!p.charthick=2.0
	ENDIF ELSE openwin,0,xsize=500,ysize=550
	ntau=n(nz[0,*,0,0])+1
	xr=[0,3.5]
	;yr=[0.8*min(nz[0,ntau-1,*,*]),max(nz[0,ntau-1,*,*])*1.2
	yr=[1.0e8,1.0e31]
	plot,[0],[0],xr=xr,yr=yr,xtit='k',ytit='n!lZ,0!n @ END',chars=1.2,/xsty,/ysty,/ylog
	color=[0,30,100,200]
	FOR i=0,nnrho-1 DO BEGIN
		oplot,k,nz[0,ntau-1,*,i],color=color[i]
		xyouts,xr[0]+0.6*(xr[1]-xr[0]),10^(alog10(yr[0])+(0.9-0.1*i)*(alog10(yr[1])-alog10(yr[0]))),'NRHO='+num2str(nrho[i],1),color=color[i]
	ENDFOR

	IF keyword_set(ps) THEN BEGIN
		device, xsize=7.0,ysize=4.5,/inches
		ls=1.0
		!p.thick=6
		!p.charthick=2.0
	ENDIF ELSE openwin,0,xsize=1400,ysize=900
	openwin,1
	!p.multi=[0,3,2]
	xr=[0.98,1.08]
	color=[0,30,100,200]
	yr=[1.0,1.0e30]
	FOR i=0,nk-1 DO BEGIN
		plot,[0],[0],xr=xr,yr=yr,xtit='PSIN',ytit='n!l0!n @ START',chars=1.4,/xsty,/ysty,/ylog,tit='k='+num2str(k[i],dp=1)
		FOR j=0,nnrho-1 DO BEGIN
			tmp=where(psin[*,i,j] NE 0)
			oplot,psin[tmp,i,j],n0[tmp,1,i,j],color=color[j],psym=-5
			xyouts,xr[0]+0.2*(xr[1]-xr[0]),10^(alog10(yr[0])+(0.9-0.1*j)*(alog10(yr[1])-alog10(yr[0]))),'NRHO='+num2str(nrho[j],1),color=color[j],chars=0.5

		ENDFOR
	ENDFOR

	openwin,2
	!p.multi=[0,3,2]
	xr=[0.0,1.1]
	color=[0,30,100,200]
	yr=[0.0,1.05]
	FOR i=0,nk-1 DO BEGIN
		plot,[0],[0],xr=xr,yr=yr,xtit='PSIN',ytit='n!lZ!n(t='+num2str(time[ntau/10.0],dp=2)+')/n!lZ!n(t='+num2str(time[ntau-1],dp=2)+')',chars=1.4,/xsty,/ysty,tit='k='+num2str(k[i],dp=1)
		FOR j=0,nnrho-1 DO BEGIN
			tmp=where(psin[*,i,j] LT 1.0 AND psin[*,i,j] NE 0)
			oplot,psin[tmp,i,j],nz[tmp,ntau/10.0,i,j]/nz[tmp,ntau-1,i,j],color=color[j]
			;oplot,psin[tmp,i,j],nz[tmp,ntau/5.0,i,j]/nz[tmp,ntau-1,i,j],color=color[j]
			xyouts,xr[0]+0.65*(xr[1]-xr[0]),yr[0]+(0.5-0.1*j)*(yr[1]-yr[0]),'NRHO='+num2str(nrho[j],1),color=color[j],chars=0.5

                ENDFOR
        ENDFOR

	!p.multi=0
	xr=[0.0,1.1]
	color=[0,30,100,200]
	yr=[0.0,1.05]
	xk=[1,3,5]
	style=[0,1,2]
	plot,[0],[0],xr=xr,yr=yr,xtit='PSIN',ytit='n!lZ!n(t='+num2str(time[ntau/10.0],dp=2)+')/n!lZ!n(t='+num2str(time[ntau-1],dp=2)+')',chars=1.2,/xsty,/ysty
	FOR ii=0,n(xk) DO BEGIN
		i=xk[ii]
		FOR j=0,nnrho-1 DO BEGIN
			tmp=where(psin[*,i,j] LT 1.0 AND psin[*,i,j] NE 0)
			oplot,psin[tmp,i,j],nz[tmp,ntau/10.0,i,j]/nz[tmp,ntau-1,i,j],color=color[j],linestyle=style[ii]
			;oplot,psin[tmp,i,j],nz[tmp,ntau/5.0,i,j]/nz[tmp,ntau-1,i,j],color=color[j]
			IF ii EQ 0 THEN xyouts,xr[0]+0.8*(xr[1]-xr[0]),yr[0]+(0.5-0.1*j)*(yr[1]-yr[0]),'NRHO='+num2str(nrho[j],1),color=color[j],chars=1.0
		ENDFOR
		oplot,xr[0]+[0.1,0.18]*(xr[1]-xr[0]),yr[0]+[1.0,1.0]*(0.8-0.1*ii)*(yr[1]-yr[0]),linestyle=style[ii]
		xyouts,xr[0]+0.2*(xr[1]-xr[0]),yr[0]+(0.79-0.1*ii)*(yr[1]-yr[0]),'k='+num2str(k[i],dp=1),chars=1.0
	ENDFOR


	openwin,3
	!p.multi=[0,3,2]
	xr=[min(time),min(time)+0.05]
	color=[0,30,100,200]
	yr=[0,1.05]
	xrho=0.5
	FOR i=0,nk-1 DO BEGIN
		plot,[0],[0],xr=xr,yr=yr,xtit='Time [sec]',ytit='n!lZ!n(t)/n!lZ!n(t='+num2str(time[ntau-1],dp=2)+') @ PSIN='+num2str(xrho,dp=2),chars=1.4,/xsty,/ysty,tit='k='+num2str(k[i],dp=1)
		FOR j=0,nnrho-1 DO BEGIN
			index=ipt(psin[*,i,j],xrho)
			oplot,time,nz[index,*,i,j]/nz[index,ntau-1,i,j],color=color[j]
			xyouts,xr[0]+0.65*(xr[1]-xr[0]),yr[0]+(0.5-0.1*j)*(yr[1]-yr[0]),'NRHO='+num2str(nrho[j],1),color=color[j],chars=0.5
		ENDFOR
	ENDFOR

	!p.multi=0
	xr=[min(time),min(time)+0.05]
	color=[0,30,100,200]
	yr=[0,1.05]
	xrho=0.5
	xk=[1,3,5]
	style=[0,1,2]
	plot,[0],[0],xr=xr,yr=yr,xtit='Time [sec]',ytit='n!lZ!n(t)/n!lZ!n(t='+num2str(time[ntau-1],dp=2)+') @ PSIN='+num2str(xrho,dp=2),chars=1.4,/xsty,/ysty
	FOR ii=0,n(xk) DO BEGIN
		i=xk[ii]
		FOR j=0,nnrho-1 DO BEGIN
			index=ipt(psin[*,i,j],xrho)
			oplot,time,nz[index,*,i,j]/nz[index,ntau-1,i,j],color=color[j],linestyle=style[ii]
			IF ii EQ 0 THEN xyouts,xr[0]+0.65*(xr[1]-xr[0]),yr[0]+(0.5-0.1*j)*(yr[1]-yr[0]),'NRHO='+num2str(nrho[j],1),color=color[j],chars=1.0
		ENDFOR
		oplot,xr[0]+[0.1,0.18]*(xr[1]-xr[0]),yr[0]+[1.0,1.0]*(0.7-0.1*ii)*(yr[1]-yr[0]),linestyle=style[ii]
		xyouts,xr[0]+0.2*(xr[1]-xr[0]),yr[0]+(0.69-0.1*ii)*(yr[1]-yr[0]),'k='+num2str(k[i],dp=1),chars=1.0
	ENDFOR


	openwin,4
	!p.multi=[0,3,2]
	xr=[0.0,1.1]
	color=[0,30,100,200]
	yr=[-3,0]
	FOR i=0,nk-1 DO BEGIN
		plot,[0],[0],xr=xr,yr=yr,xtit='PSIN',ytit='-1/n!lZ!n dn!lZ!n/dPSIN @ t='+num2str(time[ntau/8.0],dp=2)+')',chars=1.4,/xsty,/ysty,tit='k='+num2str(k[i],dp=1)
		FOR j=0,nnrho-1 DO BEGIN
			tmp=where(psin[*,i,j] LT 1.0 AND psin[*,i,j] NE 0)
			index=ntau/8.0
			oplot,psin[tmp,i,j],-1.0*deriv(psin[tmp,i,j],nz[tmp,index,i,j])/nz[tmp,index,i,j],color=color[j]
			;stop
			xyouts,xr[0]+0.50*(xr[1]-xr[0]),yr[0]+(0.5-0.1*j)*(yr[1]-yr[0]),'NRHO='+num2str(nrho[j],1),color=color[j],chars=0.5

                ENDFOR
        ENDFOR

	!p.multi=[0,3,2]
	xr=[0.0,1.1]
	color=[0,30,100,200]
	yr=[-3,0]
	FOR i=0,nk-1 DO BEGIN
		plot,[0],[0],xr=xr,yr=yr,xtit='PSIN',ytit='-1/n!lZ!n dn!lZ!n/dPSIN @ t='+num2str(time[ntau/4.0],dp=2)+')',chars=1.4,/xsty,/ysty,tit='k='+num2str(k[i],dp=1)
		FOR j=0,nnrho-1 DO BEGIN
			tmp=where(psin[*,i,j] LT 1.0 AND psin[*,i,j] NE 0)
			index=ntau/4.0
			oplot,psin[tmp,i,j],-1.0*deriv(psin[tmp,i,j],nz[tmp,index,i,j])/nz[tmp,index,i,j],color=color[j]
			;stop
			xyouts,xr[0]+0.50*(xr[1]-xr[0]),yr[0]+(0.5-0.1*j)*(yr[1]-yr[0]),'NRHO='+num2str(nrho[j],1),color=color[j],chars=0.5

                ENDFOR
        ENDFOR
	!p.multi=0
	xr=[0.0,1.1]
	color=[0,30,100,200]
	yr=[-3,0]
	xk=[1,3,5]
	style=[0,1,2]
	plot,[0],[0],xr=xr,yr=yr,xtit='PSIN',ytit='-1/n!lZ!n dn!lZ!n/dPSIN @ t='+num2str(time[ntau/8.0],dp=2)+')',chars=1.4,/xsty,/ysty
	FOR ii=0,n(xk) DO BEGIN
		i=xk[ii]
		FOR j=0,nnrho-1 DO BEGIN
			tmp=where(psin[*,i,j] LT 1.0 AND psin[*,i,j] NE 0)
			index=ntau/8.0
			oplot,psin[tmp,i,j],-1.0*deriv(psin[tmp,i,j],nz[tmp,index,i,j])/nz[tmp,index,i,j],color=color[j],linestyle=style[ii]
			IF ii EQ 0 THEN xyouts,xr[0]+0.65*(xr[1]-xr[0]),yr[0]+(0.5-0.1*j)*(yr[1]-yr[0]),'NRHO='+num2str(nrho[j],1),color=color[j],chars=1.0
		ENDFOR
		oplot,xr[0]+[0.2,0.28]*(xr[1]-xr[0]),yr[0]+[1.0,1.0]*(0.4-0.1*ii)*(yr[1]-yr[0]),linestyle=style[ii]
		xyouts,xr[0]+0.3*(xr[1]-xr[0]),yr[0]+(0.39-0.1*ii)*(yr[1]-yr[0]),'k='+num2str(k[i],dp=1),chars=1.0
	ENDFOR


	!p.multi=0
	stop
	
END

PRO run_dredge_scan
	k=[1.0,2.0,3.0]
	dredge=[0.07,0.1,0.2,0.3]
	shot=1101209004
	z=18
	t1=1.0
	t2=1.2
	nk=n(k)+1
	ndr=n(dredge)+1
	FOR i=0,nk-1 DO BEGIN
		FOR j=0,ndr-1 DO BEGIN
			path=0
			print, k[i],dredge[j]
			RUN_CMOD_STRAHL,shot,z,t1,t2,path=path,k=k[i],dr=[0.5,dredge[j]]
			savepath='/home/mlreinke/idl/impurities/strahl/dredge_scan_'+num2str(shot,1)+'_k'+num2str(i,1)+'_dr'+num2str(j,1)+'.dat'
			spawn,'cp '+path+' '+savepath			
			print, savepath
                ENDFOR
	ENDFOR
END

PRO plot_dredge_scan,load=load
	IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0
	k=[1.0,2.0,3.0]
	dredge=[0.07,0.1,0.2,0.3]
	nk=3
	ndr=4
	nmax=400
	shot=1101209004
	z=18
	t1=1.0
	t2=1.2
	savepath='/home/mlreinke/idl/impurities/strahl/dredge_scan_'+num2str(shot,1)+'.dat'
	IF NOT keyword_set(load) THEN BEGIN
		FOR i=0,nk-1 DO BEGIN
			FOR j=0,ndr-1 DO BEGIN
				path='/home/mlreinke/idl/impurities/strahl/dredge_scan_'+num2str(shot,1)+'_k'+num2str(i,1)+'_dr'+num2str(j,1)+'.dat'
				RUN_CMOD_STRAHL,shot,z,t1,t2,path=path,/nostrahl,data=idata
				IF i EQ 0 AND j EQ 0 THEN BEGIN
					psin=fltarr(nmax,nk,ndr)-1.0
					ntime=n(idata.time)
					nz=fltarr(nmax,ntime,nk,ndr)
					n0=fltarr(nmax,ntime,nk,ndr)
					etemp=fltarr(nmax,nk,ndr)
					edens=fltarr(nmax,nk,ndr)
              				time=idata.time
		                ENDIF
				irho=n(idata.rho)+1.0
				psin[0:irho-1,i,j]=idata.psin
				etemp[0:irho-1,i,j]=idata.temp[*,0]
				edens[0:irho-1,i,j]=idata.temp[*,0]
				FOR q=0,z DO nz[0:irho-1,*,i,j]+=idata.csden[*,q,*]
				FOR q=0,0 DO n0[0:irho-1,*,i,j]+=idata.csden[*,q,*]

			ENDFOR
		ENDFOR			
		save,psin,nz,n0,etemp,edens,time,filename=savepath
        ENDIF ELSE restore, savepath
	
	IF keyword_set(ps) THEN BEGIN
		device, xsize=5.0,ysize=5.5,/inches
		ls=1.0
		!p.thick=6
		!p.charthick=2.0
	ENDIF ELSE openwin,0,xsize=500,ysize=550
	ntau=n(nz[0,*,0,0])+1
	xr=[0,3.5]
	;yr=[0.8*min(nz[0,ntau-1,*,*]),max(nz[0,ntau-1,*,*])*1.2]
	yr=[1.0e8,1.0e31]
	plot,[0],[0],xr=xr,yr=yr,xtit='k',ytit='n!lZ!n @ END',chars=1.2,/xsty,/ysty,/ylog,tit='RCORE=0.5'
	color=[0,30,100,200]
	FOR i=0,ndr-1 DO BEGIN
		oplot,k,nz[0,ntau-1,*,i],color=color[i]
		xyouts,xr[0]+0.6*(xr[1]-xr[0]),10^(alog10(yr[0])+(0.9-0.1*i)*(alog10(yr[1])-alog10(yr[0]))),'REDGE='+num2str(dredge[i],dp=2),color=color[i]
	ENDFOR

	IF keyword_set(ps) THEN BEGIN
		device, xsize=7.0,ysize=2.25,/inches
		ls=1.0
		!p.thick=6
		!p.charthick=2.0
	ENDIF ELSE openwin,0,xsize=1400,ysize=900

	openwin,1
	!p.multi=[0,3,1]
	xr=[0.95,1.1]
	color=[0,30,100,200]
	yr=[1.0,1.0e30]
	FOR i=0,nk-1 DO BEGIN
		plot,[0],[0],xr=xr,yr=yr,xtit='PSIN',ytit='n!l0!n @ START',chars=1.2,/xsty,/ysty,/ylog,tit='k='+num2str(k[i],dp=1)
		FOR j=0,ndr-1 DO BEGIN
			tmp=where(psin[*,i,j] NE -1.0)
			oplot,psin[tmp,i,j],n0[tmp,1,i,j],color=color[j],psym=-2
			xyouts,xr[0]+0.2*(xr[1]-xr[0]),10^(alog10(yr[0])+(0.9-0.1*j)*(alog10(yr[1])-alog10(yr[0]))),'REDGE='+num2str(dredge[j],dp=2),color=color[j],chars=0.5

		ENDFOR
	ENDFOR

	IF keyword_set(ps) THEN BEGIN
		device, xsize=7.0,ysize=4.5,/inches
		ls=1.0
		!p.thick=6
		!p.charthick=2.0
	ENDIF ELSE openwin,0,xsize=1400,ysize=900
	!p.multi=0
	xr=[0.0,1.1]
	color=[0,30,100,200]
	yr=[-3,0]
	xk=[0,1,2]
	style=[0,1,2]
	plot,[0],[0],xr=xr,yr=yr,xtit='PSIN',ytit='-1/n!lZ!n dn!lZ!n/dPSIN @ t='+num2str(time[ntau/8.0],dp=2)+')',chars=1.4,/xsty,/ysty,tit='RCORE=0.5'
	FOR ii=0,n(xk) DO BEGIN
		i=xk[ii]
		FOR j=0,ndr-1 DO BEGIN
			tmp=where(psin[*,i,j] NE 0)
			index=ntau/8.0
			oplot,psin[tmp,i,j],-1.0*deriv(psin[tmp,i,j],nz[tmp,index,i,j])/nz[tmp,index,i,j],color=color[j],linestyle=style[ii]
			IF ii EQ 0 THEN xyouts,xr[0]+0.65*(xr[1]-xr[0]),yr[0]+(0.5-0.1*j)*(yr[1]-yr[0]),'REDGE='+num2str(dredge[j],dp=2),color=color[j],chars=1.0
		ENDFOR
		oplot,xr[0]+[0.2,0.28]*(xr[1]-xr[0]),yr[0]+[1.0,1.0]*(0.4-0.1*ii)*(yr[1]-yr[0]),linestyle=style[ii]
		xyouts,xr[0]+0.3*(xr[1]-xr[0]),yr[0]+(0.39-0.1*ii)*(yr[1]-yr[0]),'k='+num2str(k[i],dp=1),chars=1.0
	ENDFOR

	!p.multi=0
	stop
	
END

PRO run_drcore_scan
	k=[1.0,2.0,3.0]
	drcore=[0.1,0.5,1.0,2.0]
	shot=1101209004
	z=18
	t1=1.0
	t2=1.2
	nk=n(k)+1
	ndr=n(drcore)+1
	FOR i=0,nk-1 DO BEGIN
		FOR j=0,ndr-1 DO BEGIN
			path=0
			print, k[i],drcore[j]
			RUN_CMOD_STRAHL,shot,z,t1,t2,path=path,k=k[i],dr=[drcore[j],0.2]
			savepath='/home/mlreinke/idl/impurities/strahl/drcore_scan_'+num2str(shot,1)+'_k'+num2str(i,1)+'_dr'+num2str(j,1)+'.dat'
			spawn,'cp '+path+' '+savepath			
			print, savepath
                ENDFOR
	ENDFOR
END


PRO plot_drcore_scan,load=load
	IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0
	k=[1.0,2.0,3.0]
	drcore=[0.1,0.5,1.0,2.0]
	nk=3
	ndr=4
	nmax=400
	shot=1101209004
	z=18
	t1=1.0
	t2=1.2
	savepath='/home/mlreinke/idl/impurities/strahl/drcore_scan_'+num2str(shot,1)+'.dat'
	IF NOT keyword_set(load) THEN BEGIN
		FOR i=0,nk-1 DO BEGIN
			FOR j=0,ndr-1 DO BEGIN
				path='/home/mlreinke/idl/impurities/strahl/drcore_scan_'+num2str(shot,1)+'_k'+num2str(i,1)+'_dr'+num2str(j,1)+'.dat'
				RUN_CMOD_STRAHL,shot,z,t1,t2,path=path,/nostrahl,data=idata
				IF i EQ 0 AND j EQ 0 THEN BEGIN
					psin=fltarr(nmax,nk,ndr)-1.0
					ntime=n(idata.time)
					nz=fltarr(nmax,ntime,nk,ndr)
					n0=fltarr(nmax,ntime,nk,ndr)
					etemp=fltarr(nmax,nk,ndr)
					edens=fltarr(nmax,nk,ndr)
              				time=idata.time
		                ENDIF
				irho=n(idata.rho)+1.0
				psin[0:irho-1,i,j]=idata.psin
				etemp[0:irho-1,i,j]=idata.temp[*,0]
				edens[0:irho-1,i,j]=idata.temp[*,0]
				FOR q=0,z DO nz[0:irho-1,*,i,j]+=idata.csden[*,q,*]
				FOR q=0,0 DO n0[0:irho-1,*,i,j]+=idata.csden[*,q,*]

			ENDFOR
		ENDFOR			
		save,psin,nz,n0,etemp,edens,time,filename=savepath
        ENDIF ELSE restore, savepath
	
	IF keyword_set(ps) THEN BEGIN
		device, xsize=5.0,ysize=5.5,/inches
		ls=1.0
		!p.thick=6
		!p.charthick=2.0
	ENDIF ELSE openwin,0,xsize=500,ysize=550
	ntau=n(nz[0,*,0,0])+1
	xr=[0,3.5]
	;yr=[0.8*min(nz[0,ntau-1,*,*]),max(nz[0,ntau-1,*,*])*1.2]
	yr=[1.0e8,1.0e31]
	plot,[0],[0],xr=xr,yr=yr,xtit='k',ytit='n!lZ!n @ END',chars=1.2,/xsty,/ysty,/ylog,tit='REDGE=0.2'
	color=[0,30,100,200]
	FOR i=0,ndr-1 DO BEGIN
		oplot,k,nz[0,ntau-1,*,i],color=color[i]
		xyouts,xr[0]+0.6*(xr[1]-xr[0]),10^(alog10(yr[0])+(0.9-0.1*i)*(alog10(yr[1])-alog10(yr[0]))),'RCORE='+num2str(drcore[i],dp=2),color=color[i]
	ENDFOR

	IF keyword_set(ps) THEN BEGIN
		device, xsize=7.0,ysize=2.25,/inches
		ls=1.0
		!p.thick=6
		!p.charthick=2.0
	ENDIF ELSE openwin,0,xsize=1400,ysize=900

	openwin,1
	!p.multi=[0,3,1]
	xr=[0.95,1.1]
	color=[0,30,100,200]
	yr=[1.0,1.0e30]
	FOR i=0,nk-1 DO BEGIN
		plot,[0],[0],xr=xr,yr=yr,xtit='PSIN',ytit='n!l0!n @ START',chars=1.2,/xsty,/ysty,/ylog,tit='k='+num2str(k[i],dp=1)
		FOR j=0,ndr-1 DO BEGIN
			tmp=where(psin[*,i,j] NE -1.0)
			oplot,psin[tmp,i,j],n0[tmp,1,i,j],color=color[j],psym=-2
			xyouts,xr[0]+0.2*(xr[1]-xr[0]),10^(alog10(yr[0])+(0.9-0.1*j)*(alog10(yr[1])-alog10(yr[0]))),'RCORE='+num2str(drcore[j],dp=2),color=color[j],chars=0.5

		ENDFOR
	ENDFOR

	IF keyword_set(ps) THEN BEGIN
		device, xsize=7.0,ysize=4.5,/inches
		ls=1.0
		!p.thick=6
		!p.charthick=2.0
	ENDIF ELSE openwin,0,xsize=1400,ysize=900
	!p.multi=0
	xr=[0.0,1.1]
	color=[0,30,100,200]
	yr=[-3,0]
	xk=[0,1,2]
	style=[0,1,2]
	plot,[0],[0],xr=xr,yr=yr,xtit='PSIN',ytit='-1/n!lZ!n dn!lZ!n/dPSIN @ t='+num2str(time[ntau/8.0],dp=2)+')',chars=1.4,/xsty,/ysty,tit='REDGE=0.2'
	FOR ii=0,n(xk) DO BEGIN
		i=xk[ii]
		FOR j=0,ndr-1 DO BEGIN
			tmp=where(psin[*,i,j] NE 0)
			index=ntau/8.0
			oplot,psin[tmp,i,j],-1.0*deriv(psin[tmp,i,j],nz[tmp,index,i,j])/nz[tmp,index,i,j],color=color[j],linestyle=style[ii]
			IF ii EQ 0 THEN xyouts,xr[0]+0.65*(xr[1]-xr[0]),yr[0]+(0.5-0.1*j)*(yr[1]-yr[0]),'RCORE='+num2str(drcore[j],dp=2),color=color[j],chars=1.0
		ENDFOR
		oplot,xr[0]+[0.2,0.28]*(xr[1]-xr[0]),yr[0]+[1.0,1.0]*(0.4-0.1*ii)*(yr[1]-yr[0]),linestyle=style[ii]
		xyouts,xr[0]+0.3*(xr[1]-xr[0]),yr[0]+(0.39-0.1*ii)*(yr[1]-yr[0]),'k='+num2str(k[i],dp=1),chars=1.0
	ENDFOR

	!p.multi=0
	stop
	
END

PRO run_evoltemp_scan

	shot=1101209004
	z=18
	t1=1.0
	t2=1.2
	k=1.0
	dr=[0.5,0.1]
	savepath='/home/mlreinke/idl/impurities/strahl/evoltemp_'+num2str(shot,1)+'.dat'
	RUN_CMOD_STRAHL,shot,z,t1,t2,k=k,dr=dr,data=avedat

	temp={temp:avedat.temp,rho:avedat.psin,time:avedat.time}
	dens={dens:avedat.dens,rho:avedat.psin,time:avedat.time}
	time=make(1.0,1.19,20)
	dt=time*0.0+0.01
	RUN_CMOD_STRAHL,shot,z,time,dt,k=k,dr=dr,data=constdat,temp=temp,dens=dens,/nogrid

	tmp=where(avedat.time GE 0.5*(t1+t2))
	FOR i=0,n(tmp) DO temp.temp[*,tmp[i]]*=0.8
	RUN_CMOD_STRAHL,shot,z,time,dt,k=k,dr=dr,data=evoldat,temp=temp,dens=dens,/nogrid

	save,avedat,constdat,evoldat,filename=savepath
END

PRO plot_evoltemp_scan
	IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0
	shot=1101209004
	z=18
	savepath='/home/mlreinke/idl/impurities/strahl/evoltemp_'+num2str(shot,1)+'.dat'
	restore,savepath

	IF keyword_set(ps) THEN BEGIN
		device, xsize=6.0,ysize=8.5,/inches
		ls=1.0
		!p.thick=3
		!p.charthick=2.0
	ENDIF ELSE openwin,0,xsize=600,ysize=850
	xr=[1.0,1.2]
	nrho=n(avedat.psin)+1
	!p.multi=[0,0,3]

	index=ipt(avedat.rho,1.02)
	yr=[0,4.0]
	plot,[0],[0],xr=xr,yr=yr,xtit='Time [sec]',ytit='T!le!n [keV] r/a=0.0',chars=1.8,/xsty,/ysty
	oplot,avedat.time,avedat.temp[0,*]/1.0e3
	oplot,constdat.time,constdat.temp[0,*]/1.0e3,linestyle=2.0,color=30
	oplot,evoldat.time,evoldat.temp[0,*]/1.0e3,linestyle=3.0,color=200
;	oplot,avedat.time,avedat.temp[index,*]/1.0e3*10
;	oplot,constdat.time,constdat.temp[index,*]/1.0e*10,linestyle=2.0,color=30
;	oplot,evoldat.time,evoldat.temp[index,*]/1.0e3*10,linestyle=3.0,color=200


	index=ipt(avedat.rho, 0.8)
	yr=[0,max(avedat.csden[index,16,*])]*1.1e-10
	plot,[0],[0],xr=xr,yr=yr,xtit='Time [sec]',ytit='n!u16+!n @ r/a=0.8 [10!u10!n m!u-3!n]',chars=1.8,/xsty,/ysty
	oplot,avedat.time,avedat.csden[index,16,*]*1.0e-10
	oplot,constdat.time,constdat.csden[index,16,*]*1.0e-10,linestyle=2.0,color=30
	oplot,evoldat.time,evoldat.csden[index,16,*]*1.0e-10,linestyle=3.0,color=200

	index=0
	yr=[0,max(avedat.csden[index,18,*])]*1.1e-10
	plot,[0],[0],xr=xr,yr=yr,xtit='Time [sec]',ytit='n!u18+!n @ r/a=0.0 [10!u10!n m!u-3!n]',chars=1.8,/xsty,/ysty
	oplot,avedat.time,avedat.csden[index,18,*]*1.0e-10
	oplot,constdat.time,constdat.csden[index,18,*]*1.0e-10,linestyle=2.0,color=30
	oplot,evoldat.time,evoldat.csden[index,18,*]*1.0e-10,linestyle=3.0,color=200

END

PRO run_evoldens_scan

	shot=1101209004
	z=18
	t1=1.0
	t2=1.2
	k=1.0
	dr=[0.5,0.1]
	savepath='/home/mlreinke/idl/impurities/strahl/evoldens_'+num2str(shot,1)+'.dat'
	RUN_CMOD_STRAHL,shot,z,t1,t2,k=k,dr=dr,data=avedat

	temp={temp:avedat.temp,rho:avedat.psin,time:avedat.time}
	dens={dens:avedat.dens,rho:avedat.psin,time:avedat.time}
	time=make(1.0,1.19,20)
	dt=time*0.0+0.01
	RUN_CMOD_STRAHL,shot,z,time,dt,k=k,dr=dr,data=constdat,temp=temp,dens=dens,/nogrid

	tmp=where(avedat.time GE 0.5*(t1+t2))
	FOR i=0,n(tmp) DO dens.dens[*,tmp[i]]*=0.8
	RUN_CMOD_STRAHL,shot,z,time,dt,k=k,dr=dr,data=evoldat,temp=temp,dens=dens,/nogrid

	save,avedat,constdat,evoldat,filename=savepath
END

PRO plot_evoldens_scan
	IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0
	shot=1101209004
	z=18
	savepath='/home/mlreinke/idl/impurities/strahl/evoldens_'+num2str(shot,1)+'.dat'
	restore,savepath

	IF keyword_set(ps) THEN BEGIN
		device, xsize=6.0,ysize=8.5,/inches
		ls=1.0
		!p.thick=3
		!p.charthick=2.0
	ENDIF ELSE openwin,0,xsize=600,ysize=850
	xr=[1.0,1.2]
	nrho=n(avedat.psin)+1
	!p.multi=[0,0,3]

	index=ipt(avedat.rho,1.02)
	yr=[0,2.5]
	plot,[0],[0],xr=xr,yr=yr,xtit='Time [sec]',ytit='n!le!n 10!u20!n] r/a=0.0',chars=1.8,/xsty,/ysty
	oplot,avedat.time,avedat.dens[0,*]/1.0e20
	oplot,constdat.time,constdat.dens[0,*]/1.0e20,linestyle=2.0,color=30
	oplot,evoldat.time,evoldat.dens[0,*]/1.0e20,linestyle=3.0,color=200
;	oplot,avedat.time,avedat.temp[index,*]/1.0e3*10
;	oplot,constdat.time,constdat.temp[index,*]/1.0e*10,linestyle=2.0,color=30
;	oplot,evoldat.time,evoldat.temp[index,*]/1.0e3*10,linestyle=3.0,color=200


	index=ipt(avedat.rho, 0.8)
	yr=[0,max(avedat.csden[index,16,*])]*1.1e-10
	plot,[0],[0],xr=xr,yr=yr,xtit='Time [sec]',ytit='n!u16+!n @ r/a=0.8 [10!u10!n m!u-3!n]',chars=1.8,/xsty,/ysty
	oplot,avedat.time,avedat.csden[index,16,*]*1.0e-10
	oplot,constdat.time,constdat.csden[index,16,*]*1.0e-10,linestyle=2.0,color=30
	oplot,evoldat.time,evoldat.csden[index,16,*]*1.0e-10,linestyle=3.0,color=200

	index=0
	yr=[0,max(avedat.csden[index,18,*])]*1.1e-10
	plot,[0],[0],xr=xr,yr=yr,xtit='Time [sec]',ytit='n!u18+!n @ r/a=0.0 [10!u10!n m!u-3!n]',chars=1.8,/xsty,/ysty
	oplot,avedat.time,avedat.csden[index,18,*]*1.0e-10
	oplot,constdat.time,constdat.csden[index,18,*]*1.0e-10,linestyle=2.0,color=30
	oplot,evoldat.time,evoldat.csden[index,18,*]*1.0e-10,linestyle=3.0,color=200
	!p.multi=0
	stop
END

PRO run_evolcoredens_scan

	shot=1101209004
	z=18
	t1=1.0
	t2=1.2
	k=1.0
	dr=[0.5,0.1]
	savepath='/home/mlreinke/idl/impurities/strahl/evolcoredens_'+num2str(shot,1)+'.dat'
	RUN_CMOD_STRAHL,shot,z,t1,t2,k=k,dr=dr,data=avedat

	temp={temp:avedat.temp,rho:avedat.psin,time:avedat.time}
	dens={dens:avedat.dens,rho:avedat.psin,time:avedat.time}
	time=make(1.0,1.19,20)
	dt=time*0.0+0.01
	RUN_CMOD_STRAHL,shot,z,time,dt,k=k,dr=dr,data=constdat,temp=temp,dens=dens,/nogrid

	tmp=where(avedat.time GE 0.5*(t1+t2))
	index=ipt(avedat.rho,0.95)
	FOR i=0,n(tmp) DO dens.dens[0:index,tmp[i]]*=0.8
	RUN_CMOD_STRAHL,shot,z,time,dt,k=k,dr=dr,data=evoldat,temp=temp,dens=dens,/nogrid

	save,avedat,constdat,evoldat,filename=savepath
END

PRO plot_evolcoredens_scan
	IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0
	shot=1101209004
	z=18
	savepath='/home/mlreinke/idl/impurities/strahl/evolcoredens_'+num2str(shot,1)+'.dat'
	restore,savepath

	IF keyword_set(ps) THEN BEGIN
		device, xsize=6.0,ysize=8.5,/inches
		ls=1.0
		!p.thick=3
		!p.charthick=2.0
	ENDIF ELSE openwin,0,xsize=600,ysize=850
	xr=[1.0,1.2]
	nrho=n(avedat.psin)+1
	!p.multi=[0,0,3]

	index=ipt(avedat.rho,1.02)
	yr=[0,2.5]
	plot,[0],[0],xr=xr,yr=yr,xtit='Time [sec]',ytit='n!le!n 10!u20!n] r/a=0.0',chars=1.8,/xsty,/ysty
	oplot,avedat.time,avedat.dens[0,*]/1.0e20
	oplot,constdat.time,constdat.dens[0,*]/1.0e20,linestyle=2.0,color=30
	oplot,evoldat.time,evoldat.dens[0,*]/1.0e20,linestyle=3.0,color=200
;	oplot,avedat.time,avedat.temp[index,*]/1.0e3*10
;	oplot,constdat.time,constdat.temp[index,*]/1.0e*10,linestyle=2.0,color=30
;	oplot,evoldat.time,evoldat.temp[index,*]/1.0e3*10,linestyle=3.0,color=200


	index=ipt(avedat.rho, 0.8)
	yr=[0,max(avedat.csden[index,16,*])]*1.1e-10
	plot,[0],[0],xr=xr,yr=yr,xtit='Time [sec]',ytit='n!u16+!n @ r/a=0.8 [10!u10!n m!u-3!n]',chars=1.8,/xsty,/ysty
	oplot,avedat.time,avedat.csden[index,16,*]*1.0e-10
	oplot,constdat.time,constdat.csden[index,16,*]*1.0e-10,linestyle=2.0,color=30
	oplot,evoldat.time,evoldat.csden[index,16,*]*1.0e-10,linestyle=3.0,color=200

	index=0
	yr=[0,max(avedat.csden[index,18,*])]*1.1e-10
	plot,[0],[0],xr=xr,yr=yr,xtit='Time [sec]',ytit='n!u18+!n @ r/a=0.0 [10!u10!n m!u-3!n]',chars=1.8,/xsty,/ysty
	oplot,avedat.time,avedat.csden[index,18,*]*1.0e-10
	oplot,constdat.time,constdat.csden[index,18,*]*1.0e-10,linestyle=2.0,color=30
	oplot,evoldat.time,evoldat.csden[index,18,*]*1.0e-10,linestyle=3.0,color=200
	!p.multi=0
END

PRO run_evolsource_scan
	shot=1101209004
	z=18
	t1=1.0
	t2=1.2
	k=1.0
	dr=[0.5,0.1]
	savepath='/home/mlreinke/idl/impurities/strahl/evolsource_'+num2str(shot,1)+'.dat'
	RUN_CMOD_STRAHL,shot,z,t1,t2,k=k,dr=dr,source=1.05,data=allon

	RUN_CMOD_STRAHL,shot,z,t1,t2,k=k,dr=dr,source=-1.05,data=inj,/nogrid

	wave=sin(2.0*!pi*(allon.time-1.0)/0.025)
	tmp=where(wave GT 0)
	source=allon.time*0.0
	source[tmp]=1.0e17
	source={source:source,time:allon.time}
	RUN_CMOD_STRAHL,shot,z,t1,t2,k=k,dr=dr,source=source,data=square,/nogrid

	save,allon,inj,square,filename=savepath
END

PRO plot_evolsource_scan
	IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0
	shot=1101209004
	z=18
	savepath='/home/mlreinke/idl/impurities/strahl/evolsource_'+num2str(shot,1)+'.dat'
	restore,savepath

	
	IF keyword_set(ps) THEN BEGIN
		device, xsize=6.0,ysize=8.5,/inches
		ls=1.0
		!p.thick=3
		!p.charthick=2.0
	ENDIF ELSE openwin,0,xsize=600,ysize=850
	xr=[1.0,1.2]
	nrho=n(allon.psin)+1
	!p.multi=[0,0,3]
	
	yr=[0,max(allon.csden[nrho-1,0,*])]*1.1e-8
	plot,[0],[0],xr=xr,yr=yr,xtit='Time [sec]',ytit='n!u0+!n @ EDGE [10!u8!n m!u-3!n]',chars=1.8,/xsty,/ysty
	oplot, square.time,square.csden[nrho-1,0,*]*0.98e-8,color=30
	oplot, allon.time,allon.csden[nrho-1,0,*]*1.02e-8,color=100
	oplot, inj.time,inj.csden[nrho-1,0,*]*1.0e-8,color=200

	index=ipt(allon.rho, 0.8)
	yr=[0,max(allon.csden[index,16,*])]*1.1e-10
	plot,[0],[0],xr=xr,yr=yr,xtit='Time [sec]',ytit='n!u16+!n @ r/a=0.8 [10!u10!n m!u-3!n]',chars=1.8,/xsty,/ysty
	oplot, square.time,square.csden[index,16,*]*1.0e-10,color=30
	oplot, allon.time,allon.csden[index,16,*]*1.0e-10,color=100
	oplot, inj.time,10.0*inj.csden[index,16,*]*1.0e-10,color=200
	xpt=inj.time[maxloc(reform(inj.csden[index,16,*]))]
	xyouts,xpt[0],0.85*yr[1],'x10',color=200,chars=1.0

	index=0
	yr=[0,max(allon.csden[index,18,*])]*1.1e-10
	plot,[0],[0],xr=xr,yr=yr,xtit='Time [sec]',ytit='n!u18+!n @ r/a=0.0 [10!u10!n m!u-3!n]',chars=1.8,/xsty,/ysty
	oplot, square.time,square.csden[index,18,*]*1.0e-10,color=30
	oplot, allon.time,allon.csden[index,18,*]*1.0e-10,color=100
	oplot, inj.time,10.0*inj.csden[index,18,*]*1.0e-10,color=200
	xpt=inj.time[maxloc(reform(inj.csden[index,18,*]))]
	xyouts,xpt[0]-0.01,0.25*yr[1],'x10',color=200,chars=1.0
	!p.multi=0
END


PRO run_evoldv_scan
	shot=1101209004
	z=18
	t1=1.0
	t2=1.2
	k=1.0
	dr=[0.5,0.1]
	savepath='/home/mlreinke/idl/impurities/strahl/evoldv_'+num2str(shot,1)+'.dat'
	RUN_CMOD_STRAHL,shot,z,t1,t2,k=k,dr=dr,data=data

	diff={diff:data.diff,rho:data.psin,time:data.time}
	conv={conv:data.conv,rho:data.psin,time:data.time}
	tmp=where(data.time GE 0.5*(t1+t2))
	FOR i=0,n(tmp) DO diff.diff[*,tmp[i]]*=0.75
	RUN_CMOD_STRAHL,shot,z,t1,t2,k=k,dr=dr,data=datad,diff=diff,conv=conv,/nogrid

	diff={diff:data.diff,rho:data.psin,time:data.time}
	conv={conv:data.conv,rho:data.psin,time:data.time}
	tmp=where(data.time GE 0.5*(t1+t2))
	FOR i=0,n(tmp) DO conv.conv[*,tmp[i]]=-1.0*sqrt(data.psin)*2.5
	RUN_CMOD_STRAHL,shot,z,t1,t2,k=k,dr=dr,data=datav,diff=diff,conv=conv,/nogrid
	save,data,datad,datav,filename=savepath

END

PRO plot_evoldv_scan
	IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0
	shot=1101209004
	z=18
	savepath='/home/mlreinke/idl/impurities/strahl/evoldv_'+num2str(shot,1)+'.dat'
	restore,savepath

	IF keyword_set(ps) THEN BEGIN
		device, xsize=6.0,ysize=8.5,/inches
		ls=1.0
		!p.thick=3
		!p.charthick=2.0
	ENDIF ELSE openwin,0,xsize=600,ysize=850
	xr=[1.0,1.2]
	nrho=n(data.psin)+1
	!p.multi=[0,0,3]
	
	index=ipt(data.rho,0.5)
	yr=[-1,0.6]
	plot,[0],[0],xr=xr,yr=yr,xtit='Time [sec]',ytit='D [m!u2!n/s], v [m/s] at r/a=0.5',chars=1.8,/xsty,/ysty
	oplot,data.time,data.diff[index,*]
	oplot,data.time,data.conv[index,*],linestyle=2.0
	oplot,datad.time,datad.diff[index,*],color=120
	oplot,datad.time,datad.conv[index,*],linestyle=2.0,color=120
	oplot,datav.time,datav.diff[index,*],color=150
	oplot,datav.time,datav.conv[index,*],linestyle=2.0,color=150


	index=ipt(data.rho, 0.8)
	yr=[0,max(datav.csden[index,16,*])]*1.1e-10
	plot,[0],[0],xr=xr,yr=yr,xtit='Time [sec]',ytit='n!u16+!n @ r/a=0.8 [10!u10!n m!u-3!n]',chars=1.8,/xsty,/ysty
	oplot,data.time,data.csden[index,16,*]*1.0e-10
	oplot,datad.time,datad.csden[index,16,*]*1.0e-10,color=120
	oplot,datav.time,datav.csden[index,16,*]*1.0e-10,color=150

	index=0
	yr=[0,max(datav.csden[index,18,*])]*1.1e-10
	plot,[0],[0],xr=xr,yr=yr,xtit='Time [sec]',ytit='n!u18+!n @ r/a=0.0 [10!u10!n m!u-3!n]',chars=1.8,/xsty,/ysty
	oplot,data.time,data.csden[index,18,*]*1.0e-10
	oplot,datad.time,datad.csden[index,18,*]*1.0e-10,color=120
	oplot,datav.time,datav.csden[index,18,*]*1.0e-10,color=150

END

FUNCTION strahl_grid_profile,rho,tau,iptr
	nrho=n(rho)+1
	ntau=n(tau)+1
	out=fltarr(nrho,ntau)
	drho=rho[1]-rho[0]
	dtau=tau[1]-tau[0]
	FOR i=0,nrho-1 DO BEGIN
		tmpa=where(iptr.rho GE rho[i]-drho/2.0 AND iptr.rho LT rho[i]+drho/2.0)
		xem=reform(sum_array(iptr.emiss[tmpa,*],/j))/(n(tmpa)+1.0)
		FOR j=0,ntau-1 DO BEGIN
			tmpb=where(iptr.time GE tau[j]-dtau/2.0 AND iptr.time LT tau[j]+dtau/2.0)
			out[i,j]=mean(xem[tmpb])
                ENDFOR
        ENDFOR
	RETURN,out
END

FUNCTION strahl_calc_emiss,x,a,shot=shot,z=z,t1=t1,t2=t2,source=source,y=y,sigy=sigy,nrho=nrho,ntau=ntau,plot=plot
	print,a
	dvrho=[0.0,make(0.1,0.7,5),0.8,0.9,1.0]
	diff={diff:10^[a[1],a[1],a[2],a[3],a[4],a[5],a[5],a[5],a[5]],rho:dvrho^2,time:t1}	;specify in psin ~ r/a^2
	;diff={diff:10^[a[4],a[4],a[4],a[4],a[4],a[4],a[4],a[4]],rho:dvrho^2,time:t1}
	v_over_d=[0.0,a[6],a[7],a[8],a[9],a[10],0.0,0.0,0.0]
	conv={conv:v_over_d*diff.diff,rho:dvrho^2,time:t1}
	k=1.0
	dr=[0.5,0.1]
	RUN_CMOD_STRAHL,shot,z,t1,t2,k=k,dr=dr,source=source,data=data,diff=diff,conv=conv,/nogrid,/nopp

	eptr=-1	
	sptr=0		
	nobr=1		
	rho=x[0:nrho-1]
	tau=x[nrho:nrho+ntau-1]

	CASE z OF

		42 : BEGIN
			tht=0
			line=8
			hirexsr_genrad_moly,shot,t1,t2,tht=tht,line=line,eptr=eptr,efit=efit,sptr=sptr,data=data,nobr=nobr
      			iptr=*sptr[0]
			f=strahl_grid_profile(rho,tau,iptr)*10^(a[0])*1.0e-10
                END

		18 : BEGIN
			tht=0
			line=[2,3]	
			hirexsr_genrad_argon,shot,t1,t2,tht=tht,line=line,eptr=eptr,efit=efit,sptr=sptr,data=data,nobr=nobr
			iptr=*sptr[0]
			xzem=strahl_grid_profile(rho,tau,iptr)*10^(a[0])*1.0e-10
			iptr=*sptr[1]
			xlyaem=strahl_grid_profile(rho,tau,iptr)*10^(a[0])*1.0e-10
			f=[[xzem],[xlyaem]]
                END

	ENDCASE

 
	IF keyword_set(plot) THEN BEGIN

	ENDIF
	RETURN,f
END

PRO run_hirexsr_ex,sload=sload,noini=noini,plot=plot,mo=mo
	IF NOT keyword_set(plot) THEN plot=0
	
	;default to using argon
	shot=1101209004
	z=18
	t1=1.0
	t2=1.08
	source=-1.01
    	
	IF keyword_set(mo) THEN BEGIN	;alternate demonstration using moly
		shot=1120913016
		z=42
		t1=1.0
		t2=1.08
		source=-1.01
	ENDIF

	k=1.0
	dr=[0.5,0.1]
	savepath='/home/mlreinke/idl/impurities/strahl/hirexsr_ex_'+num2str(shot,1)+'.dat'
	savepath2='/home/mlreinke/idl/impurities/strahl/hirexsr_fit_'+num2str(shot,1)+'.dat'
	IF NOT keyword_set(sload) THEN BEGIN
		xrho=make(0,1.0,100)
		diff=0.8*(1-0.9*exp(-xrho^2/(2.0*0.15^2)))
		conv=-2.0*xrho*exp(-(xrho-0.3)^2/(2.0*0.1^2))
		diff={diff:diff,rho:xrho^2,time:t1}	;make diff vs. psin~rho^2
		conv={conv:conv,rho:xrho^2,time:t1}	
		
		RUN_CMOD_STRAHL,shot,z,t1,t2,k=k,dr=dr,source=source,data=data,diff=diff,conv=conv

 		eptr=-1		;avoid loading experimental emission data
		sptr=0		;collect simulated emission data based on STRAHL output
		nobr=1		;calculate brightness profiles
		CASE z OF 
			42: BEGIN
				tht=0
				line=8		;calculate emissivity for line=8 (Ne-like Mo on high-Te configuration)
				hirexsr_genrad_moly,shot,t1,t2,tht=tht,line=line,eptr=eptr,efit=efit,sptr=sptr,data=data,nobr=nobr
                        END
			18 : BEGIN
				tht=0
				line=[2,3]	;calculate emissivity for both line=2 (He-like z) and line=3 (H-like Lya1)
				hirexsr_genrad_argon,shot,t1,t2,tht=tht,line=line,eptr=eptr,efit=efit,sptr=sptr,data=data,nobr=nobr
			END
		ENDCASE
	  	save,data,sptr,filename=savepath
	ENDIF ELSE restore,savepath

	dt=0.002			;assume 4 ms time resolution
	ntau=round((t2-t1)/dt)
	tau=make(t1,t2,ntau+1)
	tau=tau[0:ntau-1]+dt/2.0

	drho=0.04			;assume r/a of 0.04
	nrho=round(1.0/drho)
	rho=make(0,1.0,nrho+1)
	rho=rho[0:nrho-1]+drho/2.0
	x=[rho,tau]

	;make fake 'measurements'
	err=0.05
	CASE z OF 
		42 : BEGIN
			iptr=*sptr[0]
			y=strahl_grid_profile(rho,tau,iptr)*1.0e-10
			y=y*(1.0+err*randomn(seed,nrho*ntau,/normal))
			yerr=err*y+0.01
                END
		18 : BEGIN
			iptr=*sptr[0]
			xzem=strahl_grid_profile(rho,tau,iptr)*1.0e-10
			xzem=xzem*(1.0+err*randomn(seed,nrho*ntau,/normal))
			xzerr=err*xzem+0.01
			iptr=*sptr[1]
			xlyaem=strahl_grid_profile(rho,tau,iptr)*1.0e-10
			xlyaem=xlyaem*(1.0+err*randomn(seed,nrho*ntau,/normal))
			xlyaerr=err*xlyaem+0.01
			y=[[xzem],[xlyaem]]
			yerr=[[xzerr],[xlyaerr]]
               END
	ENDCASE

 	;make an initial run to allow /nogrid,/nopp to be specified in the loop
	IF NOT keyword_set(noini) THEN RUN_CMOD_STRAHL,shot,z,t1,t2,k=k,dr=dr,source=source

	stop
	estimate=[0.0,-0.2,-0.2,-0.2,-0.2,-0.2,0.0,0.0,0.0,0.0,0.0]	;[scale factor, x4 alog10(diff), x4 v/d]
	yest=strahl_calc_emiss(x,estimate,shot=shot,z=z,t1=t1,t2=t2,source=source,y=y,sigy=sigy,nrho=nrho,ntau=ntau,plot=plot)

	parinfo = replicate({value:0., fixed:0, limited:[0,0],limits:[0.,0], tied:'',step:[0.0],mpside:[0]}, n_elements(estimate))
	parinfo.value=estimate
	parinfo[0].step=0.5
	parinfo[0].fixed=1
	;parinfo[1:3].fixed=1			;use only on diffusivity
	parinfo[1:5].limited[*]=1
	parinfo[1:5].limits[*]=[-2.0,1.0]
	parinfo[6:10].limited[*]=1
	parinfo[6:10].limits[*]=[-10,10.0]
	parinfo[1:5].step=0.1
	parinfo[6:10].step=0.2
	parinfo[*].mpside=2
	functargs={shot:shot,z:z,t1:t1,t2:t2,source:source,y:y,sigy:yerr,nrho:nrho,ntau:ntau,plot:plot}

	coefs = mpfitfun('strahl_calc_emiss', x,y,yerr,estimate, parinfo=parinfo,status=status,niter=niter,functargs=functargs,xtol=1.0e-4)

	yfit=strahl_calc_emiss(x,coefs,shot=shot,z=z,t1=t1,t2=t2,source=source,y=y,sigy=sigy,nrho=nrho,ntau=ntau,plot=plot)

	save,x,y,yerr,yest,yfit,coefs,estimate,rho,tau,filename=savepath2
	stop
END

PRO run_hirexsr_fit,noini=noini,plot=plot,mo=mo
	IF NOT keyword_set(plot) THEN plot=0

	;demonstration using moly
	shot=1120913016
	z=42
	line=8
	tht=0
	t1=1.04
	t2=1.16
	tinj=1.052
	ta=[1.0,1.04]
	tb=[1.22,1.24]
	source=-1.0*tinj
	

	k=1.0
	dr=[0.5,0.1]
	savepath='/home/mlreinke/idl/impurities/strahl/hirexsr_mo_'+num2str(shot,1)+'.dat'

	mdsopen,'efit20',shot
	efit_rmid=mdsvalue('\efit_rmid')
	efit_time=mdsvalue('dim_of(\efit_rmid,0)')
	efit_psin=mdsvalue('dim_of(\efit_rmid,1)')
	mdsclose,'efit20',shot
	tmp=where(efit_time GE t1 and efit_time LE t2)
	efit_rmid=reform(sum_array(efit_rmid[tmp,*],/j))/(n(tmp)+1.0)
	Ro=efit_rmid[0]
	a=last(efit_rmid)-Ro

	hirexsr_load_profile,shot,line,prof,err,rho,tau,tht=tht
	tmp=where(tau NE -1)
	IF shot LT 1120913019 THEN nrad=(n(rho[*,0])+1)/2 ELSE nrad=n(rho[*,0])+1
	em=prof[0:nrad-1,tmp,0]
	sig=err[0:nrad-1,tmp,0]
	psin=rho[0:nrad-1,tmp]		;radius in psin
	time=tau[tmp]
	time+=(time[1]-time[0])*1.5		;appears to be off by a frame+ (need to investigate)
	bl=fltarr(nrad,n(tmp)+1)
	FOR i=0,nrad-1 DO BEGIN
		tmp=[where(time GE ta[0] AND time LE ta[1]),where(time GE tb[0] AND time LE tb[1])]
		xfit=time[tmp]
		yfit=em[i,tmp]
		coefs=poly_fit(xfit,yfit,1)
		bl[i,*]=time*coefs[1]+coefs[0]
	ENDFOR
	FOR i=0,n(time) DO em[*,i]-=bl[*,i]
	tmp=where(time GE t1 AND time LE t2)
	em=em[*,tmp]
	sig=sig[*,tmp]
	tau=time[tmp]

	;prepare experimental inputs
	rho=(interpol(efit_rmid,efit_psin,psin[*,ipt(time,tinj)])-Ro)/a
	nrho=n(rho)+1
	ntau=n(tau)+1
	x=[rho,tau]
	y=em*200		;scale to roughly match experiment
	yerr=sig*200

	stop

 	;make an initial run to allow /nogrid,/nopp to be specified in the loop
	IF NOT keyword_set(noini) THEN RUN_CMOD_STRAHL,shot,z,t1,t2,k=k,dr=dr,source=source

	estimate=[-0.5,-0.2,-0.2,-0.2,-0.2,-0.2,0.0,0.0,-10.0,-10.0,-10.0] ;[scale factor, x4 alog10(diff), x4 v/d]
	yest=strahl_calc_emiss(x,estimate,shot=shot,z=z,t1=t1,t2=t2,source=source,y=y,sigy=sigy,nrho=nrho,ntau=ntau,plot=plot)

	parinfo = replicate({value:0., fixed:0, limited:[0,0],limits:[0.,0], tied:'',step:[0.0],mpside:[0]}, n_elements(estimate))
	parinfo.value=estimate
	parinfo[0].step=0.2
	parinfo[0].fixed=0
	parinfo[0].limited[*]=1
	parinfo[0].limits[*]=[-1.0,1.0]
	parinfo[1:5].limited[*]=1
	parinfo[1:5].limits[*]=[-1.9,1.0]
	parinfo[6:10].limited[*]=1
	parinfo[6:10].limits[*]=[-50,10.0]
	parinfo[1:5].step=0.1
	parinfo[6:10].step=0.2
	parinfo[*].mpside=2
	functargs={shot:shot,z:z,t1:t1,t2:t2,source:source,y:y,sigy:yerr,nrho:nrho,ntau:ntau,plot:plot}

	coefs = mpfitfun('strahl_calc_emiss', x,y,yerr,estimate, parinfo=parinfo,status=status,niter=niter,functargs=functargs,xtol=1.0e-4)

	yfit=strahl_calc_emiss(x,coefs,shot=shot,z=z,t1=t1,t2=t2,source=source,y=y,sigy=sigy,nrho=nrho,ntau=ntau,plot=plot)

	save,x,y,yerr,yest,yfit,coefs,estimate,rho,tau,filename=savepath
	stop
END
