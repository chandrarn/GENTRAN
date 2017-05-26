PRO plot_bright,u
	index=u.index
	br=*u.dat.br[index]
	brchk=*u.dat.brchk[index]
	brchka=*u.dat.brchka[index]
	brdel=u.dat.brdel
	ch=u.dat.ch
	IF u.stat.ps THEN BEGIN
		xsize=5.5
		ysize=5.0*900/700.0
		ls=0.65
		col=u.stat.pscol
		d_old=!d
		device,xsize=xsize,ysize=ysize,/in
		tit=num2str(u.shot,1)+'   t='+num2str(u.stat.time,dp=3)+' [sec]'
		pdx=0.1
		pdy=0.02
	ENDIF ELSE BEGIN
		ls=1.0
		col=u.stat.col
		widget_control,u.id.draw1,get_value=draw_win
		window,0,xsize=u.plot.bsize[0],ysize=u.plot.bsize[1],/pixmap
		tit=''
		pdx=0
		pdy=0.0
	ENDELSE

	;plot brightness profile
	position=[0.12,0.37,0.975,0.98-pdy]
	IF u.plot.auto[0] THEN BEGIN
		u.plot.low[0]=min(ch)
		widget_control,u.id.chmin,set_value=num2str(int(u.plot.low[0]),1)
		u.plot.up[0]=max(ch)
		widget_control,u.id.chmax,set_value=num2str(int(u.plot.up[0]),1)
	ENDIF
	xr=[u.plot.low[0],u.plot.up[0]]
	IF u.plot.auto[1] THEN BEGIN
		u.plot.up[1]=max(br-brdel)*1.1
		widget_control,u.id.brmax,set_value=num2str(u.plot.up[1],dp=2)
		u.plot.low[1]=(min(brchk-brdel) < min(brchka-brdel) < 0)*1.01  
		widget_control,u.id.brmin,set_value=num2str(u.plot.low[1],dp=2)
	ENDIF
	yr=[u.plot.low[1],u.plot.up[1]]
	plot, ch,br-brdel,psym=8,xr=xr,/xsty,yr=yr,/ysty,ytit='Br [kW/m!u2!n]',chars=1.2*ls,$
		pos=position,xtickname=replicate(' ',5),symsize=1.0*ls,color=col[0],tit=tit
	oplot,ch,brchk-brdel,color=col[0]
	IF u.stat.m1stat THEN oplot,ch,brchka-brdel,color=col[5],linestyle=3.0
	oplot,[ch[0],last(ch)],[0,0],linestyle=2.0

	;plot legend
	oplot,[mean(xr)]-6,[0.8*yr[1]],psym=8,symsize=1.0*ls,color=col[0]
	xyouts,mean(xr),0.795*yr[1],'EXP',color=col[0],chars=ls*1.0
	oplot,mean(xr)+[-2,2]-6,[1,1]*0.75*yr[1],color=col[0]
	xyouts,mean(xr),0.745*yr[1],'m=0 CHK',color=col[0],chars=ls*1.0
	IF u.stat.m1stat THEN BEGIN
		oplot,mean(xr)+[-2,2]-6,[1,1]*0.7*yr[1],color=col[5],linestyle=3.0
		xyouts,mean(xr),0.695*yr[1],'m=1 CHK',color=col[5],chars=ls*1.0

	ENDIF

	position=[0.12,0.06,0.975,0.37]
	dbr=br-brchk
	IF u.stat.m1stat THEN dbra=br-brchka ELSE dbra=0
	IF u.plot.auto[2] THEN BEGIN
		u.plot.up[2]=max(dbr) > max(dbra)
		widget_control,u.id.dbrmax,set_value=num2str(u.plot.up[2],dp=2)
		u.plot.low[2]=min(dbr) < min(dbra)
		widget_control,u.id.dbrmin,set_value=num2str(u.plot.low[2],dp=2)
	ENDIF
	yr=[u.plot.low[2],u.plot.up[2]]
	plot, ch,dbr,psym=-8,xr=xr,/xsty,yr=yr*1.1,/ysty,xtit='CH #', ytit=n2g('Delta')+'Br [kW/m!u2!n]',$
		chars=1.2*ls,pos=position,/noerase,symsize=1.0*ls,color=col[0]
	IF u.stat.m1stat THEN oplot,ch,dbra,color=col[5],psym=-8,symsize=1.0*ls
	oplot,[ch[0],last(ch)],[0,0],linestyle=2.0
	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.bsize[0],u.plot.bsize[1],0,0,0]
	ENDIF ELSE device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm
END

PRO plot_emiss,u
	index=u.index
	em=*u.dat.em[index]
	emdel=u.dat.emdel
	emr=*u.dat.emr[index]
	emrdel=u.dat.emrdel
	emc1=*u.dat.emc1[index]
	emc1del=u.dat.emc1del
	emc2=*u.dat.emc2[index]
	emc2del=u.dat.emc2del
	ems1=*u.dat.ems1[index]
	ems1del=u.dat.ems1del
	IF u.plot.rho THEN BEGIN
		r=u.dat.rho 
		xtit='r/a'
	ENDIF ELSE BEGIN
		r=*u.dat.r[index]
		xtit='RMID [m]'
	ENDELSE
	IF u.stat.ps THEN BEGIN
		xsize=7.0
		ysize=6.0
		ls=0.85
		col=u.stat.pscol 
		tit=num2str(u.shot,1)+'   t='+num2str(u.stat.time,dp=3)+' [sec]'
		pdx=0.1
		d_old=!d
		device,xsize=xsize,ysize=ysize,/in
	ENDIF ELSE BEGIN
		ls=1.0
		col=u.stat.col
		widget_control,u.id.draw2,get_value=draw_win
		window,1,xsize=u.plot.esize[0],ysize=u.plot.esize[1],/pixmap
		tit=''
		pdx=0
	ENDELSE

	;plot emissivity profiles
	IF u.plot.auto[3] THEN BEGIN
		u.plot.low[3]=min(r)
		widget_control,u.id.rmin,set_value=num2str(u.plot.low[3],dp=2)
		u.plot.up[3]=max(r)
		widget_control,u.id.rmax,set_value=num2str(u.plot.up[3],dp=2)
	ENDIF
	xr=[u.plot.low[3],u.plot.up[3]]
	em0=emr-emrdel
	IF u.plot.auto[4] THEN BEGIN
		u.plot.up[4]=max(em-emdel) > max(emr-emrdel) > max(emc1-emc1del) > max(ems1-ems1del) > max(emc2-emc2del)
		widget_control,u.id.emmax,set_value=num2str(u.plot.up[4],dp=2)
		u.plot.low[4]=min(em-emdel) < min(emr-emrdel) < min(emc1-emc1del) < min(ems1-ems1del) < min(emc2-emc2del) < 0
		widget_control,u.id.emmin,set_value=num2str(u.plot.low[4],dp=2)
	ENDIF
	IF u.plot.ratio THEN BEGIN
		;yr=[-0.1,0.3] 
		yr=[u.plot.low[4],u.plot.up[4]]
		ytit='EM(m=X)/EMR'
	ENDIF ELSE BEGIN
		yr=[u.plot.low[4],u.plot.up[4]]
		ytit='Emissivity [kW/m!u3!n]'
	ENDELSE
	IF u.stat.m1stat THEN linestyle=3 ELSE linestyle=0.0
	plot, [0],[0],ytit=ytit,xtit=xtit,xr=xr,/xsty,yr=yr,/ysty,chars=ls,tit=tit
	IF u.plot.ratio THEN BEGIN
		IF u.stat.m1stat THEN oplot,r,(emc1-emc1del)/em0,color=col[2]
		IF u.stat.m2stat THEN oplot,r,(emc2-emc2del)/em0,color=col[3]
		IF u.stat.m1stat THEN oplot,r,(ems1-ems1del)/em0,color=col[4]	
	ENDIF ELSE BEGIN
		oplot,r,em-emdel,color=col[0],linestyle=linestyle
		IF u.stat.m1stat THEN oplot,r,emr-emrdel,color=col[1]
		IF u.stat.m1stat THEN oplot,r,emc1-emc1del,color=col[2]
		IF u.stat.m2stat THEN oplot,r,emc2-emc2del,color=col[3]
		IF u.stat.m1stat THEN oplot,r,ems1-ems1del,color=col[4]
	ENDELSE
	oplot,xr,[0,0],color=col[0],linestyle=2.0
	
	;plot legend
	dy=yr[1]-yr[0]
	yo=yr[0]
	dx=xr[1]-xr[0]
	xo=xr[0]-pdx*dx
	IF NOT u.plot.ratio THEN oplot,xo+dx*[0.8,0.85],yo+dy*0.9*[1.0,1.0],color=col[0],linestyle=linestyle
	IF NOT u.plot.ratio THEN xyouts,xo+dx*0.9,yo+dy*0.9,'RADIAL ONLY'
	IF u.stat.m1stat THEN BEGIN
		IF NOT u.plot.ratio THEN oplot,xo+dx*[0.8,0.85],yo+dy*0.85*[1,1],color=col[1]
		IF NOT u.plot.ratio THEN xyouts,xo+dx*0.9,yo+dy*0.85,'m=0 RAD',color=col[1]
		oplot,xo+dx*[0.8,0.85],yo+dy*0.8*[1,1],color=col[2]
		xyouts,xo+dx*0.9,yo+dy*0.8,'m=1 COS',color=col[2]
		oplot,xo+dx*[0.8,0.85],yo+dy*0.75*[1,1],color=col[4]
		xyouts,xo+dx*0.9,yo+dy*0.75,'m=1 SIN',color=col[4]
	ENDIF
	IF u.stat.m2stat THEN BEGIN
		oplot,xo+dx*[0.8,0.85],yo+dy*0.7*[1,1],color=col[3]
		xyouts,xo+dx*0.9,yo+dy*0.7,'m=2 COS',color=col[3]
	ENDIF
	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.esize[0],u.plot.esize[1],0,0,1]
	ENDIF ELSE device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm
END

PRO plot_cont,u
	index=u.index
	IF u.plot.rho THEN BEGIN
		r=u.dat.rho 
		ylab='r/a'
	ENDIF ELSE BEGIN
		r=*u.dat.r[index]
		ylab='RMID [m]'
	ENDELSE
	CASE u.plot.cont OF
		0 : BEGIN
			labels={ilab:n2g('theta')+' ['+n2g('pi')+']',jlab:ylab,klab:'Emissivity [kW/m!u3!n]',ctit:'2D Emissivity Profiles',$
				itit:'Radial Profiles at Fixed '+n2g('theta'),jtit:'Polar Profiles at Fixed r/a'}
			emr=*u.dat.emr[index]-u.dat.emrdel
			emc1=*u.dat.emc1[index]-u.dat.emc1del
			emc2=*u.dat.emc2[index]-u.dat.emc2del
			ems1=*u.dat.ems1[index]-u.dat.ems1del
			th=make(0,2.0*!pi,u.stat.nth)
			array=fltarr(u.stat.nth,u.stat.nr)
			FOR i=0,u.stat.nth-1 DO BEGIN
				FOR j=0,u.stat.nr-1 DO BEGIN
					IF u.stat.m2stat THEN array[i,j]=emr[j]+emc1[j]*cos(th[i])+emc2[j]*cos(2.0*th[i])+ems1[j]*sin(th[i]) ELSE $
						array[i,j]=emr[j]+emc1[j]*cos(th[i])+ems1[j]*sin(th[i])
				ENDFOR
			ENDFOR
			ivec=th/!pi
			jvec=r
			IF NOT u.plot.auto[3] THEN jr=[u.plot.low[3],u.plot.up[3]]
		END
		1 : BEGIN
			labels={jlab:'t [sec]',ilab:ylab,klab:'Emissivity [kW/m!u3!n]',ctit:'Radial Emissivity Profiles',$
				jtit:'Radial Profiles at Fixed Time',itit:'Temporal Profiles at Fixed Radius'}
			array=*u.dat.xemr
			;FOR i=0,u.stat.ntime-1 DO array[*,i]-=u.dat.emrdel
			jvec=u.dat.t
			ivec=r
			jr=[u.plot.low[5],u.plot.up[5]]
			ir=[u.plot.low[3],u.plot.up[3]]
			
		END
		2 : BEGIN
			labels={jlab:'t [sec]',ilab:ylab,klab:'Emissivity [kW/m!u3!n]',ctit:'m=1Cosine Emissivity Profiles',$
				jtit:'Radial Profiles at Fixed Time',itit:'Temporal Profiles at Fixed Radius'}
			array=*u.dat.xemc1
			IF u.plot.ratio THEN BEGIN
				emr=*u.dat.xemr
				array=array/emr
			ENDIF
			;FOR i=0,u.stat.ntime-1 DO array-=u.dat.emc1del
			jvec=u.dat.t
			ivec=r
			jr=[u.plot.low[5],u.plot.up[5]]
			ir=[u.plot.low[3],u.plot.up[3]]
		END
		3 : BEGIN
			labels={jlab:'t [sec]',ilab:ylab,klab:'Emissivity [kW/m!u3!n]',ctit:'m=2 Cosine Emissivity Profiles',$
				jtit:'Radial Profiles at Fixed Time',itit:'Temporal Profiles at Fixed Radius'}
			array=*u.dat.xemc2
			;FOR i=0,u.stat.ntime-1 DO array[*,i]-=u.dat.emc2del
			jvec=u.dat.t
			ivec=r
			jr=[u.plot.low[5],u.plot.up[5]]
			ir=[u.plot.low[3],u.plot.up[3]]
		END
		4 : BEGIN
			labels={jlab:'t [sec]',ilab:ylab,klab:'Emissivity [kW/m!u3!n]',ctit:'m=1Sine Emissivity Profiles',$
				jtit:'Radial Profiles at Fixed Time',itit:'Temporal Profiles at Fixed Radius'}
			array=*u.dat.xems1
			;FOR i=0,u.stat.ntime-1 DO array[*,i]-=u.dat.ems1del
			jvec=u.dat.t
			ivec=r
			jr=[u.plot.low[5],u.plot.up[5]]
			ir=[u.plot.low[3],u.plot.up[3]]
		END
	ENDCASE

	widget_control,u.id.draw3,get_value=dw3
	widget_control,u.id.draw4,get_value=dw4
	widget_control,u.id.draw5,get_value=dw5
	io=*u.plot.io
	jo=*u.plot.jo
	IF NOT u.plot.auto[4] THEN BEGIN
		minpt=u.plot.low[4]
		maxpt=u.plot.up[4]
		IF u.plot.ratio THEN BEGIN
			minpt=-1.0
			maxpt=1.0
		ENDIF
	ENDIF ELSE BEGIN
		minpt=0.0
		maxpt=0.0
	ENDELSE
	genplt,array,ivec,jvec,cct=39,ncntrs=60,labels=labels,pixmap=u.plot.csize,set=set,jr=jr,ir=ir,io=io,jo=jo,minpt=minpt,maxpt=maxpt
	IF n(*u.plot.io) EQ 0 THEN BEGIN
		*u.plot.io=[io[0],io]
		iostr=num2str(io[0],dp=2)
		FOR i=1,n(io) DO iostr=iostr+', '+num2str(io[i],dp=2)
		u.plot.iostr=iostr
		widget_control,u.id.io_val,set_value=u.plot.iostr
	ENDIF
	IF n(*u.plot.jo) EQ 0 THEN BEGIN
		*u.plot.jo=[jo[0],jo]
		jostr=num2str(jo[0],dp=2)
		FOR i=1,n(jo) DO jostr=jostr+', '+num2str(jo[i],dp=2)
		u.plot.jostr=jostr
		widget_control,u.id.jo_val,set_value=u.plot.jostr
	ENDIF
	IF NOT u.stat.ps THEN BEGIN
		wset,dw3
		device,copy=[0,0,u.plot.csize[0],u.plot.csize[1],0,0,set[0]]
		wset,dw4
		device,copy=[0,0,u.plot.csize[2],u.plot.csize[3],0,0,set[1]]
		wset,dw5
		device,copy=[0,0,u.plot.csize[4],u.plot.csize[5],0,0,set[2]]
	ENDIF
END

PRO calc_delta,u
	;calculate the time-avereged profile to subtract off (good for showing injection evolution)
	tmp=where(u.dat.t GE u.plot.t1[0] AND u.dat.t LE u.plot.t1[1])
	IF tmp[0] NE -1 THEN BEGIN
		npts=n(tmp)+1.0
		brdel=u.dat.brdel*0.0
		emdel=u.dat.emdel*0.0
		emrdel=u.dat.emrdel*0.0
		emc1del=u.dat.emc1del*0.0
		emc2del=u.dat.emc2del*0.0
		ems1del=u.dat.ems1del*0.0
		FOR i=0,npts-1 DO BEGIN
			ibr=*u.dat.br[tmp[i]]
			brdel+=ibr
			iem=*u.dat.em[tmp[i]]
			emdel+=iem
			iemr=*u.dat.emr[tmp[i]]
			emrdel+=iemr
			iemc1=*u.dat.emc1[tmp[i]]
			emc1del+=iemc1
			iemc2=*u.dat.emc2[tmp[i]]
			emc2del+=iemc2
			iems1=*u.dat.ems1[tmp[i]]
			ems1del+=iems1
		ENDFOR
		u.dat.brdel=brdel/npts
		u.dat.emdel=emdel/npts
		u.dat.emrdel=emrdel/npts
		u.dat.emc1del=emc1del/npts
		u.dat.emc2del=emc2del/npts
		u.dat.ems1del=ems1del/npts
	ENDIF
END

PRO load_data,u

	;load data from the tree
	shot=u.shot
	mdsopen,'xtomo',shot
	xem=mdsvalue('_sig=\XTOMO::TOP.RESULTS.CORE:EMISS',/quiet,status=cstatus)
	rho=mdsvalue('dim_of(_sig,0)',/quiet)
	t=mdsvalue('dim_of(_sig,1)',/quiet)
	xemr=mdsvalue('dim_of(_sig,2)',/quiet,status=m1status)
	xemc1=mdsvalue('dim_of(_sig,3)',/quiet,status=m1status)
	xemc2=mdsvalue('dim_of(_sig,4)',/quiet,status=m2status)
	xems1=mdsvalue('dim_of(_sig,5)',/quiet,status=m1status)
	etree=mdsvalue('dim_of(_sig,6)',/quiet,status=tree_status)
	IF NOT tree_status THEN etree='ANALYSIS'
	xr=mdsvalue('\XTOMO::TOP.RESULTS.CORE:RMAJ',/quiet)
	xbr=mdsvalue('_sig=\XTOMO::TOP.RESULTS.CORE:BRIGHT',/quiet,status=status)
	ch=mdsvalue('dim_of(_sig,0)',/quiet)
	xbrchk=mdsvalue('_sig=\XTOMO::TOP.RESULTS.CORE:BRCHK',/quiet)
	xbrchka=mdsvalue('dim_of(_sig,2)',/quiet)
	mdsclose,'xtomo',shot

	IF NOT cstatus THEN widget_control,u.id.message,set_value='SHOT: '+num2str(shot,1)+' NO CORE EMISS DATA',/app
	IF NOT m1status THEN widget_control,u.id.message,set_value='SHOT: '+num2str(shot,1)+' NO m=1 ASYM EMISS DATA',/app
	IF NOT m2status THEN widget_control,u.id.message,set_value='SHOT: '+num2str(shot,1)+' NO m=2 ASYM EMISS DATA',/app
	IF NOT cstatus THEN RETURN
	IF NOT m1status THEN u.stat.m1stat=0 ELSE u.stat.m1stat=1
	IF NOT m2status THEN u.stat.m2stat=0 ELSE u.stat.m2stat=1

	widget_control,u.id.message,set_value='SHOT: '+num2str(shot,1)+' LOADED - '+etree,/app
	
	IF NOT u.stat.m1stat THEN BEGIN
		reset_cont_button,u
		widget_control,u.id.emr,set_button=1
		xemr=xem
		xemc1=0.0
		xems1=0.0
		u.plot.cont=1
	ENDIF
	IF u.stat.m2stat THEN BEGIN
		IF max(xemc2) EQ 0 AND min(xemc2) EQ 0 THEN u.stat.m2stat=0
	ENDIF 
	IF NOT u.stat.m2stat THEN xemc2=0.0

	;initialize delta arrays
	nrho=n(rho)+1
	brdel=fltarr(n(xbr[*,0])+1)
	emdel=fltarr(nrho)
	emrdel=fltarr(nrho)
	emc1del=fltarr(nrho)
	emc2del=fltarr(nrho)
	ems1del=fltarr(nrho)

	;reorganize into pointers
	u.stat.ntime=n(t)+1
	em=ptrarr(u.stat.ntime,/allocate)
	emr=ptrarr(u.stat.ntime,/allocate)
	ems1=ptrarr(u.stat.ntime,/allocate)
	emc1=ptrarr(u.stat.ntime,/allocate)
	emc2=ptrarr(u.stat.ntime,/allocate)
	br=ptrarr(u.stat.ntime,/allocate)
	brchk=ptrarr(u.stat.ntime,/allocate)
	brchka=ptrarr(u.stat.ntime,/allocate)
	r=ptrarr(u.stat.ntime,/allocate)
	FOR i=0L,u.stat.ntime-1 DO BEGIN
		*em[i]=xem[*,i]*1.0e-3
		*br[i]=xbr[*,i]*1.0e-3
		*brchk[i]=xbrchk[*,i]*1.0e-3
		*r[i]=xr[*,i]

		IF u.stat.m1stat THEN BEGIN
			*emr[i]=xemr[*,i]*1.0e-3
			*emc1[i]=xemc1[*,i]*1.0e-3
			*ems1[i]=xems1[*,i]*1.0e-3
			*brchka[i]=xbrchka[*,i]*1.0e-3
		ENDIF ELSE BEGIN
			*emr[i]=0.0
			*emc1[i]=0.0
			*ems1[i]=0.0
			*brchka[i]=0.0
		ENDELSE
		IF u.stat.m2stat THEN *emc2[i]=xemc2[*,i]*1.0e-3 ELSE *emc2[i]=0.0
	ENDFOR

	u.plot.low[5]=min(t)
	widget_control,u.id.tmin,set_value=num2str(u.plot.low[5],dp=4)
	u.plot.up[5]=max(t)
	widget_control,u.id.tmax,set_value=num2str(u.plot.up[5],dp=4)
	u.index=ipt(t,u.stat.time)
	widget_control,u.id.t_text,set_value=num2str(u.stat.time,dp=5)
	dat={em:em,rho:rho,t:t,emr:emr,emc1:emc1,emc2:emc2,ems1:ems1,r:r,br:br,ch:ch,brchk:brchk,brchka:brchka,xemr:ptr_new(xemr*1.0e-3,/allocate_heap,/no_copy),$
		xemc1:ptr_new(xemc1*1.0e-3,/allocate_heap,/no_copy),xemc2:ptr_new(xemc2*1.0e-3,/allocate_heap,/no_copy),xems1:ptr_new(xems1*1.0e-3,/allocate_heap,/no_copy),$
		brdel:brdel,emdel:emdel,emrdel:emrdel,emc1del:emc1del,emc2del:emc2del,ems1del:ems1del}
	u.stat.dat=1
	u.stat.nr=nrho
	u={id:u.id,shot:u.shot,index:u.index,ch:u.ch,stat:u.stat,plot:u.plot,dat:dat}
	widget_control,u.id.base, set_uvalue=u	
END


PRO cleanup_xtomo_ptr,u
	IF u.stat.dat THEN BEGIN
		heap_free,u.dat.em
		heap_free,u.dat.emr
		heap_free,u.dat.ems1
		heap_free,u.dat.emc1
		heap_free,u.dat.emc2
		heap_free,u.dat.r
		heap_free,u.dat.br
		heap_free,u.dat.brchk
		heap_free,u.dat.brchka
		heap_free,u.dat.xemr
		heap_free,u.dat.xemc1
		heap_free,u.dat.xems1
		heap_free,u.dat.xemc2
	ENDIF
END

PRO reset_cont_button,u
	id=[u.id.em2d,u.id.emr,u.id.emc1,u.id.emc2,u.id.ems1]
	widget_control,id[u.plot.cont],set_button=0
END

PRO w_xtomo_event,event
	widget_control,event.top,get_uvalue=u
	id = u.id
	tag = tag_names(event,/st)
	button=' '
	idtags=tag_names(id)
	FOR i=0,n(idtags) DO IF id.(i) EQ event.id THEN ename=idtags[i]
	CASE tag OF
		"WIDGET_BASE" : BEGIN

		END
		"WIDGET_BUTTON": BEGIN
			widget_control,event.id,get_value=button,get_uvalue=uvalue
			IF NOT keyword_set(uvalue) THEN uvalue='none'
			CASE ename OF
				"LOAD" : BEGIN
					WIDGET_CONTROL,/hourglass
					IF u.stat.dat THEN BEGIN
						cleanup_xtomo_ptr,u
						heap_gc
					ENDIF
					load_data,u
					plot_bright,u
					plot_emiss,u
					plot_cont,u
				END
				"QUIT" : BEGIN
					widget_control,event.top,/destroy
					cleanup_xtomo_ptr,u
					heap_free,u.plot.io
					heap_free,u.plot.jo
					heap_gc
					!except=1
				END
				"STOP" : BEGIN
					stop
				END
				"PRINT" : BEGIN
					u.stat.ps=1
					psplot
					plot_bright,u
					plot_emiss,u
					psc
					xwplot
					set_plot,'x'
					u.stat.ps=0
				END
				"EM2D" : BEGIN
					reset_cont_button,u
					widget_control,u.id.em2d,set_button=1
					u.plot.cont=0
					IF u.stat.m1stat THEN plot_cont,u
				END
				"EMR" : BEGIN
					WIDGET_CONTROL,/hourglass
					reset_cont_button,u
					widget_control,u.id.emr,set_button=1
					u.plot.cont=1
					IF u.stat.m1stat THEN plot_cont,u
				END
				"EMC2" : BEGIN
					WIDGET_CONTROL,/hourglass
					reset_cont_button,u
					widget_control,u.id.emc2,set_button=1
					u.plot.cont=3
					IF u.stat.m2stat THEN plot_cont,u
				END
				"EMC1" : BEGIN
					WIDGET_CONTROL,/hourglass
					reset_cont_button,u
					widget_control,u.id.emc1,set_button=1
					u.plot.cont=2
					IF u.stat.m1stat THEN plot_cont,u
				END
				"EMS1" : BEGIN
					WIDGET_CONTROL,/hourglass
					reset_cont_button,u
					widget_control,u.id.ems1,set_button=1
					u.plot.cont=4
					IF u.stat.m1stat THEN plot_cont,u
				END
				"CHAUTO" : IF event.select EQ 1 THEN BEGIN
						u.plot.auto[0]=1 
						plot_bright,u
					   ENDIF ELSE u.plot.auto[0]=0
				"BRAUTO" : IF event.select EQ 1 THEN BEGIN
						u.plot.auto[1]=1 
						plot_bright,u
					   ENDIF ELSE u.plot.auto[1]=0
				"DBRAUTO" : IF event.select EQ 1 THEN BEGIN
						u.plot.auto[2]=1 
						plot_bright,u
					    ENDIF ELSE u.plot.auto[2]=0
				"RAUTO" : IF event.select EQ 1 THEN BEGIN
				     		u.plot.auto[3]=1 
						plot_emiss,u
						plot_cont,u
					  ENDIF ELSE u.plot.auto[3]=0
				"EMAUTO" : IF event.select EQ 1 THEN BEGIN
						u.plot.auto[4]=1 
						plot_emiss,u
						plot_cont,u
					   ENDIF ELSE u.plot.auto[4]=0
				"TAUTO" : IF event.select EQ 1 THEN BEGIN
						u.plot.auto[5]=1 
						IF u.plot.cont NE 0 THEN plot_cont,u
					  ENDIF ELSE u.plot.auto[5]=0

				"XRMAJ" : BEGIN
					IF event.select EQ 1 THEN BEGIN
						widget_control,u.id.xrho,set_button=0
						u.plot.rho=0	
						plot_emiss,u
						*u.plot.io=0
						*u.plot.jo=0
						plot_cont,u
					ENDIF ELSE widget_control,u.id.xrmaj,set_button=1
				END
				"XRHO" : BEGIN
					IF event.select EQ 1 THEN BEGIN
						widget_control,u.id.xrmaj,set_button=0
						u.plot.rho=1	
						plot_emiss,u
						*u.plot.io=0
						*u.plot.jo=0
						plot_cont,u
					ENDIF ELSE widget_control,u.id.xrho,set_button=1
				END
				"PDEL" : BEGIN
					IF event.select EQ 1 THEN BEGIN
						calc_delta,u
					ENDIF ELSE BEGIN
						u.dat.brdel*=0.0
						u.dat.emdel*=0.0
						u.dat.emrdel*=0.0
						u.dat.emc1del*=0.0
						u.dat.emc2del*=0.0
						u.dat.ems1del*=0.0
					ENDELSE	
					plot_bright,u
					plot_emiss,u
					plot_cont,u	
				END
				"PRAT" : BEGIN
					IF event.select EQ 1 THEN u.plot.ratio=1 ELSE u.plot.ratio=0
					plot_emiss,u
				END				
				ELSE:
			ENDCASE
		END
		"WIDGET_SLIDER": BEGIN
			widget_control,event.id,get_value=slider
			IF NOT keyword_set(uvalue) THEN uvalue='none'
			CASE ename OF 
				'T_SLIDER' : BEGIN
					IF u.stat.dat THEN BEGIN
						index=ipt(u.dat.t,slider/5.0e4)
						tmp=where(u.dat.t GT 0)
						IF slider/5.0e4 GE max(u.dat.t[tmp]) THEN index=n(tmp)
						IF slider/5.0e4 LE min(u.dat.t[tmp]) THEN index=0
						IF u.stat.time NE u.dat.t[index] THEN BEGIN
							u.stat.time=u.dat.t[index]
							u.index=index
							plot_bright,u
							plot_emiss,u
							IF u.stat.m1stat AND u.plot.cont EQ 0 THEN plot_cont,u
							widget_control,u.id.t_text,set_value=num2str(u.stat.time,dp=5)
						ENDIF
					ENDIF
				END
				ELSE:
			ENDCASE
		END
   		"WIDGET_TEXT_CH": BEGIN
			widget_control,event.id,get_value=data
			CASE event.id OF
			u.id.shotid : u.shot=long(data)
			u.id.rmin : BEGIN
				u.plot.low[3]=float(data)
				plot_emiss,u
				plot_cont,u
			END		
			u.id.rmax : BEGIN
				u.plot.up[3]=float(data)
				plot_emiss,u
				plot_cont,u
			END
			u.id.tmin : BEGIN
				u.plot.low[5]=float(data)
				IF u.plot.cont NE 0 THEN plot_cont,u
			END		
			u.id.tmax : BEGIN
				u.plot.up[5]=float(data)
				IF u.plot.cont NE 0 THEN plot_cont,u
			END
			u.id.emmin : BEGIN
				u.plot.low[4]=float(data)
				plot_emiss,u
				IF u.plot.cont NE 0 THEN plot_cont,u
			END		
			u.id.emmax : BEGIN
				u.plot.up[4]=float(data)
				plot_emiss,u
				IF u.plot.cont NE 0 THEN plot_cont,u
			END

			u.id.t_text : BEGIN
				u.stat.time=float(data)
				u.index=ipt(u.dat.t,u.stat.time)
				widget_control,u.id.t_slider,set_value=u.stat.time*5.0e4
				plot_bright,u
				plot_emiss,u
				IF u.stat.m1stat AND u.plot.cont EQ 0 THEN plot_cont,u
			END
			u.id.io_val : BEGIN
				widget_control,u.id.io_val,get_value=iostr
				io=float(strsplit(iostr,',',/extract))
				IF u.stat.dat THEN BEGIN	
					*u.plot.io=io
					plot_cont,u
				ENDIF

			END
			u.id.jo_val : BEGIN
				widget_control,u.id.jo_val,get_value=jostr
				jo=float(strsplit(jostr,',',/extract))
				IF u.stat.dat THEN BEGIN	
					*u.plot.jo=jo
					plot_cont,u
				ENDIF

			END
			u.id.t1_val : BEGIN
				widget_control,u.id.t1_val,get_value=t1str
				t1=float(strsplit(t1str,',',/extract))
				IF u.stat.dat THEN BEGIN	
					u.plot.t1=t1
				ENDIF
			END
			u.id.t2_val : BEGIN
				widget_control,u.id.t2_val,get_value=t2str
				t2=float(strsplit(t2str,',',/extract))
				IF u.stat.dat THEN BEGIN	
					u.plot.t2=t2
				ENDIF
			END
			ELSE: 
			ENDCASE
		END
		ELSE:
	ENDCASE
	IF button NE 'QUIT' THEN widget_control,event.top,set_uvalue=u		
END

;+
;NAME:
;	W_XTOMO
;
;PURPOSE:
;	This procedure launches the W_XTOMO widget which is used to
;	visualize the results of the GENPOS inversions of the x-ray
;	tomography data.  Note this is different then the tomographic
;	analysis performed by Granetz.
;
;CALLING SEQUENCE:
;	@run_wxtomo.pro 	will compile the proper files, set the color
;				tables and launch the widget
;PROCEUDRE:
;	See the C-Mod Wiki page on the W_XTOMO widget for instructions
;	(listed under Shot Analysis Tools)
;
;RESTRICTIONS:
;	Due to the limited namespace, many of the GENIE widgets have
;	overlapping function/procedure names so these widgets may not run
;	well in contiguous IDL sessions.
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 2010
;
;-

PRO w_xtomo,shot=shot,time=time
	IF NOT keyword_set(shot) THEN shot=1110105031
	IF NOT keyword_set(time) THEN time=1.0
	user=logname()
	loadct,12,/silent
	mdim=get_screen_size()
	IF mdim[0] LT 1600 AND mdim[1] LT 1000 THEN base=widget_base(title='XTOMO ASYMMETRY ANALYSIS',/row,tlb_size_events=1,/scroll,$
		x_scroll_size=mdim[0]*0.95,y_scroll_size=mdim[1]*0.85) ELSE base=widget_base(title='XTOMO ASYMMETRY ANALYSIS',/row,tlb_size_events=1)
	A=widget_base(base,/column)
	B=widget_base(base,/column)

	xsA=600
	dum=widget_label(A,value='BRIGHTNESS PROFILES')
	A1=widget_base(A,frame=2)
	draw1=widget_draw(A1,xsize=xsA,ysize=570)

	AX=widget_base(A,/row)
	A2=widget_base(AX,/column,xsize=xsA*0.55,ysize=360,/frame)
	dum=widget_label(A2,value='SETUP')
	A2p1=widget_base(A2,/row)
	dum = widget_label(A2p1,value='SHOT: ')
	shotid = widget_text(A2p1,xsize=10,ysize=1,/edit)
	dum = widget_label(A2p1,value='')
	load= widget_button(A2p1,value='LOAD')
	dum = widget_label(A2p1,value='')
	quit= widget_button(A2p1,value='QUIT')
	dum = widget_label(A2p1,value='')
	print= widget_button(A2p1,value='PRINT')
	dum = widget_label(A2p1,value='')
	stop= widget_button(A2p1,value='STOP')
	A2p2=widget_base(A2,/row)
	message = widget_text(A2p2,xsize=40,ysize=5,/scroll)

	dum=widget_label(A2,value='CONTOUR OPTIONS')
	A2p3=widget_base(A2,/row)
	dum=widget_label(A2p3,value='CONTOUR: ')
	A2p3x=widget_base(A2p3,/row, /nonexcl)
	em2d=widget_button(A2p3x,value='EM2D')
	emr=widget_button(A2p3x,value='EMR')
	emc1=widget_button(A2p3x,value='EMC1')
	emc2=widget_button(A2p3x,value='EMC2')
	ems1=widget_button(A2p3x,value='EMS1')
	A2p4=widget_base(A2,/row)
	dum=widget_label(A2p4,value='io: ')
	io_val=widget_text(A2p4,xsize=35,ysize=1,/edit)
	A2p5=widget_base(A2,/row)
	dum=widget_label(A2p5,value='jo: ')
	jo_val=widget_text(A2p5,xsize=35,ysize=1,/edit)
	A2p6=widget_base(A2,/row)
	dum=widget_label(A2p6,value='t1: ')
	t1_val=widget_text(A2p6,xsize=10,ysize=1,/edit)
	dum=widget_label(A2p6,value=' t2: ')
	t2_val=widget_text(A2p6,xsize=10,ysize=1,/edit)


	
	;plotting options
	Px=widget_base(AX,/column,xsize=xsA*0.4,ysize=320,/frame)
	dum=widget_label(Px, value='PLOTTING OPTIONS')
	P1=widget_base(Px, /row)
	P1a=widget_base(P1,/row,/nonexcl)
	chauto=widget_button(P1a,value='AUTO ')
	chmin=widget_text(P1,xsize=5,ysize=1,/edit)
	dum=widget_label(P1,value=' < CH < ')
	chmax=widget_text(P1,xsize=5,ysize=1,/edit)
	P2=widget_base(Px, /row)
	P2a=widget_base(P2,/row,/nonexcl)
	brauto=widget_button(P2a,value='AUTO ')
	brmin=widget_text(P2,xsize=5,ysize=1,/edit)
	dum=widget_label(P2,value=' < BR < ')
	brmax=widget_text(P2,xsize=5,ysize=1,/edit)
	P3=widget_base(Px, /row)
	P3a=widget_base(P3,/row,/nonexcl)
	dbrauto=widget_button(P3a,value='AUTO ')
	dbrmin=widget_text(P3,xsize=5,ysize=1,/edit)
	dum=widget_label(P3,value=' < DBR < ')
	dbrmax=widget_text(P3,xsize=5,ysize=1,/edit)
	P4=widget_base(Px, /row)
	P4a=widget_base(P4,/row,/nonexcl)
	rauto=widget_button(P4a,value='AUTO ')
	rmin=widget_text(P4,xsize=5,ysize=1,/edit)
	dum=widget_label(P4,value=' < r < ')
	rmax=widget_text(P4,xsize=5,ysize=1,/edit)
	P5=widget_base(Px, /row)
	P5a=widget_base(P5,/row,/nonexcl)
	emauto=widget_button(P5a,value='AUTO ')
	emmin=widget_text(P5,xsize=5,ysize=1,/edit)
	dum=widget_label(P5,value=' < EM < ')
	emmax=widget_text(P5,xsize=5,ysize=1,/edit)
	P6=widget_base(Px, /row)
	P6a=widget_base(P6,/row,/nonexcl)
	tauto=widget_button(P6a,value='AUTO ')
	tmin=widget_text(P6,xsize=5,ysize=1,/edit)
	dum=widget_label(P6,value=' < t < ')
	tmax=widget_text(P6,xsize=5,ysize=1,/edit)
	P7=widget_base(Px, /row)
	dum=widget_label(P7,value='RADIAL:  ')
	P7a=widget_base(P7,/row,/nonexcl)
	xrmaj=widget_button(P7a,value='RMAJ ')
	xrho=widget_button(P7a,value='r/a ')
	P8=widget_base(Px, /row)
	dum=widget_label(P8,value='DELTA:  ')
	P8a=widget_base(P8,/row,/nonexcl)
	pdel=widget_button(P8a,value=' ')
	dum=widget_label(P8,value='RATIO:  ')
	P8b=widget_base(P8,/row,/nonexcl)
	prat=widget_button(P8b,value=' ')


	xsB=900
	dum=widget_label(B,value='RADIAL PROFILES')
	B1=widget_base(B,frame=2)
	draw2=widget_draw(B1,xsize=xsB+50,ysize=550)
	Bt=widget_base(B,/row)
	dum=widget_label(Bt, value='TIME: ')
	t_slider=widget_slider(Bt,xsize=xsB-75,min=0,max=100000,value=50000,/drag,/suppress)
	t_text=widget_text(Bt,xsize=8,ysize=1,/edit)
	B2=widget_base(B,/row)
	draw3=widget_draw(B2,xsize=550,ysize=350)
	B2a=widget_base(B2,/column)
	draw4=widget_draw(B2a,xsize=xsb-510,ysize=175)
	draw5=widget_draw(B2a,xsize=xsb-510,ysize=175)
	
	id={base:base,draw1:draw1,draw2:draw2,draw3:draw3,draw4:draw4,draw5:draw5,t_slider:t_slider,t_text:t_text,$
		shotid:shotid,load:load,quit:quit,print:print,stop:stop,message:message,$
		em2d:em2d,emr:emr,emc1:emc1,emc2:emc2,ems1:ems1,io_val:io_val,jo_val:jo_val,t1_val:t1_val,t2_val:t2_val,$
		chauto:chauto,brauto:brauto,dbrauto:dbrauto,rauto:rauto,emauto:emauto,tauto:tauto,$
		chmax:chmax,brmax:brmax,dbrmax:dbrmax,rmax:rmax,emmax:emmax,tmax:tmax,$
		chmin:chmin,brmin:brmin,dbrmin:dbrmin,rmin:rmin,emmin:emmin,tmin:tmin,$
		xrmaj:xrmaj,xrho:xrho,pdel:pdel,prat:prat}


	;plots ch,br,dbr,r,em,t
	auto=[1,1,1,1,1,1]
	up=float([0,0,0,0,0,0])
	low=float([0,0,0,0,0,0])
	csize=[550,400,390,175,390,175]
	bsize=[600,570]
	esize=[950,550]
	plot={auto:auto,up:up,low:low,rho:0,cont:0,csize:csize,esize:esize,bsize:bsize,ratio:0,io:ptr_new([0],/allocate_heap),jo:ptr_new([0],/allocate_heap),t1:[0.0,0.0],t2:[0.0,0.0],$
		iostr:'',jostr:'',t1str:'0.0,0.0',t2str:'0.0,0.0'}
	stat={ps:0,dat:0,m1stat:0,m2stat:0,ntime:0L,time:time,col:[255,255,50,80,200,150],pscol:[0,0,30,100,200,150],nth:100,nr:0}
	u={id:id,shot:shot,index:0L,ch:0,stat:stat,plot:plot}
	widget_control,base,set_uvalue=u
	widget_control,u.id.shotid,set_value=num2str(u.shot,1)
	
	autoid=[u.id.chauto,u.id.brauto,u.id.dbrauto,u.id.rauto,u.id.emauto,u.id.tauto]
	upid=[u.id.chmax,u.id.brmax,u.id.dbrmax,u.id.rmax,u.id.emmax,u.id.tmax]
	lowid=[u.id.chmin,u.id.brmin,u.id.dbrmin,u.id.rmin,u.id.emmin,u.id.tmin]
	FOR i=0,n(autoid) DO BEGIN
		IF u.plot.auto[i] THEN widget_control,autoid[i],set_button=1 ELSE BEGIN
			widget_control,upid[i],set_value=num2str(u.plot.up[i],dp=2)
			widget_control,lowid[i],set_value=num2str(u.plot.low[i],dp=2)
		ENDELSE
	ENDFOR
	
	IF u.plot.rho THEN widget_control,u.id.xrho,set_button=1 ELSE widget_control,u.id.xrmaj,set_button=1
	CASE u.plot.cont OF
		0 : widget_control,u.id.em2d,set_button=1
		1 : widget_control,u.id.emr,set_button=1
		2 : widget_control,u.id.emc1,set_button=1
		3 : widget_control,u.id.emc2,set_button=1
		4 : widget_control,u.id.ems1,set_button=1
	ENDCASE
	widget_control,u.id.t1_val,set_value=u.plot.t1str		
	widget_control,u.id.t2_val,set_value=u.plot.t2str		


	!except=0
	widget_control,base,/realize
	xmanager,'w_xtomo',base
END
