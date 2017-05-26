PRO wspec_load_linestr,line
	get_genie_env, VUV_LINE_LIST=VUV_LINE_LIST, VIS_LINE_LIST=VIS_LINE_LIST

	print,' loading lines from '+VUV_LINE_LIST
	tsv_read,VUV_LINE_LIST
		lam0=wave
		elem0=elem
		cs0=cs
		note0=note
		iso0=iso
		priority0=priority
		source0=source

	combine=1
	IF combine THEN BEGIN
		print,' loading lines from '+VIS_LINE_LIST
		tsv_read,VIS_LINE_LIST
		lam=[lam0,wave*10] ; [A]
		elem=[elem0,elem]
		cs=[cs0,cs]
		note=[note0,note]
		iso=[iso0,iso]
		priority=[priority0,priority]
		source=[source0,source]
        ENDIF ELSE BEGIN
		lam=lam0
		elem=elem0
		cs=cs0
		note=note0
		iso=iso0
		priority=priority0
		source=source0
	ENDELSE

	; make a list of the unique elements
	tmp=elem[sort(elem)]
	list=strjoin(tmp[uniq(tmp)],',')

	line={lam:lam,elem:elem,cs:cs,note:note,iso:iso,priority:priority,source:source,list:list,usrlist:''}
END

PRO wspec_calc_br,u,line
	specfound=0
	n_lines=n((*u).mach.line.nam)+1

	for j=0,n((*u).mach.inst) do begin
	 i=n((*u).mach.inst)-j
	 if (*u).inst.(i).s then begin
		tmp=where((*u).inst.(i).l GE (*u).mach.line.lam0[line] $
		      AND (*u).inst.(i).l LE (*u).mach.line.lam1[line])
		if tmp[0] ne -1 then begin
	 		specfound=1
			ntime=n((*u).inst.(i).t)+1
			nch=n((*u).inst.(i).ch_nam)+1
			br=fltarr(ntime,nch)
			fwhm=fltarr(ntime,nch)
			FOR k=0,nch-1 DO BEGIN
				br[*,k]=(sum_array((*u).inst.(i).dch[tmp,*,k],/j) $ ; sum over the range, use total instead of sum_array?
			    -0.5*((*u).inst.(i).dch[tmp[0],*,k]+(*u).inst.(i).dch[last(tmp),*,k])*(n(tmp)+1) $ ; subtract offset below line
			   )*( (*u).mach.line.mult[line] )
			ENDFOR
			tbr=(*u).inst.(i).t
			spec=(*u).mach.inst[i]
		endif
	 endif
	endfor

	if specfound eq 0 then begin
		; st=execute('widget_control,(*u).id.p'+(*u).mach.inst[i]+',set_button=0')
		;(*u).mach.line.plot[line]=0
		br=[0,0]
		tbr=br
		spec='none'
		print, 'wspec_calc_br error: cant find loaded instrument for line '+num2str(line,1) $
		      +' in spectral region '+num2str((*u).mach.line.lam0[line],1)+'-'+num2str((*u).mach.line.lam1[line],1)
	endif

	*(*u).dat.br[line]=br
	*(*u).dat.br[line+n_lines]=tbr
	(*u).mach.line.spec[line]=spec
END

PRO wspec_calc_timetr,u
	xr=[(*u).plot.x0[0],(*u).plot.x1[0]]
	tau=make(xr[0],xr[1],(*u).plot.ntau)

	for i=0,n_tags((*u).timetr)-1 do begin
	 if (*u).timetr.(i).s $
	  then (*u).dat.timetr[i,*]=interpol((*u).timetr.(i).d,(*u).timetr.(i).t,tau) $
	  else (*u).dat.timetr[i,*]=fltarr((*u).plot.ntau)

	 ; set yrange for time traces
	 if (*u).timetr.(i).oplot eq 0 then begin
	 	loc_max=0.0
		if (*u).timetr.(i).s then begin
			tmp=where((*u).timetr.(i).t GE xr[0] $
			      AND (*u).timetr.(i).t LE xr[1])
			IF tmp[0] NE -1 THEN loc_max=max((*u).timetr.(i).d[tmp])
		endif
		j=i+1
		while j le n_tags((*u).timetr)-1 do begin
		 if (*u).timetr.(j).oplot eq 0 then break
		 if (*u).timetr.(j).s then begin
			tmp=where((*u).timetr.(j).t GE xr[0] $
			      AND (*u).timetr.(j).t LE xr[1])
			IF tmp[0] NE -1 THEN BEGIN
			 loc_max1=max((*u).timetr.(j).d[tmp])
			 if loc_max1 gt loc_max then loc_max=loc_max1
			ENDIF
		 endif
		 j=j+1
		endwhile

		yr=[-0.05,1.05]*loc_max
		(*u).timetr.(i).yr=yr
	 endif
	endfor

	; set yrange for tracked lines
	ymax=0.1
	n_lines=n((*u).mach.line.nam)+1
	FOR i=0,n_lines-1 DO $ 
	 IF (*u).mach.line.plot[i] THEN BEGIN
		IF (*u).mach.line.spec[i] EQ 'chromex' THEN widget_control,(*u).id.slchromex,get_value=index ELSE index=0
		x=*(*u).dat.br[i+n_lines]
		y=*(*u).dat.br[i]
		y=y[*,index]
		tmp=where(x GE xr[0] AND x LE xr[1])
		IF tmp[0] NE -1 THEN ymax=ymax > max(y[tmp])
         ENDIF
	(*u).line.yr[1]=ymax*1.05
END

PRO wspec_plot_time,u
	IF (*u).plot.isbr THEN wspec_plot_br,u
	IF (*u).plot.issp THEN wspec_plot_sp,u
	IF (*u).stat.ps THEN BEGIN
		xsize=6.0
		ysize=9.5
		ls=1.1
		tit=num2str((*u).shot,1)+' t='+num2str((*u).stat.time,dp=3)
		d_old=!d
		lcol=(*u).line.pscol
		device,xsize=xsize,ysize=ysize,/in
	ENDIF ELSE BEGIN
		ls=1.2
		widget_control,(*u).id.draw1,get_value=draw_win
		window,0,xsize=(*u).plot.asize[0],ysize=(*u).plot.asize[1],/pixmap
		lcol=(*u).line.col
		tit=''
	ENDELSE

	xr=[(*u).plot.x0[0],(*u).plot.x1[0]]
	tau=make(xr[0],xr[1],(*u).plot.ntau)
	nt=(*u).plot.nt
	nm=(*u).plot.nm
	low=0.07
	del=(0.96-low)/5
	tau0=[(*u).stat.time,(*u).stat.time]
	tlinecol=0

	for i=0,n_tags((*u).timetr)-1 do begin
	 if (*u).timetr.(i).oplot eq 0 then begin
	  yr=(*u).timetr.(i).yr
	  im=n_tags((*u).timetr)-1 - i
	  pos=[0.13,low+im*del,0.9,low+(im+1)*del]
	  plot,[0],[0],xr=xr,yr=yr,pos=pos,/xsty,chars=1.0*ls,xticks=xt,xminor=nm,xtickname=replicate(' ',nt+2),/noerase,/ysty,ytit=(*u).timetr.(i).dlabel
	  xyouts,1.05*(xr[1]-xr[0])+xr[0],0.6*(yr[1]-yr[0])+yr[0],(*u).timetr.(i).title,chars=1.0*ls,orient=90
	 endif else $
	  xyouts,1.05*(xr[1]-xr[0])+xr[0],0.2*(yr[1]-yr[0])+yr[0],(*u).timetr.(i).title,chars=1.0*ls,orient=90,color=128
	 IF (*u).timetr.(i).s ne 0 $
	  THEN IF (*u).timetr.(i).oplot eq 0 $
	  THEN oplot,tau,(*u).dat.timetr[i,*] $
	  ELSE oplot,tau,(*u).dat.timetr[i,*],color=128
	  ;THEN oplot,tau,interpol((*u).timetr.(i).d,(*u).timetr.(i).t,tau) $
	  ;ELSE oplot,tau,interpol((*u).timetr.(i).d,(*u).timetr.(i).t,tau),color=128
	 IF (*u).stat.tline THEN oplot,tau0,yr,linestyle=1.0,color=tlinecol
	 IF i eq 0 && (*u).stat.ps THEN xyouts,0.8*(xr[1]-xr[0])+xr[0],1.05*(yr[1]-yr[0])+yr[0],'shot '+num2str((*u).shot,1),chars=0.8*ls
	endfor

	i=0
	yr=(*u).line.yr
	pos=[0.13,low+i*del,0.9,low+(i+1)*del]
	plot,[0],[0],xr=xr,pos=pos,/xsty,chars=1.0*ls,xticks=xt,xminor=nm,xtit='Time [sec]',yr=yr,/ysty,/noerase,ytit='[AU]'
	xyouts,1.05*(xr[1]-xr[0])+xr[0],0.15*(yr[1]-yr[0])+yr[0],'Line Brightness',chars=1.0*ls,orient=90
	ylab=0.9*(yr[1]-yr[0])+yr[0]
	xlab=0.85*(xr[1]-xr[0])+xr[0]
	n_lines=n((*u).mach.line.nam)+1

	for i=0,n_lines-1 do begin
	 IF (*u).mach.line.plot[i] THEN BEGIN
		IF (*u).mach.line.spec[i] EQ 'chromex' THEN widget_control,(*u).id.slchromex,get_value=index ELSE index=0
		loadct,(*u).mach.line.ct[i],/silent
		ibr=*(*u).dat.br[i]
		IF index EQ -1 THEN ibr=sum_array(ibr,/i)/(n(ibr[0,*])+1.0) ELSE ibr=ibr[*,index]
		oplot,*(*u).dat.br[n_lines+i],ibr,color=(*u).mach.line.col[i]
		xyouts,xlab,ylab,(*u).mach.line.nam[i],color=(*u).mach.line.col[i],chars=0.8*ls
		ylab-=0.1*(yr[1]-yr[0])
	 ENDIF
      endfor
	loadct,12,/silent
	IF (*u).stat.tline THEN oplot,tau0,yr,linestyle=1.0,color=tlinecol,thick=2.0

	IF NOT (*u).stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,(*u).plot.asize[0],(*u).plot.asize[1],0,0,0]
		(*u).plot.str1={X:!X,Y:!Y,Z:!Z,P:!P}
	ENDIF ELSE device, xsize=float(d_old.x_size)/d_old.x_px_cm, $
			   ysize=float(d_old.y_size)/d_old.y_px_cm


END

PRO wspec_plot_br,u
	
	IF (*u).stat.ps THEN BEGIN
		xsize=6.0
		ysize=6.0*6.0/8.0
		ls=1.1
		tit=num2str((*u).shot,1)+' t='+num2str((*u).stat.time,dp=3)
		d_old=!d
		lcol=(*u).line.pscol
		device,xsize=xsize,ysize=ysize,/in
	ENDIF ELSE BEGIN
		ls=1.2
		widget_control,(*u).id.br_draw,get_value=draw_win
		window,0,xsize=(*u).plot.brsize[0],ysize=(*u).plot.brsize[1],/pixmap
		lcol=(*u).line.col
		tit=''
	ENDELSE

	xr=[(*u).plot.x0[0],(*u).plot.x1[0]]
	tau=make(xr[0],xr[1],(*u).plot.ntau)
	nt=(*u).plot.nt
	nm=(*u).plot.nm
	low=0.07
	del=(0.96-low)/5
	tau0=[(*u).stat.time,(*u).stat.time]
	tlinecol=0
	chk=where((*u).mach.line.spec EQ 'chromex')
	IF chk[0] NE -1 THEN BEGIN
		widget_control,(*u).id.slchromex,get_value=index
		inst=where((*u).mach.inst EQ 'chromex')
		IF index EQ -1 THEN tit='AVE' ELSE tit=(*u).inst.(inst).ch_nam[index]
        ENDIF ELSE tit=''

	i=0
	yr=(*u).line.yr
	plot,[0],[0],xr=xr,/xsty,chars=1.0*ls,xticks=xt,xminor=nm,xtit='Time [sec]',yr=yr,/ysty,ytit='[AU]',tit=tit
	xyouts,1.05*(xr[1]-xr[0])+xr[0],0.15*(yr[1]-yr[0])+yr[0],'Line Brightness',chars=1.0*ls,orient=90
	ylab=0.9*(yr[1]-yr[0])+yr[0]
	xlab=0.85*(xr[1]-xr[0])+xr[0]
	n_lines=n((*u).mach.line.nam)+1

	for i=0,n_lines-1 do begin
	IF (*u).mach.line.plot[i] THEN BEGIN
		IF (*u).mach.line.spec[i] EQ 'chromex' THEN widget_control,(*u).id.slchromex,get_value=index ELSE index=0
		loadct,(*u).mach.line.ct[i],/silent
		ibr=*(*u).dat.br[i]
		IF index EQ -1 THEN ibr=sum_array(ibr,/i)/(n(ibr[0,*])+1.0) ELSE ibr=ibr[*,index]
		oplot,*(*u).dat.br[n_lines+i],ibr,color=(*u).mach.line.col[i]
		xyouts,xlab,ylab,(*u).mach.line.nam[i],color=(*u).mach.line.col[i],chars=0.8*ls
		ylab-=0.1*(yr[1]-yr[0])
	 ENDIF
      endfor
	loadct,12,/silent
	IF (*u).stat.tline THEN oplot,tau0,yr,linestyle=1.0,color=tlinecol,thick=2.0

	IF NOT (*u).stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,(*u).plot.brsize[0],(*u).plot.brsize[1],0,0,0]
		(*u).plot.str1={X:!X,Y:!Y,Z:!Z,P:!P}
	ENDIF ELSE device, xsize=float(d_old.x_size)/d_old.x_px_cm, $
			   ysize=float(d_old.y_size)/d_old.y_px_cm

END

PRO wspec_plot_sp,u
	IF (*u).stat.ps THEN BEGIN
		xsize=8.0*600/1000.0
		ysize=8.0
		ls=1.1
		tit=num2str((*u).shot,1)+' t='+num2str((*u).stat.time,dp=3)
		d_old=!d
		lcol=(*u).line.pscol
		device,xsize=xsize,ysize=ysize,/in
	ENDIF ELSE BEGIN
		ls=1.2
		widget_control,(*u).id.sp_draw,get_value=draw_win
		window,0,xsize=(*u).plot.spsize[0],ysize=(*u).plot.spsize[1],/pixmap
		lcol=(*u).line.col
		tit=''
	ENDELSE
	k=where((*u).mach.inst EQ 'chromex')
	widget_control,(*u).id.slchromex,get_value=sl
	IF sl NE -1 THEN BEGIN
		slperi=(*u).inst.(k).c.periscop[sl]
		slfib=(*u).inst.(k).c.per_fibr[sl]
        ENDIF 
	peri=['A_BOTTOM','K_BOTTOM','K_TOP']
	x0=[5,0,0]
	x1=[18,19,13]
	yr=(*u).line.yr
	widget_control, (*u).id.sp_slider,get_value=j
	nlines=n((*u).mach.line.nam)+1
	time=*(*u).dat.br[nlines]	;time base,stored in weird way
	tau=ipt(time,(*u).stat.time)	;index in timebase
	!p.multi=[0,0,n(peri)+1.0]
	FOR i=0,n(peri) DO BEGIN
		xr=[x0[i],x1[i]]
		xlab=xr[0]+0.1*(xr[1]-xr[0])
		ylab=yr[0]+0.8*(yr[1]-yr[0])
		plot,[0],[0],xr=xr,/xsty,chars=1.5*ls,xticks=xt,xminor=nm,xtit='CH#',yr=yr,/ysty,ytit='[AU]',tit=peri[i]
		index=where((*u).inst.(k).c.periscop EQ peri[i])
		IF index[0] NE -1 THEN BEGIN
			ch=(*u).inst.(k).c.per_fibr[index]
			ibr=*(*u).dat.br[j]
			br=ibr[tau,index]
			loadct,(*u).mach.line.ct[j],/silent
			order=sort(ch)
			makesym,9
			oplot,ch[order],br[order],psym=-8,color=(*u).mach.line.col[j]
			makesym,10
			IF sl NE -1 THEN BEGIN
				IF slperi EQ peri[i] THEN BEGIN
					oplot,[slfib],br[where(ch EQ slfib)],psym=8,color=(*u).mach.line.col[j]
					oplot,slfib*[1,1],yr,linestyle=1.0,color=0
				ENDIF
                        ENDIF ELSE BEGIN
				oplot,ch[order],br[order],psym=-8,color=(*u).mach.line.col[j]
			ENDELSE
			xyouts,xlab,ylab,(*u).mach.line.nam[j],color=(*u).mach.line.col[j],chars=1.5*ls
			xyouts,xlab,ylab-0.1*(yr[1]-yr[0]),'t='+num2str(time[tau],dp=2),chars=1.5*ls
                ENDIF
	ENDFOR
	!p.multi=0
	IF NOT (*u).stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,(*u).plot.spsize[0],(*u).plot.spsize[1],0,0,0]
		(*u).plot.str1={X:!X,Y:!Y,Z:!Z,P:!P}
	ENDIF ELSE device, xsize=float(d_old.x_size)/d_old.x_px_cm, $
			   ysize=float(d_old.y_size)/d_old.y_px_cm

END

FUNCTION wspec_maxsbr,u,min=min
	xr=[(*u).plot.x0[1],(*u).plot.x1[1]]
	ymax=0.0
	for i=0,n((*u).mach.inst) do begin
		IF (*u).inst.(i).s && (*u).mach.plot[i] THEN BEGIN
			index=ipt((*u).stat.time,(*u).inst.(i).t)
			IF index[0] ne -1 THEN BEGIN
			 IF keyword_set(min) THEN BEGIN
				tmp=where((*u).inst.(i).l GE xr[0] $
				      AND (*u).inst.(i).l LE xr[1] $
				      AND (*u).inst.(i).d[*,index] GT 0.0)
				if tmp[0] ne -1 then ymax=ymax > min((*u).inst.(i).d[tmp,index])
			 ENDIF ELSE BEGIN
				tmp=where((*u).inst.(i).l GE xr[0] $
				      AND (*u).inst.(i).l LE xr[1])
				if tmp[0] ne -1 then ymax=ymax > max((*u).inst.(i).d[tmp,index])
			 ENDELSE
			ENDIF
       		ENDIF
	endfor
	RETURN,ymax
END

PRO wspec_plot_spec,u
	; setup plot window
	IF (*u).stat.ps THEN BEGIN
		xsize=7.0
		ysize=4.5
		ls=0.9
		col=(*u).plot.pscol 
		lcol=(*u).line.pscol
		tit=num2str((*u).shot,1)+' t='+num2str((*u).stat.time,dp=3)
		d_old=!d
		device,xsize=xsize,ysize=ysize,/in
	ENDIF ELSE BEGIN
		ls=1.0
		col=(*u).plot.col
		lcol=(*u).line.col
		widget_control,(*u).id.draw2,get_value=draw_win
		window,0,xsize=(*u).plot.bsize[0],ysize=(*u).plot.bsize[1],/pixmap ; draw into a hidden pixmap buffer offscreen, then copy into the window
		tit=''
        ENDELSE
	xr=[(*u).plot.x0[1],(*u).plot.x1[1]]
	IF (*u).plot.auto THEN BEGIN
		if (*u).plot.log $
		 then ymn=wspec_maxsbr(u,/min) $
		 else ymn=0.0
		yr=[ymn,wspec_maxsbr(u)*1.05]
		(*u).plot.yr=yr
		widget_control, (*u).id.sbr0,set_value=num2str((*u).plot.yr[0],dp=2)
		widget_control, (*u).id.sbr1,set_value=num2str((*u).plot.yr[1],dp=2)
        ENDIF ELSE BEGIN 
		yr=(*u).plot.yr
	ENDELSE
        IF (*u).plot.log THEN ylog=1
	for i=0,n((*u).mach.inst) do $
	 if (*u).inst.(i).s then begin
	 	llabel=(*u).inst.(i).llabel
		dlabel=(*u).inst.(i).dlabel
		break
	 endif
	plot,[0],[0],xr=xr,yr=yr,ylog=ylog,/xsty,/ysty,xtit=llabel,ytit=dlabel,chars=1.2*ls

	; plot nearest time slice for loaded and selected instruments
	for i=0,n((*u).mach.inst) do begin
		IF (*u).inst.(i).s && (*u).mach.plot[i] THEN BEGIN
			index=ipt((*u).stat.time,(*u).inst.(i).t)
			IF index NE -1 THEN oplot,(*u).inst.(i).l,(*u).inst.(i).d[*,index],col=col[0]
			(*u).inst.(i).index=index
	        ENDIF 
	endfor
 
	ylab=1.02*(yr[1]-yr[0])+yr[0]
	xlab=0.1*(xr[1]-xr[0])+xr[0]

	; overplot time history tracked line regions
	FOR i=0,n((*u).mach.line.plot) DO BEGIN
		spec=(*u).mach.line.spec[i]
		IF (*u).mach.line.plot[i] THEN BEGIN
			specfound=0
			for j=0,n((*u).mach.inst) do begin
			 if spec eq (*u).mach.inst[j] && (*u).inst.(j).s then begin
			 	specfound=1
				loadct,(*u).mach.line.ct[i],/silent
				tmp=where((*u).inst.(j).l GE (*u).mach.line.lam0[i] $
				      AND (*u).inst.(j).l LE (*u).mach.line.lam1[i])
				index=ipt((*u).stat.time,(*u).inst.(j).t)
				IF index NE -1 THEN oplot,(*u).inst.(j).l[tmp],(*u).inst.(j).d[tmp,index],col=(*u).mach.line.col[i]
			 endif
			endfor
			;if specfound eq 0 then print,'spec not found:'+spec

			IF (*u).stat.ps THEN xyouts,xlab,ylab,(*u).mach.line.nam[i],color=lcol[i]
			xlab+=0.2*(xr[1]-xr[0])
               ENDIF
	ENDFOR
	loadct,12,/silent
	; overplot labels for known lines
	IF (*u).plot.label THEN BEGIN
		widget_control,(*u).id.lelem,get_value=lelem
		eplot=strsplit(lelem,',',/extract)
		FOR order=1,(*u).plot.order DO BEGIN
		 tmp=where( order*(*(*u).dat.line).lam GT xr[0] $
		        AND order*(*(*u).dat.line).lam LT xr[1] $
			AND (*(*u).dat.line).priority LE (*u).plot.priority )
		 IF tmp[0] NE -1 THEN BEGIN
			FOR i=0,n(tmp) DO BEGIN
				IF total(where( eplot EQ strtrim((*(*u).dat.line).elem[tmp[i]]) )) NE -1 AND (*(*u).dat.line).priority[tmp[i]] LE (*u).plot.priority THEN BEGIN
				 oplot, order*(*(*u).dat.line).lam[tmp[i]]*[1,1],yr,linestyle=1.0,col=col[order-1]
				 xyouts,order*(*(*u).dat.line).lam[tmp[i]]-0.005*(xr[1]-xr[0]),yr[0]+(yr[1]-yr[0])*0.99, $
					      (*(*u).dat.line).elem[tmp[i]]+' '+(*(*u).dat.line).cs[tmp[i]], $
					      orient=90,align=1,col=col[order-1]
				ENDIF
			ENDFOR
		 ENDIF
		ENDFOR
	ENDIF

	; overplot nist line labels
;	IF (*u).nist.p && (*u).nist.s THEN BEGIN
;		FOR i=0,n(*(*u).nist.ln) DO BEGIN
;		 oplot, (*(*u).nist.ln)[i].wl*[1,1],yr,linestyle=1.0,col=col[0]
;		 xyouts,(*(*u).nist.ln)[i].wl-0.005*(xr[1]-xr[0]),yr[0]+(yr[1]-yr[0])*0.99, $
;			(*(*u).nist.ln)[i].spec,orient=90,align=1,col=col[0]
;		ENDFOR
;	ENDIF
;	IF (*u).stat.ps THEN xyouts,1.03*(xr[1]-xr[0])+xr[0],0.1*(yr[1]-yr[0])+yr[0], $
;				 num2str((*u).shot,1)+'     t='+num2str((*u).stat.time,dp=3), $
;				 chars=0.8*ls,orient=90

	; copy the buffer into the selected device
	IF NOT (*u).stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,(*u).plot.bsize[0],(*u).plot.bsize[1],0,0,0]
		(*u).plot.str2={X:!X,Y:!Y,Z:!Z,P:!P}
	ENDIF ELSE device, xsize=float(d_old.x_size)/d_old.x_px_cm, $
			   ysize=float(d_old.y_size)/d_old.y_px_cm

END


PRO wspec_plot,u
	wspec_plot_time,u
	wspec_plot_spec,u
	wspec_plot_views,u
END

PRO wspec_plot_views,u
	IF (*u).stat.ps THEN BEGIN
		xsize=7.0
		ysize=7.25
		ls=1.1
		tit=num2str((*u).shot,1)+' t='+num2str((*u).stat.time,dp=3)
		d_old=!d
		lcol=(*u).line.pscol
		device,xsize=xsize,ysize=ysize,/in
	ENDIF ELSE BEGIN
		ls=1.2
		widget_control,(*u).id.draw3,get_value=draw_win
		window,0,xsize=(*u).plot.csize[0],ysize=(*u).plot.csize[1],/pixmap
		lcol=(*u).line.col
		tit=''
	ENDELSE

	;hard coded for post 2007 structure
	yrange=[-0.67,-0.15]
	xrange=[0.38,0.85]
	titstr=num2str((*u).shot,1)+' t='+num2str((*u).stat.time,dp=2)
	plot, [0],[0],title=titstr,chars=0.9,xrange=xrange,yrange=yrange,/xsty, /ysty
	IF (*u).shot LT 1070101000 THEN widget_control,(*u).id.message,set_value='Warning: vessel structure out of date',/append
	restore, "/home/labombard/minicad/vv_tiles_cryo_2007_s.vctr"
	FOR i=0,nvctr-1 DO oplot,xvctr(0:lvctr(i)-1,i),yvctr(0:lvctr(i)-1,i)

	index=ipt((*u).dat.fs.t,(*u).stat.time)
	IF index[0] NE -1 THEN BEGIN
		z_xpt=(*u).dat.fs.zxpt[index]
		rho=[make(min((*u).dat.fs.rhofine),-.01,8),-0.004,0.0,0.004,0.008]
        	FOR i=0,n(rho) DO BEGIN
        		contour,(*u).dat.fs.rhofine[*,*,index],(*u).dat.fs.rfine,(*u).dat.fs.zfine,levels=[rho[i]],$
				/overplot,path_xy=path_fs,path_info=info_fs,/path_data
	        	rfs=reform(path_fs[0,*])
        		zfs=reform(path_fs[1,*])
                	good=inside(rfs,zfs,(*u).dat.fs.rw,(*u).dat.fs.zw)
	                IF z_xpt LT 0 THEN BEGIN
                		IF rho[i] LT 0 THEN BEGIN
                        		tmp=where(good EQ 1 AND zfs GT z_xpt) 
		                        color=100
        		        ENDIF ELSE BEGIN
                		        tmp=where(good EQ 1)
                        		IF rho[i] EQ 0 THEN color=200 ELSE color=30
		                ENDELSE
        			IF tmp[0] NE -1 THEN oplot,rfs[tmp],zfs[tmp],color=color
                		IF rho[i] LT 0 THEN BEGIN
              				tmp=where(good EQ 1 AND zfs LT z_xpt) 
                                	color=30
		                        IF tmp[0] NE -1 THEN oplot,rfs[tmp],zfs[tmp],color=color
        	                ENDIF
                	ENDIF ELSE BEGIN
                		IF rho[i] LT 0 THEN BEGIN
                        		tmp=where(good EQ 1 AND zfs LT z_xpt) 
		                        color=100
        		        ENDIF ELSE BEGIN
                		        tmp=where(good EQ 1 AND zfs GT -0.4)
                        		IF rho[i] EQ 0 THEN color=200 ELSE color=30
		                ENDELSE
        			IF tmp[0] NE -1 THEN oplot,rfs[tmp],zfs[tmp],color=color
                		IF rho[i] LT 0 THEN BEGIN
              				tmp=where(good EQ 1 AND zfs GT z_xpt)
	                        	IF tmp[0] NE -1 THEN oplot,rfs[tmp],zfs[tmp],color=color
                        	ENDIF
                
                        ENDELSE	
		ENDFOR
	ENDIF
	FOR i=0,n_elements((*u).mach.inst)-1 DO BEGIN
		IF strlowcase((*u).mach.inst[i]) EQ 'chromex' AND (*u).inst.(i).s THEN BEGIN
			st=execute('widget_control,(*u).id.sl'+(*u).mach.inst[i]+',get_value=slider')
			rzview=(*u).inst.(i).rzview
			FOR k=0,n_elements(rzview[0,*])-1 DO BEGIN
				oplot,[rzview[0,k],rzview[2,k]],[rzview[1,k],rzview[3,k]],color=120
				IF k EQ slider OR slider EQ -1 THEN oplot,[rzview[0,k],rzview[2,k]],[rzview[1,k],rzview[3,k]],color=120,thick=5.0
			ENDFOR	
		ENDIF
	ENDFOR
	IF NOT (*u).stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,(*u).plot.csize[0],(*u).plot.csize[1],0,0,0]
		(*u).plot.str2={X:!X,Y:!Y,Z:!Z,P:!P}
	ENDIF ELSE device, xsize=float(d_old.x_size)/d_old.x_px_cm, $
			   ysize=float(d_old.y_size)/d_old.y_px_cm
			

END

PRO wspec_msg,u,str
	widget_control,(*u).id.message,set_value=str,/append
	print,str	
END

PRO wspec_load_data,u
	widget_control,/hourglass
	shot=(*u).shot
	mach=(*u).mach.name
	
	wspec_msg,u,'loading shot '+num2str(shot,1)+'...'


	for i=0,n((*u).mach.timetr) do begin
	 ttr=(*u).mach.timetr[i]
	 st=execute('load_'+mach+'_'+ttr+',shot,d'+num2str(i,1))
	endfor
	str='timetr={'
	for i=0,n((*u).mach.timetr) do begin
	 if i lt n((*u).mach.timetr) then sep=',' else sep=''
	 str=str+'d'+num2str(i,1)+':d'+num2str(i,1)+sep
	endfor
	str=str+'}'
	st=execute(str)
	wspec_msg,u,' time traces loaded'
	
	; calculate the timetrace values
	; just initialize the timetrace array for now, it will be filled later in wspec_calc_timetr
	dattimetr=fltarr(n((*u).mach.timetr)+1,(*u).plot.ntau)

	;load spectroscopic instruments
	t_all=0.0
	for i=0,n((*u).mach.inst) do begin
		IF (*u).mach.load[i] $
		 THEN st=execute('load_'+(*u).mach.name+'_'+(*u).mach.inst[i]+',shot,d') $ ; st=1 on success
		 ELSE d={s:0}
		IF d.s THEN BEGIN
			d.index=ipt(d.t,(*u).stat.time)
			st=execute('widget_control,(*u).id.p'+(*u).mach.inst[i]+',set_button=1')
			IF (*u).mach.nch[i] GT 1 THEN BEGIN
				st=execute('widget_control,(*u).id.sl'+(*u).mach.inst[i]+',set_slider_min=-1,set_slider_max=' $
					+num2str( (size(d.dch,/dim))[2]-1 ,1) ) ; just in case it changed
				st=execute('widget_control,(*u).id.sl'+(*u).mach.inst[i]+',get_value=slider')
				;stop
				IF slider EQ -1 THEN BEGIN
					dim=size(d.dch,/dim)
					nch=dim[2]
					d.d=total(d.dch,3)/nch
					d.l=total(d.lch,2)/nch
					st=execute('widget_control,(*u).id.lch'+(*u).mach.inst[i]+',set_value=''SUM''')
				ENDIF ELSE BEGIN
					d.d=d.dch(*,*,slider)
					d.l=d.lch(*,slider)
					st=execute('widget_control,(*u).id.lch'+(*u).mach.inst[i]+',set_value=(*u).inst.(i).ch_nam['+num2str(slider,1)+']')
				ENDELSE
			ENDIF
			t_all=[t_all,d.t]
			wspec_msg,u,' '+(*u).mach.inst[i]+' loaded'+d.msg
	        ENDIF ELSE BEGIN
			; Don't change any settings, just let plotting routines work around missing data.
			; st=execute('widget_control,(*u).id.p'+(*u).mach.inst[i]+',set_button=0')
			; (*u).mach.line.plot[line]=0

			IF (*u).mach.load[i] THEN wspec_msg,u,' no '+(*u).mach.inst[i]+' data'
	        ENDELSE
		st=execute('d'+num2str(i,1)+'=d')
	endfor
	str='inst={'
	for i=0,n((*u).mach.inst) do begin
	 if i lt n((*u).mach.inst) then sep=',' else sep=''
	 is=num2str(i,1)
	 str=str+'d'+is+':d'+is+sep
	endfor
	str=str+'}'
	st=execute(str) & if not st then stop

	(*u).plot.x0[0]=min(t_all)
	(*u).plot.x1[0]=max(t_all)
	widget_control,(*u).id.t_slider,set_slider_max=(*u).plot.x1[0]*1000,set_slider_min=(*u).plot.x0[0]*1000

	wspec_msg,u,'done loading data'

	;initialize brightness ptrarr
	nlines=n((*u).mach.line.nam)+1
	IF (*u).stat.dat THEN BEGIN
		heap_free,(*u).dat.br 		;release br
		heap_free,(*u).dat.sbr		;release sbr
	ENDIF
	br =ptrarr(nlines*2,/allocate_heap)
	sbr=ptrarr(nlines*2,/allocate_heap)

	;load LINE table
	line=ptr_new(1)
	wspec_load_linestr,*line

	;load FS structure
	fs=make_fs_struc(shot)

	(*u).stat.dat=1
	dat={line:line,br:br,sbr:sbr,timetr:dattimetr,fs:fs}
	*u={mach:(*u).mach,inst:inst,timetr:timetr,id:(*u).id,shot:(*u).shot,stat:(*u).stat,plot:(*u).plot,line:(*u).line,dat:dat,nist:(*u).nist}
	widget_control,(*u).id.base, set_uvalue=u
	
	FOR i=0,n((*u).mach.line.nam) DO wspec_calc_br,u,i
	wspec_calc_timetr,u
	wspec_set_tzoomout,u
	wspec_set_lzoomout,u
	wspec_plot_time,u
	wspec_plot_spec,u
END

PRO wspec_load_nistasd,u
	if (*u).nist.spec ne '' && (*u).nist.order ge 1 then begin
	  nist_asd,(*u).nist.spec,(*u).plot.x0[1],(*u).plot.x1[1],d
		*(*u).nist.ln=d.lines
	  wspec_msg,u,d.msg
	  if (*(*u).nist.ln)[0].spec eq 'none' $
	   then (*u).nist.s=0 $
	   else (*u).nist.s=1
	endif

	if (*u).nist.s eq 0 then wspec_msg,u,'could not load nist asd...'
END

PRO wspec_set_lzoom,u,lr
	(*u).plot.x0[1]=lr[0]
	(*u).plot.x1[1]=lr[1]
	widget_control, (*u).id.lam0,set_value=num2str((*u).plot.x0[1],dp=3)
	widget_control, (*u).id.lam1,set_value=num2str((*u).plot.x1[1],dp=3)
	wspec_plot_spec,u
END

PRO wspec_set_lzoomout,u
	l_all=!Values.F_NAN
	for i=0,n((*u).mach.inst) do $
	 if (*u).inst.(i).s && (*u).mach.plot[i] then l_all=[(*u).inst.(i).l,l_all]
	tmp=where(finite(l_all))
	if tmp[0] ne -1 then l_all = l_all[tmp]
	(*u).plot.x0[1]=min(l_all)
	(*u).plot.x1[1]=max(l_all)
	widget_control, (*u).id.lam0,set_value=num2str((*u).plot.x0[1],dp=3)
	widget_control, (*u).id.lam1,set_value=num2str((*u).plot.x1[1],dp=3)

	wspec_set_lzoom,u,[(*u).plot.x0[1],(*u).plot.x1[1]]
END

PRO wspec_set_tzoom,u,tr
	(*u).plot.x0[0]=tr[0]
	(*u).plot.x1[0]=tr[1]
	widget_control, (*u).id.tau0,set_value=num2str((*u).plot.x0[0],dp=3)
	widget_control, (*u).id.tau1,set_value=num2str((*u).plot.x1[0],dp=3)

	wspec_calc_timetr,u
	wspec_plot_time,u
END

PRO wspec_set_tzoomout,u
	t_all=!Values.F_NAN
	for i=0,n((*u).mach.inst) do $
	 if (*u).inst.(i).s && (*u).mach.plot[i] then t_all=[(*u).inst.(i).t,t_all]

	(*u).plot.x0[0]=min(t_all)
	(*u).plot.x1[0]=max(t_all)
	widget_control, (*u).id.tau0,set_value=num2str((*u).plot.x0[0],dp=3)
	widget_control, (*u).id.tau1,set_value=num2str((*u).plot.x1[0],dp=3)

	wspec_calc_timetr,u
	wspec_plot_time,u
END

PRO wspec_br_resize,u,event

END


PRO wspec_base_resize,u,event
	; following http://www.idlcoyote.com/widget_tips/resize_draw.html
	x=event.x
	y=event.y
	xmn=(*u).id.x0
	ymn=(*u).id.y0
	if x lt 2*xmn then x=2*xmn
	if y lt   ymn then y=ymn
	tfrac=0.3 ; fraction of window width for time plots
	
	widget_control,(*u).id.base, xsize=x, ysize=y ; also draw_xsize, scr_xsize, and xsize
	x1=(x-xmn)*tfrac-6
	x2=(x-xmn)*(1-tfrac)-6
	geo_B0=widget_info((*u).id.B0,/geometry) & yB0=geo_B0.ysize+2*geo_B0.margin+2*geo_B0.space
	geo_B2=widget_info((*u).id.B2,/geometry) & yB2=geo_B2.ysize+2*geo_B2.margin+2*geo_B2.space
	geo_draw1=widget_info((*u).id.draw1,/geometry) & yxtra=2*geo_draw1.margin+2*geo_draw1.space
	ydraw=y-yB0-yB2-yxtra-6
	(*u).plot.asize=[x1,ydraw]
	(*u).plot.bsize=[x2,ydraw]

	widget_control,(*u).id.A, scr_xsize=x1, scr_ysize=y
	 widget_control,(*u).id.A1,    xsize=(*u).plot.asize[0],ysize=(*u).plot.asize[1]
	 widget_control,(*u).id.draw1, scr_xsize=(*u).plot.asize[0],scr_ysize=(*u).plot.asize[1]
	widget_control,(*u).id.A, scr_xsize=0 ; I don't know why this is necessary, but in frustration I tried it and it works! Without this line the viewport size does not update

	widget_control,(*u).id.B, scr_xsize=x2, scr_ysize=y
	 widget_control,(*u).id.B1,    xsize=(*u).plot.bsize[0],ysize=(*u).plot.bsize[1]
	 widget_control,(*u).id.draw2, draw_xsize=(*u).plot.bsize[0],draw_ysize=(*u).plot.bsize[1]
	 widget_control,(*u).id.t_slider,xsize=x2-80
	widget_control,(*u).id.B, scr_xsize=0

	wspec_plot,u

	;id=(*u).id.A1 & help,/str,widget_info(id,/geometry),widget_info(id,/update)
	;id=(*u).id.draw1 & help,/str,widget_info(id,/geometry),widget_info(id,/update)
	;stop
END

PRO w_spec_event,event
	widget_control,event.top,get_uvalue=u
	id = (*u).id
	tag = tag_names(event,/st)
	button=' '
	idtags=tag_names(id)
	ename=''
	FOR i=0,n(idtags) DO IF id.(i) EQ event.id THEN ename=idtags[i]
	;print,tag+' '+ename
	CASE tag OF
		"WIDGET_BASE" : BEGIN
			IF event.id EQ (*u).id.base THEN wspec_base_resize,u,event ; resize events
			IF event.id EQ (*u).id.br_base THEN wspec_br_resize,u,event
		END
		"WIDGET_BUTTON": BEGIN
			widget_control,event.id,get_value=button,get_uvalue=uvalue
			IF NOT keyword_set(uvalue) THEN uvalue='none'
			CASE ename OF
				"LOAD" : wspec_load_data,u
				"QUIT" : BEGIN
					IF (*u).stat.dat THEN BEGIN
						heap_free,(*u).dat.br
						heap_free,(*u).dat.sbr
						heap_free,(*u).nist
						heap_free,u
					ENDIF
					widget_control,event.top,/destroy
					!except=1
					heap_gc
				END
				"PRINT" : BEGIN
					psplot
					(*u).stat.ps=1
					wspec_plot_time,u
					wspec_plot_spec,u
					psc
					(*u).stat.ps=0
					xwplot
					set_plot,'x'	
				END
				"STOP" : BEGIN
					stop
                                END
				"HELP" : BEGIN
					wspec_msg,u,'Launching wspec_ng help...'
					spawn,'firefox http://www.psfc.mit.edu/research/alcator/cmodwiki/index.php/W_spec_ng &'
                                END
        			"SETCURR": BEGIN
            				(*u).shot = mdsvalue('current_shot("cmod")')
            				widget_control,(*u).id.shotid,set_value=num2str((*u).shot,1)
					wspec_load_data,u
                                END
				'OUTZ' : BEGIN
					wspec_set_lzoomout,u
                                END
				'BRWIN' : BEGIN
					IF event.select EQ 1 THEN BEGIN
						(*u).plot.isbr=1
						widget_control,(*u).id.br_base,map=1 
						wspec_plot_br,u
                                        ENDIF  ELSE  BEGIN
						(*u).plot.isbr=0
						widget_control,(*u).id.br_base,map=0
					ENDELSE
				END
				'SPWIN' : BEGIN
					IF event.select EQ 1 THEN BEGIN
						(*u).plot.issp=1
						widget_control,(*u).id.sp_base,map=1 
						wspec_plot_sp,u
                                        ENDIF  ELSE  BEGIN
						(*u).plot.issp=0
						widget_control,(*u).id.sp_base,map=0
					ENDELSE
				END
				'AELEM' : BEGIN
					IF event.select EQ 1 THEN BEGIN
					 widget_control,(*u).id.lelem,get_value=data
					 (*(*u).dat.line).usrlist=data
					 widget_control,(*u).id.lelem,set_value=(*(*u).dat.line).list 
					ENDIF ELSE $
					 widget_control,(*u).id.lelem,set_value=(*(*u).dat.line).usrlist

					wspec_plot_spec,u
				END
				'OELEM' : BEGIN
					IF event.select EQ 1 THEN (*u).plot.order=2 ELSE (*u).plot.order=1
					wspec_plot_spec,u
				END
				'AUTOSBR' : BEGIN
					IF event.select EQ 1 THEN (*u).plot.auto=1 ELSE (*u).plot.auto=0
					wspec_plot_spec,u
				END
                                'LOGSBR' : BEGIN
                                        IF event.select EQ 1 THEN (*u).plot.log=1 ELSE (*u).plot.log=0
                                        wspec_plot_spec,u
                                END
				'LABEL' : BEGIN
					IF event.select EQ 1 THEN (*u).plot.label=1 ELSE (*u).plot.label=0
					wspec_load_linestr,(*(*u).dat.line)
					wspec_plot_spec,u
				END
				'TLINE' : BEGIN
					IF event.select EQ 1 THEN (*u).stat.tline=1 ELSE (*u).stat.tline=0
					wspec_plot_time,u
                                END
				'PLINE' : BEGIN
					IF event.select EQ 1 THEN (*u).stat.pline=1 ELSE (*u).stat.pline=0
					wspec_plot_views,u
                                END
				'LZOOM' : BEGIN
					wspec_plot_spec,u
					cursor,x1,y1,/up
					wait,0.1
					cursor,x2,y2,/up
					(*u).plot.x0[1]=x1
					(*u).plot.x1[1]=x2
					widget_control, (*u).id.lam0,set_value=num2str((*u).plot.x0[1],dp=1)
					widget_control, (*u).id.lam1,set_value=num2str((*u).plot.x1[1],dp=1)
					wspec_plot_spec,u
                                 END
				'TZOOM' : BEGIN
					wspec_plot_time,u
					cursor,x1,y1,/up
					wait,0.1
					cursor,x2,y2,/up
					(*u).plot.x0[0]=x1
					(*u).plot.x1[0]=x2
					widget_control, (*u).id.tau0,set_value=num2str((*u).plot.x0[0],dp=1)
					widget_control, (*u).id.tau1,set_value=num2str((*u).plot.x1[0],dp=1)
					wspec_plot_time,u
				END
				'NISTPLT' : BEGIN
				 (*u).nist.p=event.select
				 if (*u).nist.s eq 1 then wspec_plot_spec,u
				END
				'NISTGET' : BEGIN
				 wspec_load_nistasd,u
				 if (*u).nist.s eq 1 && (*u).nist.p eq 1 then wspec_plot_spec,u
				END
				ELSE : BEGIN
				 print,ename,(*u).stat.dat,event.select
				 IF (*u).stat.dat THEN BEGIN
				  for i=0,n((*u).mach.inst) do begin
				  	IF ename eq strupcase('l'+(*u).mach.inst[i]) $ ; instrument LOAD buttons
						THEN IF event.select EQ 1 $
						 THEN (*u).mach.load[i]=1 $
						 ELSE (*u).mach.load[i]=0
				  	IF ename eq strupcase('p'+(*u).mach.inst[i]) THEN BEGIN ; instrument PLOT buttons
						IF (*u).inst.(i).s && event.select eq 1 $
						 THEN (*u).mach.plot[i]=1 $
						 ELSE (*u).mach.plot[i]=0
						wspec_plot_spec,u
					ENDIF
				  	IF ename eq strupcase('z'+(*u).mach.inst[i]) && (*u).inst.(i).s $ ; instrument ZOOM buttons
					 THEN wspec_set_lzoom,u,[min((*u).inst.(i).l),max((*u).inst.(i).l)]
				  endfor
                                 ENDIF
				 for i=0,n((*u).mach.line.nam) do begin
				  is=num2str(i,1)
				  if ename eq strupcase('l'+is+'plot') then begin
				   if event.select eq 1 $
				    then begin
					(*u).mach.line.plot[i]=1
					wspec_calc_br,u,i
				    endif $
				    else (*u).mach.line.plot[i]=0
				   wspec_plot,u
				  endif
				  if ename eq strupcase('l'+is+'find') then begin
					wspec_plot_spec,u
					widget_control,(*u).id.message,set_value='Select Lower Wavelength',/append
					cursor,x1,y1,/up
					wait,0.15
					widget_control,(*u).id.message,set_value='Select Upper Wavelength',/append
					cursor,x2,y2,/up
					(*u).mach.line.lam0[i]=x1
					(*u).mach.line.lam1[i]=x2
					st=execute('widget_control, (*u).id.l'+is+'lam0,set_value=num2str((*u).mach.line.lam0['+is+'],dp=2)')
					st=execute('widget_control, (*u).id.l'+is+'lam1,set_value=num2str((*u).mach.line.lam1['+is+'],dp=2)')
					wspec_plot_time,u
					wspec_calc_br,u,i
					(*u).mach.line.plot[i]=1
					st=execute('widget_control,(*u).id.l'+is+'plot,set_button=1')
					wspec_plot,u
				  endif
				 endfor

                                END
			ENDCASE
		END
		"WIDGET_SLIDER": BEGIN
			widget_control,event.id,get_value=slider
			IF NOT keyword_set(uvalue) THEN uvalue='none'
			CASE ename OF 
				'T_SLIDER' : BEGIN
					IF (*u).stat.dat THEN BEGIN
						(*u).stat.time=slider/1.0e3
						widget_control, (*u).id.t_text,set_value=num2str((*u).stat.time,dp=3)
						if (*u).stat.tline then wspec_plot_time,u
						IF (*u).stat.pline THEN wspec_plot_views,u
						IF (*u).plot.isbr THEN wspec_plot_br,u
						wspec_plot_spec,u
					ENDIF
				END
				'SP_SLIDER' : BEGIN
					IF (*u).plot.issp THEN wspec_plot_sp,u		
				END
				ELSE : BEGIN
					for i=0,n((*u).mach.inst) do begin
					 if ename eq strupcase('sl'+(*u).mach.inst[i]) then begin
					  if slider eq -1 then begin
					   dim=size((*u).inst.(i).dch,/dim)
					    nch=dim[2]
					   (*u).inst.(i).d=total((*u).inst.(i).dch,3)/nch
					   (*u).inst.(i).l=total((*u).inst.(i).lch,2)/nch
					   st=execute('widget_control,(*u).id.lch'+(*u).mach.inst[i]+',set_value=''AVE''')
					  endif else begin
					   (*u).inst.(i).d=(*u).inst.(i).dch(*,*,slider)
					   (*u).inst.(i).l=(*u).inst.(i).lch(*,slider)
					   st=execute('widget_control,(*u).id.lch'+(*u).mach.inst[i]+',set_value=(*u).inst.(i).ch_nam['+num2str(slider,1)+']')
					  endelse
					  tmp=where((*u).mach.line.spec eq (*u).mach.inst[i])
					  ;if tmp[0] ne -1 then for j=0,n(tmp) do wspec_calc_br,u,tmp[j]
					  wspec_plot,u
					 endif
					endfor
				END
			ENDCASE
		END
		'WIDGET_DRAW' : begin ;process a button event in the draw window
		 redraw=0
		 IF event.id EQ (*u).id.draw1 OR event.id EQ (*u).id.br_draw THEN redraw=1
		 IF event.id EQ (*u).id.draw2 THEN redraw=2
		 case redraw of
		  1 : begin ; in time pane
	
		   widget_control,event.id,get_value=draw_win 
		   wset,draw_win & !X=(*u).plot.str1.X & !Y=(*u).plot.str1.Y & !Z=(*u).plot.str1.Z & !P=(*u).plot.str1.P
		   if (event.press ne 0) or (event.release ne 0) then begin ; ie not a motion event
		    click_loc = convert_coord(event.x, event.y, /device, /to_data)
		    ;print,click_loc
	       	    if event.press eq 1 then begin ; LMB click,click - zoom time L,R
		
		    ;also event.release
					if (*u).plot.zoom[0] eq 0 then begin
					 (*u).plot.zoom[0]=1
					 (*u).plot.x0[0]=click_loc[0]
					 widget_control, (*u).id.tau0,set_value=num2str((*u).plot.x0[0],dp=3)
					endif else begin
					 (*u).plot.zoom[0]=0
					 if click_loc[0] lt (*u).plot.x0[0] then begin
					  (*u).plot.x1[0]=(*u).plot.x0[0]
					  (*u).plot.x0[0]=click_loc[0]
					  widget_control, (*u).id.tau0,set_value=num2str((*u).plot.x0[0],dp=3)
					 endif else (*u).plot.x1[0]=click_loc[0]
					 widget_control, (*u).id.tau1,set_value=num2str((*u).plot.x1[0],dp=3)
					 wspec_calc_timetr,u
					 wspec_plot_time,u
					endelse
		    endif
		    if event.press eq 2 then begin ; MMB click - set spectrum time
					(*u).stat.time=click_loc[0]
					widget_control, (*u).id.t_slider,set_value=(*u).stat.time*1000
					widget_control, (*u).id.t_text,set_value=num2str((*u).stat.time,dp=3)
					wspec_plot,u
		    endif
		    if event.press eq 4 then begin ; RMB click - zoom out (ie zoom 0,2)
		    			 (*u).plot.zoom[0]=0
					 wspec_set_tzoomout,u
		    endif
		   endif else begin
		   	ptr_loc = convert_coord(event.x, event.y, /device,/to_data)
			widget_control, (*u).id.tim_ind,set_value='t='+num2str(ptr_loc[0],dp=3)+' y='+num2str(ptr_loc[1],dp=3)
		   endelse
		  end
		  2 : begin ; in wl pane
		   widget_control,(*u).id.draw2,get_value=draw_win
		   wset,draw_win & !X=(*u).plot.str2.X & !Y=(*u).plot.str2.Y & !Z=(*u).plot.str2.Z & !P=(*u).plot.str2.P
		   if (event.press ne 0) or (event.release ne 0) then begin ; ie not a motion event
		    click_loc = convert_coord(event.x, event.y, /device, /to_data)
	       	    if event.press eq 1 then begin ; LMB click,click - zoom time L,R
		    ;also event.release
					if (*u).plot.zoom[1] eq 0 then begin
					 (*u).plot.zoom[1]=1
					 (*u).plot.x0[1]=click_loc[0]
					 widget_control, (*u).id.lam0,set_value=num2str((*u).plot.x0[1],dp=3)
					endif else begin
					 (*u).plot.zoom[1]=0
					 if click_loc[0] lt (*u).plot.x0[1] then begin
					  (*u).plot.x1[1]=(*u).plot.x0[1]
					  (*u).plot.x0[1]=click_loc[0]
					  widget_control, (*u).id.lam0,set_value=num2str((*u).plot.x0[1],dp=3)
					 endif else (*u).plot.x1[1]=click_loc[0]
					 widget_control, (*u).id.lam1,set_value=num2str((*u).plot.x1[1],dp=3)
					 wspec_set_lzoom,u,[(*u).plot.x0[1],(*u).plot.x1[1]]
					endelse
		    endif
		    if event.press eq 4 then begin ; RMB click - zoom out (ie zoom 1,300)
		    			 (*u).plot.zoom[1]=0
					 wspec_set_lzoomout,u
		    endif
		   endif else begin ; pointer motion in pane
		   	ptr_loc = convert_coord(event.x, event.y, /device,/to_data)
			widget_control, (*u).id.br_ind,set_value='l='+num2str(ptr_loc[0],dp=3)+' y='+num2str(ptr_loc[1],dp=3)
		   endelse
		  end
		 endcase
		end
   		"WIDGET_TEXT_CH": BEGIN
			widget_control,event.id,get_value=data
			CASE event.id OF
				(*u).id.shotid : BEGIN
					(*u).shot=long(data)
					wspec_load_data,u
				END
				(*u).id.lam0 : BEGIN
					(*u).plot.x0[1]=float(data)
					wspec_plot_spec,u
				END
				(*u).id.lam1 : BEGIN
					(*u).plot.x1[1]=float(data)
					wspec_plot_spec,u
				END			
				(*u).id.sbr0 : BEGIN
					(*u).plot.yr[0]=float(data)
					wspec_plot_spec,u
                                END
				(*u).id.sbr1 : BEGIN
					(*u).plot.yr[1]=float(data)
					wspec_plot_spec,u
				END
				(*u).id.priore : BEGIN
					(*u).plot.priority=fix(data)
					wspec_plot_spec,u
				END
				(*u).id.tau0 : BEGIN
					(*u).plot.x0[0]=float(data)
					wspec_plot_time,u
				END
				(*u).id.tau1 : BEGIN
					(*u).plot.x1[0]=float(data)
					wspec_plot_time,u
				END	
				(*u).id.lelem : BEGIN
					(*(*u).dat.line).usrlist=data
					wspec_plot_spec,u
				END
				(*u).id.t_text : BEGIN
					(*u).stat.time=float(data)
					widget_control,(*u).id.t_slider,set_value=(*u).stat.time*1.0e3
					wspec_plot_time,u
					wspec_plot_spec,u
				END
				(*u).id.nistord : (*u).nist.order=float(data)
				(*u).id.nistspec : (*u).nist.spec=data
			ELSE: BEGIN
			 for i=0,n((*u).mach.line.nam) do begin
			  st=execute('iid=(*u).id.l'+num2str(i,1)+'lab')
			  if event.id eq iid then begin
				(*u).mach.line.nam[i]=data
				wspec_plot_time,u
			  endif
			  st=execute('iid=(*u).id.l'+num2str(i,1)+'mult')
			  if event.id eq iid then begin
				(*u).mach.line.mult[i]=float(data)
				wspec_calc_br,u,i
				wspec_plot,u
			  endif
			  st=execute('iid=(*u).id.l'+num2str(i,1)+'lam0')
			  if event.id eq iid then begin
				(*u).mach.line.lam0[i]=float(data)
				wspec_calc_br,u,i
				wspec_plot,u
			  endif
			  st=execute('iid=(*u).id.l'+num2str(i,1)+'lam1')
			  if event.id eq iid then begin
				(*u).mach.line.lam1[i]=float(data)
				wspec_calc_br,u,i
				wspec_plot,u
			  endif
			 endfor
			END
			ENDCASE
		END
		ELSE:
	ENDCASE
	IF button NE 'QUIT' THEN widget_control,event.top,set_uvalue=u		
END

;+
;NAME:
;	W_SPEC
;
;PURPOSE:
;	This procedure launches the W_SPEC widget which is used to
;	analyze the SXR/VUV spectroscopy data from HIREXSR, the
;	McPherson, XEUS and LoWEUS
;
;CALLING SEQUENCE:
;	@run_wspec.pro 		will compile the proper files, set the color
;				tables and launch the widget
;PROCEDURE:
;	See the C-Mod Wiki page on the W_SPEC widget for instructions
;	(listed under Shot Analysis Tools)
;	set keyword block to prevent returning to IDL prompt while running
;
;RESTRICTIONS:
;	Due to the limited namespace, many of the GENIE widgets have
;	overlapping function/procedure names so these widgets may not run
;	well in contiguous IDL sessions.
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 1/2012
;       Ported to general tokamak environment by A.N. James - 07/2012
;-

PRO w_spec,shot=shot,time=time, block=block
	; load machine configuration
	get_genie_env,MACHINE=mach
	st=execute('load_'+mach+',mach')
	
	user=logname()
	stdef=1150625023
	IF NOT keyword_set(shot) THEN st=execute('shot=load_'+mach.name+'_shot()')
	IF NOT keyword_set(time) THEN time=1.0
	IF shot EQ -1 THEN shot=stdef

	stat={time:time,dat:0,ps:0,col:!p.color,pscol:0,tline:1,pline:1}
	ysz=800
	str={X:!X,Y:!Y,Z:!Z,P:!P}
	plot={	asize:[400,ysz],$
		bsize:[ysz,ysz],$
		csize:[300,300*740.0/700],$
		brsize:[800,600],$
		isbr:0,$
		spsize:[600,1000],$
		issp:0,$
		yr:[0.0,1.0],$
		x0:[0.0,90.0],$ ;time, wavelength
		x1:[5.0,1200.0],$
		str1:str,str2:str, $ ; structures to store the coordinate values in a window
		zoom:[0,0],$
		nt:9,	nm:2,	auto:1,	label:1,priority:0,order:1,ntau:1000,log:0,$
		col:[!p.color,196,96],$
		pscol:[0,0,0]$
		}
	n_lines=float(n(mach.line.nam)+1)
	line={	yr:[0.0,1.0],$
		col:  1 + (findgen(n_lines))/n_lines*253,  $
		pscol:findgen(n_lines)/n_lines*255 $
		}


	base=widget_base(title='GENIE SPECTROSCOPY VIEWER',/row,/tlb_size_events)

	A=widget_base(base,/column)
	B=widget_base(base,/column)
	C=widget_base(base,/column)

	A0=widget_base(A,/row)
	dum=widget_label(A0,value='TIME EVOLUTION')
	tim_ind=widget_label(A0,value='                          ')
	A1=widget_base(A,frame=0)
	draw1=widget_draw(A1,xsize=plot.asize[0],ysize=plot.asize[1], /motion_events, /button_events)

	B0=widget_base(B,/row)
	dum=widget_label(B0,value='SPECTRAL BRIGHTNESS')
	br_ind=widget_label(B0,value='                          ')
	B1=widget_base(B,frame=0)
	draw2=widget_draw(B1,xsize=plot.bsize[0],ysize=plot.bsize[1], /motion_events, /button_events)
	B2=widget_base(B,/row)
	t_text=widget_text(B2,xsize=8,ysize=1,/edit)
	 widget_control, t_text,set_value=num2str(stat.time,dp=3)
	t_slider=widget_slider(B2,xsize=ysz-80,min=0,max=5000,value=1000,/drag,/suppress)

	C0=widget_base(C,/row)
	 dum=widget_label(C0,value='SETUP                     ')
	 stop= widget_button(C0,value='STOP')
	 print= widget_button(C0,value='PRINT')
	 help= widget_button(C0,value='HELP')
	 quit= widget_button(C0,value='QUIT')
	C1=widget_base(C,/column,xsize=320,ysize=ysz-plot.csize[1],frame=0)
	C1p0=widget_base(C1,/row)
	 setcurr= widget_button(C1p0,value='GET CURRENT')
	 dum = widget_label(C1p0,value='SHOT:')
	 shotid = widget_text(C1p0,xsize=10,ysize=1,/edit) & widget_control,shotid,set_value=num2str(shot,1)
	 load= widget_button(C1p0,value='LOAD')

	C1TA=widget_tab(C1,location=0) ; loc=0/1 -> tabs on top/bottom
	C1TA4=widget_base(C1TA,title='Message',/column)
	C1TA1=widget_base(C1TA,title='View',/column)
	C1TA2=widget_base(C1TA,title='Spectrometer',/column)
	C1TA3=widget_base(C1TA,title='Lines',/column)

	;C2=widget_base(C,/column,xsize=320,ysize=0.8*ysz/3,/frame)
	C1p3=widget_base(C1TA4,/column)
	 message=widget_text(C1p3,xsize=45,ysize=12,/scroll)
	  C1p3a=widget_base(C1p3,/row,/nonexcl)
          brwin=widget_button(C1p3a, value='EXPANDED BRIGHTNESS')
	  C1p3b=widget_base(C1p3,/row,/nonexcl)
          spwin=widget_button(C1p3b, value='SPATIAL PROFILES')
	  C1p3c=widget_base(C1p3,/row)
   	  dum=widget_label(C1p3c,value='LINE #')
	  sp_slider=widget_slider(C1p3c,xsize=100,/drag)
	  widget_control,sp_slider,set_slider_min=0,set_slider_max=n(mach.line.nam)

	C1p5=widget_base(C1TA1,/row)
	 dum=widget_label(C1p5,value='LAM RANGE:')
 	 lam0=widget_text(C1p5,xsize=6,ysize=1,/edit) & widget_control, lam0,set_value=num2str(plot.x0[1],dp=1)
	 dum=widget_label(C1p5,value='to')
 	 lam1=widget_text(C1p5,xsize=6,ysize=1,/edit) & widget_control, lam1,set_value=num2str(plot.x1[1],dp=1)
	 dum=widget_label(C1p5,value='[Ang]')
	 lzoom=widget_button(C1p5,value='ZOOM')
	 outz=widget_button(C1p5, value='OUT')

	C1p6=widget_base(C1TA1,/row)
	 dum=widget_label(C1p6,value='SPECBR RANGE:')
 	 sbr0=widget_text(C1p6,xsize=5,ysize=1,/edit) & widget_control, sbr0,set_value=num2str(plot.yr[0],dp=2)
	 dum=widget_label(C1p6,value='to')
 	 sbr1=widget_text(C1p6,xsize=5,ysize=1,/edit) & widget_control, sbr1,set_value=num2str(plot.yr[1],dp=2)
	 dum=widget_label(C1p6,value='[ph/s/m^2]')
	 C1p6=widget_base(C1TA1,/row)
	  dum=widget_label(C1p6,value='')
	  C1p6a=widget_base(C1p6,/row,/nonexcl)
          logsbr=widget_button(C1p6a, value='LOG Y')
	  autosbr=widget_button(C1p6a, value='AUTOSCALE Y') & widget_control, autosbr,set_button=plot.auto

	C1p7=widget_base(C1TA1,/row)
	 dum=widget_label(C1p7,value='TIME RANGE:')
 	 tau0=widget_text(C1p7,xsize=5,ysize=1,/edit) & widget_control, tau0,set_value=num2str(plot.x0[0],dp=1)
	 dum=widget_label(C1p7,value='to')
 	 tau1=widget_text(C1p7,xsize=5,ysize=1,/edit) & widget_control, tau1,set_value=num2str(plot.x1[0],dp=1)
	 dum=widget_label(C1p7,value='[sec]')
	 tzoom=widget_button(C1p7,value='ZOOM')
	 C1p7=widget_base(C1TA1,/row)
	  dum=widget_label(C1p7,value='')
	  C1p7a=widget_base(C1p7,/row,/nonexcl)
	  tline=widget_button(C1p7a,value='OPLOT TIME LINE') & widget_control, tline,set_button=stat.tline
	  pline=widget_button(C1p7a,value='OPLOT SPEC VIEWS') & widget_control, pline,set_button=stat.pline

	space=widget_base(C1TA1,/row)
	 dum=widget_label(space,value='')
	C1p8=widget_base(C1TA1,/row)
	 C1p8a=widget_base(C1p8,/row,/nonexcl)
	 label=widget_button(C1p8a,value='ELEM:') & widget_control, label,set_button=plot.label
	 lelem=widget_text(C1p8,xsize=38,ysize=1,/edit) & widget_control, lelem,set_value=mach.elem
	 C1p8=widget_base(C1TA1,/row)
	  C1p8a=widget_base(C1p8,/row,/nonexcl)
	  aelem=widget_button(C1p8a, value='ALL') ; buttons in nonexcl rows are checkboxes
	  oelem=widget_button(C1p8a, value='+2nd order')
	  dum=widget_label(C1p8,value='PRIORITY<=')
	  priore=widget_text(C1p8, xsize=6,ysize=1,/edit) & widget_control, priore,set_value=num2str(plot.priority,1)

	space=widget_base(C1TA1,/row)
	 dum=widget_label(space,value='')
	C1p9=widget_base(C1TA1,/row)
	 dum=widget_label(C1p9,value='NIST ASD:')
	 C1p9a=widget_base(C1p9,/row,/nonexcl)
	 nistplt=widget_button(C1p9a, value='plot')
	 dum=widget_label(C1p9,value='order:')
	 nistord=widget_text(C1p9,xsize=1,ysize=1,/edit) & widget_control, nistord,set_value='1'
	C1p9=widget_base(C1TA1,/row)
	 dum=widget_label(C1p9,value='spec:')
	 nistspec=widget_text(C1p9,xsize=5,ysize=1,/edit) & widget_control, nistspec,set_value=''
	 nistget=widget_button(C1p9, value='GET')

	
	C1p1=widget_base(C1TA2,/row)
	 C1p1a=widget_base(C1p1,/row)
	 dum=widget_label(C1p1a,value='INSTRUMENT    LOAD   PLOT')
	 for i=0,n(mach.inst) do begin
	  C1p1a=widget_base(C1TA2,/row)
	  lbl='            ' & strput,lbl,mach.inst[i],0
	  dum=widget_label(C1p1a,value=lbl)
	  C1p1b=widget_base(C1p1a,/row,/nonexcl)
	  st=execute('l'+mach.inst[i]+'=widget_button(C1p1b,value='''')')
	  st=execute('widget_control, l'+mach.inst[i]+',set_button=mach.load['+num2str(i,1)+']')
	  C1p1a=widget_base(C1p1a,/row)
	  dum=widget_label(C1p1a,value=' ')
	  C1p1b=widget_base(C1p1a,/row,/nonexcl)
	  st=execute('p'+mach.inst[i]+'=widget_button(C1p1b,value='''')')
	  st=execute('widget_control, p'+mach.inst[i]+',set_button=mach.plot['+num2str(i,1)+']')
	  st=execute('z'+mach.inst[i]+'=widget_button(C1p1a,value=''ZOOM'')')

	  if mach.nch[i] gt 1 then begin
	   C1p1c=widget_base(C1p1a,/column)
	   st=execute('sl'+mach.inst[i]+'=widget_slider(C1p1c,xsize=100,/drag,/suppress_value)')
	   st=execute('lch'+mach.inst[i]+'=widget_label(C1p1c,value=''                  '')')
	   st=execute('widget_control,sl'+mach.inst[i]+',set_slider_min=-1,set_slider_max=mach.nch[i]-1')
	  endif
	 endfor


	dum=widget_label(C1TA3,value='LINE BRIGHTNESS')
	 C2p=widget_base(C1TA3,/row)
	 dum=widget_label(C2p,value='     LINE      LAM0        LAM1         PLOT  MULT')
	 for i=0,n(mach.line.nam) do begin
	  C2px=widget_base(C1TA3,/row)
	  is=num2str(i,1)
	  st=execute('dum=widget_label(C2px,value='''+is+''')')
	  st=execute('l'+is+'lab=widget_text(C2px,xsize=8,ysize=1,/edit)')
	  st=execute('widget_control,l'+is+'lab,set_value=mach.line.nam['+is+']')
	  st=execute('l'+is+'lam0=widget_text(C2px,xsize=5,ysize=1,/edit)')
	  st=execute('widget_control,l'+is+'lam0,set_value=num2str(mach.line.lam0['+is+'],dp=2)')
	  dum=widget_label(C2px,value='to')
	  st=execute('l'+is+'lam1=widget_text(C2px,xsize=5,ysize=1,/edit)')
	  st=execute('widget_control,l'+is+'lam1,set_value=num2str(mach.line.lam1['+is+'],dp=2)')
	  st=execute('l'+is+'find=widget_button(C2px,value=''FIND'')')
	  C2pxx=widget_base(C2px,/row,/nonexcl)
	  st=execute('l'+is+'plot=widget_button(C2pxx,value='''')')
	  st=execute('widget_control,l'+is+'plot,set_button=mach.line.plot['+is+']')
	  st=execute('l'+is+'mult=widget_text(C2px,xsize=3,value=num2str(mach.line.mult['+is+'],1),/edit)')
	 endfor

	draw3=widget_draw(C,xsize=plot.csize[0],ysize=plot.csize[1])

	geo_C=widget_info(C, /geometry)
		x0=geo_C.xsize ; are these client size or screen size?
		y0=geo_C.ysize

	dstr=''
	for i=0,n(mach.inst) do begin
		dstr=dstr+'l'+mach.inst[i]+':l'+mach.inst[i] $
			+',p'+mach.inst[i]+':p'+mach.inst[i] $
			+',z'+mach.inst[i]+':z'+mach.inst[i]+','
		if mach.nch[i] gt 1 then begin
		 dstr=dstr+'lch'+mach.inst[i]+':lch'+mach.inst[i] $
			  +',sl'+mach.inst[i]+':sl' +mach.inst[i]+','
		endif
	endfor
	lstr=''
	for i=0,n(mach.line.nam) do begin
		li='l'+num2str(i,1)
		lstr=lstr +li+'lab:'+li+'lab,' +li+'lam0:'+li+'lam0,' +li+'lam1:'+li+'lam1,' $
		          +li+'find:'+li+'find,' +li+'mult:'+li+'mult,' +li+'plot:'+li+'plot,'
	endfor
	br_base=widget_base(title='BRIGHT',/column,group_leader=base,/tlb_size_events)
	br_draw=widget_draw(br_base,xsize=800,ysize=600,/button_events,/motion_events)
	sp_base=widget_base(title='SPATIAL PROFILES',/column,group_leader=base,/tlb_size_events)
	sp_draw=widget_draw(sp_base,xsize=600,ysize=1000,/button_events,/motion_events)
	str='id={base:base,t_slider:t_slider,t_text:t_text,' $
		+'A:A,A1:A1,B:B,B0:B0,B1:B1,B2:B2,C:C,x0:x0,y0:y0,draw1:draw1,draw2:draw2,draw3:draw3,' $ ; ids needed for resize
		+'shotid:shotid,load:load,quit:quit,print:print,stop:stop,help:help,message:message,setcurr:setcurr,' $
		+'br_ind:br_ind,tim_ind:tim_ind,' $
		+dstr +lstr $
		+'lam0:lam0,lam1:lam1,sbr0:sbr0, sbr1:sbr1,tau0:tau0, tau1:tau1,' $
		+'nistplt:nistplt,nistord:nistord,nistspec:nistspec,nistget:nistget,' $
		+'logsbr:logsbr,autosbr:autosbr,label:label,lelem:lelem,aelem:aelem,oelem:oelem,priore:priore,' $
		+'lzoom:lzoom,outz:outz,' $
		+'tline:tline,tzoom:tzoom,pline:pline,'$
		+'sp_base:sp_base,sp_draw:sp_draw,spwin:spwin,sp_slider:sp_slider,'$
		+'br_base:br_base,br_draw:br_draw,brwin:brwin}'
	st=execute(str)
	
	u=ptr_new(1,/allocate_heap)
	nist={spec:'',lam0:'',lam1:'',order:1,ln:ptr_new(1,/allocate_heap),s:0,p:0}
	*u={mach:mach,id:id,shot:shot,stat:stat,plot:plot,line:line,nist:nist}

	widget_control,base,/realize
	widget_control,br_base,/realize,map=0
	widget_control,sp_base,/realize,map=0
	wspec_load_data,u
	widget_control,base,set_uvalue=u
	widget_control,br_base,set_uvalue=u
	widget_control,sp_base,set_uvalue=u

	if not keyword_set(block) then block=0
	xmanager,'w_spec',base, event_handler='w_spec_event',no_block=(block eq 0)
	xmanager,'w_spec',br_base, event_handler='w_spec_event',no_block=(block eq 0)
END


;Main program executed when this file is compiled and run with .run
w_spec
END

