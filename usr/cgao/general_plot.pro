;;;;;;;;;;;;;;;;;;;;;;
;GENERAL 2D PLOT
;;;;;;;;;;;;;;;;;;;;;;

;+
;NAME:
;	GENPLT
;
;PURPOSE:
;	This procedure provides a generic interface for creating x-window and postscript plots
;	of contours and contour slices.
;
;CALLING SEQUENCE
;	GENPLT,array,ivec,jvec
;
;INPUTS:
;	array:	FLTARR 2D (m x n) array of data that you want contour plotted
;	ivec:	FLTARR of m values that is the i scaling for array
;	jvec:	FLTARR of n values that is the j scaling for array
;
;OPTIONAL INPUTS:
;	io: 	FLTARR of values in ivec at which slices are taken (array[io,*]) [DEFAULT: 4 evenly spaced]
;	jo:	FLTARR of values in jvec at which slices are taken (array[*,jo]) [DEFAULT: 4 evenly spaced]
;	ir:	FLTARR [min_i,max_i] setting plot range for ivec [DEFAULT: full scale]
;	jr:	FLTARR [min_j,max_j] setting plot range for jvec [DEFAULT: full scale]
;	maxpt:	FLT maximum of plot range [DEFAULT: 1.1*max in array[ir,jr]]
;	minpt:	FLT minimum of plot range [DEFAULT: 0]
;	ncntrs:	FLT number of contours [DEFAULT: 18]
;	cct:	INT color table for the contour plot [DEFAULT: 12]
;	pct:	INT color table for the slice plots [DEFAULT: 12]
;	win:	INT plot to window group 0 or 1 [DEFAULT: 0]
;	dp:	INT number of decimal places in io/jo labels [DEFAULT: 2]
;	path:	STR path to put the plots if /ps is used [DEFAULT: /home/username/idl/plots]
;	suffix:	STR suffix to add to plot names [DEFAULT: empty]
;	prefix:	STR prefix to add to plot names [DEFAULT: empty]
;	labels:	STRUC structure of strings used to label the plots
;			*.ilab  - label for all ivec axes
;			*.jlab  - label for all jvec axes
;			*.klab  - label for all array axes
;			*.ctit  - title for contour plot
;			*.itit  - title for slices at io plot
;			*.jtit  - title for slices at jo plot
;
;KEYWORD PARAMETERS:
;	ps:	/ps allows plotting to PS device
;	debug:	/debug stops in places
;
;OPTIONAL OUTPUTS:
;	If a /ps is used three postscripts will be placed in the directory defined by path.
;		prefix+'genplt_cont'+suffix+'.ps'
;		prefix+'genplt_io'+suffix+'.ps'
;		prefix+'genplt_jo'+suffix+'.ps'
;
;	NOTE: The PS device must be properly configured be fore calling GENPLT
;
;SIDE EFFECTS:
;	Color tables will be hijacked for your IDL session if cct and pct aren't used.
;
;RESTRICTIONS:
;	This code uses lots of functions in MLR_FUNCTIONS
;	(/home/mlreinke/idl/general/mlr_functions)	
;	
;MODIFICATION HISTORY:
;	Written by:	ML Reinke 5-19-05 (update on my previous GENPLT which sucked)
;	6-26-06:	ML Reinke: Set io/jo plots to use inclusive boundries when determining maxplot
;				   Tweaked the contour plot parameters for a proper fit when psplot 
;				   is used to setup the PS device and /ps is sent to GENPLT
;	6-20-06:	ML Reinke: Set the io/jo plots to use /xsty and /ysty and change ymax to 1.1 the maxplot
;	11-1-06:	ML Reinke: Modified so that numbers less than zero are shown on the contour/slice plots
;	1-26-11:	ML Reinke: adapted to allow pixmap features to be used for better widget plotting
;-      9-17-14:        C Gao: enabled user to specify windows set, instead of [13,14,15] or [16,17,18]

PRO genplt,array,ivec,jvec,io=io,jo=jo,ir=ir,jr=jr,maxpt=maxpt,minpt=minpt,ncntrs=ncntrs,labels=labels,$
		path=path,suffix=suffix,prefix=prefix,dp=dp,chars=chars,debug=debug,cct=cct,pct=pct,ps=ps,$
		win=win,pixmap=pixmap,set=set
	x=size(array)
	IF NOT keyword_set(ivec) THEN ivec=indgen(x[1])
	IF NOT keyword_set(jvec) THEN jvec=indgen(x[2])
	
	;setup windowing
        IF not keyword_set(set) THEN BEGIN
           IF keyword_set(win) THEN BEGIN
              IF win[0] EQ 0 THEN BEGIN
                 reset=8
                 set=[13,14,15]
              ENDIF
              IF win[0] EQ 1 THEN BEGIN
                 reset=9
                 set=[16,17,18]
              ENDIF
              IF keyword_set(pixmap) THEN BEGIN
                 reset=0
                 set=[13,14,15]
              ENDIF
           ENDIF ELSE BEGIN
              reset=8
              set=[13,14,15]
           ENDELSE
        ENDIF ELSE BEGIN
           reset=8
        ENDELSE
    

	;setup path variables
	IF NOT keyword_set(path) THEN path='/home/'+logname()+'/idl/plots/'
	IF NOT keyword_set(suffix) THEN suffix=''
	IF NOT keyword_set(prefix) THEN prefix=''	

	;set plotting variables
	IF NOT keyword_set(ir) THEN ir=[min(ivec), max(ivec)]
	IF NOT keyword_set(jr) THEN jr=[min(jvec), max(jvec)]
	IF NOT keyword_set(ncntrs) THEN  ncntrs=18.0
	extremepts=localmaxmin(array,ivec,jvec,ir,jr)
	IF NOT keyword_set(maxpt) THEN  maxplot=extremepts[1]*1.1 ELSE maxplot=maxpt
	IF NOT keyword_set(minpt) THEN minplot=0.0 < extremepts[0]*1.1 ELSE minplot=minpt
	IF NOT keyword_set(io) THEN  BEGIN
		io=make(ir[0],ir[1],6)
		io=io[1:4]
	ENDIF
	IF NOT keyword_set(jo) THEN BEGIN
		jo=make(jr[0],jr[1],6)
		jo=jo[1:4]
	ENDIF
	IF NOT keyword_set(dp) THEN dp=2
	IF NOT keyword_set(chars) THEN chars=1.0
	IF NOT keyword_set(cct) THEN cct=12
	IF NOT keyword_set(pct) THEN pct=12

	;correnct selection of i,j ranges outside of ivec,jvec
	IF (ir[0] LT min(ivec)) THEN ir[0]=min(ivec)
	IF (ir[1] GT max(ivec)) THEN ir[1]=max(ivec)
	IF (jr[0] LT min(jvec)) THEN jr[0]=min(jvec)
	IF (jr[1] GT max(jvec)) THEN jr[1]=max(jvec)

	IF NOT keyword_set(ps) THEN BEGIN
		IF reset NE 0 THEN wreset,reset
		IF keyword_set(pixmap) THEN window,set[0],xsize=pixmap[0],ysize=pixmap[1],/pixmap ELSE wset,set[0]
	ENDIF

	;this is weird, but it seems to be necessary, to get good contours, must truncate the emiss array to the plot size
	loadct,cct,/silent
	tmp1=where(jvec GE jr[0] AND jvec LE jr[1])
	tmp2=where(ivec GE ir[0] AND ivec LE ir[1])
	array_temp1=array[*,tmp1]
	array_plot=array_temp1[tmp2,*]

	;setup titles for contour
	IF keyword_set(labels) THEN BEGIN
		xtit=labels.ilab
		ytit=labels.jlab
		tit=labels.ctit
		ztit=labels.klab
	ENDIF ELSE BEGIN
		xtit='i-vec'
		ytit='j-vec'
		tit='contour plot'
		ztit='array'
	ENDELSE
	;levels=(maxplot-minplot)*findgen(ncntrs)/ncntrs
	levels=make(minplot,maxplot,ncntrs)
	IF keyword_set(ps) THEN position=[0.1,0.08,0.85,0.95] ELSE position=[0.085,0.08,0.825,0.95]	

	contour,array_plot,ivec[tmp2],jvec[tmp1],min=minplot,max=maxplot,fill=1,levels=levels,$
		xtit=xtit,ytit=ytit,tit=tit,xrange=ir,yrange=jr,xstyle=1,ystyle=1,position=position,chars=chars

	cbar=fltarr(2,n(levels)+1)
	cbar[0,*]=levels
	cbar[1,*]=levels
	IF keyword_set(ps) THEN position=[0.95,0.08,0.99,0.95] ELSE position=[0.92,0.08,0.97,0.95]
	contour,cbar,[0,1],levels,fill=1,levels=levels,ytit=ztit,position=position,/noerase,chars=chars,$
		/xsty,/ysty,xticklayout=0,xticks=1,xtickname=[' ',' ']

	;close ps plot if chosen
	IF keyword_set(ps) THEN BEGIN
		device,/close
		filename=prefix+'genplt_cont'+suffix+'.ps'
		spawn, 'cp idl.ps '+path+filename
	ENDIF

	IF keyword_set(debug) THEN stop
	
	
	;plot the chosen or deafult i-coords
	loadct,pct,/silent					
	iplots=fltarr(n(io)+1,n(jvec)+1)
	FOR k=0,n(io) DO iplots[k,*]=array_slice(array,ivec,jvec,i=io[k])
	colormap=colormap(io)
	IF NOT keyword_set(ps) THEN BEGIN
		IF keyword_set(pixmap) THEN window,set[1],/pixmap,xsize=pixmap[2],ysize=pixmap[3] ELSE wset,set[1]
	ENDIF
	IF keyword_set(maxpt) THEN maxplot=maxpt ELSE maxplot=max(iplots[*,where(jvec GE jr[0] AND jvec LE jr[1])])
	IF keyword_set(labels) THEN BEGIN
		xtit=labels.jlab
		ytit=labels.klab
		tit=labels.itit
		note=labels.ilab
	ENDIF ELSE BEGIN
		xtit='j-vec'
		ytit='array'
		tit='array @ various io'
		note='io'
	ENDELSE
	
	IF keyword_set(debug) THEN stop
	plot,[0],[0],xtit=xtit,ytit=ytit,tit=tit,xrange=jr,yrange=[minplot,maxplot*1.1],chars=chars,/xsty,/ysty
	xyouts,0.17,0.85,note, /norm
	FOR i=0,n(io) DO BEGIN
		oplot,jvec, iplots[i,*],color=colormap[i]
		xyouts,0.2,0.80-0.05*i,num2str(io[i],dp=dp),color=colormap[i], /norm
	ENDFOR
	IF minplot LT 0 THEN oplot,jr,[0,0],linestyle=2

	;close ps plot if chosen
	IF keyword_set(ps) THEN BEGIN
		device,/close
		filename=prefix+'genplt_io'+suffix+'.ps'
		spawn, 'cp idl.ps '+path+filename
	ENDIF

	;plot the chosen or deafult j-coords					
	jplots=fltarr(n(jo)+1,n(ivec)+1)
	FOR k=0,n(jo) DO jplots[k,*]=array_slice(array,ivec,jvec,j=jo[k])
	colormap=colormap(jo)
	IF NOT keyword_set(ps) THEN BEGIN
			IF keyword_set(pixmap) THEN window,set[2],/pixmap,xsize=pixmap[4],ysize=pixmap[5] ELSE wset,set[2]
	ENDIF
	IF keyword_set(maxpt) THEN maxplot=maxpt ELSE maxplot=max(jplots[*,where(ivec GE ir[0] AND ivec LE ir[1])])
	IF keyword_set(labels) THEN BEGIN
		xtit=labels.ilab
		ytit=labels.klab
		tit=labels.jtit
		note=labels.jlab
	ENDIF ELSE BEGIN
		xtit='i-vec'
		ytit='array'
		tit='array @ various jo'
		note='jo'
	ENDELSE
	

	plot,[0],[0],xtit=xtit,ytit=ytit,tit=tit,xrange=ir,yrange=[minplot,maxplot*1.1],chars=chars,/xsty,/ysty
	xyouts,0.17,0.85,note, /norm
	FOR i=0,n(jo) DO BEGIN
		oplot,ivec, jplots[i,*],color=colormap[i]
		xyouts,0.2,0.80-0.05*i,num2str(jo[i],dp=dp),color=colormap[i], /norm
	ENDFOR
	IF minplot LT 0 THEN oplot,ir,[0,0],linestyle=2

	;close ps plot if chosen
	IF keyword_set(ps) THEN BEGIN
		device,/close
		filename=prefix+'genplt_jo'+suffix+'.ps'
		spawn, 'cp idl.ps '+path+filename
	ENDIF

	IF keyword_set(debug) THEN stop

END




;taken from Martin's W_FLUCT Widget so I can use arbitrary signals
PRO transform,data,time,ft,t,f,samples=samples,nf=nf,nt=nt
	nraw=n(data)+1
	dt=time[1]-time[0]
	IF NOT keyword_set(samples) THEN num = 1024 ELSE num=samples
	nspec = nraw/num
	f = findgen(num)/dt/(num-1)
	t = fltarr(nspec)
	ft = complexarr(nspec,num)

	FOR i = 0L,nspec-1 do begin
  		temp = fft(data[long(i)*num:long(i+1)*num-1],-1)
  		ft[i,*] = temp
 		t[i] = time[long(i+.5)*num]
	ENDFOR
	IF NOT keyword_set(nf) THEN nf=1
	IF NOT keyword_set(nt) THEN nt=1
	IF nt*nf GT 1 THEN ft=smooth(abs(ft),[nt,nf]) ELSE ft=abs(ft)	

END
