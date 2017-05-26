;+
;NAME:
;	GENSPEC_MOMENT_RANGE
;	
;PURPOSE:
;	This procedure is used to help find the best wavelength interval over which to perform the GENSPEC_LI_MOMENTS.
;	The moments
;
;CALLING SEQUENCE:
;	GENSPEC_MOMENT_RANGE,lam,int,lam_o,i_time
;
;INPUTS:
;	lam:	FLTARR [n_pts,n_ch,n_time] of the wavelength values [in lam_o units] for each channel at each time
;			Can also be [n_pts] too.
;	int:	FLTARR [n_pts,n_ch,n_time] of the intensity values [arbitrary "power" units] at those wavelength values
;	lam_o	FLT	of the unshifted wavelength of the line
;	i_time:	INT	of the point in n_time for which this procedure should be run
;
;KEYWORD PARAMETERS:
;	bsub		/bsub sends /backsub to GENSPEC_LI_MOMENTS
;
;OPTIONAL INPUTS:
;	n_r:	INT	number of points to for which to plot spectra
;	low	FLT	lam_o-low
;	up:	FLT	lam_o+up
;
;OUTPUTS:
;	Outputs are graphical to the currently selected plotting device.
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke - 7/19/07
;	3/6/08:		ML Reinke - added a bsub keyword and allowed lam to be 1D
;-

PRO genspec_moment_range,lam,int,lam_o,i_time,n_r=n_r,low=low,up=up,bsub=bsub
	
	IF NOT keyword_set(n_r) THEN n_r=8
	IF NOT keyword_set(low) THEN low=1.0
	IF NOT keyword_set(up) THEN up=5.0

	lr=[transpose(make(lam_o-up,lam_o-low,n_r)),transpose(make(lam_o+up,lam_o+low,n_r))]

        lam_back=lam
        int_back=int

        x=size(lam)
        IF x[0] EQ 3 THEN lam=reform(lam[*,*,i_time])
	int=reform(int[*,*,i_time])

	mom=genspec_li_moments(lam,int,lam_o,lr=lr[*,0])
	x=size(mom)
	moms=fltarr(3,x[2],n_r)
        chpt=x[2]/2.0
	moms[*,*,0]=mom
	FOR i=1,n_r-1 DO BEGIN
		mom=genspec_li_moments(lam,int,lam_o,debug=mom_debug,plot=mom_plot,lr=lr[*,i],backsub=bsub)
		moms[*,*,i]=mom
	ENDFOR

	colors=colormap(lr[0,*])
	openwin,0
	plot,moms[0,*,0],xtit='CH #',ytit='0th MOMENT',tit='Time Frame: '+num2str(i_time,1),chars=1.3,yr=[min(moms[0,*,*]),max(moms[0,*,*])]
	FOR i=1,n_r-1 DO BEGIN
		oplot, moms[0,*,i],color=colors[i]
		xyouts,chpt,max(moms[0,*,*])*(1.0-0.1*i),'['+num2str(lr[0,i],dp=2)+','+num2str(lr[1,i],dp=2)+']',chars=1.2,color=colors[i]
	ENDFOR
	xyouts,1,0.9*max(moms[0,*,*]),'line '+n2g('lambda')+'!lo!n='+num2str(lam_o,dp=2),chars=1.2
	openwin,1
	plot,moms[1,*,0],xtit='CH #',ytit='1st MOMENT',tit='Time Frame: '+num2str(i_time,1),chars=1.3,yr=[min(moms[1,*,*]),max(moms[1,*,*])]
	FOR i=1,n_r-1 DO BEGIN
		oplot, moms[1,*,i],color=colors[i]
		xyouts,chpt,max(moms[1,*,*])*(1.0-0.1*i),'['+num2str(lr[0,i],dp=2)+','+num2str(lr[1,i],dp=2)+']',chars=1.2,color=colors[i]
	ENDFOR
	xyouts,1,0.9*max(moms[1,*,*]),'line '+n2g('lambda')+'!lo!n='+num2str(lam_o,dp=2),chars=1.2

	openwin,2		
	plot,moms[2,*,0],xtit='CH #',ytit='2nd MOMENT',tit='Time Frame: '+num2str(i_time,1),chars=1.3,yr=[min(moms[2,*,*]),max(moms[2,*,*])]
	FOR i=1,n_r-1 DO BEGIN
		oplot, moms[2,*,i],color=colors[i]
		xyouts,chpt,max(moms[2,*,*])*(1.0-0.1*i),'['+num2str(lr[0,i],dp=2)+','+num2str(lr[1,i],dp=2)+']',chars=1.2,color=colors[i]
	ENDFOR
	xyouts,1,0.9*max(moms[2,*,*]),'line '+n2g('lambda')+'!lo!n='+num2str(lam_o,dp=2),chars=1.2

        lam=lam_back
        int=int_back
END
