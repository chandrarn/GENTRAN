
;+
;NAME:
;	MLR_FUNCTIONS
;
;VERSION:
;	Version 2.5
;	1-24-07
;-



;+
;NAME:
;	N
;
;PURPOSE:
;	This function is simply a space saver for FOR loops
;
;OUTPUT:
;	n_elements(x)-1
;-

FUNCTION  n,x					;just a space saver for invoking (n_elements-1) when using FOR statements
	output=n_elements(x)-1
	RETURN, output
END

FUNCTION last,vec
	output=vec[n_elements(vec)-1]
	RETURN,output
END

;+
;NAME:
;	SAMEAN
;
;PURPOSE:
;	This function calculates the mean of 'y' in between the xvalues (x0,x1)
;
;CALLING SEQUENCE:
;	result=SAMEAN(y,x,x0,x1)
;
;-

FUNCTION samean,y,x,a,b				;computes and returns the mean of a sub-array of 'y' between the 'x' values of (a,b)
	del_t=abs(x[1]-x[0])
	base_array=extrac(y, (a-x[0])/del_t, abs(b-a)/del_t)
	RETURN, mean(base_array)
END

FUNCTION sasum,y,x,a,b				;computes and returns the sum of sub-array
	del_t=abs(x[1]-x[0])
	base_array=extrac(y, (a-x[0])/del_t, abs(b-a)/del_t)
	output=total(base_array)
	RETURN, output
END

FUNCTION saminmax,y,x,a,b			;computes and returns both the max and min values in sub-array
	del_t=abs(x[1]-x[0])
	base_array=extrac(y, (a-x[0])/del_t, abs(b-a)/del_t)
	min=min(base_array)
	max=max(base_array)
	output=[min,max]
	RETURN, output
END

FUNCTION saint,y,x,a,b				;computes and returns both the the integral over the sub_array
	del_t=abs(x[5]-x[4])
	base_y=extrac(y, (a-x[0])/del_t, abs(b-a)/del_t)
	base_x=extrac(x, (a-x[0])/del_t, abs(b-a)/del_t)
	print, base_x
	output = int_tabulated(base_x,base_y)
	RETURN, output
END

FUNCTION sum,y					;same as total(y) but kept incase of needed redundancies
	output=0
	for i=0,n(y)do  output=output+y[i]
	RETURN, output
END


;+
;NAME:
;	FIXARRAY
;
;PURPOSE:
;	This function reforms the input array including only non-zero elements (rows/columns for 2D array)
;	of the vector good.  If array is 2D then the length of good should match one dimension of array.
;
;CALLING SEQUENCE:
;	result=FIXARRAY(array,good)
;
;RESTRICTIONS
;	array cannot be 2D and (n x n)
;
;-

FUNCTION fixarray,array,good			;will reform the input ARRAY so that only the values where GOOD not equal to zero	
	x=size(array)
	
	IF (x[0] EQ 1) THEN BEGIN
		fixedarray=fltarr(total(good))
		cntr=0
		FOR i=0,n(array) DO BEGIN
			IF (good[i] NE 0) THEN BEGIN
				fixedarray[cntr]=array[i]
				cntr+=1
			ENDIF
		ENDFOR
		RETURN, fixedarray
	ENDIF

	IF (x[0] EQ 2) THEN BEGIN
		IF (x[1] EQ x[2]) THEN BEGIN
			print, 'Cannot determine which dimension to reform'
			RETURN, -1
		ENDIF

		IF (x[1] EQ n(good)+1) THEN BEGIN
			fixedarray=fltarr(total(good), x[2])
			cntr=0
			FOR i=0,x[1]-1 DO BEGIN
				IF (good[i] NE 0) THEN BEGIN
					fixedarray[cntr,*]=array[i,*]
					cntr+=1
				ENDIF
			ENDFOR	

		ENDIF ELSE BEGIN
			fixedarray=fltarr(x[1],total(good))
			cntr=0
			FOR i=0,x[2]-1 DO BEGIN
				IF (good[i] NE 0) THEN BEGIN
					fixedarray[*,cntr]=array[*,i]
					cntr+=1
				ENDIF
			ENDFOR	

		ENDELSE
		
		RETURN, fixedarray
	ENDIF
	
	RETURN, -1
	
END

;+
;NAME:
;	WIN_SAVE
;
;PURPOSE:
;	This procedure is for times when you need to quickly output a plot that cannot be
;	output again to the PS device.  The quaility is poor since it TVRD's to a JPEG.
;
;CALLING SEQUENCE:
;	WIN_SAVE,filename_string
;
;RESTRICTIONS
;	Operates on currently set window
;
;-

PRO win_save, name				;saves the current set window as a JPEG using the string 'name'			
	image=tvrd(true=3)
	write_jpeg, name+'.jpg', image, true=3,quality=100
;	command='jpeg2ps '+name+'.jpg > '+name+'.ps'
;	spawn, command, result,err_result
;	print, err_result
;	command='rm '+name+'.jpg'
;	spawn, command,result,err_result
;	print, err_result
END

;+
;NAME:
;	SUB_BL
;
;PURPOSE:
;	This procedure calculates the mean,m, of 'y' in between the xvalues (x0,x1) and returns
;	writes over the 'y' with y - m.
;
;CALLING SEQUENCE:
;	SUB_BL,y,x,x0,x1
;
;-

PRO sub_bl,y,x,a,b				;performs a DC subtraction on 'y' from a level computed between values (a,b) of 'x'
	del_t=abs(x[1]-x[0])			
	base_array=extrac(y, (a-x[0])/del_t, abs(b-a)/del_t)
	y-=mean(base_array)
END

PRO saextrac,y,x,a,b				;extracts subarrays of 'y' and 'x' between the values of (a,b) of 'x'
	del_t=abs(x[5]-x[4])
	base_y=extrac(y, (a-x[0])/del_t, abs(b-a)/del_t)
	base_x=extrac(x, (a-x[0])/del_t, abs(b-a)/del_t)
	y=base_y
	x=base_x
END

;+
;NAME:
;	NUM2STR
;
;PURPOSE:
;	This function converts an input FLT,INT,LONINT to a string.  Really this program sucks, but I want to fucking
;	strangle RSI for not having a native fully functional num2str() function.  They're bastard coated
;	bastards with bastard filling.  Really, I mean whats their fuckin' problem.  You think I'm over reacting,
;	but just try and put string(n) into something and notice how fucked up it is.
;
;CALLING SEQUENCE:
;	result=NUM2STR(n)
;
;OPTIONAL INPUTS
;	raw:	Use num2str(n,raw) to convert 'n' directly w/o leading zeros (use for shot numbers)
;	dp:	dp = INT to set the number of decimal points 
;
;KEYWORD PARMAMETERS:
;	verb:	/verb to print status messages (used for debugging)
;	sn:	/sn forces scientific notation for numbers below IDL's limit (1.0e6)
;
;RESTRICTIONS
;	Can only have 8 total digits in a FLT
;
;MODIFICATION HISTORY
;	Written by: 	ML Reinke
;	8-04-05		ML Reinke - finally fixed trailing zeros issue and added a dp flag to choose the # decimal points
;	6-13-06		ML Reinke - added the ability to handle floats above 1.0e6 and force scientific notation
;-

FUNCTION num2str,n,raw,verb=verb,dp=dp,sn=sn
	
	;limitations only 8 total numbers (not including a leading zero on #<1, ie 0.013
	;unless raw is envoked
	n_old=n
	
	tmp=size(n,/type)
	IF tmp[0] EQ 7 THEN RETURN, n ; already a string!
	IF (tmp[0] EQ 4 OR tmp[0] EQ 5) AND max(n) GE 1.0e6 THEN BEGIN
		str=string(n)
		tmp=strsplit(str, '+',/extrac)
		expon=int(tmp[1])
		tmp=strsplit(tmp[0],'e',/extrac)
		n=tmp[0]
	ENDIF ELSE expon=0

	IF keyword_set(sn) AND expon EQ 0 THEN BEGIN
		i=1
		WHILE n/10.0^i GE 1.0 AND i LT 6 DO i+=1
		expon=i-1
		n=n/10.0^(i-1)
	ENDIF
		
	
	;find the number of decimal points,
	IF (n(raw) GE 0) THEN dec=0 ELSE dec=-1
	i=0
	WHILE (dec LT 0 AND i LT 10) DO BEGIN
		IF (floor(abs(n*10.0^i)) EQ abs(n*10.0^i)) THEN dec=i
		i+=1
		IF keyword_set(verb) THEN print, floor(abs(n*10.0^i)), abs(n*10.0^i)
	ENDWHILE
	
	;find the number of figures
	IF (n(raw) GE 0) THEN fig=10 ELSE fig=-1
	i=0
	WHILE (fig LT 0 AND i LT 10) DO BEGIN
		IF (floor(abs(n/10.0^i)) EQ 0) THEN fig=i
		i+=1
		IF keyword_set(verb) THEN print, floor(abs(n/10.0^i))
	ENDWHILE

	IF keyword_set(dp) THEN dec=dp
	intype=size(n,/type)
	IF intype EQ 2 THEN raw=1

	IF keyword_set(verb) THEN BEGIN
		print, 'FIG: ',fig 
		print, 'DEC: ',dec
	ENDIF
	IF (dec GE 0 AND n(raw) LT 0) THEN BEGIN
		IF dec EQ 0 THEN dec+=1
		IF (fig EQ 0) THEN fig+=1 
		IF (abs(n) NE n) THEN fig+=1
		format_string='(F'+strtrim(string(fig+dec+1),1)+'.'+strtrim(string(dec),1)+')'
		tmp=strtrim(string(n,format=format_string),1)
		IF keyword_set(verb) THEN print, format_string
	ENDIF ELSE BEGIN	
		tmp=strtrim(string(n),1)
	ENDELSE
	
	IF expon LT 10 THEN extra='0' ELSE extra=''
	IF expon GT 0 THEN tmp+='e+'+extra+strtrim(expon,1)
	
	;tmp=strtrim(string(n*10^dec),1)
	;if (abs(n) NE n) THEN fig+=1
	;fig_string=strmid(tmp,0, fig)
	;dec_string=strmid(tmp,fig,dec)
	;IF (dec EQ 0) THEN num_string=fig_string ELSE num_string=fig_string+'.'+dec_string
	
	n=n_old
	num_string=tmp
	RETURN,num_string
END


;+
;NAME:
;	WRESET
;
;PURPOSE:
;	This procedure sets up various window sets used for DATAPLOTS and other MLR plotting routines.
;	If the current window set is already open, it does nothing.  Give it a try to see what happens!
;	Windows are labeled with what they are generally used for as well as the window number for WSET refernce.
;
;CALLING SEQUENCE:
;	WRESET,n
;
;KEYWORD PARAMETERS:
;	FORCE:		Use a /force to force a reset of the window set if they are already open
;
;-

PRO wreset,n,force=force,pixmap=pixmap
	device, window_state=var
	IF keyword_set(force) THEN var=intarr(n(var)+1)

	IF (n EQ 0 AND total(var[0:3]) NE 4) THEN BEGIN
		xo=2500
		window,0,xsize=650,ysize=512,xpos=xo,ypos=670, title='foil em,0'
		window,1,xsize=650,ysize=512,xpos=xo,ypos=100, title='foil br,1'
		window,2,xsize=650,ysize=512,xpos=xo-700,ypos=100, title='axuv br,2'
		window,3,xsize=650,ysize=512,xpos=xo-700,ypos=670, title='axuv em,3'
	ENDIF

	IF (n EQ 1 AND total(var[0:3]) NE 4) THEN BEGIN
		xo=900
		window,0,xsize=650,ysize=512,xpos=xo,ypos=670
		window,1,xsize=650,ysize=512,xpos=xo,ypos=100
		window,2,xsize=650,ysize=512,xpos=xo-700,ypos=100
		window,3,xsize=650,ysize=512,xpos=xo-700,ypos=670
	ENDIF

	IF (n EQ 2 AND total(var[4:5]) NE 2) THEN BEGIN
		window,4,xsize=1025,ysize=830,xpos=2150,ypos=85, title='em time,4'
		window,5,xsize=750,ysize=750,xpos=1610,ypos=410, title='em contour,5'
	ENDIF

	IF (n EQ 3 AND total(var[7:8]) NE 2) THEN BEGIN
		window,7,xsize=1025,ysize=830,xpos=1650,ypos=85, title='emZ time,7'
		window,8,xsize=750,ysize=750,xpos=2375,ypos=410, title='emZ controur,8'
	ENDIF

	IF (n EQ 4 AND  total(var[9:10]) NE 2) THEN BEGIN
		window,9,xsize=1025,ysize=830,xpos=1650,ypos=85, title='emf time,10'
		window,10,xsize=750,ysize=750,xpos=2375,ypos=410, title='emf controur,11'
	ENDIF

	IF (n EQ 5 AND  total(var[4:6]) NE 3) THEN BEGIN
		window,4,xsize=500,ysize=250,xpos=2140,ypos=370, title='em time,4'
		window,5,xsize=500,ysize=550,xpos=2140,ypos=760, title='em contour,5'
		window,6,xsize=500,ysize=250,xpos=2140,ypos=85, title='em radii, 6'
	ENDIF
	
	IF (n EQ 6 AND total(var[7:9]) NE 3) THEN BEGIN
		window,7,xsize=500,ysize=250,xpos=2675,ypos=370, title='emz time,7'
		window,8,xsize=500,ysize=550,xpos=2675,ypos=760, title='emz contour,8'
		window,9,xsize=500,ysize=250,xpos=2675,ypos=85, title='emz radii, 9'
	ENDIF

	IF (n EQ 7 AND total(var[10:12]) NE 3) THEN BEGIN
		window,10,xsize=500,ysize=250,xpos=1610,ypos=370, title='emf time,10'
		window,11,xsize=500,ysize=550,xpos=1610,ypos=760, title='emf contour,11'
		window,12,xsize=500,ysize=250,xpos=1610,ypos=85, title='emf radii, 12'
	ENDIF

	IF (n EQ 8 AND total(var[13:15]) NE 3) THEN BEGIN
		window,13,xsize=700,ysize=540,xpos=1610,ypos=760, title='contour,13',pixmap=pixmap
		window,15,xsize=600,ysize=255,xpos=2340,ypos=1000, title='slices@jo,15',pixmap=pixmap
		window,14,xsize=600,ysize=255,xpos=2340,ypos=660, title='slices@io 14',pixmap=pixmap
	ENDIF
	IF (n EQ 9 AND total(var[16:18]) NE 3) THEN BEGIN
		window,16,xsize=700,ysize=540,xpos=1610,ypos=85, title='contour,16'
		window,18,xsize=600,ysize=255,xpos=2340,ypos=370, title='slices@jo,18'
		window,17,xsize=600,ysize=255,xpos=2340,ypos=85, title='slices@io 17'
	ENDIF

	IF (n EQ 10 AND  total(var[20]) NE 1) THEN BEGIN
		window,20,xsize=1025,ysize=830,xpos=550,ypos=285, title='CoolCurve,20'
	ENDIF

	IF (n EQ 19 AND total(var[19]) NE 1) THEN BEGIN	
		window,19,xsize=500,ysize=690,xpos=0,ypos=670,title='vessel cx,19'
	ENDIF

	IF (n EQ 21 AND total(var[21]) NE 1) THEN BEGIN	
		window,21,xsize=500,ysize=690,xpos=0,ypos=670,title='vessel cx,21'
	ENDIF
END

;+
;NAME:
;	WCLOSE
;
;PURPOSE:
;	This closes the given window set.
;
;CALLING SEQUENCE:
;	WCLOSE,n
;
;-

PRO wclose,n
	IF (n EQ -1) THEN BEGIN
		device,window_state=var
		tmp=where(var NE 0)
		FOR i=0,n(tmp) DO err=execute('wdelete,'+num2str(tmp[i],1))
	ENDIF

	IF (n EQ 0) THEN BEGIN
		wreset,0
		wdelete,0
		wdelete,1
		wdelete,2
		wdelete,3
	ENDIF

	IF (n EQ 1) THEN BEGIN
		wreset,1
		wdelete,0
		wdelete,1
		wdelete,2
		wdelete,3
	ENDIF

	IF (n EQ 2) THEN BEGIN
		wreset,2
		wdelete,4
		wdelete,5
	ENDIF

	IF (n EQ 3) THEN BEGIN
		wreset,3
		wdelete,7
		wdelete,8
	ENDIF

	IF (n EQ 4) THEN BEGIN
		wreset,4
		wdelete,9
		wdelete,10
	ENDIF

	IF (n EQ 5) THEN BEGIN
		wreset,5
		wdelete,4
		wdelete,5
		wdelete,6
	ENDIF

	IF (n EQ 6) THEN BEGIN
		wreset,6
		wdelete,7
		wdelete,8
		wdelete,9	
	ENDIF

	IF (n EQ 7) THEN BEGIN
		wreset,7
		wdelete,10
		wdelete,11		
		wdelete,12
	ENDIF

	IF (n EQ 8) THEN BEGIN
		wreset,8
		wdelete,13
		wdelete,15		
		wdelete,14
	ENDIF

	IF (n EQ 9) THEN BEGIN
		wreset,9
		wdelete,16
		wdelete,18		
		wdelete,17
	ENDIF

	IF (n EQ 10) THEN BEGIN
		wreset,10
		wdelete,20
	ENDIF

	IF (n EQ 19) THEN BEGIN
		wreset,19
		wdelete,19
	ENDIF

	IF (n EQ 21) THEN BEGIN
		wreset,21
		wdelete,21
	ENDIF

END



;+
;NAME:
;	N2G
;	
;PURPOSE:
;	This is a handy little function for easily inserting greek letters into plots without having
;	to memorize formatting codes that are device dependent.
;
;CALLING SEQUENCE
;	result=N2G(str)
;	
;INPUTS:
;	str:		STR of greek letter name (ie 'alpha', 'beta')
;			Capitalize the first letter for capital greek letters
;
;OPTIONAL INPUTS:
;	font:	INT of the font number being used that's not DEFAULT = 3
;
;OUPUTS:
;	result: 	STR of the format code to put the greek letter in a string that
;			will be sent to the postscript or x-windows device	
;
;EXAMPLE:
;	The proper way to insert a greek letter in a plot title would be:
;		title='!7t!3 vs. Radiated Power' 
;	where !7 selects the symbol font, tau=t and !3 returns to Helvetica.
;	The problem comes in when !9 is needed when using the PS device and the mapping of roman letter
;	to greek letter changes.  So N2G keeps track of all this and let's the user just say what they mean.
;		title=n2g('tau')+' vs. Radiated Power'
;	The formatting will change if the PS device is active and include the font changing commands.
;
;MODIFICATION HISTORY:
;	Written by: 	ML Reinke  - 4/16/07
;
;-

FUNCTION n2g,str,ps=ps,font=font
	IF NOT keyword_set(font) THEN f_old='3' ELSE f_old=num2str(font,1)
	IF !d.name EQ 'PS' THEN ps=1 ELSE ps = 0
	var=''
	IF str EQ strlowcase(str) THEN BEGIN
		IF keyword_set(ps) THEN BEGIN
			CASE str OF
				'alpha': 	var='a'
				'beta': 	var='b'
				'gamma':	var='g'
				'delta':	var='d'	
				'epsilon':	var='e'
				'zeta':		var='z'
				'eta':		var='h'
				'theta':	var='q'
				'iota':		var='i'
				'kappa':	var='k'
				'lambda':	var='l'
				'mu':		var='m'
				'nu':		var='n'
				'xi':		var='x'
				'omicron':	var='o'
				'pi':		var='p'
				'rho':		var='r'
				'sigma':	var='s'
				'tau':		var='t'
				'upsilon':	var='u'
				'phi':		var='f'
				'chi':		var='c'
				'psi':		var='y'
				'omega':	var='w'
				else:		var=''
			ENDCASE
		ENDIF ELSE BEGIN
			CASE str OF
				'alpha': 	var='a'
				'beta': 	var='b'
				'gamma':	var='c'
				'delta':	var='d'
				'epsilon':	var='e'
				'zeta':		var='f'
				'eta':		var='g'
				'theta':	var='h'
				'iota':		var='i'
				'kappa':	var='j'
				'lambda':	var='k'
				'mu':		var='l'
				'nu':		var='m'
				'xi':		var='n'
				'omicron':	var='o'
				'pi':		var='p'
				'rho':		var='q'
				'sigma':	var='r'
				'tau':		var='s'
				'upsilon':	var='t'
				'phi':		var='u'
				'chi':		var='v'
				'psi':		var='w'
				'omega':	var='x'
				else:		var=''
			ENDCASE		
		ENDELSE
	ENDIF ELSE BEGIN
		str=strlowcase(str)
		IF keyword_set(ps) THEN BEGIN
			CASE str OF
				'alpha': 	var='A'
				'beta': 	var='B'
				'gamma':	var='G'
				'delta':	var='D'	
				'epsilon':	var='E'
				'zeta':		var='Z'
				'eta':		var='H'
				'theta':	var='Q'
				'iota':		var='I'
				'kappa':	var='K'
				'lambda':	var='L'
				'mu':		var='M'
				'nu':		var='N'
				'xi':		var='X'
				'omicron':	var='O'
				'pi':		var='P'
				'rho':		var='R'
				'sigma':	var='S'
				'tau':		var='T'
				'upsilon':	var='U'
				'phi':		var='F'
				'chi':		var='C'
				'psi':		var='Y'
				'omega':	var='W'
				else:		var=''
			ENDCASE
		ENDIF ELSE BEGIN
			CASE str OF
				'alpha': 	var='A'
				'beta': 	var='B'
				'gamma':	var='C'
				'delta':	var='D'
				'epsilon':	var='E'
				'zeta':		var='F'
				'eta':		var='G'
				'theta':	var='H'
				'iota':		var='I'
				'kappa':	var='J'
				'lambda':	var='K'
				'mu':		var='L'
				'nu':		var='M'
				'xi':		var='N'
				'omicron':	var='O'
				'pi':		var='P'
				'rho':		var='Q'
				'sigma':	var='R'
				'tau':		var='S'
				'upsilon':	var='T'
				'phi':		var='U'
				'chi':		var='V'
				'psi':		var='W'
				'omega':	var='X'
				else:		var=''
			ENDCASE
		ENDELSE
	ENDELSE
	
	IF keyword_set(ps) THEN f='!9' ELSE f='!7'	

	RETURN,f+var+'!'+f_old
END

FUNCTION  num2rom,n		;converts input to roman numeral (0-99)
	
	nvar = n
	ones=nvar mod 10
	tens=((nvar-ones) mod 100)/10
	
	IF tens EQ 0 THEN romstr=''
	IF tens EQ 1 THEN romstr='X'
	IF tens EQ 2 THEN romstr='XX'
	IF tens EQ 3 THEN romstr='XXX'
	IF tens EQ 4 THEN romstr='XL'
	IF tens EQ 5 THEN romstr='L'
	IF tens EQ 6 THEN romstr='LX'
	IF tens EQ 7 THEN romstr='LXX'
	IF tens EQ 8 THEN romstr='LXX'
	IF tens EQ 9 THEN romstr='XC'

	IF ones EQ 0 THEN romstr+=''
	IF ones EQ 1 THEN romstr+='I'
	IF ones EQ 2 THEN romstr+='II'
	IF ones EQ 3 THEN romstr+='III'
	IF ones EQ 4 THEN romstr+='IV'
	IF ones EQ 5 THEN romstr+='V'
	IF ones EQ 6 THEN romstr+='VI'
	IF ones EQ 7 THEN romstr+='VII'
	IF ones EQ 8 THEN romstr+='VIII'
	IF ones EQ 9 THEN romstr+='IX'


	RETURN,romstr

END

FUNCTION rom2num,romstr,verb=verb	;converts roman numeral string to number (0-99)

	x=strlen(romstr)
	romarr=intarr(x)
	matcharr=['I','V','X','L']
	numarr=[1,5,10,50]
	FOR i=0,(x-1) DO romarr[i]=numarr[where(matcharr EQ strmid(romstr,i,1))]
	IF keyword_set(verb) THEN print, romarr
	n=0
	cntr=0
	WHILE (cntr LE n(romarr)) DO BEGIN
		IF cntr EQ n(romarr) THEN BEGIN
			n+=romarr[cntr]
			cntr+=1
		ENDIF ELSE BEGIN
			a = romarr[cntr]
			b= romarr[cntr+1]
			IF a GE b THEN BEGIN
				n+=a
				cntr+=1
			ENDIF ELSE BEGIN
				n+=(b-a)
				cntr+=2
			ENDELSE		
		ENDELSE	

	ENDWHILE
	RETURN, n
END

;+
;NAME:
;	LOCALMAXMIN
;
;PURPOSE:
;	This function calculates the max and min of array within a defined x,y range.
;
;CALLING SEQUENCE:
;	result=LOCALMAXMIN(array,x,y,xrange,yrange)
;
;INPUTS:
;	array:		2D array of values
;	x:		i-scaling of array
;	y:		j-scaling of array
;	xrange:		[x0,x1] x bounds for sub array
;	yrange:		[y0,y1] y bounds for sub array
;
;OUTPUTS:
;	result:		[min,max]
;-

FUNCTION localmaxmin,array,x,y,xrange,yrange

	;check correct scaling
	l=size(array)
	IF (l[1] NE n(x)+1 OR l[2] NE n(y)+1) THEN BEGIN
		print, 'Array dimensions do not match scalings'
		RETURN, -1
	ENDIF

	temp1=where(x GE xrange[0] AND x LE xrange[1])
	temp2=where(y GE yrange[0] and y LE yrange[1])

	IF (temp1[0] EQ -1 OR temp2[0] EQ -1) THEN BEGIN
		print, 'no values found in xy range'
		RETURN, -1
	ENDIF

	arr_temp1=array[temp1,*]
	arr_temp2=arr_temp1[*,temp2]
	
	output=[min(arr_temp2),max(arr_temp2)]

	RETURN,OUTPUT
END
;+
;NAME:
;	BINARRAY
;
;PURPOSE:
;	This program bins the input array in defined del_i,del_j bins and reforms the array and i,j scalings
;
;CALLING SEQUENCE:
;	BINARRAY,array,ivec,jvec,del_i,del_j,new_array,new_ivec,new_jvec
;
;INPUTS:
;	array:		2D array of values
;	ivec:		i-scaling of array
;	jvec:		j-scaling of array
;	del_i:		number of i points per bin
;	del_j:		number of j points per bin
;
;OUTPUTS:
;	new_array:	Binned 2D array
;	new_ivec:	Binned 1D vector
;	new_jvec:	Binned 1D vector
;-

PRO bin_array,array,ivec,jvec,del_i,del_j,new_array, new_ivec, new_jvec
	IF (del_i NE 1 OR del_j NE 1) THEN BEGIN
		a=size(array)
		del_i=float(del_i)
		del_j=float(del_j)
		extra_i=a[1]-floor(a[1]/del_i)*del_i
		extra_j=a[2]-floor(a[2]/del_j)*del_j
		new_i=floor(a[1]/del_i)+(extra_i NE 0)
		new_j=floor(a[2]/del_j)+(extra_j NE 0)	
		new_array=fltarr(new_i,new_j)
		
	
		;bin the core of the array
		FOR i=0,new_i-1-(extra_i NE 0) DO BEGIN
			FOR j=0,new_j-1-(extra_j NE 0) DO BEGIN
				new_array[i,j]=total(array[i*del_i:(i+1)*del_i-1,j*del_j:(j+1)*del_j-1])/(del_i*del_j)
			ENDFOR
		ENDFOR
	
		;bin the extra i-slice
		IF extra_i NE 0 THEN FOR j=0,new_j-2 DO new_array[new_i-1,j]=total(array[a[1]-extra_i:(a[1]-1),j*del_j:(j+1)*del_j-1])/(extra_i*del_j)
		;bin the extra j-slice
		IF extra_j NE 0 THEN FOR i=0,new_i-2 DO  new_array[i,new_j-1]=total(array[i*del_i:(i+1)*del_i-1, a[2]-extra_j:a[2]-1])/(del_i*extra_j)
		;bin the extra i/j-slice
		IF (extra_i NE 0 AND extra_j NE 0) THEN new_array[new_i-1,new_j-1]=total(array[a[1]-extra_i:(a[1]-1),a[2]-extra_j:(a[2]-1)])/(extra_i*extra_j)
		
		;bin the ivec and jvec
		new_ivec=fltarr(new_i)
		FOR i=0,new_i-1-(extra_i NE 0) DO new_ivec[i]=total(ivec[i*del_i:(i+1)*del_i-1])/del_i
		IF extra_i NE 0 THEN new_ivec[new_i-1]=total(ivec[a[1]-extra_i:a[1]-1])/extra_i
		new_jvec=fltarr(new_j)
		FOR j=0,new_j-1-(extra_j NE 0) DO new_jvec[j]=total(jvec[j*del_j:(j+1)*del_j-1])/del_j
		IF extra_j NE 0 THEN new_jvec[new_j-1]=total(jvec[a[2]-extra_j:a[2]-1])/extra_j

	ENDIF ELSE BEGIN
		new_ivec=ivec
		new_jvec=jvec
		new_array=array		
	ENDELSE

END

;+
;NAME:
;	MINLOC
;
;PURPOSE:
;	This function calculates the i-point of the minimum in the input data
;
;CALLING SEQUENCE:
;	result=MINLOC(data)
;
;INPUTS:
;	data:	Can be 1D or 2D array of any values
;
;OUTPUTS:
;	result:		i, for 1D data and [i,j] for 2D data
;-

FUNCTION minloc,data
	;finds the array position for the minimum value
	a=size(data)
	IF a[0] EQ 1 OR a[0] EQ 0 THEN BEGIN
		tmp=min(data)
		output=!C
		RETURN, OUTPUT
	ENDIF ELSE BEGIN
		tmp=min(data)
		output=[!C-floor(!C/float(a[1]))*float(a[1]), floor(!C/float(a[1]))]
		RETURN, int(OUTPUT)
	ENDELSE
END

;+
;NAME:
;	MAXLOC
;
;PURPOSE:
;	This function calculates the i-point of the maximum  in the input data
;
;CALLING SEQUENCE:
;	result=MAXLOC(data)
;
;INPUTS:
;	data:	Can be 1D or 2D array of any values
;
;OUTPUTS:
;	result:		i, for 1D data and [i,j] for 2D data
;-

FUNCTION maxloc,data
	;finds the array position for the minimum value
	a=size(data)
	IF a[0] EQ 1 THEN BEGIN
		tmp=max(data)
		output=!C
		RETURN, OUTPUT
	ENDIF ELSE BEGIN
		tmp=max(data)
		output=[!C-floor(!C/float(a[1]))*a[1], floor(!C/float(a[1]))]
		RETURN, int(OUTPUT)
	ENDELSE
END

;+
;NAME:
;	ARRAY_SLICE
;
;PURPOSE:
;	This function returns the slice along io or jo of array
;
;CALLING SEQUENCE:
;	result=ARRAY_SLICE(array, ivec,jvec)
;
;INPUTS:
;	array:		2D array of any format
;	ivec:		1D i scaling of array
;	jvec:		1D j scaling of array
;
;OPTIONAL INPUTS:
;	io:		Use io=# (FLT) as a number in ivec to take the slice
;	jo:		Use jo=# (FLT) as a number in jvec to take the slice
;
;OUTPUTS:
;	result:		1D vector of array interpolated at io or jo
;-

FUNCTION array_slice,array,ivec,jvec,i=io,j=jo

	IF NOT keyword_set(io) AND NOT keyword_set(jo) THEN RETURN,-1 ;no value specified
	x=size(array)

	;interpolate at ivec=io
	IF keyword_set(io) THEN BEGIN
		IF io LT min(ivec) OR io GT max(ivec) THEN RETURN, -2	;value out of range		
		array_slice=fltarr(x[2])
		IF ivec[0] LT last(ivec) THEN BEGIN
			tmp1=where(ivec GT io)
			tmp2=where(ivec LT io)
			j1=tmp1[0]			;upper bound
			j2=tmp2[n(tmp2)]		;lower bound
		ENDIF ELSE BEGIN
			tmp1=where(ivec GT io)
			tmp2=where(ivec LT io)
			j1=tmp1[n(tmp1)]		;upper bound
			j2=tmp2[0      ]		;lower bound
		ENDELSE
		array_slice=array[j2,*]+(array[j1,*]-array[j2,*])/(ivec[j1]-ivec[j2])*(io-ivec[j2])
	ENDIF
	
	;interpolate at jvec=jo
	IF keyword_set(jo) THEN BEGIN
		IF jo LT min(jvec) OR jo GT max(jvec) THEN RETURN, -2	;value out of range		
		array_slice=fltarr(x[1])
		IF jvec[0] LT last(jvec) THEN BEGIN
			tmp1=where(jvec GT jo)
			tmp2=where(jvec LT jo)
			j1=tmp1[0]			;upper bound
			j2=tmp2[n(tmp2)]		;lower bound
		ENDIF ELSE BEGIN
			tmp1=where(jvec GT jo)
			tmp2=where(jvec LT jo)
			j1=tmp1[n(tmp1)]		;upper bound
			j2=tmp2[0]			;lower bound
		ENDELSE

		array_slice=array[*,j2]+(array[*,j1]-array[*,j2])/(jvec[j1]-jvec[j2])*(jo-jvec[j2])
	ENDIF	

	RETURN,array_slice
END

;+
;NAME:
;	LOGNAME
;
;PURPOSE:
;	This function is just a time save for calling the linux shell command logname to
;	check and see who's using the code
;
;CALLING SEQUENCE:
;	result=LOGNAME()
;
;OUTPUTS:
;	result:		STR of the current logname
;-

FUNCTION logname
	spawn, 'logname',output, error
	RETURN, output
END

;+
;NAME:
;	IS_FILE
;
;PURPOSE:
;	This function uses a CHKSUM call to Linux and feeds back on the error to see if a
;	file exists or not.
;
;CALLING SEQUENCE:
;	result=IS_FILE(filename)
;
;INPUTS:
;	filename:	STR full path from root to where file is located
;
;KEYWORD_PARAMETERS
;	size:		/size will make the output the size of the file in bytes
;
;OUTPUTS:
;	result:		INT 0 or 1 depending on the file existing (see /size)
;
;-

FUNCTION is_file,path,size=size
	spawn, 'cksum '+path, output,error
	err_match='cksum: '+path+': No such file or directory
	IF (strmatch(error, '*'+err_match+'*') EQ 1) THEN file=0 ELSE file=1
	
	IF keyword_set(size) AND file EQ 1 THEN BEGIN
		tmp=strsplit(output,' ',/extrac)
		output=long(tmp[1])
	ENDIF ELSE output=file

	RETURN,output
END



;+
;NAME:
;	IPT
;
;PURPOSE:
;	This function calculates the point, i, in the input vector that val is closest to.
;
;CALLING SEQUENCE:
;	result=IPT(vec,val)
;
;INPUTS:
;	vec:		1D vector
;	val:		just a number folks
;
;OUTPUTS:
;	result:		index, i, of vec
;-

FUNCTION ipt,vec, val
	;computes and returns the index value of vec that is closest to val, checking for exactness
	;-1 is returned if val is out of range of vec	
	check_exact=where(vec EQ val)
	IF check_exact[0] EQ -1 THEN BEGIN
		up=where(vec GT val)
		low=where(vec LT val)
		IF up[0] EQ -1 OR low[0] EQ -1 THEN RETURN, -1
		;i_h=up[0]
		;i_l=low[n(low)]
		;IF abs(vec[i_h] - val) LT abs(vec[i_l]-val) THEN output=i_h ELSE output=i_l
		tmp=min(abs(vec-val))
		output=!C
	ENDIF ELSE BEGIN
		output=check_exact[0]
	ENDELSE

	RETURN, output

END

;+
;NAME:
;	IVEC
;
;PURPOSE:
;	This function makes repeated calls to IPT so that a vector of values can be input
;
;-

FUNCTION ivec,vec,val_vec
	n_pts=n(val_vec)+1
	i_vec=intarr(n_pts)
	FOR i=0,n_pts-1 DO i_vec[i]=ipt(vec,val_vec[i])
	output=i_vec
	RETURN,output
END

;+
;NAME:
;	IBOUND
;
;PURPOSE:
;	This function calculates the bounding incides of vec that surround val
;
;CALLING SEQUENCE:
;	result=IBOUND(vec,val)
;
;INPUTS:
;	vec:		1D vector
;	val:		just a number folks
;
;OUTPUTS:
;	result:		bounding indices [i0,i1] of val. -1 if val is outside vec
;-

FUNCTION ibound,vec,val,rev=rev
	;computes the bounding index values of vec that surround val
	;-1 is returned if val is out of range of vec
	;if val EQ an element in vec then output=[i,i] instead of [i1,i2]

	tmp1=where(vec GE val)
	tmp2=where(vec LE val)
	IF (tmp1[0] EQ -1 OR tmp2[0] EQ -1) THEN BEGIN
		output=-1
	ENDIF ELSE BEGIN
		IF keyword_set(rev) THEN BEGIN
			low_i=tmp1[n(tmp1)]
			high_i=tmp2[0]
		ENDIF ELSE BEGIN
			high_i=tmp1[0]
			low_i=tmp2[n(tmp2)]
		ENDELSE
		output=[low_i, high_i]
	ENDELSE

	RETURN, output	
END

;+
;NAME:
;	INTERP_VEC_REFORM
;
;PURPOSE:
;	This function reforms a vector of values according to the indices of a
;	given scaling vector, which usually be for a 2D array.  This is a preprocessing
;	utility for INTERPOLATE.  See the help file on INTERPOLATE to understand why this
;	is useful
;
;CALLING SEQUENCE:
;	result=INTERP_VEC_REFROM(scale_vec,val_vec)
;
;INPUTS:
;	scale_vec:	scaling vector of an array
;	val_vec:	vector of points of interest for use in INTERPOLATE
;
;OUTPUTS:
;	result:		interpolated index values of scale_vec for the values in val_vec
;-

FUNCTION interp_vec_reform,scale_vec, val_vec,rev=rev
	;this function is a preprocessor for interpolate it takes
	;a vector of values (val_vec) and an assumed array scaling 
	;vector (scale_vec) and outputs an index vector according to
	;scale_vec just how interpolate() likes it.  Values out of
	;range are filled as -1 and /missing in interpolate is used to assign them
	
	ind_vec=fltarr(n(val_vec)+1)
	FOR i=0L,n(val_vec) DO BEGIN
		bnd=ibound(scale_vec,val_vec[i],rev=rev)
		IF bnd[0] NE -1 THEN BEGIN
			IF bnd[0] NE bnd[1] THEN BEGIN
				ind_vec[i]=(val_vec[i]-scale_vec[bnd[0]])/(scale_vec[bnd[1]]-scale_vec[bnd[0]])+bnd[0]
			ENDIF ELSE BEGIN
				ind_vec[i]=bnd[0]
			ENDELSE
		ENDIF ELSE ind_vec[i]=bnd[0]
	ENDFOR
	output=ind_vec	
	RETURN, output
END


;+
;NAME:
;	 EQ_STR
;
;PURPOSE: 
;	This function determines whether 2 structurs are completely equal from header tags
;	right down to type and value.
;
;CALLING SEQUENCE:
;	result=EQ_STR(str1,str2)
;
;INPUTS:
;	str1,str2:	two generically formatted structurs
;
;OUTPUTS:
;	result:		0 (false) or 1 (true) depending on the comparison
;
;RESTRICTIONS:
;	This is untested on complicated things like pointers and object references but should
;	work for arrays in structures and complex value numbers w/ more than one element 
;
;-


FUNCTION eq_str,str1,str2
	
	;verify that str1 and str2 really are structures
	output=0
	IF size(str1,/type) EQ 8 and size(str2,/type) EQ 8 THEN BEGIN
		n1=n_tags(str1)
		n2=n_tags(str2)
		IF n1 EQ n2 THEN BEGIN
			names1=tag_names(str1)
			names2=tag_names(str2)
			IF total(names1 EQ names2) EQ n1 THEN BEGIN
				types1=intarr(n1)
				types2=intarr(n2)
				FOR i=0,n1-1 DO BEGIN
					types1[i] = size(str1.(i),/type)
					types2[i] = size(str2.(i),/type)
				ENDFOR
				IF total(types1 EQ types2) EQ n1 THEN BEGIN
					dataeq=intarr(n1)
					FOR i=0,n1-1 DO dataeq[i]=total(str1.(i) EQ str2.(i))/n_elements(str1.(i))
					IF total(dataeq) EQ n1 THEN output=1
				ENDIF	
			ENDIF
		
		ENDIF

	ENDIF
	RETURN, output
END


;+
;NAME:
;	RUNS
;
;PURPOSE:
;	A handy little function if I can be a little arrogant.  This function takes a list of shots
;	and returns a sorted list of run numbers.  This is useful for plotting routines which
;	compile data by run number from a non sorted list of shot numbers
;
;CALLLING SEQUENCE:
;	result=RUNS(shotlist)
;
;INPUTS:
;	shotlist:	LONARR shots
;	
;OUTPUTS:
;	result:		LONARR of runs sorted from lowest to highest
;
;-

FUNCTION runs,shotlist
	l=n(shotlist)
	runlist=fix(shotlist/1000.0,type=3)
	runs_arr=lonarr(l+1)
	FOR i=0,l DO IF where(runs_arr EQ runlist[i]) EQ -1 THEN runs_arr[i]=runlist[i]
	output=runs_arr[where(runs_arr NE 0)]
	output=output[sort(output)]
	RETURN, output
END


;+
;NAME: 
;	DAY2RUN
;
;PURPOSE:
;	This function checks the system clock and turns the date into a run number 
;	based on a SYSTIME call.  Run number is Alcator type.
;
;CALLING SEQUENCE:
;	result=DAY2RUN()
;
;OUTPUTS:
;	result: INT run number (ie ... 1050426 for April 26,2005) 
;
;-

FUNCTION day2run
	date=systime(0,/utc)
	tmp=strsplit(date,' ',/extrac)
	cal=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
	run=long(1000000+10000*(int(tmp[4])-2000)+(where(cal EQ tmp[1])+1)*100+int(tmp[2]))
	RETURN,run
END

;+
;NAME:
;	SUM_ARRAY	
;
;PURPOSE:
;	This function returns the sum of a 2D array summed along the i or j direction
;
;CALLING SEQUENCE:
;	result=SUM_ARRAY(array, (/i or /j))
;
;INPUTS:
;	array:	Any type of 2D array
;
;OPTIONAL INTPUTS:
;	i:	/i to sum in the i-direction 
;	j:	/j to sum in the j-direction
;	k:	/k to sum in the k-direction
;	
;RESTRICTIONS:
;	The summing in i and j directions only work on 2D arrays.
;-

FUNCTION sum_array,array,i=i,j=j,k=k
	x=size(array)
	IF x[0] NE 2 AND x[0] NE 3  THEN RETURN, -1
	IF NOT keyword_set(i) AND NOT keyword_set(j) AND NOT keyword_set(k) THEN RETURN, -1
	output=-1
	IF keyword_set(i) AND x[0] NE 3 THEN BEGIN
		i_sum=fltarr(x[1])
		FOR n=0L,x[1]-1 DO i_sum[n]=total(array[n,*])
		output=i_sum
	ENDIF
	IF keyword_set(j) AND x[0] NE 3 THEN BEGIN
		j_sum=fltarr(x[2])
		FOR n=0L,x[2]-1 DO j_sum[n]=total(array[*,n])
		output=j_sum
	ENDIF
	IF keyword_set(k) AND x[0] EQ 3 THEN BEGIN
		ij_sum=fltarr(x[1],x[2])
		FOR n=0L,x[3]-1 DO ij_sum[*,*]+=array[*,*,n]
		output=ij_sum
	ENDIF
	RETURN, output
END

;+
;NAME:
;	AVEARRAY
;
;PURPOSE:
;	This function allows you to take average cuts out of 1D and 2D arrays using values of scaling
;	vectors (such as time or space for emissivity profiles).
;
;CALLING SEQUENCE:
;	result=AVEARRAY(array,x,x1,x2)
;
;
;-

FUNCTION avearray,array,x,x1,x2
	;array should be [space,time]
	i_low=ipt(x,x1)
	i_high=ipt(x,x2)
	IF i_low[0] EQ -1 OR i_high[0] EQ -1 THEN RETURN,-1
	size=size(array)
	IF size[0] EQ 1 THEN ave_array=mean(array[i_low:i_high])

	IF size[0] EQ 2 THEN BEGIN
		n_x=n(x)+1
		IF size[1] NE n_x AND size[2] NE n_x THEN RETURN,-1
		IF size[1] EQ n_x THEN ave_array=reform(sum_array(array[i_low:i_high,*],/j))/(i_high-i_low+1)
		IF size[2] EQ n_x THEN ave_array=reform(sum_array(array[*,i_low:i_high],/i))/(i_high-i_low+1)
	ENDIF

	output=ave_array
	IF keyword_set(debug) THEN stop
	RETURN,output
END

;+
;NAME:
;	VALIDSHOTS
;
;PURPOSE:
;	This function will find all the valid shot numbers between a range of run days.
;	The main purpose is for running tree stuffing scripts on backlogged data
;
;CALLING SEQUENCE:
;	result=VALIDSHOTS(run0,run1)
;
;INPTUS:
;	run0:	LONGINT of the starting run
;	run1:	LONGINT of ending run
;
;OPTIONAL INPUTS:
;	ip:	FLT of the lower level discriminator for the plasma current [kA] to be counted as valid
;
;OUTPUTS:
;	result:	LONARR of shot numbers where the spectroscopy tree exists
;		No non-shot number values are included
;
;RESTRICTIONS:
;	This can't wrap around a year so you have to use the optional input
;	and break the job into multiple parts.  I don't think it works for
;	years less than 2000 yet, though.
;
;MODIFICATION HISTORY:
;	Written by MLReinke (I forget when)
;	2-09-06 	MLR - fixed a crashing bug if there were no valid shots
;			and added an autodetect to the year and removed optional input
;			for the year.
;	1-23-08		ML Reinke - added the ip optional input
;-

FUNCTION validshots,run0,run1,ip=ip
	
	year=1900+floor(run0/10000)

	IF year GE 2000 THEN base_year=1000000000L+(year-2000)*10000000L ELSE $
		base_year=(year-1900)*10000000L
	days=[$
		[1,31],$
		[2,28],$
		[3,31],$
		[4,30],$
		[5,31],$
		[6,30],$
		[7,31],$
		[8,31],$
		[9,30],$
		[10,31],$
		[11,30],$
		[12,31]]

	month0=int(strmid(num2str(run0),3,2))
	day0=int(strmid(num2str(run0),5,2))
	month1=int(strmid(num2str(run1),3,2))
	day1=int(strmid(num2str(run1),5,2))

	shotlist=lonarr(5000)

	month=month0
	day=day0
	cntr=0
	WHILE month LE month1 DO BEGIN 
		base_run=base_year+month*100000L+day*1000L
		FOR i=0,50 DO BEGIN
			mdsopen, 'spectroscopy',base_run+i,/quiet,status=status
			IF status THEN BEGIN
				mdsclose
				IF keyword_set(ip) THEN BEGIN
					mdsopen,'magnetics',base_run+i,/quiet,status=ip_status
					IF ip_status THEN BEGIN
						ip_trace=mdsvalue('\magnetics::ip',/quiet,status=ip_status)			
						IF ip_status THEN BEGIN
							ip_max=max(abs(ip_trace)*1.0e-3)
							IF ip_max GT ip THEN ip_status=1 ELSE ip_status=0
							mdsclose
						ENDIF
					ENDIF
				ENDIF ELSE ip_status=1
				IF ip_status THEN BEGIN
					shotlist[cntr]=base_run+i
					cntr+=1
				ENDIF
			ENDIF
		ENDFOR		
		day+=1
		IF day GT days[1,where(days[0,*] EQ month)] THEN BEGIN
			month+=1
			day=1
		ENDIF
		IF month EQ month1 AND day EQ day1+1 THEN month = month1+1
	ENDWHILE

	tmp=where(shotlist NE 0)
	
	IF tmp[0] EQ -1 THEN output=-1 ELSE output=shotlist[tmp]
	RETURN,output
END

FUNCTION is_shot,shot,ip=ip
	IF NOT keyword_set(ip) THEN ip=100.0
	mdsopen,'magnetics',shot,/quiet,status=ip_status
	IF ip_status THEN BEGIN
		ip_trace=mdsvalue('\magnetics::ip',/quiet,status=ip_status)			
		IF ip_status THEN BEGIN
			ip_max=max(abs(ip_trace)*1.0e-3)
			IF ip_max GT ip THEN ip_status=1 ELSE ip_status=0
			
		ENDIF
		mdsclose
	ENDIF
	IF ip_status THEN output=1 ELSE output=0
	RETURN,output
END

;+
;NAME:
;	WAVESTATS
;
;PURPOSE:
;	This procedure is duplication of an Igor Pro command line call of the same name
;	and is used to print to terminal or a store in a structure a variety of useful
;	data for a 1D or 2D input array.  Note that you can type wavestats entirely with
;	your left hand.  Neat, huh.
;
;CALLING SEQUENCE:
;	WAVESTATS,input
;
;INPTUS:
;	input:	ANY NUMERIC array of 1 or 2 dimentions
;
;OPTIONAL OUTPUTS:
;	out:	STRUC of all the data that is computed
;
;KEYWORD PARAMETERS:
;	quiet:	/quiet to supress terminal outputs
;
;RESTRICTIONS:
;	This program used MOMENT to calculate the more significant statistics and 
;	calls other MLR_FUNCTIONS routines like MINLOC and MAXLOC (hence the 2D 
;	restriction).
;
;EXAMPLE
;	WAVESTATS, backround,out=out
;
;	Number of Dim = 2
;	Size = [101,20]
;	Number of Points = 2020
;	Maximum = 1017
;	Minimum = 985
;	Maxloc = [49,2]
;	Minloc = [0,3]
;	Average = 1004.6367
;	StdDev = 2.8164318
;	Var = 7.9322877
;	Skew = -2.522629
;	Kurt = 14.520065
;
;-

PRO wavestats,input, nan=nan,output=output,quiet=quiet
	x=size(input)
	IF NOT keyword_set(quiet) THEN BEGIN
		print, ' Number of Dim = '+num2str(x[0],1)
		IF x[0] EQ 2 THEN print, ' Size = ['+num2str(x[1],1)+','+num2str(x[2],1)+']'
		print, ' Number of Points = '+num2str(x[n(x)],1)
		print, ' Maximum = '+num2str(max(input))
		print, ' Minimum = '+num2str(min(input))
		IF x[0] EQ 1 THEN BEGIN
			print, ' Maxloc = '+num2str(maxloc(input),1)
			print, ' Minloc = '+num2str(minloc(input),1)
		ENDIF
		IF x[0] EQ 2 THEN BEGIN
			tmp=maxloc(input)
			print, ' Maxloc = ['+num2str(int(tmp[0]),1)+','+num2str(int(tmp[1]),1)+']'
			tmp=minloc(input)	
			print, ' Minloc = ['+num2str(int(tmp[0]),1)+','+num2str(int(tmp[1]),1)+']'
		ENDIF
		out=moment(input,sdev=sdevout, nan=nan)
		print, ' Average = '+num2str(out[0])
		print, ' StdDev = '+num2str(sdevout)		
		print, ' Var = '+num2str(out[1])
		print, ' Skew = '+num2str(out[2])
		print, ' Kurt = '+num2str(out[3])
	ENDIF

	out=moment(input,sdev=sdevout, nan=nan)
	IF x[0] EQ 1 THEN size=[x[n(x)],1] ELSE size=[x[1],x[2]]
	output={dim:x[0],npts:x[n(x)],size:size, min:min(input),max:max(input),minl:minloc(input),$
		maxl:maxloc(input),ave:out[0],sdev:sdevout,var:out[1], skew:out[2], kurt:out[3]}
	
END

;+
;NAME:
;	REMOUT
;
;PURPOSE:
;	The purpose of this function is to remove outliers from a dataset by repeatedly
;	truncating a dataset until all it's points lie within a user-definable number of
;	standard deviations.
;
;CALLING SEQUENCE:
;	result=REMOUT(input)
;
;INPUTS:
;	input:	1D array of datapoints
;
;OPTIONAL INPUTS:
;	nsig:	FLT number of standard deviations to which to truncate dataset [DEFAULT = 3.0]
;
;KEYWORD PARAMETERS:
;	mean:	/mean uses the deviation from the mean instead of the DEFAULT median
;
;OUTPUTS:
;	resut:	1D truncated dataset of INPUT where all values lie within +/-NSIG*STDEV 
;		of the median (or mean if selected)
;
;OPTIONAL OUTPUTS:
;	good:	good=good will be a 1D of size of INPUT which will have 1's for the values kept
;		and 0's for the values that were removed.  This can then be used to truncate
;		other datasets that depend on INPUT
;
;MODIFICATION HISTORY:
;	Written by: 	ML Reinke, 2/21/06
;
;-


FUNCTION remout,input,nsig=nsig,good=good,mean=mean
	
	
	IF NOT keyword_set(nsig) THEN nsig=3.0

	goodpts=intarr(n(input)+1)+1
	tmp=where(goodpts GT 0)
	allin=0

	WHILE n_elements(tmp) GT 1 AND allin NE 1  DO BEGIN
		datavec=input[tmp]
		IF keyword_set(mean) THEN pt=mean(datavec) ELSE pt=median(datavec)
		sig=nsig*stdev(datavec)
		out=where(datavec GT pt+sig OR datavec LT pt-sig)
		IF out[0] NE -1 THEN goodpts[tmp[out]]=0 ELSE allin=1
		tmp=where(goodpts GT 0)
	ENDWHILE	

	IF total(goodpts) EQ 0 OR total(goodpts) EQ 1 THEN output = -1
	IF allin EQ 1 THEN output=input[where(goodpts GT 0)]
	good=goodpts

	IF keyword_set(debug) THEN stop
	RETURN,output
END

PRO clean_heap,i=i,j=j
	IF NOT keyword_set(i) THEN i=0L
	IF NOT keyword_set(j) THEN j=65535
	FOR k=i,j DO ptr_free,ptr_valid(k,/cast)
END

;+
;NAME:
;	MAKE 
;	
;PURPOSE:
;	This function creates a FLTARR of a specified number of points
;	inclusively between two specificed endpoints.  Really it just
;	saves me from having to remember the one line way to do this
;	using FINDGEN.  It's also a port of the IGOR PRO syntax for
;	doing the same thing.
;
;CALLING SEQUENCE:
;	output=make(a,b,n)
;
;INPUTS:
;	a:	FLT lower bound
;	b:	FLT upper bound
;	n:	INT number of data points
;
;OUTPUTS:
;	output:	FLTARR of length n with evenly spaced points between (a,b) including
;		both a and b.
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke, Nov, 2005
;	
;-


FUNCTION make,a,b,n
	output=findgen(n)/(n-1)*(b[0]-a[0])+a[0]
	RETURN,output
END

FUNCTION makepos,n,xmin=xmin,ymin=ymin,ymax=ymax,xmax=xmax
	IF NOT keyword_set(xmin) THEN xmin=0.1
	IF NOT keyword_set(xmax) THEN xmax=0.975
	IF NOT keyword_set(ymin) THEN ymin=0.1
	IF NOT keyword_set(ymax) THEN ymax=0.975
	
	dy=(ymax-ymin)/n
	pos=fltarr(4,n)
	FOR i=0,n-1 DO BEGIN
		pos[0,i]=xmin
		pos[2,i]=xmax
		pos[1,i]=ymin+dy*i
		pos[3,i]=ymin+dy*(i+1)
	ENDFOR
	output=pos
	RETURN,output
END
	
	


; Taken from http://www.dfanning.com/tips/set_operations.html

FUNCTION SetIntersection, a, b
	minab = Min(a, Max=maxa) > Min(b, Max=maxb) ;Only need intersection of ranges
	maxab = maxa < maxb

  	 ; If either set is empty, or their ranges don't intersect: result = NULL.

	IF maxab LT minab OR maxab LT 0 THEN RETURN, -1
	r = Where((Histogram(a, Min=minab, Max=maxab) NE 0) AND  $
   	       (Histogram(b, Min=minab, Max=maxab) NE 0), count)

	IF count EQ 0 THEN RETURN, -1 ELSE RETURN, r + minab
END

FUNCTION dist,x,y
	IF n(x) NE n(y) THEN RETURN,-1
	dist=0
	FOR i=0,n(x) DO dist+=(x[i]-y[i])^2
	dist=sqrt(dist)
	
	RETURN,dist
END

FUNCTION str_dbl_ext,string,del_1, del_2
	tmp=strsplit(string, del_1,/extract)
	tmp=strsplit(tmp[1],del_2,/extract)
	output=tmp[0]
	RETURN, output
END


FUNCTION stack3d,a,b
	x_a=size(a)
	x_b=size(b)
	IF total(x_a[0:2]-x_b[0:2]) NE 0 THEN RETURN,-1
	
	new=fltarr(x_a[1],x_a[2],x_a[3]+x_b[3])
	new[*,*,0:x_a[3]-1]=a
	FOR i=0,n(b[0,0,*]) DO new[*,*,x_a[3]+i]=b[*,*,i]

	RETURN,new
END

FUNCTION make_fs_struc,shot,npts=npts,tree=tree
	common PSI,Psi_Grid,Rho_Grid,B_Grid,Time,R,Z,Psi_Surf,Z_midplane,current_shot,current_option,dPsidRho
	common PSI_ext,RSurf,ZSurf,RXpt,ZXpt,Topology,bs,order,Rknot,Zknot,PSI_BSCoef,RHO_BSCoef
	read_flux,shot,option,error,/bspline,efit_tree=tree

	IF NOT keyword_set(npts) THEN npts=101
	ntime=n(time)+1
	R_fine0=R(0) & R_fine1=R(n_elements(R)-1)
	R_fine=R_fine0+(R_fine1-R_Fine0)*findgen(npts)/100
	Z_fine0=Z(0) & Z_fine1=Z(n_elements(Z)-1)
	Z_fine=Z_fine0+(Z_fine1-Z_Fine0)*findgen(npts)/100
	rho_fine=fltarr(npts,npts,ntime)
	FOR i=0,ntime-1 DO Rho_fine[*,*,i]=BS2GD(0,0,R_fine,Z_fine,order,order,Rknot,Zknot,RHO_BSCoef(*,i))
	firstwall_nopump,rw,zw,shot=shot

	out={rfine:r_fine,zfine:z_fine,rhofine:rho_fine,t:time,rxpt:rxpt,zxpt:zxpt,rw:rw,zw:zw,shot:shot}
	RETURN,out
END

;+
;NAME:
;	VESSEL_PLOT
;
;PURPOSE:
;	This procedures sends a properly scaled plot of a poloidal C-Mod vacuum vessel/tile cross-section 
;	to either the postscript or x-windows device.
;
;CALLING SEQUENCE:
;	VESSEL_PLOT
;
;OPTIONAL INPUTS:
;	n:		INT:	The window number to plot to
;	color:		INT:	Color of the vessel plot lines (not color of plot) DEFAULT: 0
;	title:		STR:	String of the plot title DEFAULT: "Alcator C-Mod"
;	thick:		INT:	Thickness of the vessel plot lines DEFAULT: 0
;	shot:		LON:	Shot number to determine which vessel cross-section to use DEFAULT: 1070510001
;
;KEYWORD PARAMETERS:
;	force:		/force will force window number n to reopen and resize it accordingly
;	oplot:		/oplot will oplot the vessel structure on the currently selected plotting device
;	div:		/div will create a properly scaled zoomed plot of the lower divertor
;	edge:		/edge will create a properly scaled zoomed plot of the outer edge
;
;OUTPUTS:
;	All outputs are sent to the currently selected plotting device.  VESSEL_PLOTS will do nothing if a device other than
;	the X or PS is selected.
;
;PROCEDURE:
;	The vessel cross-section plots are taken from Brian Labombard's data files:
;			/home/labombard/minicad/vv_tiles_cryo_2007_s.vctr	2007-current
;			/home/labombard/minicad/vv_and_tiles_2002_s.vctr   	2002-2007
;			/home/labombard/minicad/vacuumvessel.vctr
;			/home/labombard/minicad/tiles_2002_s.vctr		1993-2002
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke, Fall 2006
;	5-10-07:	ML Reinke - added the new cyropump vessel cross-section and made it the default
;
;-

PRO vessel_plot,ps=ps,n=n, force=force,title=title,x_offset=x_offset,y_offset=y_offset,color=color,thick=thick,oplot=oplot,$
		shot=shot,d_old=d_old,div=div,edge=edge,pedestal=pedestal,wall=wall
	IF NOT keyword_set(shot) THEN shot=1070510001
	IF NOT keyword_set(title) THEN title='Alcator C-Mod'
	IF NOT keyword_set(x_offset) THEN x_offset=0.0
	IF NOT keyword_set(y_offset) THEN y_offset=0.0
	
	CASE !d.name OF
		'PS': 	ps=1
		'X':	ps=0
		ELSE:	RETURN
	ENDCASE


	IF NOT keyword_set(div) AND NOT keyword_set(edge) AND NOT keyword_set(pedestal) AND NOT keyword_set(wall)  THEN BEGIN
		yrange=[-0.67,0.67]
		xrange=[0.38,1.1]
		IF keyword_set(ps) THEN BEGIN
			xsize=6.0
			ysize=6.0*810.0/500.0
		ENDIF ELSE BEGIN
			xsize=500.0
			ysize=818.0
		ENDELSE
	ENDIF ELSE BEGIN
		IF keyword_set(div) THEN BEGIN
			yrange=[-0.67,-0.15]
			xrange=[0.38,0.85]
			IF keyword_set(ps) THEN BEGIN
				xsize=7.0
				ysize=7.25
			ENDIF ELSE BEGIN
				xsize=700.0
				ysize=740.0
			ENDELSE
		ENDIF
		IF keyword_set(edge) THEN BEGIN
			yrange=[-0.35,0.35]
			xrange=[0.574,0.95]
			IF keyword_set(ps) THEN BEGIN
				xsize=6.0
				ysize=6.0*810.0/500.0
			ENDIF ELSE BEGIN
				xsize=500
				ysize=818
			ENDELSE	
		ENDIF
		IF keyword_set(wall) THEN BEGIN
			yrange=[-0.35,0.35]
			xrange=[0.574,0.95]-0.174
			IF keyword_set(ps) THEN BEGIN
				xsize=6.0
				ysize=6.0*810.0/500.0
			ENDIF ELSE BEGIN
				xsize=500
				ysize=818
			ENDELSE	
		ENDIF
		IF keyword_set(pedestal) THEN BEGIN
			yrange=[-0.1,0.1]
			xrange=[0.822,0.93]
			IF keyword_set(ps) THEN BEGIN
				xsize=6.0
				ysize=6.0*810.0/500.0
			ENDIF ELSE BEGIN
				xsize=500
				ysize=818
			ENDELSE	
		ENDIF
	ENDELSE

	IF NOT keyword_set(oplot) THEN BEGIN
		IF NOT keyword_set(ps) THEN BEGIN
			IF NOT keyword_set(n) THEN n=19
			device, window_state=var
			IF var[n] EQ 0 OR keyword_set(force) EQ 1 THEN window,n,xsize=xsize,ysize=ysize,xpos=0,ypos=670,title='vessel cx,'+num2str(n) $
				ELSE wset,n
		ENDIF ELSE BEGIN
			d_old=!d
			device, xsize=xsize, ysize=ysize, /inches
		ENDELSE
		plot, [0],[0],title=title,chars=1.3,xrange=xrange,yrange=yrange,xtit='R (m)',ytit='Z (m)',/xsty, /ysty
	ENDIF
	
	IF shot GT 1070101000 THEN BEGIN
		restore, "/home/labombard/minicad/vv_tiles_cryo_2007_s.vctr"
		for i=0,nvctr-1 do oplot,xvctr(0:lvctr(i)-1,i)+x_offset,yvctr(0:lvctr(i)-1,i)+y_offset,color=color,thick=thick
	ENDIF		

	IF shot GT 1020101000 AND shot LT 1070101000 THEN BEGIN
		restore, "/home/labombard/minicad/vv_and_tiles_2002_s.vctr"
		for i=0,nvctr-1 do oplot,xvctr(0:lvctr(i)-1,i)+x_offset,yvctr(0:lvctr(i)-1,i)+y_offset,color=color,thick=thick
	ENDIF
	
	IF shot LT 1020101000 THEN BEGIN
		restore, '/home/labombard/minicad/vacuumvessel.vctr'
		for i=0,nvctr-1 do oplot,xvctr(0:lvctr(i)-1,i)+x_offset,yvctr(0:lvctr(i)-1,i)+y_offset,color=color,thick=thick
		restore, '/home/labombard/minicad/tiles_2002_s.vctr'
		for i=0,nvctr-1 do oplot,xvctr(0:lvctr(i)-1,i)+x_offset,yvctr(0:lvctr(i)-1,i)+y_offset,color=color,thick=thick
	ENDIF

	;IF NOT keyword_set(oplot) THEN IF keyword_set(ps) THEN device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm
	
END

;+
;NAME:
;	FS_PLOT
;
;PURPOSE:
;	This procedure OPLOTS flux surfaces, presumabley onto a window where VESSEL_PLOT has been called.
;
;CALLING SEQUENCE
;	FS_PLOT,shot,time
;
;INPUTS:
;	shot:	LONG	shot number
;	time:	FLT	time point [sec]
;
;MODIFICATION HISTORY:
;	6/18/11		M.L. Reinke - modified to allow non-ANALYSIS EFIT's to be used
;-

PRO fs_plot,shot,time,num=num,color=color,style=style,thick=thick,sol=sol,nsol=nsol,rad_sol=rad_sol,invssel=invessel,tree=tree
	IF NOT keyword_set(color) THEN color=100
	IF NOT keyword_set(style) THEN style=2 ELSE IF style EQ -1 THEN style=0
	IF NOT keyword_set(thick) THEN thick=1.0
	IF NOT keyword_set(num) THEN num=3
	IF NOT keyword_set(nsol) THEN nsol=4
	IF NOT keyword_set(rad_sol) THEN rad_sol=0.905
	IF NOT keyword_set(tree) THEN tree='analysis'
	IF tree EQ 'analysis' THEN tpath='\ANALYSIS::TOP:EFIT.RESULTS.' ELSE tpath='\'+strupcase(tree)+'::TOP.RESULTS.'
		
	
	efit_psicont,psi_r,psi_z,psi_n,psi_t,shot=shot,time=time,tree=tree
	FOR i=0,floor(n(psi_r[0,*])/num)-1 DO BEGIN
		out_r=reform(psi_r[*,num*i+1])
		out_z=reform(psi_z[*,num*i+1])
		tmp=where(out_r NE 0)
		oplot, out_r[tmp], out_z[tmp], color=color,linestyle=style,thick=thick
	ENDFOR
	help, psi_r

	mdsopen, tree, (shot)
	bdry=mdsvalue(tpath+'G_EQDSK:BDRY',/quiet,status=status)
	efit_time=mdsvalue('dim_of('+tpath+'G_EQDSK:RMAXIS,0)',/quiet,status=status)
	IF n(efit_time) EQ 0 THEN efit_i = 0 ELSE efit_i=ipt(efit_time,time)
	IF NOT status THEN BEGIN
		rbbs=mdsvalue(tpath+'G_EQDSK:RBBBS',/quiet, status=rstat)
		zbbs=mdsvalue(tpath+'G_EQDSK:ZBBBS',/quiet,status=zstat)
		IF rstat AND zstat THEN BEGIN
			nbbs=n(rbbs[0,*])+1
			ntime=n(rbbs[*,0])+1
			bdry=fltarr(2,nbbs,ntime)
			FOR i=0,ntime-1 DO BEGIN
				bdry[0,*,i] = rbbs[i,*]
				bdry[1,*,i] = zbbs[i,*]
			ENDFOR
			output=bdry
		ENDIF ELSE output = -1			
	ENDIF ELSE output=bdry
	mdsclose,tree, (shot)
	efit_lcfs=output

	rbdry=efit_lcfs[0,where(efit_lcfs[0,*,efit_i] NE 0 AND efit_lcfs[1,*,efit_i] NE 0),efit_i]
	zbdry=efit_lcfs[1,where(efit_lcfs[0,*,efit_i] NE 0 AND efit_lcfs[1,*,efit_i] NE 0),efit_i]
	rbdry=rbdry[0:n(rbdry)-1]
	zbdry=zbdry[0:n(zbdry)-1]
	rb_plt=[rbdry,rbdry[0]]
	zb_plt=[zbdry,zbdry[0]]
	oplot,rb_plt,zb_plt,color=color,thick=2.0


	IF keyword_set(sol) THEN BEGIN
		z_sol=zbdry[maxloc(rbdry)]
		sol_pts=make(max(rbdry),rad_sol,nsol+2)
		sol_pts=sol_pts[1:nsol]
		FOR i=0,nsol-1 DO BEGIN
			sol_tracet,shot,efit_time[efit_i],sol_pts[i],zbdry,xcontr,ycontr,ncontr,psivl,direction=1,/accurate,/nolim
			IF keyword_set(invessel) THEN BEGIN
				inv=line_invessel(xcontr,ycontr)
				tmp=where(inv EQ 1)
				xcontr=xcontr[tmp]
				ycontr=ycontr[tmp]		
			ENDIF
			oplot, xcontr,ycontr,color=30
			sol_tracet,shot,efit_time[efit_i],sol_pts[i],zbdry,xcontr,ycontr,ncontr,psivl,direction=2, /accurate,/nolim
			IF keyword_set(invessel) THEN BEGIN
				inv=line_invessel(xcontr,ycontr)
				tmp=where(inv EQ 1)
				xcontr=xcontr[tmp]
				ycontr=ycontr[tmp]		
			ENDIF
			oplot, xcontr,ycontr,color=30
		ENDFOR
	ENDIF

END

FUNCTION Inside, x, y, px, py, Index=index

   ; Purpose: see if point is inside polygon
   ; Category: maths
   ; Input: x, y - [vector of] points
   ;        px,py - points defining polygon (will be closed automatically)
   ; Output: vector of 1's and 0's
   ;         OR
   ;         indicies of points inside (if /Index is set)
   ; Author: "Bård Krane" 
   ; Mods: wmc - make it work with x, y as vectors
   ; More-info: posted to comp.lang.idl-pvwave on Wed, 01 Apr 1998 12:26:38 +0200
   ; See-also: http://www.ecse.rpi.edu/Homepages/wrf/geom/pnpoly.html for another method, possibly better
   ; This is better than my routine "is_inside"
   ; Note: reduce test from 1e-8 to 1e-4, since we are usually in single precision. In fact, "0.1"
   ;       would do as well.

       On_Error, 1
    
       sx = Size(px)
       sy = Size(py)
       IF (sx[0] EQ 1) THEN NX = sx[1] ELSE Message,'Variable px is not a vector'
       IF (sy[0] EQ 1) THEN NY = sy[1] ELSE Message,'Varialbe py is not a vector'
       IF (NX EQ NY)   THEN N = NX ELSE Message,'Incompatible vector dimensions'

       tmp_px =  [px, px[0]]                           ; Close Polygon in x
       tmp_py =  [py, py[0]]                           ; Close Polygon in y
    
       i  = Indgen(N,/Long)                            ; indices 0...N-1
       ip = Indgen(N,/Long) + 1                        ; indices 1...N
   
       nn = N_Elements(x) 
       X1 = tmp_px(i)  # Replicate(1,nn) - Replicate(1,n) # Reform([x],nn)
       Y1 = tmp_py(i)  # Replicate(1,nn) - Replicate(1,n) # Reform([y],nn)
       X2 = tmp_px(ip) # Replicate(1,nn) - Replicate(1,n) # Reform([x],nn)
       Y2 = tmp_py(ip) # Replicate(1,nn) - Replicate(1,n) # Reform([y],nn)
 
       dp = X2*X1 + Y1*Y2                               ; Dot-product
       cp = X1*Y2 - Y1*X2                               ; Cross-product
       theta = Atan(cp,dp)

       ret = Replicate(0L, N_Elements(x))
       i = Where(Abs(Total(theta,1)) GT 0.01,count)
       IF (count GT 0) THEN ret(i)=1
       IF (N_Elements(ret) EQ 1) THEN ret=ret[0]

       IF (Keyword_Set(index)) THEN ret=(Indgen(/Long, N_Elements(x)))(Where(ret eq 1))

       RETURN, ret

END


;+
;NAME:
;	IONGYRO
;
;PURPOSE:
;	This function returns the ion gyro radius rho=v_th/Omega
;
;CALLING SEQUENCE:
;	result=IONGYRO(z,Ti,B)
;
;INPUTS:
;	z	INT charge of ion
;	Ti	FLT ion temperature [in eV]
;	B	FLT magnetic field [in Tesla]
;
;OPTIONAL INPUTS:
;	zstr	STR of the element symbol.  Use this to specify the element when not fully stripped ion
;		    or to choose different hydrogenic species.
;
;OUTPUTS:
;	result:	FLT of the ion gyro radius in meters	
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke, April 2nd, 2007
;	
;-	

FUNCTION iongyro,z,Ti,B,zstr=zstr
	IF NOT keyword_set(zstr) THEN m=read_atomic_mass(read_atomic_name(z)) ELSE m=read_atomic_mass(zstr)
	m_in_kg=m*1.6605e-27
	Ti_in_J=Ti*1.6022e-19
	rho=sqrt(2.0*Ti_in_J*m_in_kg)/(Z*1.602e-19*B)
	
	RETURN,rho
END


;+
;NAME:
;	FREQSPEC
;	
;PURPOSE:
;	This program is a quick and dirty way to plot frequency spectra to look for things like noise.  I would
;	advise NOT using it for quantitative data analysis.
;
;CALLING SEQUENCE
;	FREQSPEC,v,t
;
;INPUTS:
;	v:	FLTARR	of the the signal
;	t:	FLTARR	of the time 
;
;OPTIONAL INPUTS:
;	w:	INT	number of points to smooth frequency spectra before plotting DEFAULT: 0
;	tr:	FLTARR	[t1,t2] of the time region to use DEFAULT: all time points
;	fr:	FLTARR	[f1,f2]	of the frequence range to plot DEFAULT: [1.0, Nyquist]
;	pr:	FLTARR	[p1,p2] of the fractional power range to plot
;
;KEYWORD PARAMETERS:
;	debug	/debug will stop the code before the return statement
;
;OUTPUTS:
;	Outputs are plotted to x-windows.
;		WINDOW-0: Time history of the signal trace
;		WINDOW-1: Frequency spectrum of the signal trace
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke - 12/07
;	3/8/13		M.L. Reinke - added the data optional output
;-

PRO freqspec,v,t,w=w,debug=debug,tr=tr,fr=fr,data=data,noplot=noplot
	IF keyword_set(tr) THEN BEGIN
		i_low=ipt(t,tr[0])
		i_high=ipt(t,tr[1])
		v_copy=v
		t_copy=t
		v=v[i_low:i_high]
		t=t[i_low:i_high]
	ENDIF

	;plot raw signal for all time
	IF NOT keyword_set(noplot) THEN BEGIN
		openwin,0
		plot, t,v, xtitle='Time [s]', ytitle='Signal', title='Signal vs. Time',chars=1.2
	ENDIF
	freq_max=n_elements(t)/(t[n_elements(t)-1]-t[0])				;find maximum frequency from time history
	freq=fltarr(n_elements(t)/2.0)							;define frequency scale from 0->1/2*freq_max
	freq=make(0,n(freq)+1,n(freq)+1)*0.5*freq_max/n_elements(freq)

	sum=0
	fft=fltarr(n_elements(t))
	
	temp=fft(v)
	fft_v=abs(temp)
	sum=total(temp)
	fft_v/=sum						;normalize to total power 
	IF keyword_set(w) THEN fft_v=smooth(fft_v,w)		;smooth, cause it's rough like that...
	fft_v/=max(fft_v)					;normalize on 0-1 scale (now plot should be of fractional total power)
	

	IF keyword_set(fr) THEN freqrange=fr ELSE freqrange=[1.0,freq_max/2.0]
	IF keyword_set(pr) THEN powrange=pr ELSE powrange=[fft_v[n_elements(freq)-1],1]
	IF NOT keyword_set(noplot) THEN BEGIN
		openwin,1
		plot,freq,fft_v,yrange=powrange,xtitle='Freq [Hz]',$
			ytitle='Fractional Spectral Power (%)',/xlog,title='Freq. Spectrum', /ylog,chars=1.2,xr=freqrange,/xsty
	ENDIF
	data={f:freq/1.0e3,npow:fft_v}
	IF keyword_set(debug) THEN stop
	IF keyword_set(tr) THEN BEGIN
		t=t_copy
		v=v_copy
	ENDIF
END

;+
;NAME:
;	TVZOOM
;	
;PURPOSE:
;	This procedure interpolates a given image and redisplays it zoomed in.
;
;CALLING SEQUENCE:
;	TVZOOM,image
;
;INPUTS:
;	image		FLTARR [n_x,n_y,n_fr] of an image you'd normally TV
;
;OPTIONAL INPUTS:
;	zoom:		FLT of the zoom DEFAULT=2
;	frame:		INT of the frame number (if n_fr > 0)  DEFAULT = 0
;	win:		INT of the window number to plot to.  Will automattically resize for the new image.
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke 12/13/07
;	3/17/10:	ML Reinke - added ability to send /true and /norm to TV
;-

PRO tvzoom,image,zoom=zoom,frame=frame,win=win,true=true,norm=norm
	IF NOT keyword_set(frame) THEN frame=0
	IF NOT keyword_set(zoom) THEN zoom=2.0
	IF NOT keyword_set(win) THEN win=10
	x=size(image)
	nx=x[1]
	ny=x[2]
	IF x[0] EQ 3 THEN pic=image[*,*,frame] ELSE pic=image
	
	i_pt=indgen(nx*zoom)/zoom
	j_pt=indgen(ny*zoom)/zoom
	new_pic=interpolate(pic,i_pt,j_pt,/grid)
	openwin,win,xsize=nx*zoom,ysize=ny*zoom
	IF !d.x_size NE nx*zoom AND !d.y_size NE ny*zoom THEN BEGIN
		wdelete,win
		openwin,win,xsize=nx*zoom,ysize=ny*zoom		
	ENDIF

	tv,new_pic,true=true,norm=norm

END


FUNCTION pts_in_poly,rpoly,zpoly,rpt,zpt,debug=debug,cent=cent,norm=norm

	n_pts=n(rpt)+1
	n_poly=n(rpoly)+1
	IF n(rpt) NE n(zpt) THEN RETURN,-1
	output=intarr(n_pts)
	norm=fltarr(2,n_poly-1)
	cent=fltarr(2,n_poly-1)

	FOR i=0,n_poly-2 DO BEGIN
		pt_a=[rpoly[i], zpoly[i]]
		pt_b=[rpoly[i+1], zpoly[i+1]]
		pt_c=[0.5*(pt_a[0]+pt_b[0]),0.5*(pt_a[1]+pt_b[1])]
		sl=(pt_b[1]-pt_a[1])/(pt_b[0]-pt_a[0])
		IF finite(sl) EQ 0 THEN invsl=0.0 ELSE invsl=-1.0/sl
		vec_ab=pt_b-pt_a
		dx=0.5*(pt_b[0]-pt_a[0])
		IF finite(invsl) EQ 0 THEN BEGIN
			n1pt=pt_c+[0,(dx)]
			n2pt=pt_c+[0,-1.0*dx]
		ENDIF ELSE BEGIN
			n1pt=pt_c+[dx,invsl*(dx)]
			n2pt=pt_c+[-1.0*dx,invsl*(-1.0*dx)]
		ENDELSE
		vec_n1=n1pt-pt_c
		vec_n2=n2pt-pt_c
		cross1=crossp([vec_ab,0],[vec_n1,0])
		cross2=crossp([vec_ab,0],[vec_n2,0])
		IF cross1[2] LT 0 THEN vec_n=vec_n1 ELSE vec_n=vec_n2

		cent[*,i]=pt_c
		norm[*,i]=vec_n
	ENDFOR
	rmax=max(rpoly)
	rmin=min(rpoly)
	zmax=max(zpoly)
	zmin=min(zpoly)

	FOR i=0,n_pts-1 DO BEGIN
		IF rpt[i] LT rmax AND rpt[i] GT rmin AND zpt[i] GT zmin AND zpt[i] LT zmax THEN BEGIN
			pt=[rpt[i],zpt[i]]
			check=intarr(n_poly-1)
			FOR j=0,n_poly-2 DO BEGIN	
				vec_pt=cent[*,j]-pt
				IF total(norm[*,j]*vec_pt) GE 0 THEN check[j]=1 ELSE check[j]=0
			ENDFOR
			IF total(check) EQ n_poly-1 THEN output[i]=1 ELSE output[i]=0
			IF keyword_set(debug) THEN BEGIN
				print, 'Outward Normal', vec_n
				print, 'Vector to Point', vec_pt
			ENDIF	
		ENDIF
	ENDFOR

	RETURN, output
END

FUNCTION Inside, x, y, px, py, Index=index

   ; Purpose: see if point is inside polygon
   ; Category: maths
   ; Input: x, y - [vector of] points
   ;        px,py - points defining polygon (will be closed automatically)
   ; Output: vector of 1's and 0's
   ;         OR
   ;         indicies of points inside (if /Index is set)
   ; Author: "Bård Krane" 
   ; Mods: wmc - make it work with x, y as vectors
   ; More-info: posted to comp.lang.idl-pvwave on Wed, 01 Apr 1998 12:26:38 +0200
   ; See-also: http://www.ecse.rpi.edu/Homepages/wrf/geom/pnpoly.html for another method, possibly better
   ; This is better than my routine "is_inside"
   ; Note: reduce test from 1e-8 to 1e-4, since we are usually in single precision. In fact, "0.1"
   ;       would do as well.

       On_Error, 1
    
       sx = Size(px)
       sy = Size(py)
       IF (sx[0] EQ 1) THEN NX = sx[1] ELSE Message,'Variable px is not a vector'
       IF (sy[0] EQ 1) THEN NY = sy[1] ELSE Message,'Varialbe py is not a vector'
       IF (NX EQ NY)   THEN N = NX ELSE Message,'Incompatible vector dimensions'

       tmp_px =  [px, px[0]]                           ; Close Polygon in x
       tmp_py =  [py, py[0]]                           ; Close Polygon in y
    
       i  = Indgen(N,/Long)                            ; indices 0...N-1
       ip = Indgen(N,/Long) + 1                        ; indices 1...N
   
       nn = N_Elements(x) 
       X1 = tmp_px(i)  # Replicate(1,nn) - Replicate(1,n) # Reform([x],nn)
       Y1 = tmp_py(i)  # Replicate(1,nn) - Replicate(1,n) # Reform([y],nn)
       X2 = tmp_px(ip) # Replicate(1,nn) - Replicate(1,n) # Reform([x],nn)
       Y2 = tmp_py(ip) # Replicate(1,nn) - Replicate(1,n) # Reform([y],nn)
 
       dp = X2*X1 + Y1*Y2                               ; Dot-product
       cp = X1*Y2 - Y1*X2                               ; Cross-product
       theta = Atan(cp,dp)

       ret = Replicate(0L, N_Elements(x))
       i = Where(Abs(Total(theta,1)) GT 0.01,count)
       IF (count GT 0) THEN ret(i)=1
       IF (N_Elements(ret) EQ 1) THEN ret=ret[0]

       IF (Keyword_Set(index)) THEN ret=(Indgen(/Long, N_Elements(x)))(Where(ret eq 1))

       RETURN, ret

END



;----------------------------------------------------------------------------------
;MLR FUNCTIONS/PROCEDURES FOR PLOTTING
;	these aren't documented because they're of a 'personal nature'
;----------------------------------------------------------------------------------

PRO helpsym,n=n
	IF NOT keyword_set(n) THEN n=10
	openwin,n
	plot, [0],[0], xr=[0,24],yr=[0,24],/xsty,/ysty
	FOR i=0,7 DO oplot, [i],[i],psym=i,symsize=3.5,color=100
	FOR i=9,23 DO BEGIN
		makesym,i
		oplot, [i],[i],psym=8,symsize=3.5,color=200
	ENDFOR
END

FUNCTION colormap,vec

	ntimes=n(vec)
	;set colors
	colormap=intarr(ntimes+1)
	n_colors=!d.table_size-20	;those last few suck
	FOR i=0,ntimes DO colormap[i]=i*(n_colors)/(ntimes+1)
	
	RETURN, colormap
END

PRO log_bolo,shot 
	wset,5
	win_save,'em_cont_'+num2str(shot,1)
	wset,4
	win_save,'em_time_'+num2str(shot,1)
	wset,6
	win_save,'em_rad_'+num2str(shot,1)

	run=runs(shot)
	run=run[0]

	basepath='/usr/local/cmod/logbook/'+logname()+'/bolo/'
	spawn, 'cd '+basepath, output,error
	IF (strmatch(error, '*No such file or directory*') EQ 1) THEN $
		spawn, 'mkdir '+basepath, output,error
	runpath=basepath+num2str(run,1)+'/'
	spawn, 'cd '+runpath, output,error
	IF (strmatch(error, '*No such file or directory*') EQ 1) THEN $
		spawn, 'mkdir '+runpath, output,error
	
	spawn,'mv em_*_'+num2str(shot,1)+'.jpg '+runpath,output,error
	print, error
	print, ''
	print, '<img src="logbook/'+logname()+'/bolo/'+num2str(run,1)+'/em_rad_'+num2str(shot,1)+'.jpg">
	print, '<img src="logbook/'+logname()+'/bolo/'+num2str(run,1)+'/em_time_'+num2str(shot,1)+'.jpg">
	print, '<img src="logbook/'+logname()+'/bolo/'+num2str(run,1)+'/em_cont_'+num2str(shot,1)+'.jpg">
END

PRO log_lymid,shot 
	wset,8
	win_save,'ly_cont_'+num2str(shot,1)
	wset,7
	win_save,'ly_time_'+num2str(shot,1)
	wset,9
	win_save,'ly_rad_'+num2str(shot,1)

	run=runs(shot)
	run=run[0]

	basepath='/usr/local/cmod/logbook/'+logname()+'/edge/'
	spawn, 'cd '+basepath, output,error
	IF (strmatch(error, '*No such file or directory*') EQ 1) THEN $
		spawn, 'mkdir '+basepath, output,error
	runpath=basepath+num2str(run,1)+'/'
	spawn, 'cd '+runpath, output,error
	IF (strmatch(error, '*No such file or directory*') EQ 1) THEN $
		spawn, 'mkdir '+runpath, output,error
	
	spawn,'mv ly_*_'+num2str(shot,1)+'.jpg '+runpath,output,error
	print, error
	print, ''
	print, '<img src="logbook/'+logname()+'/edge/'+num2str(run,1)+'/ly_rad_'+num2str(shot,1)+'.jpg">
	print, '<img src="logbook/'+logname()+'/edge/'+num2str(run,1)+'/ly_time_'+num2str(shot,1)+'.jpg">
	print, '<img src="logbook/'+logname()+'/edge/'+num2str(run,1)+'/ly_cont_'+num2str(shot,1)+'.jpg">
END

;+
;NAME:
;	OPENWIN
;
;PURPOSE:
;	This procedure opens new x-windows or sets them to active if they are open and the
;	right size.  If the PS device is active then nothing happens.  This procedure should
;	be used for all window management since it is device independent (for X and PS)
;
;CALLING SEQUENCE:
;	openwin,n
;
;INPUTS:
;	n:	INT	window nubmer (0-31)
;
;OPTIONAL INPUTS:
;	xsize,ysize,xpos,ypos,title are carred through to the window command
;
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke - 2007
;	1-24-08:	ML Reinke - added the ability to force a new window with xsize and ysize
;
;-

PRO openwin,n,var=var,xsize=xsize,ysize=ysize,xpos=xpos,ypos=ypos,title=title
	IF !d.name EQ 'X' THEN BEGIN
		IF NOT keyword_set(var) THEN device, window_state=var
		IF var[n] EQ 1 THEN BEGIN
			newwin=0
			x_ok=1
			IF keyword_set(xsize) THEN IF !d.x_size NE xsize THEN x_ok=0
			y_ok=1
			IF keyword_set(ysize) THEN IF !d.y_size NE ysize THEN y_ok=0
			IF x_ok AND y_ok THEN wset,n ELSE newwin=1
		ENDIF ELSE newwin=1
		IF newwin THEN window,n,xsize=xsize,ysize=ysize,xpos=xpos,ypos=ypos,title=title
	ENDIF
END

PRO psplot_old
	set_plot,'ps'
	!p.font=0
	device, /school
	device, /color
	device, /portrait
	device, xsi=9, ysi=7, /in
;	device, yoffset=24.5
	device, encapsulated=1
	device, preview=1
	!p.thick=3.0
	!x.thick=3.0
	!y.thick=3.0
END

PRO psplot_2
	set_plot,'ps'
	!p.font=0
	device, /school
	device, /color
	device, /portrait
	device, encapsulated=0
	device, preview=0
	!p.thick=3.0
	!x.thick=3.0
	!y.thick=3.0
END

PRO psplot
	set_plot,'ps'
	!p.font=00
	device, /color
	device, /portrait
	device, encapsulated=0
	device, preview=0
	device, yoffset=2.0
	!p.thick=3.0
	!x.thick=3.0
	!y.thick=3.0
END

PRO psplot_3
	psplot
END

;mlr 9/10/12 - replaced the loadct,12 functionality of xwplot.
;              Whoever messed with this please come to be for an
;              explanation of how MLR_FUNCTIONS is not yours to fuck with.
PRO xwplot,ct=ct
	IF NOT keyword_set(ct) THEN ct=12
	set_plot,'x'
	device,true_color=24,retain=2,decomposed=0  
	if(!d.table_size gt 16) then loadct,ct,/silent
	!p.background=!d.table_size-1
	!p.color=0
	!p.font=-1
	!p.thick=1.0
	!x.thick=1.0
	!y.thick=1.0
END

PRO psc
	device,/close
END

PRO helpr,string
	help, /rout,names=string
END

FUNCTION ev2ang,energy
	h=6.626068e-34	 ;m^2 kg/s
	e=1.60217646e-19 ;J/eV
	c=2.99792458e8	 ;m/s
	;E=hc/lam
	ang=h*c/energy/e*1.0e10
	RETURN, ang
END

FUNCTION ang2ev,wave
	h=6.626068e-34	 ;m^2 kg/s
	e=1.60217646e-19 ;J/eV
	c=2.99792458e8	 ;m/s
	;E=hc/lam
	energy=h*c/(wave*1.0e-10)/e
	RETURN, energy
END

FUNCTION cm2ang,cm
	ang=1.0/cm*1.0e8
	RETURN,ang
END

FUNCTION cm2ev,cm
	ang=1.0/cm*1.0e8
	energy=ang2ev(ang)
	RETURN,energy
END

FUNCTION dcos,th
	output=cos(float(th)/180.0*!pi)
	RETURN,output
END

FUNCTION dsin,th
	output=sin(float(th)/180.0*!pi)
	RETURN,output
END

FUNCTION dtan,th
	output=tan(float(th)/180.0*!pi)
	RETURN,output
END

PRO lorentz_fit,x,A,F
	gamma=A[0]
	xo=A[1]
	h=A[2]
	L=gamma/2.0/((x-xo)^2+(0.5*gamma)^2)
	F=L*h/max(L)
	
END

PRO two_gauss_fit,x,A,F
	z1=(x-A[1])/A[2]
	z2=(x-A[4])/A[5]
	F=A[0]*exp(-z1^2/2)+A[3]*exp(-z2^2/2)+A[6]+A[7]*x
END

PRO voigt_fit,x,A,F
	F=voigt(A[0],x-A[1])*A[2]
END
	

FUNCTION bessel_zeros,file=file
	
	IF keyword_set(file) THEN BEGIN
		path='/home/mlreinke/idl/genie/data/bessel_zeros'
		openr,lun,path, /get_lun
		line=strarr(1)
		readf,lun,line
		tmp=strsplit(line,'=',/extract)
		bessel_zeros=fltarr(tmp[1]+1)
		FOR i=0,tmp[1] DO BEGIN
			readf,lun,line
			tmp=strsplit(line,',',/extract)
			bessel_zeros[i]=tmp[1]
		ENDFOR
		close,lun
		free_lun,lun
	ENDIF ELSE BEGIN
		path='/home/mlreinke/idl/genie/data/bessel_zeros.dat'
		restore, path
	ENDELSE
	
	output=bessel_zeros
	RETURN,output
END

PRO greeklabels, str,ps=ps
	IF keyword_set(ps) THEN title='!9'+str+'!3' ELSE title='!7'+str+'!3'
	plot,[0],[0],title=title,chars=1.5
END

	
;HOW TO CONCATENATE IN ONE LINE
;a=intarr(4,4)+1
;b=intarr(4)
;c=[a,transpose(b)]	;adds as another i-index
;print, size(c,/dim)
;           5           4
;
;c=[[a],[b]]		;adds as another j-index
;print, size(c, /dim)
;           4           5
