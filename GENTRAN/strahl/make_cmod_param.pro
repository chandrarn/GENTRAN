;+
;NAME:
;	MAKE_CMOD_PARAM
;
;PURPOSE:
;	This procedure to generate the setup file for use
;	with the STRAHL impurity transport code
;
;CALLING SEQUENCE:
;	MAKE_CMOD_PARAM,shot,z,t1,t2
;	
;INPUTS:
;	shot	LONG	shot number
;	z	INT	atomic number of impurity 
;	t1	FLOAT	time to start simulation
;	t2	FLOAT 	time to end simulation
;
;OPTIONAL INPUTS:
;	dt	FLOAT	time step DEFAULT: 0.1 ms
;	ncyc	INT	number of time steps per cycle DEFAULT: 10
;	ion	STRING	of the background ion DEAFULT 'D' (deuterium)
;	tneut	FLOAT	temperature at which the neutrals enter the domain DEFAULT 1 eV
;	niter	INT	number of iterations DEFAULT=2000
;	fiter	FLOAT	stop condition for frac. change DEFAULT=0.02
;	start	STRUC	of shot/time to use as a starting point for simulation DEFAULT not used
;		*.shot	shot number of reference
;		*.time	time at which to start
;		
;	******* RADIAL GRID OPTIONS *********
;	k	INT	equally spaced points in rho^k DEFAULT=1.5
;	nrho	INT	number of points in radial grid DEFAULT=100
;	dr	FLTARR	[dr_cent,dr_edge] if > 0 enables alternate
;			gridding as described in (1.37) in STRAHL Manual (IPP-10/30) DEFAULT [-1.0,-1.0]
;
;	******* PLASMA TRANSPORT OPTIONS *********	
;	### See GENTRAN_DIFFCONV_PROFILES for details ###
;	exp	INT	of the pre-determined transport "experiments" to use for diff/conv profiles DEFAULT = 1
;	dff	FLOAT	"fudge-factor" used to multiply diffusion profiles DEFAULT=0.5
;	cff	FLOAT	"fudge-factor" used to multiply convection profiles DEFAULT=1.0
;	doff	FLOAT	DC offset in diffusion profile DEFAULT=0.01 [m^2/s]
;	
;	### Directly specify the profiles ###
;	diff	STRUC	of the diffusion profile
;		*.diff	[nrho,ntransp]	of the diffusion profile [m^2/s]
;		*.rho	[nrho] of the radial grid [norm. pol. flux]
;		*.time	[ntransp] of time grid [sec]
;	conv	STRUC	of the convection profile
;		*.conv	[nrho,ntransp]	of the convection profile [m/s]
;		*.rho	[nrho] of the radial grid [norm. pol. flux]
;		*.time	[ntransp] of time grid [sec]
;
;	If these profiles are fixed in time *.time can be ommitted and *.diff (*.conv) need only be [nrho]
;
;	### Sawtooth Modeling ###
;	saw	STRUC	of the 
;
;	******* IMPURITY SOURCE OPTIONS *********	
;	fz	FLOAT	scaling of the source rate DEFAULT=1.0
;	source	FLOAT/STRING/STRUC that sets the time evolving impurity source	
;	
;	source	FLOAT	source=0.0 means fz*10^17 1/s is on continuously from start to finish DEFAULT
;			source > 0 means fz*10^17 1/s is on and stays on after t=source
;			source < 0 means fz*10^17 1/s is on for 1 ms, simulating an injection
;	source	STRING	assumed to be a path to file properly formatted for STRAHL, and is merely copied to spath
;	source	STRUC	structure of the time-evolving source
;		*.(0)	[ntau] time evolving source [1/s]
;		*.(1)	[ntau] time vector of the source [sec]
;
;OUTPUTS:
;	All data are written to an ASCII file formatted as per the instructions in the STRAHL manual.  This file written to the
;	users home directory in ~/strahl/cmod/run_SHOTNUMBER.0 and, if necessary, the source file is also written to the same
;	directory as ZzflxSHOTNUMBER.dat
;
;OPTIONAL OUTPUTS:
;	filepath	STRING	of the path to the parameter file
;	spath		STRING	of the path to the source file (if generated )
;	diff		STRUC 	(see above) available as an output if not set as input
;	conv		STRUC 	(see above) available as an output if not set as input
;		
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - January 2013 (adapted from N.T. Howards's WRITE_STRAHL_FILES.PRO)
;	1/14/13		M.L. Reinke - modified the transport profiles to accept time-evolving diff/conv data
;	1/16/13		M.L. Reinke - added the fz optional input, changed PSIN max to 0.99 for
;                                     neoclassical computation and made NEOART=2 the default option
;	10/12/14	M.L. Reinke - modified the steady-state source to allow for lower rates
;       11/11/14	M.L. Reinke - added a check to interpolate diff/conv onto a 100-pt radial grid to avoid a bug in STRAHL                             
;
;-

PRO make_cmod_param,shot,z,t1,t2,dt=dt,ncyc=ncyc,ion=ion,tneut=tneut,start=start,niter=niter,fiter=fiter,k=k,nrho=nrho,dr=dr,source=source,tau=tau,fz=fz,exp=exp,$
	dff=dff,cff=cff,doff=doff,diff=diff,conv=conv,neo=neo,saw=saw,debug=debug,filepath=filepath,spath=spath

	IF NOT keyword_set(dt) THEN dt=1.0e-4			;time step is 0.1 ms
	IF NOT keyword_set(ncyc) THEN ncyc=10			;number of time steps per cycle
	IF NOT keyword_set(tneut) THEN tneut=1.0		;neutrals enter at 1 eV
	IF NOT keyword_set(ion) THEN ion='D'			;set main ion to deuterium as default
	elemstr=num2elem(z)
	IF strlen(elemstr) EQ 1 THEN elemstr+='_'
	zmass=read_atomic_mass(z)
	imass=read_atomic_mass(ion)
	ichrg=read_atomic_charge(ion)
	IF NOT keyword_set(start) THEN start={shot:0 ,time:0.0}	;set optional shot/time to start simulation at
	IF start.shot EQ 0 THEN istart='0' ELSE istart='1'
	zrecycl=[2,10,18,36,54]
	IF total(where(zrecycl EQ z)) EQ -1 THEN irecycl='0' ELSE irecycl='1'
	IF NOT keyword_set(niter) THEN niter=200		;number of iterations	(even on transient, looks like only 
	IF NOT keyword_set(fiter) THEN fiter=0.001		;iteration stop condition (%change)
	IF NOT keyword_set(k) THEN k=1.5			;non-uniform radial grid (density increases toward edge for k > 1)
	IF NOT keyword_set(nrho) THEN nrho=100			;number of radial points
	IF NOT keyword_set(dr) THEN dr=[-1.0,-1.0]		;[dr_center,dr_edge] alternate way of determining grid	
	IF keyword_set(neo) THEN ineo='100' ELSE ineo='-100'	;use neoclassical transport if /neo

	IF NOT keyword_set(source) THEN source=0.0		;if source not specified, assume steady-state edge source
	IF NOT keyword_set(fz) THEN fz=1.0			;keep 10^17 1/s as the injection rate
	x=size(source,/type)
	spath=strcompress('/home/'+logname()+'/strahl/cmod/'+elemstr+'flx'+string(shot)+'.dat',/remove_all)
	CASE x[0] OF 
		4 : BEGIN	;specify injection/turn-on
			IF source[0] NE 0.0 THEN BEGIN
				IF NOT keyword_set(tau) THEN tau=1.0e-3				;length of transient phase
				time=[t1,abs(source[0])-dt,make(abs(source[0]),abs(source[0])+tau,tau/dt),abs(source[0])+tau+dt,t2]
				inj=time*0.0
				npts=n(time)+1
				IF source[0] GT 0 THEN BEGIN		;source turns on and stays on after t=source[0]
					tmp=where(time GE source[0])
					inj[tmp]=1.0
				ENDIF
				IF source[0] LT 0 THEN BEGIN		;source turns on for 1 ms (injection) after t=abs(source[0])
					tmp=where(time GE abs(source[0]) AND time LE abs(source[0])+tau)
					inj[tmp]=1.0
                                ENDIF
				inj*=1.0e17*fz				;scale the source

				openw,5,spath
				printf,5,num2str(int(npts),1)
				FOR i=0,npts-1 DO printf,5,'     '+strtrim(time[i],2)+'   '+strtrim(inj[i],2)
				close,5
				free_lun,5
				isource='1'
                        ENDIF ELSE BEGIN
				isource='0'
				spath=-1
			ENDELSE
                END
		7 : BEGIN	;specified in file, copy to staging location
			path=spath
			spawn, 'cp '+source+' '+path
			isource='1'	
                END
		8 : BEGIN	;source input as a structure written to file
			openw,5,spath
			npts=n(source.(0))+1
			printf,5,num2str(npts,1)
			FOR i=0,npts-1 DO printf,5,'     '+strtrim(source.(1)[i],2)+'   '+strtrim(source.(0)[i],2)
			close,5
			free_lun,5
			isource='1'
                END
		ELSE : BEGIN
			print, 'ERROR - specification of impurity source not recognized'
			RETURN
		END
	ENDCASE
	
       ;determine #of time points in transport profile, generate from GENTRAN if not specified
	IF NOT keyword_set(diff) AND NOT keyword_set(conv) THEN BEGIN
		ntransp=1
		ttransp=[0.0]
		IF NOT keyword_set(exp) THEN exp=1
		IF NOT keyword_set(dff) THEN dff=0.5
		IF NOT keyword_set(cff) THEN cff=1.0
		IF NOT keyword_set(doff) THEN doff=0.01
		gentran_diffconv_profiles,exp,diff,conv,dff=dff,cff=cff,doff=doff,shot=shot,time=0.5*(t1+t2),/psin
        ENDIF ELSE BEGIN
		x=size(diff.(0))
		IF x[0] EQ 2 THEN ntransp=x[2] ELSE ntransp=1 
		y=size(conv.(0))
		IF total((x[0:2]-y[0:2])) NE 0 THEN BEGIN
			print, 'ERROR - size mismatch between diff and conv structures'
			RETURN
                ENDIF
		IF n_tags(diff) EQ 3 THEN ttransp=diff.(2) ELSE ttransp=[0.0]
	ENDELSE	
	dvmax=100
	IF n(diff.rho) GT dvmax-1 THEN BEGIN		;apparently a hard limit on the size of the radial grid
		x=size(diff.(0))
		drho=make(0.0,last(diff.(1)),dvmax)
		crho=make(0.0,last(conv.(1)),dvmax)
		CASE x[0] OF 
			1 : BEGIN
				odiff={diff:interpol(diff.(0),diff.(1),drho),rho:drho}
				oconv={conv:interpol(conv.(0),conv.(1),drho),rho:crho}
                        END
			2 : BEGIN
				xdiff=fltarr(dvmax,x[2])
				xconv=fltarr(dvmax,x[2])
				FOR i=0,x[2]-1 DO BEGIN
					xdiff[*,i]=interpol(diff.(0)[*,i],diff.(1),drho)
					xconv[*,i]=interpol(conv.(0)[*,i],conv.(1),crho)
				ENDFOR
				odiff={diff:xdiff,rho:drho,time:diff.(2)}
				oconv={conv:xconv,rho:crho,time:conv.(2)}
			END
		ENDCASE
        ENDIF ELSE BEGIN
		odiff=diff
		oconv=conv
	ENDELSE

	paraloss=ttransp*0.0+2.5	;hard code 2.5 ms parallel loss at each time point
	ndrho=n(odiff.(0)[*,0])+1
	ncrho=n(oconv.(0)[*,0])+1

	IF keyword_set(saw) THEN BEGIN
		sawt=saw.(0)
		sawrad=saw.(1)
		nsaw=n(sawt)+1
        ENDIF ELSE BEGIN
		nsaw=0
		sawt=0.0
		sawrad=1.0
	END
	
	;stop before opening and writing the param file
	IF keyword_set(debug) THEN stop		
	filepath=strcompress('/home/'+logname()+'/strahl/cmod/run_'+string(shot)+'.0',/remove_all)
	openw,4,filepath
	printf,4,'                    E L E M E N T'
	printf,4,'cv    element   atomic weight(amu)   energy of neutrals(eV)'
	printf,4,"      '"+elemstr+"'            "+num2str(zmass,dp=2)+"                   "+num2str(tneut,dp=2)
	printf,4,''
	printf,4,'cv    main ion:  atomic weight(amu)    charge'
	printf,4,'                       '+num2str(imass,dp=3)+'               '+num2str(ichrg,1)   
	printf,4,''
	printf,4,'                    G R I D - F I L E'
	printf,4,'cv   shot         index'
	printf,4,'    '+strtrim(shot,2)+'      0'
	printf,4,''
	printf,4,'                    G R I D  P O I N T S  A N D  I T E R A T I O N'
	printf,4,'cv    K    number of grid points  dr_center(cm)  dr_edge(cm)'
	printf,4,'     '+num2str(k,dp=1)+'            '+num2str(nrho,1)+'               '+num2str(dr[0],dp=2)+'            '+num2str(dr[1],dp=2)
	printf,4,''
	printf,4,'cv         max iterations at fixed time      stop iteration if change below (%)'
	printf,4,'                      '+num2str(niter,1)+'                                  '+num2str(fiter,dp=3) 
	printf,4,''
	printf,4,'                    S T A R T  C O N D I T I O N S'
	printf,4,'cv    start new=0/from old calc=1   take distr. from shot   at time'
	printf,4,'                 '+istart+'               '+num2str(start.shot,1)+'        '+num2str(start.time,dp=3)
	printf,4,''
	printf,4,''
	printf,4,'                    O U T P U T'
	printf,4,'cv    save all cycles = 1, save final and start distribution = 0'
	printf,4,'            1'
	printf,4,''
	printf,4,'                    T I M E S T E P S'
	printf,4,'cv    number of changes(start-time+....+stop-time)'
	printf,4,'                  2'
	printf,4,''
	printf,4,'cv    time   dt at start   increase of dt after cycle   steps per cycle'
	printf,4,'    '+strtrim(t1,2)+'     '+num2str(dt,dp=5)+'               1.001                      '+num2str(ncyc,1)
	printf,4,'    '+strtrim(t2,2)+'     '+num2str(dt,dp=5)+'               1.001                      '+num2str(ncyc,1)
	printf,4,''
	printf,4,'                    S O U R C E'
	printf,4,'cv    position(cm)    constant rate (1/s)    time dependent rate from file'
	printf,4,'         90.5            '+num2str(1.0e17*fz,dp=6)+'                         '+isource
	printf,4,''
	printf,4,'cv    divertor puff   source width in(cm)    source width out(cm)'
	printf,4,'           0                  -1                       -1
	printf,4,''
	printf,4,'                    E D G E ,  R E C Y C L I N G'
	printf,4,'cv    decay length of impurity outside last grid point (cm)'
	printf,4,'                              1.0 '
	printf,4,''
	printf,4,'cv    Rec.:ON=1/OFF=0   wall-rec   Tau-div->SOL(ms)   Tau-pump(ms)'
	printf,4,'             0             '+irecycl+'           1.0             1000.'
	printf,4,''
	printf,4,'cv    SOL=width(cm)'
	printf,4,'          1'
	printf,4,''
	printf,4,'                    D E N S I T Y,  T E M P E R A T U R E,  A N D  N E U T R A L  H Y D R O G R E N  F O R  C X'
	printf,4,'cv    take from file with:    shot      index'
	printf,4,'                           '+strtrim(shot,2)+'     0'
	printf,4,''
	printf,4,'                    N E O C L A S S I C A L  T R A N S P O R T'
	printf,4,'                    method'
	printf,4,'    0 = off,    >0 = % of Drift,   1 = approx'
	printf,4,'cv  <0 = figure out, but dont use  2/3 = NEOART   neoclassics for rho_pol <'
	printf,4,'                '+ineo+'                      2                   0.99'
	printf,4,''
	printf,4,'                    A N A M A L O U S  T R A N S P O R T'
	printf,4,'cv   # of changes for transport'
	printf,4,'      '+num2str(ntransp,1)
	printf,4,''
	printf,4,'cv   time-vector'
	printf,4,ttransp
	printf,4,''
	printf,4,'cv   parallel loss times (ms)'
	printf,4,paraloss
	printf,4,''
	printf,4,'cv   Diffusion [m^2/s]'
	printf,4,"     'interp'"
	printf,4,''
	printf,4,'cv   # of interpolation points'
	printf,4,'        '+num2str(ndrho,1)
	printf,4,''
	printf,4,'cv   rho_pol grid'
	printf,4,sqrt(odiff.(1)[*,0])
	printf,4,''
	printf,4,'cv   Diffusion Coefficient Grid'
	printf,4,odiff.(0)
	printf,4,''
	printf,4,'cv   Drift function        only for drift'
	printf,4,"     'interp'             'velocity'"
	printf,4,''
	printf,4,'cv   # of interpolation points'
	printf,4,'            '+num2str(ncrho,1)
	printf,4,''
	printf,4,'cv   rho_pol grid'
	printf,4,sqrt(oconv.(1)[*,0])
	printf,4,''
	printf,4,'cv   Velocity Coefficient Grid'
	printf,4,''
	printf,4,oconv.(0)
	printf,4,''
	printf,4,'cv   # of sawteeth       inversion radius (cm)'
	printf,4,'          '+num2str(nsaw,1)+'                       '+num2str(sawrad,dp=2)
	printf,4,''
	printf,4,'cv   times of sawteeth'
	printf,4,sawt
	close,4
	free_lun,4

END

 
