;6/6/12 - modified to use version controlled atomic physics directory
PRO jphysb_cvs2dat
	;wd='/home/mlreinke/idl/impurities/data/jphysb_39_4457/'
	wd='/usr/local/cmod/idl/atomic_physics/jphysb_39_4457/'
	csv_convert,wd+'kr_ion.csv',wd+'kr_ion.dat'
	csv_convert,wd+'kr_rec.csv',wd+'kr_rec.dat'
	csv_convert,wd+'mo_ion1.csv',wd+'mo_ion1.dat'
	csv_convert,wd+'mo_ion2.csv',wd+'mo_ion2.dat'
	csv_convert,wd+'mo_rec.csv',wd+'mo_rec.dat'

	;KRYPTON IONIZATION DATA
	restore,wd+'kr_ion.dat'
	ion_c=fltarr(37,5)
	FOR i=0,4 DO out=execute('ion_c[*,i]=c'+num2str(i+1,1))
	ion_e=fltarr(37,5)
	FOR i=0,4 DO out=execute('ion_e[*,i]=e'+num2str(i+1,1))
	;KRYPTON RECOMBINATION DATA
	restore,wd+'kr_rec.dat'
	rr=fltarr(37,4)
	rr[*,0]=a
	rr[*,1]=b
	rr[*,2]=t0
	rr[*,3]=t1
	dr_c=fltarr(37,7)
	FOR i=0,6 DO out=execute('dr_c[*,i]=c'+num2str(i+1,1))
	dr_e=fltarr(37,7)
	FOR i=0,6 DO out=execute('dr_e[*,i]=e'+num2str(i+1,1))
	save,ion_c,ion_e,rr,dr_c,dr_e,filename=wd+'kr.dat'

	;MOLYBDENUM IONIZATION DATA	
	restore,wd+'mo_ion1.dat'
	ion_c=fltarr(43,4)
	FOR i=0,3 DO out=execute('ion_c[*,i]=c'+num2str(i+1,1))
	ion_e=fltarr(43,4)
	FOR i=0,3 DO out=execute('ion_e[*,i]=e'+num2str(i+1,1))
	restore,wd+'mo_ion2.dat'
	n_sup=n(q)+1
	ion_supp=fltarr(n_sup,6)
	ion_supp[*,0]=q
	ion_supp[*,1]=ip
	ion_supp[*,2]=a
	ion_supp[*,3]=b
	ion_supp[*,4]=c
	ion_supp[*,5]=d

	;MOLYBDENUM RECOMBINATION DATA
	restore,wd+'mo_rec.dat'
	rr=fltarr(43,4)
	rr[*,0]=a
	rr[*,1]=b
	rr[*,2]=t0
	rr[*,3]=t1
	dr_c=fltarr(43,7)
	FOR i=0,6 DO out=execute('dr_c[*,i]=c'+num2str(i+1,1))
	dr_e=fltarr(43,7)
	FOR i=0,6 DO out=execute('dr_e[*,i]=e'+num2str(i+1,1))

	save,ion_c,ion_e,ion_supp,rr,dr_c,dr_e,filename=wd+'mo.dat'
END

FUNCTION jphysb_f1,x
	f1=exp(x)*expint(1,x)
	RETURN,f1
END

FUNCTION jphysb_f2,x
	;pj=[1.0,2.1658e2,2.0336e4,1.0911e6,3.7114e7,8.3963e8,1.2889e10,1.3449e11,9.4002e11,4.2571e12,1.1743e13,1.7549e13,1.0806e13,4.9776e11]
	;qj=[1.0,2.1958e2,2.0984e4,1.1517e6,4.0349e7,9.4900e8,1.5345e10,1.7182e11,1.3249e12,6.9071e12,2.3531e13,4.9432e13,5.7760e13,3.0225e13,3.3641e12]
	;p=dblarr(n(x)+1)
	;FOR j=0,n(pj) DO p+=pj[j]/x^j
	;q=dblarr(n(x)+1)
	;FOR j=0,n(qj) DO q+=qj[j]/x^j
	;f2=1.0/x^2*p/q
	;f2=float(f2)
	f2=fltarr(n(x)+1)
	amat=[[0.30492e6,0.10672e4,-0.86781,-0.32275e6,0.97598e6,0.35238e5],$
	     [0.60050e6,0.16560e4,-0.89271,0.63419e6,0.15300e7,0.98814e5],$
	     [0.68406e6,0.45186e4,-0.46548e1,0.68665e6,0.19197e7,-0.3703e5]]	
	tmp=where(x LT 0.2)
	IF tmp[0] NE -1 THEN f2[tmp]=1.0/x[tmp]^2*(x[tmp]^3+amat[0,0]*x[tmp]^2+amat[1,0]*x[tmp]+amat[2,0])/(x[tmp]^3+amat[3,0]*x[tmp]^2+amat[4,0]*x[tmp]+amat[5,0])
	tmp=where(x GE 0.2 AND x LT 1.0)
	IF tmp[0] NE -1 THEN f2[tmp]=1.0/x[tmp]^2*(x[tmp]^3+amat[0,1]*x[tmp]^2+amat[1,1]*x[tmp]+amat[2,1])/(x[tmp]^3+amat[3,1]*x[tmp]^2+amat[4,1]*x[tmp]+amat[5,1])
	tmp=where(x GE 1.0)
	IF tmp[0] NE -1 THEN f2[tmp]=1.0/x[tmp]^2*(x[tmp]^3+amat[0,2]*x[tmp]^2+amat[1,2]*x[tmp]+amat[2,2])/(x[tmp]^3+amat[3,2]*x[tmp]^2+amat[4,2]*x[tmp]+amat[5,2])
	;f2=fltarr(n(x)+1)
	;a=10.0^(make(0,3,3000))
	;FOR i=0,n(x) DO f2[i]=exp(x[i])*int_tabulated(a,alog(a)/a*exp(-1.0*x[i]*a))
	;f2=(-0.0005725+0.01345/x+0.8691/x^2+0.03404/x^3)/(1.0+2.197/x+0.2475/x^2+0.002053/x^3)
	;stop
	RETURN,f2
END

;6/6/12 - modified to use version controlled atomic physics directory
FUNCTION jphysb_dat2ionrec,z,temp,debug=debug
	IF z NE 42 AND z NE 36 THEN RETURN,-1
	ntemp=n(temp)+1
	ion=dblarr(z+1,ntemp)
	rr_rec=dblarr(z+1,ntemp)
	dr_rec=dblarr(z+1,ntemp)
	;wd='/home/mlreinke/idl/impurities/data/jphysb_39_4457/'
	wd='/usr/local/cmod/idl/atomic_physics/jphysb_39_4457/'
	CASE z OF 
		36 : path=wd+'kr.dat'
		42 : path=wd+'mo.dat'
	ENDCASE
	restore,path

	FOR i=0,z DO BEGIN
		;ionization
		nterms=n(ion_c[0,*])+1
		CASE total(ion_c[i,*]) OF 
			0 : BEGIN
				;do nothing
			END
	
			-1.0*nterms : BEGIN
				j=min(where(ion_supp[*,0] EQ i))
				x=double(ion_supp[j,1]/temp)
				f1=jphysb_f1(x)
				f2=jphysb_f2(x)
				F=ion_supp[j,2]*(1.0-x*f1)+ion_supp[j,3]*(1.0+x-x*(2.0+x)*f1)+ion_supp[j,4]*f1+ion_supp[j,5]*x*f2 ;(eq3)
				F*=1.0e-14
				ion[i,*]=6.69e7*exp(-1.0*x)/temp^(1.5)*F/x
			
			END
		
			-2.0*nterms : BEGIN
				a=4.5e-14
				j=min(where(ion_supp[*,0] EQ i))
				nelec=z-i
				IF nelec LE 2 THEN nq=nelec
				IF nelec GT 2 AND nelec LE 4 THEN nq=nelec-2
				IF nelec GT 4 AND nelec LE 10 THEN nq=nelec-4
				x=ion_supp[j,1]/temp
				ion[i,*]=6.69e7*a*nq/temp^(1.5)*expint(1,x)/x		;ionization of the outer electron (eq8)
				IF nelec GT 4 THEN BEGIN
					be_like=min(where(ion_supp[*,0] EQ z-4))
					x=ion_supp[be_like,1]/temp					
					ion[i,*]+=6.69e7*a*2.0/temp^(3/2.0)*expint(1,x)/x	;inner shell ionization of the 2s electron (eq8)
				ENDIF
			END

			ELSE : BEGIN
				FOR j=0,nterms-1 DO ion[i,*]+=1.0e-11/sqrt(temp)*ion_c[i,j]*exp(-ion_e[i,j]/temp)	;(eq9)
			END
		ENDCASE
	
		;dielectronic recombination
		nterms=n(dr_c[0,*])+1
		CASE total(dr_c[i,*]) OF 
			0 : BEGIN
				;do nothing
			END
	
			ELSE : BEGIN
				FOR j=0,nterms-1 DO dr_rec[i,*]+=1.0/(temp)^(1.5)*dr_c[i,j]*exp(-dr_e[i,j]/temp) ;(eq11)
			END
		ENDCASE
		
		;radiative recombination
		CASE total(rr[i,*]) OF 
			0 : BEGIN
				;do nothing
			END
	
			ELSE : BEGIN
				rr_rec[i,*]=rr[i,0]/(sqrt(temp/rr[i,2])*(1.0+sqrt(temp/rr[i,2]))^(1.0-rr[i,1])*(1.0+sqrt(temp/rr[i,3]))^(1.0+rr[i,1])) ;(eq10)
			END
		ENDCASE
	ENDFOR 
	
	output={z:z, temp:temp,ion:ion*1.0e-6,rr:rr_rec*1.0e-6,dr:dr_rec*1.0e-6}	;all rates in m^3/s
	IF keyword_set(debug) THEN stop
	RETURN,output

END

;6/6/12 - modified to use version controlled atomic physics directory
FUNCTION read_matt06_frac_data,z
	;path='/home/mlreinke/idl/impurities/data/jphysb_39_4457/'
	path='/usr/local/cmod/idl/atomic_physics/jphysb_39_4457/'
	CASE z OF 
		36 : BEGIN
			path=path+'kr_frac.dat'
			str='Kr'
		END
		42 : BEGIN
			path=path+'mo_frac.dat'
			str='Mo'
		END
	ENDCASE
	restore,path
	
	temp=10^te
	frac=fltarr(z+1,n(temp)+1)
	FOR i=0,z DO err=execute('frac[i,*]='+str+num2str(i+1,1))
	output={frac:frac,temp:temp}
	RETURN,output
END
	


FUNCTION read_matt06_rec_data,z,nt=nt
	IF NOT keyword_set(nt) THEN nt=400
	temp=10^(make(0.0,4.0,nt))
	out=jphysb_dat2ionrec(z,temp)
	dens=[1.0e17,1.0e22]
	path='jphysb_dat2ionrec'
	rec=out.rr+out.dr
	output={rates:[[[rec]],[[rec]]],temp:temp,dens:dens,z:z,path:path}
	RETURN,output
END

FUNCTION read_matt06_ion_data,z,nt=nt
	IF NOT keyword_set(nt) THEN nt=400
	temp=10^(make(0.0,4.0,nt))
	out=jphysb_dat2ionrec(z,temp)
	dens=[1.0e17,1.0e22]
	path='jphysb_dat2ionrec'
	ion=out.ion
	output={rates:[[[ion]],[[ion]]],temp:temp,dens:dens,z:z,path:path}
	RETURN,output
END
