; modified from /usr/local/cmod/idl/GENIE/GENRAD/genrad_hirexsr_thaco.pro

;+
;NAME:
;	HIREXSR_GENRAD_ARGON
;
;PURPOSE:
;	This procedure loads experimental Argon data from HIREXSR and computes
;	the equivilent using data from an impurity transport simulation
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 1/18/13
;       Modified by:    C. Gao        8/25/14
;                       C. Gao       10/11/14: add the tangency radius of the chord norm minor radius as output (brRho)
;                                              fixed a typo on index (i1,i2) "FOR i=0,i2-i1 DO BEGIN"

;Example:
;result=hirexsr_genrad_argon_cgao(1140221016,0,0.7,1.4,csden,cserr,temp,temperr,dens,denserr,rhovec,filename='Ar_1140221016_1.sav')
;result=hirexsr_genrad_argon_cgao(1140221016,[0,2],0.7,1.4,csden,cserr,temp,temperr,dens,denserr,rhovec,filename='Ar_1140221016_1.sav')


FUNCTION hirexsr_genrad_argon_cgao,shot,line,t1,t2,csden,cserr,temp,temperr,dens,denserr,rhovec,tht=tht,eptr=eptr,efit=efit,sptr=sptr,data=data,filename=filename,nobr=nobr
	IF keyword_set(data) THEN BEGIN
           print,'loading data from var:data......'
        ENDIF ELSE IF keyword_set(filename) THEN BEGIN
           path='~/strahl/result/'
           print,'loading data from '+path+filename+'......'
           restore,path+filename,/v
        ENDIF ELSE BEGIN
           MESSAGE,'No data or file specified. Exit program.'
        ENDELSE
        index=where(data.time ge t1 and data.time le t2)
        time=data.time[index]
        csden=data.csden[*,*,index]
        cserr=data.cserr[*,*,index]
        temp=data.temp[*,index]
        temperr=data.terr[*,index]
        dens=data.dens[*,index]
        denserr=data.derr[*,index]
        rhovec=data.rho

        IF NOT keyword_set(tht) THEN tht=0
	nlines=n(line)+1

	;get efit_data if not included as optional input
	IF NOT keyword_set(efit) THEN BEGIN
		mdsopen,'analysis',shot
		rmid=mdsvalue('\efit_rmid')
		efit_t=mdsvalue('dim_of(\efit_rmid)')
		ipsin=mdsvalue('dim_of(\efit_rmid,1)')
		mdsclose,'analysis',shot
		i1=ipt(efit_t,t1)
		i2=ipt(efit_t,t2)
		ro=min(sum_array(rmid[i1:i2,*],/j)/(i2-i1+1.0))
		a=max(sum_array(rmid[i1:i2,*],/j)/(i2-i1+1.0))-ro
		irmid=reform(sum_array(rmid[i1:i2,*],/j))/(i2-i1+1.0)
		irho=(irmid-ro)/a
		efit={rmid:irmid,rho:irho,psin:ipsin,ro:ro,a:a}
        ENDIF 
	
	
        ;load emissivity and brightness profiles from experiment if not included as optional inputs
	IF NOT keyword_set(eptr) THEN BEGIN
		eptr=ptrarr(nlines,/allocate_heap)
		FOR j=0,nlines-1 DO BEGIN
			hirexsr_load_line_pos,shot,line[j],pos,tpos=tpos,tht=tht		;load POS vector
			tmp=where(pos[1,*] NE -1)
			pos=pos[*,tmp]
			CASE line[j] OF 
				0 : BEGIN
					tit='XICS: Ar!u16+!n'
					label='w'
				END
				2 : BEGIN
					tit='XICS: Ar!u16+!n'
					label='z'
				END
				3 : BEGIN
					tit='XICS: Ar!u17+!n'
					label='lya'
				END
			ENDCASE


			hirexsr_load_mlintptr,shot,line[j],lint,tau,tht=tht		;load brightness profile
			i1=ipt(tau,t1)
			i2=ipt(tau,t2)
			ilint=*lint[i1]
                        bright=fltarr(n(ilint[*,8])+1,i2-i1+1)
                        brerr=bright
                        brRho = bright
			;bright=ilint[*,8]
			;brerr=ilint[*,9]^2

			;FOR i=0,i2-i1 DO BEGIN
			FOR i=i1,i2 DO BEGIN
                                ilint=*lint[i]
				bright[*,i-i1]=ilint[*,8]
				brerr[*,i-i1]=ilint[*,9]^2
                                brRho[*,i-i1]=ilint[*,10]
                	ENDFOR
			heap_free,lint
			;bright/=(i2-i1+1.0)
			;brerr=sqrt(brerr)/(i2-i1+1.0)

			hirexsr_load_profile,shot,line[j],prof,err,psin,tau,tht=tht	;load emissivity profile
			i1=ipt(tau,t1)
			i2=ipt(tau,t2)
		
			FOR i=i1,i2 DO BEGIN		;check to make sure radial grid is constant over time
				chk=total(psin[*,i]-psin[*,i1])
				IF chk NE 0 THEN BEGIN
					message, 'time-evolving radial grid for THT='+num2str(tht,1)+', - incompatible'
		                	ENDIF
			ENDFOR
	
			;check for m=1 term, isolate m=0 only
			npsi=n(psin[*,i1])+1
			IF psin[0] EQ psin[npsi/2] THEN ism=1 ELSE ism=0
			CASE ism OF
				0 : BEGIN
					emiss=prof[*,i1:i2,0]
					emerr=err[*,i1:i2,0]
					emrho=interpol(irho,ipsin,psin[*,i1])		;convert to r/a grid
        	        	END
				1 : BEGIN
					emiss=prof[0:npsi/2-1,i1:i2,0]
					emerr=err[0:npsi/2-1,i1:i2,0]
					emrho=interpol(irho,ipsin,psin[0:npsi/2-1,i1])	;convert to r/a grid	
		                END
			ENDCASE
			;emiss=sum_array(emiss,/i)/(i2-i1+1.0)
			;emerr=sqrt(sum_array(emerr^2,/i))/(i2-i1+1.0)
			ied={emiss:emiss,emerr:emerr,rho:emrho,time:tau[i1:i2],bright:bright,brerr:brerr,brRho:brRho,pos:pos,shot:shot,t1:t1,t2:t2,line:line,tht:tht,tit:tit,label:label,enorm:0.0,bnorm:0.0}
			*eptr[j]=ied
		ENDFOR
             ENDIF

	;calculate emissivity and brightness profiles from simulations
	IF NOT keyword_set(sptr) THEN sptr=ptrarr(nlines,/allocate_heap)
	FOR j=0,nlines-1 DO BEGIN
           IF size(*sptr[j],/type) NE 0 THEN BEGIN ;restore from sdata
              sdata=*sptr[j]
              pos=sdata.pos
           ENDIF ELSE BEGIN
              hirexsr_load_line_pos,shot,line[j],pos,tpos=tpos,tht=tht ;load POS vector
              tmp=where(pos[1,*] NE -1)
              pos=pos[*,tmp]
           ENDELSE
           nt=n(time)+1
           nr=n(rhovec)+1
           emsim=fltarr(nr,nt)
           print, 'Calculating Simulated Emissivity Profile: line = '+num2str(line[j],1)
           for it=0,nt-1 do begin
              csden_tmp=csden[*,*,it]
              temp_tmp=temp[*,it]
              dens_tmp=dens[*,it]
              x=size(rhovec)
              y=size(csden_tmp)
              IF y[1] EQ x[1] THEN BEGIN
                 csden_tmp=rotate(csden_tmp,4)
              ENDIF
              CASE line[j] OF 
                 0 : BEGIN
                    calc_ar_line_rates,'w',temp_tmp/1.0e3,q,rates
                    emsim[*,it]=dens_tmp*csden_tmp[q[0],*]*rates[0,*]
                    emsim[*,it]+=dens_tmp*csden_tmp[q[1],*]*rates[1,*]
                    emsim[*,it]+=dens_tmp*csden_tmp[q[2],*]*rates[2,*]
                 END
                 2 : BEGIN
                    calc_ar_line_rates,'z',temp_tmp/1.0e3,q,rates
                    emsim[*,it]=dens_tmp*csden_tmp[q[0],*]*rates[0,*]
                    emsim[*,it]+=dens_tmp*csden_tmp[q[1],*]*rates[1,*]
                    emsim[*,it]+=dens_tmp*csden_tmp[q[2],*]*rates[2,*]
                 END
                 3 : BEGIN	
                    rates=reform_ark_data([3.73105,3.73105],/load) ;hard coded for the lya1 transition in the ArK database 
                    emsim[*,it]=dens_tmp*csden_tmp[16,*]*interpol(rates.ion,alog10(rates.temp),alog10(temp/1.0e3))
                    emsim[*,it]+=dens_tmp*csden_tmp[17,*]*interpol(rates.exc,alog10(rates.temp),alog10(temp/1.0e3))
                    emsim[*,it]+=dens_tmp*csden_tmp[18,*]*interpol(rates.rec,alog10(rates.temp),alog10(temp/1.0e3))
                 END
              ENDCASE
              eserr=emsim*0.0   ;no error propogation for now
           endfor
           print, 'Simulated Emissivity Profile Computed: line = '+num2str(line[j],1)
           
           
           IF keyword_set(nobr) THEN BEGIN
              nch=n(pos[0,*])+1
              brsim=fltarr(nch,nt)+1.0
              bserr=fltarr(nch,nt)
              print, 'Simulated Brightness Profile Skipped: line = '+num2str(line[j],1)
           ENDIF ELSE BEGIN
              print,'Calculating Simulated Brightness Profile......: line = '+num2str(line[j],1)
              print,'This may take a while... Be patient...'
              xpos=pos
              genpos_pos_reform,xpos,[0.44,1.0,-0.6,0.6]	
              brsim=genpos_line_br(xpos,emsim,rhovec,time,shot,time,/rho,brerr=bserr)
              print, 'Simulated Brightness Profile Computed: line = '+num2str(line[j],1)
           ENDELSE
           isd={emiss:emsim,emerr:eserr,rho:rhovec,time:time,bright:brsim,brerr:bserr,brRho:brRho,pos:pos,shot:shot,t1:t1,t2:t2,enorm:0.0,bnorm:0.0}
           *sptr[j]=isd
        ENDFOR
        result={exp:eptr,sim:sptr}
        
        return,result
     END
