PRO bspline_fs,shot,t,rho=rho,div=div,edge=edge

	common PSI,Psi_Grid,Rho_Grid,B_Grid,Time,R,Z,Psi_Surf,Z_midplane,current_shot,current_option,dPsidRho
	common PSI_ext,RSurf,ZSurf,RXpt,ZXpt,Topology,bs,order,Rknot,Zknot,PSI_BSCoef,RHO_BSCoef
	read_flux,shot,option,error,/bspline
	i=ipt(time,t)
        z_xpt=zxpt[i]
	rho_data=reform(Rho_Grid(i,*,*))
	R_fine0=R(0) & R_fine1=R(n_elements(R)-1)
	R_fine=R_fine0+(R_fine1-R_Fine0)*findgen(101)/100
	Z_fine0=Z(0) & Z_fine1=Z(n_elements(Z)-1)
	Z_fine=Z_fine0+(Z_fine1-Z_Fine0)*findgen(101)/100
	Rho_fine=BS2GD(0,0,R_fine,Z_fine,order,order,Rknot,Zknot,RHO_BSCoef(*,i))
	vessel_plot,shot=shot,tit=num2str(shot,1)+' Time: '+num2str(time[i],dp=2),div=div,edge=edge
        IF NOT keyword_set(rho) THEN rho=[make(min(rho_fine),0,15),-0.002,0.0015,0.0039,0.006]
        IF z_xpt GT 0.0 THEN rho=rho[0:n(rho)-1]
        firstwall_nopump,rw,zw,shot=shot
        FOR i=0,n(rho) DO BEGIN
        	contour,Rho_fine,R_fine,Z_fine,levels=[rho[i]],/overplot,path_xy=path_fs,path_info=info_fs,/path_data
        	rfs=reform(path_fs[0,*])
        	zfs=reform(path_fs[1,*])
                good=inside(rfs,zfs,rw,zw)
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
END

PRO bspline_lcfs,shot,t,div=div,edge=edge,color=color,linestyle=linestyle,oplot=oplot

	IF NOT keyword_set(color) THEN color=[0,30,100,200]
	IF NOT keyword_set(linestyle) THEN linestyle=[0,0,0,0]
        IF n(shot) GT 3 THEN RETURN
	common PSI,Psi_Grid,Rho_Grid,B_Grid,Time,R,Z,Psi_Surf,Z_midplane,current_shot,current_option,dPsidRho
	common PSI_ext,RSurf,ZSurf,RXpt,ZXpt,Topology,bs,order,Rknot,Zknot,PSI_BSCoef,RHO_BSCoef

	IF keyword_set(oplot) THEN vessel_plot,shot=shot[0],tit=' ',div=div,edge=edge

        FOR i=0,n(shot) DO BEGIN
            	firstwall_nopump,rw,zw,shot=shot[i]
        	read_flux,shot[i],option,error,/bspline
                j=ipt(time,t[i])
		rho_data=reform(Rho_Grid(j,*,*))
		R_fine0=R(0) & R_fine1=R(n_elements(R)-1)
		R_fine=R_fine0+(R_fine1-R_Fine0)*findgen(101)/100
                Z_fine0=Z(0) & Z_fine1=Z(n_elements(Z)-1)
                Z_fine=Z_fine0+(Z_fine1-Z_Fine0)*findgen(101)/100
                Rho_fine=BS2GD(0,0,R_fine,Z_fine,order,order,Rknot,Zknot,RHO_BSCoef(*,j))
        	contour,Rho_fine,R_fine,Z_fine,levels=[0.00],/overplot,path_xy=path_fs,path_info=info_fs,/path_data
        	rfs=reform(path_fs[0,*])
        	zfs=reform(path_fs[1,*])
                good=inside(rfs,zfs,rw,zw)
                tmp=where(good EQ 1)
                oplot,rfs[tmp],zfs[tmp],color=color[i],linestyle=linestyle[i]
        ENDFOR  
END
