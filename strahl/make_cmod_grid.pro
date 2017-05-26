;+
;NAME:
;	MAKE_CMOD_GRID
;	
;PURPOSE:
;	This procedure constructs the grid input file for STRAHL and
;	produces graphical output.  This works for steady-state equilibria
;
;CALLING SEQUENCE:
;	MAKE_CMOD_GRID,shot,grid_time
;	
;INPUTS:
;	shot		LONG	shot number
;	grid_time	FLOAT	time point
;
;KEYWORD PARAMETERS
;	plot	/plot will graphically output equilibrium data
;	ps	/ps will force graphical output to the PS device (if /plot)
;
;OUTPUTS:
;	ASCII file (/home/username/strahl/cmod/grid_shotnumber.0)
;	XWIN/PS of the EFIT equilibrium, the interpolated equilibrium
;		and the neoclassical transport parameters
;
;OPTIONAL OUTPUTS:
;	filepath	STRING gives the path to the ASCII file
;
;MODIFICATION HISTORY:
;	Written by:	T. Putterich 5/2010
;	3/5/12		M.L. Reinke - modified the default ASCII file path
;	1/14/12		M.L. Reinke - modified path to match Mercurial strahl version
;
;-

pro make_cmod_grid,shot,grid_time,ps=ps,filepath=filepath,plot=plot


;grid properties - this is only about the output! Take care that the input EFIT data is sufficiently detailed!
n_theta = 301           ; must be odd

n_rho_sol = 21l
n_rho_ped= 45l
n_rho_conf= 55l         ; must be odd

n_lambda = 101l         ;must be odd, nlambda is used for the integral during evaluating the fraction of circulating particles


;set-up colors and ps-properties
IF keyword_set(plot) THEN BEGIN
	if keyword_set(ps) then begin
	    set_plot,'ps'
	    psportrait
	    set_color_right
	    !P.font=0
	endif else begin
	    set_color_right
	    cw
	    window,0,xsize=500,ysize=900
	    !P.font=-1
	endelse
	plot,[0.4,1.0],[-0.45,0.45],/nodata,xstyle=4,ystyle=4,/isotropic
	plot_tiles
ENDIF



;get contours and |B|
efit_psicont, xc,yc,nc,tc, shot = shot
bmod = bmod_grid(grid_time,shot=shot,tree=tree)
;GET CORRESPONDING R AND Z
mdsopen,'analysis',shot
    psirz = mdsvalue('_psi=\efit_psirz')
    rgrid = mdsvalue('dim_of(_psi,0)')
    zgrid = mdsvalue('dim_of(_psi,1)')
mdsclose,'analysis',shot



;get volume and v_loop

mdsopen,'analysis',shot

volume_lcfs=mdsvalue('_volume=\analysis::top.efit.results.a_eqdsk:vout')/1e6
time_volume_lcfs=mdsvalue('dim_of(_volume,0)')

vloop=mdsvalue('_vloop=\analysis::top.efit.results.a_eqdsk:vloopt')
time_vloop=mdsvalue('dim_of(_vloop,0)')

q=mdsvalue('_q=\analysis::efit_fitout:qpsi')
time_q=mdsvalue('dim_of(_q,0)')
rho_q=sqrt(mdsvalue('dim_of(_q,1)'))

mdsclose,'analysis',shot

v_loop = interpol(vloop,time_vloop,grid_time)
total_volume = interpol(volume_lcfs,time_volume_lcfs,grid_time)



;find the right time_index
t_index = where(abs(abs(grid_time-tc)-min(abs(grid_time-tc))) le 0.0001)
t_index=t_index[0]



;plot the flux surface as given by EFIT


; magnetic axis:
r_axis = xc(0,0,t_index)
z_axis = yc(0,0,t_index)

IF keyword_set(plot) THEN BEGIN
	oplot,[r_axis],[z_axis],/psym,thick=3

	; all internal flux surfaces
	for j = 1., n_elements(xc(0,*,t_index))-2 do begin
	    usefull_inds = where(xc(*,j,t_index) ne 0 or yc(*,j,t_index) ne 0)
	    oplot,xc(usefull_inds,j,t_index),yc(usefull_inds,j,t_index),thick=1
	endfor

	;LCFS
	j=n_elements(xc(0,*,t_index))-1
	usefull_inds = where(xc(*,j,t_index) ne 0 or yc(*,j,t_index) ne 0)
	oplot,xc(usefull_inds,j,t_index),yc(usefull_inds,j,t_index),thick=3

	;determine the closest structure at LFS - coordinates of limiters were read of a plot with cursor
	;lim_rhos=fltarr(5)

	oplot,[0.83766, 0.874, 0.9050, 0.874, 0.83766],[0.192,0.13458,0.0, -0.13562, -0.192] , psym = 1,symsize=2 ,thick =2
	oplot,[0.44],[z_axis] , psym = 1 , thick =2,symsize=2 
ENDIF
lim_rhos=efit_rz2rho([0.83766, 0.874, 0.9050, 0.874, 0.83766],  $
                    [0.192,0.13458,0.0, -0.13562, -0.192],  $
                    grid_time,shot=shot,/PSINORM,/sqrt)

lim_rhos = min(lim_rhos)

inner_lim_rho=efit_rz2rho(0.44,z_axis,grid_time,shot=shot,/PSINORM,/sqrt)

if inner_lim_rho lt lim_rhos then begin
    print,'This plasma is limited by the inner limiter!'
endif
print,'rho(inner_limiter) =' ,inner_lim_rho
print,'rho(outer_limiter) =' ,lim_rhos

rho_pols_SOL = (lim_rhos -1.0) * findgen(n_rho_sol)/(n_rho_sol-1)+1.0


; now do the theta angle evaluation - needs magnetic axis
theta = fltarr(n_elements(xc(*,0,t_index)),n_elements(xc(0,*,t_index)))
theta(*,*) = atan((yc(*,*,t_index)-z_axis),(xc(*,*,t_index)-r_axis))
index_phase=where(yc(*,*,t_index)-z_axis le 0.0)
theta(index_phase) = theta(index_phase)+2.0*!Pi


; interpolate R and z to the new theta grid
; later R and z are rescaled using the theta grid

theta_new=findgen(n_theta)/(n_theta-1)*2*!Pi
r_new=fltarr(n_theta,n_elements(nc(*,0)))
z_new=fltarr(n_theta,n_elements(nc(*,0)))
bmod_new=fltarr(n_theta,n_elements(nc(*,0)))
r_new_tmp=fltarr(n_theta)
z_new_tmp=fltarr(n_theta)
bmod_new_tmp=fltarr(n_theta)

for i_flux = 0, n_elements(nc(*,0))-1 do begin

        if nc(i_flux,t_index) ge 2 then begin

            for i_theta = 0, n_theta-1 do begin
            
                r_new_tmp(i_theta) = interpol([xc(0:nc(i_flux,t_index)-2,i_flux,t_index),xc(0,i_flux,t_index)]-xc(0,0,t_index), $
                                     [0.,theta(1:nc(i_flux,t_index)-2,i_flux),( 2.0* !Pi -1e-4) ],theta_new(i_theta),/lsquadratic)+xc(0,0,t_index)
                z_new_tmp(i_theta) = interpol([yc(0:nc(i_flux,t_index)-2,i_flux,t_index),yc(0,i_flux,t_index)]-yc(0,0,t_index), $
                                     [0.,theta(1:nc(i_flux,t_index)-2,i_flux),(2.0 * !Pi -1e-4) ],theta_new(i_theta),/lsquadratic)+yc(0,0,t_index)
                
                bmod_new_tmp(i_theta) = bilin_interpol(bmod,rgrid,zgrid,r_new_tmp(i_theta),z_new_tmp(i_theta))

            endfor


    ; in the center there are too few points, which leads to discontinuities, which leads to problems for the neoclassical calculations....
    ; workaround: smooth in the center! not good!
;stop
                if i_flux lt n_elements(nc(*,0))/3.0 then begin
                    r_new(*,i_flux) = smooth(smooth(r_new_tmp(*),n_theta/25.),n_theta/10.)
                    z_new(*,i_flux) = smooth(smooth(z_new_tmp(*),n_theta/25.),n_theta/10.)
                    bmod_new(*,i_flux) = smooth(smooth(bmod_new_tmp(*),n_theta/25.,/edge),n_theta/10.,/edge)
                endif else begin
                    r_new(*,i_flux) = r_new_tmp(*)
                    z_new(*,i_flux) = z_new_tmp(*)
                    bmod_new(*,i_flux) = bmod_new_tmp(*)
                endelse        


        endif else begin
        
            r_new(*,i_flux)=theta_new(*)*0.0+r_axis
            z_new(*,i_flux)=theta_new(*)*0.0+z_axis
            bmod_new(*,i_flux)=theta_new(*)*0.0+bilin_interpol(bmod,rgrid,zgrid,r_axis,z_axis)
         
        endelse

endfor




; evaluate the scaling for flux surfaces in the SOL


r_cut= findgen(5000)/4999. * 0.5 + 0.42
r_inner = fltarr(n_rho_sol)
r_outer = fltarr(n_rho_sol)
r_mid = fltarr(n_rho_sol)
scale_factor = fltarr(n_rho_sol)

rho_cut = efit_rz2rho(r_cut,  replicate(z_axis,n_elements(r_cut)), grid_time,shot=shot,/PSINORM,/sqrt)

    inds_smaller=where(r_cut lt r_axis)
    inds_larger=where(r_cut gt r_axis)

for i_rho_sol = 0.,n_rho_sol-1 do begin
    ind_inner = where(abs(abs(rho_cut(inds_smaller)-rho_pols_sol(i_rho_sol))-min(abs(rho_cut(inds_smaller)-rho_pols_sol(i_rho_sol)))) lt 1e-4)
    ind_outer = where(abs(abs(rho_cut(inds_larger)-rho_pols_sol(i_rho_sol))-min(abs(rho_cut(inds_larger)-rho_pols_sol(i_rho_sol)))) lt 1e-4)
    r_inner_tmp = r_cut(inds_smaller(ind_inner(0)))
    r_outer_tmp = r_cut(inds_larger(ind_outer(0)))
    r_inner(i_rho_sol) = r_inner_tmp
    r_outer(i_rho_sol) = r_outer_tmp

    r_mid(i_rho_sol) = 0.5 * (r_inner_tmp+r_outer_tmp)
    scale_factor(i_rho_sol) = (r_outer(i_rho_sol)-r_mid(i_rho_sol))/(r_outer(0)-r_mid(0))
endfor 


; make final rho_pol-grid
; rescale R and z for each theta
;inside LCFS: rescale on the above obtained R and Z grid
;outside LCFS: use the obtained scaling factors

rho_final = fltarr(n_rho_conf+n_rho_sol)
rho_final(0:n_rho_conf-1-n_rho_ped)=(findgen(n_rho_conf-n_rho_ped)/(n_rho_conf-n_rho_ped-1))*0.95
rho_final(n_rho_conf-n_rho_ped:n_rho_conf-1)=(findgen(n_rho_ped)/(n_rho_ped))*(0.05)+0.95+1d-3
rho_final(n_rho_conf:n_rho_conf+n_rho_sol-1)=rho_pols_sol(0:n_rho_sol-1)
rho_pol=rho_final

; r and z for the final grid
r_final=fltarr(n_theta,n_rho_conf+n_rho_sol)
z_final=fltarr(n_theta,n_rho_conf+n_rho_sol)
bmod_final=fltarr(n_theta,n_rho_conf)



for i_theta = 0, n_theta-1 do begin

    rho_tmp = efit_rz2rho(r_new(i_theta,*),z_new(i_theta,*),grid_time,shot=shot,/PSINORM,/sqrt)


    r_final(i_theta,0:n_rho_conf)= interpol(r_new(i_theta,*), rho_tmp , rho_final(0:n_rho_conf))
    z_final(i_theta,0:n_rho_conf)= interpol(z_new(i_theta,*), rho_tmp , rho_final(0:n_rho_conf))
    bmod_final(i_theta,0:n_rho_conf-1) = abs(interpol(bmod_new(i_theta,*), rho_tmp , rho_final(0:n_rho_conf-1)))

    r_final(i_theta,n_rho_conf:n_rho_conf+n_rho_sol-1)=reform((r_final(i_theta,n_rho_conf)-r_axis)*scale_factor(0:n_rho_sol-1),1,n_rho_sol)+r_axis
    z_final(i_theta,n_rho_conf:n_rho_conf+n_rho_sol-1)=reform((z_final(i_theta,n_rho_conf)-z_axis)*scale_factor(0:n_rho_sol-1),1,n_rho_sol)+z_axis


endfor
                    

;endfor


IF keyword_set(plot) THEN BEGIN
	plot,[0.4,1.0],[-0.45,0.45],/nodata,xstyle=4,ystyle=4,/isotropic
	plot_tiles

	for i = 0.,n_rho_conf-1 do oplot,r_final(*,i),z_final(*,i),color=!mycol.black
	oplot,r_final(*,n_rho_conf),z_final(*,n_rho_conf),color=!mycol.black,thick=2
	for i = n_rho_conf+1,n_rho_conf+n_rho_sol-1 do oplot,r_final(*,i),z_final(*,i),color=!mycol.blue

	oplot,[0.83766, 0.874, 0.9050, 0.874, 0.83766],[0.192,0.13458,0.0, -0.13562, -0.192] , psym = 1,symsize=2 ,thick =2
	oplot,[0.44],[z_axis] , psym = 1 , thick =2,symsize=2 
ENDIF


;get rho_vol

rho_vol = fltarr(n_rho_conf+n_rho_sol)


;first inside LCFS (easy!)

rho_vol(0:n_rho_conf-1)=sqrt(efit_rz2rho(r_final(0,0:n_rho_conf-1),z_final(0,0:n_rho_conf-1),grid_time,shot=shot,/VOLNORM)/2.0/!Pi^2/r_axis*total_volume)

; calculate rho_vol outside of the the SOL  

total_volume_tmp = total_volume

for i = 0., n_rho_sol-1 do begin 

    delta_vol = 0 
    for i_theta = 1., n_theta-1 do begin


         vector_along_fs = [r_final(i_theta-1,n_rho_conf+i)-r_final(i_theta,n_rho_conf+i),z_final(i_theta-1,n_rho_conf+i)-z_final(i_theta,n_rho_conf+i)]
         vector_to_point = [r_final(i_theta,n_rho_conf+i)-r_axis,z_final(i_theta,n_rho_conf+i)-z_axis]
         delta_vector_to_point = [r_final(i_theta,n_rho_conf+i-1)-r_final(i_theta,n_rho_conf+i),z_final(i_theta,n_rho_conf+i-1)-z_final(i_theta,n_rho_conf+i)]

         angle = acos(transpose(vector_along_fs)#vector_to_point/(sqrt(transpose(vector_along_fs)#vector_along_fs)*sqrt(transpose(vector_to_point)#vector_to_point)))

         area = sqrt(transpose(vector_along_fs)#vector_along_fs) * sqrt(transpose(delta_vector_to_point)#delta_vector_to_point) * sin(angle)


         average_large_radius_of_volume_element = ( r_final(i_theta,n_rho_conf+i)+r_final(i_theta,n_rho_conf+i-1) ) *0.5

         delta_vol=delta_vol+ 2*!Pi * average_large_radius_of_volume_element * area
         if i_theta eq 1 then delta_vol = 2.0 * delta_vol        

    endfor
    total_volume_tmp =  total_volume_tmp + delta_vol 
    rho_vol(n_rho_conf+i) = sqrt(total_volume_tmp / 2.0/!Pi^2/r_axis)

endfor

; get r_outer and r_inner for final rho_pol grid
r_in = fltarr(n_rho_sol+n_rho_conf)
r_out = fltarr(n_rho_sol+n_rho_conf)

r_in(0)  = r_axis
r_out(0) = r_axis

for i_rho = 1.,n_rho_conf+n_rho_sol-1 do begin

    ind_inner = where(abs(abs(rho_cut(inds_smaller)-rho_pol(i_rho))-min(abs(rho_cut(inds_smaller)-rho_pol(i_rho)))) lt 1e-4)
    ind_outer = where(abs(abs(rho_cut(inds_larger)-rho_pol(i_rho))-min(abs(rho_cut(inds_larger)-rho_pol(i_rho)))) lt 1e-4)
    r_inner_tmp = r_cut(inds_smaller(ind_inner(0)))
    r_outer_tmp = r_cut(inds_larger(ind_outer(0)))
    r_in(i_rho) = r_inner_tmp
    r_out(i_rho) = r_outer_tmp

endfor 



; interpolate q-profile at grid_time

q_out = bilin_interpol(q,time_q,rho_q,replicate(grid_time,n_rho_conf),rho_pol(0:n_rho_conf-1))

;-----------------------------------------------------------------------
; now flux surface averages.... only in core possible!
;-----------------------------------------------------------------------

;get magnetic field for all gridpoints

bz=fltarr(n_theta,n_rho_conf)
br=fltarr(n_theta,n_rho_conf)
bpol=fltarr(n_theta,n_rho_conf)

for i_rho = 0.0, n_rho_conf-1 do begin

    psi=efit_rz2psi(r_final(*,i_rho),z_final(*,i_rho),grid_time,bz_tmp,br_tmp,debug=debug,shot=shot,tree=tree)
    bz(*,i_rho) = bz_tmp(*)
    br(*,i_rho) = br_tmp(*)
    bpol(*,i_rho) = sqrt(br_tmp(*)^2+bz_tmp(*)^2)

endfor


; first normalization integral  < 1/ B_poloidal > 
norm_int = dblarr(n_rho_conf)
length_vector_along_fs_array =fltarr(n_rho_conf,n_theta) 
length_around=fltarr(n_rho_conf) 


for i_rho = 1.0, n_rho_conf-1 do begin

    for i_theta = 1., n_theta-1 do begin

         vector_along_fs = double([r_final(i_theta-1,i_rho)-r_final(i_theta,i_rho),z_final(i_theta-1,i_rho)-z_final(i_theta,i_rho)])
         length_vector_along_fs_array(i_rho,i_theta)=sqrt(transpose(vector_along_fs)#vector_along_fs)

       ;  print,length_vector_along_fs_array(i_rho,i_theta),double([r_final(i_theta-1,i_rho)-r_final(i_theta,i_rho),z_final(i_theta-1,i_rho)-z_final(i_theta,i_rho)])

    endfor

    length_vector_along_fs_array(i_rho,*) = median(reform(length_vector_along_fs_array(i_rho,*)),5)

    
endfor


for i_rho = 1.0, n_rho_conf-1 do begin

    for i_theta = 1., n_theta-1 do begin


         norm_int(i_rho) = norm_int(i_rho) + length_vector_along_fs_array(i_rho,i_theta) / bpol(i_theta,i_rho)
         length_around(i_rho) = length_around(i_rho) + length_vector_along_fs_array(i_rho,i_theta)
         if i_theta eq 1 then norm_int(i_rho) = 2.0 * norm_int(i_rho)  

    endfor
    
endfor
norm_int(0) = norm_int(1)
length_around(0) = 0.0








; now <B_total>

b_total = dblarr(n_rho_conf)
for i_rho = 1.0, n_rho_conf-1 do begin

    for i_theta = 1., n_theta-1 do begin

         b_total(i_rho) = b_total(i_rho) + length_vector_along_fs_array(i_rho,i_theta) / bpol(i_theta,i_rho) * bmod_final(i_theta,i_rho)
         if i_theta eq 1 then b_total(i_rho) = 2.0 * b_total(i_rho) 

    endfor
         b_total(i_rho) = b_total(i_rho)/norm_int(i_rho) 
    
endfor
b_total(0) = bmod_final(0,0)


; now <B_total**2>

b_totalsq = dblarr(n_rho_conf)
for i_rho = 1.0, n_rho_conf-1 do begin

    for i_theta = 1., n_theta-1 do begin

         b_totalsq(i_rho) = b_totalsq(i_rho) + length_vector_along_fs_array(i_rho,i_theta) / bpol(i_theta,i_rho) * bmod_final(i_theta,i_rho)^2
         if i_theta eq 1 then b_totalsq(i_rho) = 2.0 * b_totalsq(i_rho) 

    endfor
         b_totalsq(i_rho) = b_totalsq(i_rho)/norm_int(i_rho) 
    
endfor
b_totalsq(0) = bmod_final(0,0)^2


; now <1/B_total**2>

inv_b_totalsq = dblarr(n_rho_conf)
for i_rho = 1.0, n_rho_conf-1 do begin

    for i_theta = 1., n_theta-1 do begin

         inv_b_totalsq(i_rho) = inv_b_totalsq(i_rho) + length_vector_along_fs_array(i_rho,i_theta) / bpol(i_theta,i_rho) / bmod_final(i_theta,i_rho)^2
         if i_theta eq 1 then inv_b_totalsq(i_rho) = 2.0 * inv_b_totalsq(i_rho) 

    endfor
         inv_b_totalsq(i_rho) = inv_b_totalsq(i_rho)/norm_int(i_rho) 
    
endfor
inv_b_totalsq(0) = 1./bmod_final(0,0)^2



; now <R B_tor>

r_bt_total = dblarr(n_rho_conf)
for i_rho = 1.0, n_rho_conf-1 do begin

    for i_theta = 1., n_theta-1 do begin

         r_bt_total(i_rho) = r_bt_total(i_rho) + length_vector_along_fs_array(i_rho,i_theta) / bpol(i_theta,i_rho)  $
                             *r_final(i_theta,i_rho) * sqrt(bmod_final(i_theta,i_rho)^2-bpol(i_theta,i_rho)^2)
                             
         if i_theta eq 1 then r_bt_total(i_rho) = 2.0 * r_bt_total(i_rho) 

    endfor
         r_bt_total(i_rho) = r_bt_total(i_rho) / norm_int(i_rho) 
         
endfor
r_bt_total(0) = sqrt(bmod_final(0,0)^2-bpol(0,0)^2) * r_axis



; now <R^2 B_pol^2 / B^2 >

r2_bpol2_inv_b2_total = dblarr(n_rho_conf)
for i_rho = 1.0, n_rho_conf-1 do begin

    for i_theta = 1., n_theta-1 do begin

         r2_bpol2_inv_b2_total(i_rho) = r2_bpol2_inv_b2_total(i_rho) + length_vector_along_fs_array(i_rho,i_theta) / bpol(i_theta,i_rho)  $
                             *r_final(i_theta,i_rho)^2 * bpol(i_theta,i_rho)^2 / bmod_final(i_theta,i_rho)^2
                             
         if i_theta eq 1 then r2_bpol2_inv_b2_total(i_rho) = 2.0 * r2_bpol2_inv_b2_total(i_rho) 

    endfor
         r2_bpol2_inv_b2_total(i_rho) = r2_bpol2_inv_b2_total(i_rho) / norm_int(i_rho) 
         
endfor
r2_bpol2_inv_b2_total(0) = 0.0



; now < 1 / R^2 >

inv_r2_total = dblarr(n_rho_conf)
for i_rho = 1.0, n_rho_conf-1 do begin

    for i_theta = 1., n_theta-1 do begin

         inv_r2_total(i_rho) = inv_r2_total(i_rho) + length_vector_along_fs_array(i_rho,i_theta) / bpol(i_theta,i_rho)  $
                             / r_final(i_theta,i_rho)^2 
                             
         if i_theta eq 1 then inv_r2_total(i_rho) = 2.0 * inv_r2_total(i_rho) 

    endfor
         inv_r2_total(i_rho) = inv_r2_total(i_rho) / norm_int(i_rho) 
         
endfor
inv_r2_total(0) = 1/r_axis^2



; now the fourier coefficients for NEOART........

; fourier NEOART, normalized theta

norm_theta = dblarr(n_theta,n_rho_conf)
norm_theta_normalization = dblarr(n_rho_conf)

for i_rho = 1.0, n_rho_conf-1 do begin

    for i_theta = 1., n_theta-1 do begin

        
        for j_theta = 1., i_theta do begin

         norm_theta(i_theta, i_rho) = norm_theta(i_theta,i_rho) + length_vector_along_fs_array(i_rho,j_theta) / bpol(j_theta,i_rho) * bmod_final(j_theta,i_rho)
                             

        endfor
         if i_theta eq 1 then norm_theta(i_theta, i_rho) = 2.0 * norm_theta(i_theta, i_rho)

    endfor
    norm_theta_normalization(i_rho) = norm_theta(n_theta-1, i_rho)
    norm_theta(*, i_rho) = 2 * !Pi * norm_theta(*, i_rho) / norm_theta_normalization(i_rho)
         
endfor
norm_theta(*, 0) = 0.0

;fourier coefficient < cos( m theta ) B^2 >
n_m = 3
fcoefficients_cos_b2 = dblarr(n_rho_conf,n_m)


for i_m = 0., n_m -1 do begin

    for i_rho = 1.0, n_rho_conf-1 do begin

        for i_theta = 1., n_theta-1 do begin

             fcoefficients_cos_b2(i_rho,i_m) = fcoefficients_cos_b2(i_rho,i_m) + length_vector_along_fs_array(i_rho,i_theta) / bpol(i_theta,i_rho)  $
                                 * cos((i_m+1) * norm_theta(i_theta,i_rho) ) * bmod_final(i_theta,i_rho)^2  

             if i_theta eq 1 then fcoefficients_cos_b2(i_rho,i_m) = 2.0 * fcoefficients_cos_b2(i_rho,i_m) 

        endfor
             fcoefficients_cos_b2(i_rho,i_m) = fcoefficients_cos_b2(i_rho,i_m) / norm_int(i_rho) 

    endfor
    fcoefficients_cos_b2(0,i_m) = 0

endfor



;fourier coefficient < sin( m theta ) B^2 >
n_m = 3
fcoefficients_sin_b2 = dblarr(n_rho_conf,n_m)


for i_m = 0., n_m -1 do begin

    for i_rho = 1.0, n_rho_conf-1 do begin

        for i_theta = 1., n_theta-1 do begin

             fcoefficients_sin_b2(i_rho,i_m) = fcoefficients_sin_b2(i_rho,i_m) + length_vector_along_fs_array(i_rho,i_theta) / bpol(i_theta,i_rho)  $
                                 * sin((i_m+1) * norm_theta(i_theta,i_rho) ) * bmod_final(i_theta,i_rho)^2  

             if i_theta eq 1 then fcoefficients_sin_b2(i_rho,i_m) = 2.0 * fcoefficients_sin_b2(i_rho,i_m) 

        endfor
             fcoefficients_sin_b2(i_rho,i_m) = fcoefficients_sin_b2(i_rho,i_m) / norm_int(i_rho) 

    endfor
    fcoefficients_sin_b2(0,i_m) = 0

endfor



;fourier coefficient < cos( m theta ) B lnB >
n_m = 3
fcoefficients_cos_blnb = dblarr(n_rho_conf,n_m)


for i_m = 0., n_m -1 do begin

    for i_rho = 1.0, n_rho_conf-1 do begin

        for i_theta = 1., n_theta-1 do begin

             fcoefficients_cos_blnb(i_rho,i_m) = fcoefficients_cos_blnb(i_rho,i_m) + length_vector_along_fs_array(i_rho,i_theta) / bpol(i_theta,i_rho)  $
                                 * cos((i_m+1) * norm_theta(i_theta,i_rho) ) * bmod_final(i_theta,i_rho) * alog(bmod_final(i_theta,i_rho)) 

             if i_theta eq 1 then fcoefficients_cos_blnb(i_rho,i_m) = 2.0 * fcoefficients_cos_blnb(i_rho,i_m) 

        endfor
             fcoefficients_cos_blnb(i_rho,i_m) = fcoefficients_cos_blnb(i_rho,i_m) / norm_int(i_rho) 

    endfor
    fcoefficients_cos_blnb(0,i_m) = 0

endfor

;fourier coefficient < sin( m theta ) B lnB >
n_m = 3
fcoefficients_sin_blnb = dblarr(n_rho_conf,n_m)


for i_m = 0., n_m -1 do begin

    for i_rho = 1.0, n_rho_conf-1 do begin

        for i_theta = 1., n_theta-1 do begin

             fcoefficients_sin_blnb(i_rho,i_m) = fcoefficients_sin_blnb(i_rho,i_m) + length_vector_along_fs_array(i_rho,i_theta) / bpol(i_theta,i_rho)  $
                                 * sin((i_m+1) * norm_theta(i_theta,i_rho) ) * bmod_final(i_theta,i_rho) * alog(bmod_final(i_theta,i_rho)) 

             if i_theta eq 1 then fcoefficients_sin_blnb(i_rho,i_m) = 2.0 * fcoefficients_sin_blnb(i_rho,i_m) 

        endfor
             fcoefficients_sin_blnb(i_rho,i_m) = fcoefficients_sin_blnb(i_rho,i_m) / norm_int(i_rho) 

    endfor
    fcoefficients_sin_blnb(0,i_m) = 0

endfor

; calculate <grad rho_vol> = 
;c
;c    < grad(rho_vol) > = <R B_p> dPsi/drho_vol
;c                      = <R B_p> * Int(dl_p/B_p) / (2 pi R0 rho_vol)      

;first  <(R B_p)>

r_bp_total = dblarr(n_rho_conf)
for i_rho = 1.0, n_rho_conf-1 do begin

    for i_theta = 1., n_theta-1 do begin

         r_bp_total(i_rho) = r_bp_total(i_rho) + length_vector_along_fs_array(i_rho,i_theta) / bpol(i_theta,i_rho)  $
                             * r_final(i_theta,i_rho) * bpol(i_theta,i_rho)
                             
         if i_theta eq 1 then r_bp_total(i_rho) = 2.0 * r_bp_total(i_rho) 

    endfor
         r_bp_total(i_rho) = r_bp_total(i_rho) / norm_int(i_rho) 
         
endfor
r_bp_total(0) = 0.0


;now <grad rho_vol>
grad_rho_vol = r_bp_total * norm_int / ( 2.0 * !Pi * r_axis * rho_vol(0:n_rho_conf-1))
grad_rho_vol(0) = interpol(grad_rho_vol(1:n_rho_conf-1),rho_pol(1:n_rho_conf-1),0.0)

;now <grad rho_vol**2>
;c
;c    < grad(rho_vol)**2 > 
;c       = <(R B_p)**2> (dPsi/drho_vol)**2
;c       = <(R B_p)**2> * (Int(dl_p/B_p) / (2 pi R0 rho_vol) )**2 

;first  <(R B_p)**2>

r_bp2_total = dblarr(n_rho_conf)
for i_rho = 1.0, n_rho_conf-1 do begin

    for i_theta = 1., n_theta-1 do begin

         r_bp2_total(i_rho) = r_bp2_total(i_rho) + length_vector_along_fs_array(i_rho,i_theta) / bpol(i_theta,i_rho)  $
                             *(r_final(i_theta,i_rho) * bpol(i_theta,i_rho))^2
                             
         if i_theta eq 1 then r_bp2_total(i_rho) = 2.0 * r_bp2_total(i_rho) 

    endfor
         r_bp2_total(i_rho) = r_bp2_total(i_rho) / norm_int(i_rho) 
         
endfor
r_bp2_total(0) = 0.0

grad_rho_vol2 = r_bp2_total * norm_int / ( 2.0 * !Pi * r_axis * rho_vol(0:n_rho_conf-1))^2
grad_rho_vol2(0) = interpol(grad_rho_vol2(1:n_rho_conf-1),rho_pol(1:n_rho_conf-1),0.0)


;calculate fraction of circulating particles
;
;  Int_0^1/Bmax (lambda dlambda / < sqrt(1 - lambda B) >) 
;
frac_circ = dblarr(n_rho_conf)
integrand=dblarr(n_lambda)
dum1 = dblarr(n_theta)
lambda = dblarr(n_lambda)
denomin = dblarr(n_theta)
denom_integral=dblarr(n_lambda,n_rho_conf)
integral_all=dblarr(n_rho_conf)

frac_circ(0) = 1.
integrand(0) = 0.

  for i_rho = 1 , n_rho_conf-1 do begin

     b_max = max(bmod_final(*,i_rho))  ;! find maximum field

     dlambda = 1. / b_max / (n_lambda-1.)

     for k = 1,N_LAMBDA-1 do begin

        lambda(k) =dlambda * (k-1)
;c
;c     < sqrt(1 - lambda B) > 
;c            
        for j_theta = 1,N_theta-1 do begin

           denomin(j_theta) = sqrt(max([0.,1.-bmod_final(j_theta,i_rho)*lambda(k)]))
           denom_integral(k,i_rho) = denom_integral(k,i_rho)+ denomin(j_theta) * length_vector_along_fs_array(i_rho,j_theta) / bpol(j_theta,i_rho)
           if i_theta eq 1 then denom_integral(k,i_rho) = 2.0 * denom_integral(k,i_rho) 

        endfor
        denom_integral(k,i_rho)=denom_integral(k,i_rho)/norm_int(i_rho)
        integrand(k) = lambda(k) / denom_integral(k,i_rho)
        integral_all(i_rho) = integral_all(i_rho) + integrand(k) * dlambda

     endfor

     frac_circ(i_rho) = 3./4. * b_totalsq(i_rho) * integral_all(i_rho)

  endfor



; output data into file/ on screen
IF keyword_set(plot) THEN !P.multi=[0.,3.,6.]
filepath=strcompress('/home/'+logname()+'/strahl/cmod/grid_'+string(shot)+'.0',/remove_all)
openw,lun,filepath,/get_lun


printf,lun,' '
printf,lun,'cv  rho volume(LCFS)[cm]   R_axis[cm]   U_loop[V]    time[s]' 
printf,lun,format='(2F16.1,2F13.4)',rho_vol(n_rho_conf)*100.,r_axis*100.,abs(v_loop),grid_time
printf,lun,' '
printf,lun,' '
printf,lun,'cv  number of grid points  points up to separtrix  fourier coefficients'
printf,lun,format='(1I15,2I22)',n_rho_conf+n_rho_sol,n_rho_conf,n_m
printf,lun,' '
printf,lun,' '







stdf='(7f10.5)'
stde='(7e10.2)'

print,'V_loop',v_loop
print,'R_axis',r_axis
print,'minor radius a',rho_vol(n_rho_conf)

print,'rho_pol'
printf,lun,'cv  rho_pol'
print,format=stdf,rho_pol
printf,lun,format=stdf,rho_pol
printf,lun,' '
printf,lun,' '

print,'rho_vol/rho_vol_LCFS'
printf,lun,'cv  rho_vol/rho_vol_LCFS'
print,format=stdf,rho_vol/rho_vol(n_rho_conf)
printf,lun,format=stdf,rho_vol/rho_vol(n_rho_conf)
IF keyword_set(plot) THEN plot,rho_pol,rho_vol/rho_vol(n_rho_conf),title = 'rho_vol/rho_vol_LCFS',xtitle='rho_pol'
printf,lun,' '
printf,lun,' '

print,'rmid_out/r_axis'
printf,lun,'cv  rmid_out/r_axis'
print,format=stdf,r_out/r_axis
printf,lun,format=stdf,r_out/r_axis
IF keyword_set(plot) THEN plot,rho_pol,r_out/r_axis,title = 'rmid_out/r_axis',xtitle='rho_pol'
printf,lun,' '
printf,lun,' '

print,'rmid_in/r_axis'
printf,lun,'cv  rmid_in/r_axis'
print,format=stdf,r_in/r_axis
printf,lun,format=stdf,r_in/r_axis
IF keyword_set(plot) THEN plot,rho_pol,r_in/r_axis,title = 'rmid_in/r_axis',xtitle='rho_pol'
printf,lun,' '
printf,lun,' '

print,'safety factor'
printf,lun,'cv  safety factor'
print,format=stdf,q_out
printf,lun,format=stdf,q_out
IF keyword_set(plot) THEN plot,rho_pol,q_out,title = 'safety factor q',xtitle='rho_pol',xr=[0.,1.]
printf,lun,' '
printf,lun,' '

print,'fraction of circulating particles'
printf,lun,'cv  fraction of circulating particles'
print,format=stdf,frac_circ
printf,lun,format=stdf,frac_circ
IF keyword_set(plot) THEN plot,rho_pol,frac_circ,title = 'fraction of circulating particles',xtitle='rho_pol',xr=[0.,1.]
printf,lun,' '
printf,lun,' '


print,'Integral(dl_p/Bp) [m/T]'
printf,lun,'cv  Integral(dl_p/Bp) [m/T]'
print,format=stdf,norm_int
printf,lun,format=stdf,norm_int
IF keyword_set(plot) THEN plot,rho_pol,norm_int,title = '<dl_p/Bp>',xtitle='rho_pol',xr=[0.,1.]
printf,lun,' '
printf,lun,' '

print,'<B_total> [T]'
printf,lun,'cv  <B_total> [T]'
print,format=stdf,b_total
printf,lun,format=stdf,b_total
IF keyword_set(plot) THEN plot,rho_pol,b_total,title = '<B_total> [T]',xtitle='rho_pol',xr=[0.,1.]
printf,lun,' '
printf,lun,' '

print,'<B_total**2> [T**2]'
printf,lun,'cv  <B_total**2> [T**2]'
print,format=stdf,b_totalsq
printf,lun,format=stdf,b_totalsq
IF keyword_set(plot) THEN plot,rho_pol,b_totalsq,title = '<B_total**2> [T]',xtitle='rho_pol',xr=[0.,1.]
printf,lun,' '
printf,lun,' '

print,'<1/B_total**2> [1/T**2]'
printf,lun,'cv  <1/B_total**2> [1/T**2]'
print,format=stdf,inv_b_totalsq
printf,lun,format=stdf,inv_b_totalsq
IF keyword_set(plot) THEN plot,rho_pol,inv_b_totalsq,title = '<1/B_total**2> [T]',xtitle='rho_pol',xr=[0.,1.]
printf,lun,' '
printf,lun,' '

print,'<R B_T> [mT]'
printf,lun,'cv  <R B_T> [mT]'
print,format=stdf,r_bt_total
printf,lun,format=stdf,r_bt_total
IF keyword_set(plot) THEN plot,rho_pol,r_bt_total,title = '<R B_tor> [T]',xtitle='rho_pol',xr=[0.,1.]
printf,lun,' '
printf,lun,' '

print,'<R^2 B_pol^2 / B^2> [m^2]'
printf,lun,'cv  <R^2 B_pol^2 / B^2> [m^2]'
print,format=stde,r2_bpol2_inv_b2_total
printf,lun,format=stde,r2_bpol2_inv_b2_total
IF keyword_set(plot) THEN plot,rho_pol,r2_bpol2_inv_b2_total,title = '<R^2 B_pol^2 / B^2> [m^2]',xtitle='rho_pol',xr=[0.,1.]
printf,lun,' '
printf,lun,' '

print,'< 1 / R^2 > [1/m^2]'
printf,lun,'cv  < 1 / R^2 > [1/m^2]'
print,format=stdf,inv_r2_total
printf,lun,format=stdf,inv_r2_total
IF keyword_set(plot) THEN plot,rho_pol,inv_r2_total,title = '< 1 / R^2 > [1/m^2]',xtitle='rho_pol',xr=[0.,1.]
printf,lun,' '
printf,lun,' '


IF keyword_set(plot) THEN plot,[0.,1.],[0.,1.],/nodata,xr=[0.,1.],yr=[min(fcoefficients_cos_b2),max(fcoefficients_cos_b2)],/ystyle,title='< cos(m theta_norm) B^2 > [T^2]',xtitle='rho_pol'
print,'< cos(m theta_norm) B^2 > [T^2]'
printf,lun,'cv  < cos(m theta_norm) B^2 > [T^2]'
for i_m = 0., n_m-1 do begin
    print,format=stde,fcoefficients_cos_b2(*,i_m)
    printf,lun,format=stde,fcoefficients_cos_b2(*,i_m)
    IF keyword_set(plot) THEN oplot,rho_pol,fcoefficients_cos_b2(*,i_m)
endfor
printf,lun,' '
printf,lun,' '

IF keyword_set(plot) THEN plot,[0.,1.],[0.,1.],/nodata,xr=[0.,1.],yr=[min(fcoefficients_sin_b2),max(fcoefficients_sin_b2)],/ystyle,title='< sin(m theta_norm) B^2 > [T^2]',xtitle='rho_pol'
print,'< sin(m theta_norm) B^2 > [T^2]'
printf,lun,'cv  < sin(m theta_norm) B^2 > [T^2]'
for i_m = 0., n_m-1 do begin
    print,format=stde,fcoefficients_sin_b2(*,i_m)
    printf,lun,format=stde,fcoefficients_sin_b2(*,i_m)
    IF keyword_set(plot) THEN oplot,rho_pol,fcoefficients_sin_b2(*,i_m)
endfor
printf,lun,' '
printf,lun,' '


IF keyword_set(plot) THEN plot,[0.,1.],[0.,1.],/nodata,xr=[0.,1.],yr=[min(fcoefficients_cos_blnb),max(fcoefficients_cos_blnb)],/ystyle,title='< cos(m theta_norm) B lnB > [T lnT]',xtitle='rho_pol'
print,'< cos(m theta_norm) B lnB > [T lnT]'
printf,lun,'cv  < cos(m theta_norm) B lnB > [T lnT]'
for i_m = 0., n_m-1 do begin
    print,format=stde,fcoefficients_cos_blnb(*,i_m)
    printf,lun,format=stde,fcoefficients_cos_blnb(*,i_m)
    IF keyword_set(plot) THEN oplot,rho_pol,fcoefficients_cos_blnb(*,i_m)
endfor
printf,lun,' '
printf,lun,' '

IF keyword_set(plot) THEN plot,[0.,1.],[0.,1.],/nodata,xr=[0.,1.],yr=[min(fcoefficients_sin_blnb),max(fcoefficients_sin_blnb)],/ystyle,title='< sin(m theta_norm) B lnB > [T lnT]',xtitle='rho_pol'
print,'< sin(m theta_norm) B lnB > [T lnT]'
printf,lun,'cv  < sin(m theta_norm) B lnB > [T lnT]'
for i_m = 0., n_m-1 do begin
    print,format=stde,fcoefficients_sin_blnb(*,i_m)
    printf,lun,format=stde,fcoefficients_sin_blnb(*,i_m)
    IF keyword_set(plot) THEN oplot,rho_pol,fcoefficients_sin_blnb(*,i_m)
endfor
printf,lun,' '
printf,lun,' '

print,'<grad(rho_vol)> '
printf,lun,'cv  <grad(rho_vol)> '
print,format=stdf,grad_rho_vol
printf,lun,format=stdf,grad_rho_vol
IF keyword_set(plot) THEN plot,rho_pol,grad_rho_vol,title = '<grad(rho_vol)>',xtitle='rho_pol',xr=[0.,1.]
printf,lun,' '
printf,lun,' '

print,'<grad(rho_vol)**2> '
printf,lun,'cv  <grad(rho_vol)**2> '
print,format=stdf,grad_rho_vol2
printf,lun,format=stdf,grad_rho_vol2
IF keyword_set(plot) THEN plot,rho_pol,grad_rho_vol2,title = '<grad(rho_vol)**2>',xtitle='rho_pol',xr=[0.,1.]
printf,lun,' '
printf,lun,' '


close,lun
free_lun,lun

IF keyword_set(plot) THEN !P.multi=[0.,1.,1.]




;set-up colors and ps-properties
IF keyword_set(plot) THEN BEGIN
	if keyword_set(ps) then begin
	    device,/close
	    set_plot,'x'
	    set_color_right
	    !P.font=-1
	endif 
ENDIF

end
