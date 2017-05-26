;+
; NAME:
;
; AUTHOR:
;   Chi Gao, PSFC/MIT
;   cgao@mit.edu
;
; PURPOSE:
;   This routine is written to create a widget to run strahl
;   simulation with the ability to interactively specify:
;       a. simulation times
;       b. impurity puffing level and period
;       c. diffusivity 2d profile: D(r,t)
;       d. convectivity 2d profile: V(r,t)
;
;   The widget will calculate the charge states of impurity.
;   The widget can calculate the emissivity and brightness (line integral
;   of emissivity) of lines (more emission lines will be added):
;       1. Ar16+  w line
;       2. Ar16+  z line
;
;   The widget can load experimental results:
;       1. Ar16+ w line from HiReX-Sr (load from mdsplus, THACO analysis required)    
;       2. Ar16+ z line from HiReX-Sr (load from mdsplus, THACO analysis required)    
;
; DESCRIPTION:
;   
; INPUT:
;
; OUTPUT:
;         strahl result is stored in ~/strahl/result
;         figures are stored in ~/strahl/widget/result
;  
; HISTORY:
;   08/22/2014: Created by cgao.
;            
;   04/16/2015: use the new code to run strahl
;               /usr/local/cmod/idl/GENIE/GENTRAN/strahl/run_cmod_strahl.pro             

 


;**************************************************************************************
;@genie_ini.bat
;xwplot
@/usr/local/cmod/idl/GENIE/usr/cgao/idl/legend.pro
@/usr/local/cmod/idl/GENIE/usr/cgao/idl/hirexsr_genrad_argon_cgao.pro
@/usr/local/cmod/idl/GENIE/usr/cgao/idl/hirexsr_load_result.pro

;@idl/den2emiss_cgao.pro
;@/home/cgao/GENIE/general_plot.pro
;**************************************************************************************


;@idl/run_cmod_strahl.pro
; As of 11/04/15, Matt has updated run_cmod_strahl to run strahl with
; time-dependence density and temperature profiles.
; It is located in
; /usr/local/cmod/idl/GENIE/GENTRAN/strahl/run_cmod_strahl.pro
; I am updating my code to use this code, instead of my own code in idl/run_cmod_strahl.pro





pro wstrahl_copy,ufrom,uto
  ;keep the wset id
  ufrom.stat.strahl.set_dens_abs = uto.stat.strahl.set_dens_abs

  ufrom.stat.strahl.set_emiss_sim_w_abs = uto.stat.strahl.set_emiss_sim_w_abs
  ufrom.stat.strahl.set_emiss_sim_z_abs = uto.stat.strahl.set_emiss_sim_z_abs
  ufrom.stat.strahl.set_bright_sim_w_abs = uto.stat.strahl.set_bright_sim_w_abs
  ufrom.stat.strahl.set_bright_sim_z_abs = uto.stat.strahl.set_bright_sim_z_abs

  ufrom.stat.hirex.set_emiss_exp_w_abs = uto.stat.hirex.set_emiss_exp_w_abs
  ufrom.stat.hirex.set_emiss_exp_z_abs = uto.stat.hirex.set_emiss_exp_z_abs
  ufrom.stat.hirex.set_bright_exp_w_abs = uto.stat.hirex.set_bright_exp_w_abs
  ufrom.stat.hirex.set_bright_exp_z_abs = uto.stat.hirex.set_bright_exp_z_abs


  ufrom.stat.strahl.set_dens_norm = uto.stat.strahl.set_dens_norm

  ufrom.stat.strahl.set_emiss_sim_w_norm = uto.stat.strahl.set_emiss_sim_w_norm
  ufrom.stat.strahl.set_emiss_sim_z_norm = uto.stat.strahl.set_emiss_sim_z_norm
  ufrom.stat.strahl.set_bright_sim_w_norm = uto.stat.strahl.set_bright_sim_w_norm
  ufrom.stat.strahl.set_bright_sim_z_norm = uto.stat.strahl.set_bright_sim_z_norm

  ufrom.stat.hirex.set_emiss_exp_w_norm = uto.stat.hirex.set_emiss_exp_w_norm
  ufrom.stat.hirex.set_emiss_exp_z_norm = uto.stat.hirex.set_emiss_exp_z_norm
  ufrom.stat.hirex.set_bright_exp_w_norm = uto.stat.hirex.set_bright_exp_w_norm
  ufrom.stat.hirex.set_bright_exp_z_norm = uto.stat.hirex.set_bright_exp_z_norm


  ;make changes 
  uto.stat = ufrom.stat
  uto.shot = ufrom.shot
  uto.Z = ufrom.Z
  uto.Zplot = ufrom.Zplot
  uto.t1 = ufrom.t1
  uto.t2 = ufrom.t2
  uto.Nt = ufrom.Nt
  uto.fz_text = ufrom.fz_text
  *uto.fz = *ufrom.fz
  uto.fztime_text = ufrom.fztime_text
  *uto.fztime = *ufrom.fztime
end



pro wstrahl_load_info,u
  widget_control,u.id.shotid,set_value=num2str(u.shot,1)
  widget_control,u.id.t1id,set_value=num2str(u.t1,dp=2)
  widget_control,u.id.t2id,set_value=num2str(u.t2,dp=2)
  widget_control,u.id.Ntid,set_value=num2str(u.Nt,1)
  widget_control,u.id.Zid,set_value=num2str(u.Z,1)
  widget_control,u.id.Zplotid,set_value=num2str(u.Zplot,1)
  widget_control,u.id.fzid,set_value=u.fz_text
  widget_control,u.id.fztimeid,set_value=u.fztime_text
  widget_control,u.id.fits,set_button=u.stat.profile.fits
  widget_control,u.id.qfits,set_button=u.stat.profile.qfits
  widget_control,u.id.nogrid,set_button=u.stat.optional.nogrid
  widget_control,u.id.nopp,set_button=u.stat.optional.nopp
  widget_control,u.id.noparam,set_button=u.stat.optional.noparam
  widget_control,u.id.nostrahl,set_button=u.stat.optional.nostrahl
  widget_control,u.id.quiet,set_button=u.stat.optional.quiet
  widget_control,u.id.savepath,set_value=u.savepath
  widget_control,u.id.loadpath,set_value=u.loadpath
  widget_control,u.id.strahltmin,set_value=num2str(u.stat.strahl.tr[0],dp=2)
  widget_control,u.id.strahltmax,set_value=num2str(u.stat.strahl.tr[1],dp=2)
  widget_control,u.id.strahlrmin,set_value=num2str(u.stat.strahl.rr[0],dp=2)
  widget_control,u.id.strahlrmax,set_value=num2str(u.stat.strahl.rr[1],dp=2)
  widget_control,u.id.hirextmin,set_value=num2str(u.stat.hirex.tr[0],dp=2)
  widget_control,u.id.hirextmax,set_value=num2str(u.stat.hirex.tr[1],dp=2)
  widget_control,u.id.hirexrmin,set_value=num2str(u.stat.hirex.rr[0],dp=2)
  widget_control,u.id.hirexrmax,set_value=num2str(u.stat.hirex.rr[1],dp=2)


  widget_control,u.id.coef_Nr,set_value=num2str(u.stat.coef.coef_Nr)
  widget_control,u.id.coef_Nt,set_value=num2str(u.stat.coef.coef_Nt)
 
  widget_control,u.id.CHIGmodel,set_droplist_select=u.stat.coef.chiGdroplist
  widget_control,u.id.CHIc0,set_value=num2str(u.stat.coef.chiGparam.c0,dp=2)
  widget_control,u.id.CHIa0,set_value=num2str(u.stat.coef.chiGparam.a0,dp=2)
  widget_control,u.id.CHIa1,set_value=num2str(u.stat.coef.chiGparam.a1,dp=2)
  widget_control,u.id.CHIa2,set_value=num2str(u.stat.coef.chiGparam.a2,dp=2)
  widget_control,u.id.CHIa3,set_value=num2str(u.stat.coef.chiGparam.a3,dp=2)
  widget_control,u.id.CHIx0,set_value=num2str(u.stat.coef.chiGparam.x0,dp=2)

  widget_control,u.id.CHIHmodel,set_droplist_select=u.stat.coef.chiHdroplist
  widget_control,u.id.CHIc,set_value=num2str(u.stat.coef.chiHparam.c,dp=2)
  widget_control,u.id.CHIt0,set_value=num2str(u.stat.coef.chiHparam.t0,dp=2)
  widget_control,u.id.CHIb0,set_value=num2str(u.stat.coef.chiHparam.b0,dp=2)
  widget_control,u.id.CHIb1,set_value=num2str(u.stat.coef.chiHparam.b1,dp=2)
  widget_control,u.id.CHIb2,set_value=num2str(u.stat.coef.chiHparam.b2,dp=2)
  widget_control,u.id.CHIb3,set_value=num2str(u.stat.coef.chiHparam.b3,dp=2)
  
  ;widget_control,u.id.CHI_time,set_value=num2str(u.stat.times.CHI_time,dp=2)

  widget_control,u.id.VGmodel,set_droplist_select=u.stat.coef.VGdroplist
  widget_control,u.id.Vc0,set_value=num2str(u.stat.coef.VGparam.c0,dp=2)
  widget_control,u.id.Va0,set_value=num2str(u.stat.coef.VGparam.a0,dp=2)
  widget_control,u.id.Va1,set_value=num2str(u.stat.coef.VGparam.a1,dp=2)
  widget_control,u.id.Va2,set_value=num2str(u.stat.coef.VGparam.a2,dp=2)
  widget_control,u.id.Va3,set_value=num2str(u.stat.coef.VGparam.a3,dp=2)
  widget_control,u.id.Vx0,set_value=num2str(u.stat.coef.VGparam.x0,dp=2)

  widget_control,u.id.VHmodel,set_droplist_select=u.stat.coef.VHdroplist
  widget_control,u.id.Vc,set_value=num2str(u.stat.coef.VHparam.c,dp=2)
  widget_control,u.id.Vt0,set_value=num2str(u.stat.coef.VHparam.t0,dp=2)
  widget_control,u.id.Vb0,set_value=num2str(u.stat.coef.VHparam.b0,dp=2)
  widget_control,u.id.Vb1,set_value=num2str(u.stat.coef.VHparam.b1,dp=2)
  widget_control,u.id.Vb2,set_value=num2str(u.stat.coef.VHparam.b2,dp=2)
  widget_control,u.id.Vb3,set_value=num2str(u.stat.coef.VHparam.b3,dp=2)

  ;widget_control,u.id.V_time,set_value=num2str(u.stat.times.V_time,dp=2)

  widget_control,u.id.density_time,set_value=num2str(u.stat.times.density_time,dp=2)
  widget_control,u.id.density_time_slider,$
                 set_value= int(2000*(u.stat.times.density_time - u.t1)/(u.t2-u.t1))
  ;widget_control,u.id.emissivity_time,set_value=num2str(u.stat.times.emissivity_time,dp=2)
end







pro wstrahl_plot_binning,u,ps=ps
  if keyword_set(ps) then begin
     set_plot,'ps'
     savepath = '~/strahl/widget/result/'+u.savepath
     device,/encapsulated,/color,filename=savepath+'/binning.eps'
  endif else begin
     widget_control,u.id.draw_binning,get_value=draw_win
     wset,draw_win
  endelse

  t1=u.t1
  t2=u.t2
  Nt=u.nt
  dt=(t2-t1)/Nt
  plot,[t1,t2],[0,0],yrange=[0,nt+1],xtit='time [sec]',charsize=1,xstyle=1
  for in=0,Nt-1 do begin
     oplot,[t1+in*dt,t1+(in+1)*dt],[in,in]+1
     oplot,[t1+in*dt,t1+in*dt],[in,in+1]     
  endfor
  oplot,[t2,t2],[0,nt]
  
  if keyword_set(ps) then begin
     device,/close
     set_plot,'x'
  endif
end



pro wstrahl_plot_source,u,ps=ps
  if keyword_set(ps) then begin
     set_plot,'ps'
     savepath = '~/strahl/widget/result/'+u.savepath
     device,/encapsulated,/color,filename=savepath+'/source.eps'
  endif else begin
     widget_control,u.id.draw_source,get_value=draw_win
     wset,draw_win
  endelse

  t1=u.t1
  t2=u.t2
  fz=*u.fz
  fztime=*u.fztime
  nfz=n(fz)+1
  nfztime=n(fztime)+1
  if (nfz ne nfztime) then begin
     message,'ERROR: fz and fztime lengths are unequal'
  endif else begin
     plot,fztime,fz,xrange=[t1,t2],yrange=[0,max(fz)*1.1],xtit='time [sec]',ytit='fz(1e17)',charsize=1,xstyle=1
  endelse

  ;if (nfz ne nfztime) then begin
  ;   message,'ERROR: fz and fztime lengths are unequal'
  ;endif else if (nfz eq 1) then begin
  ;  plot,[fztime[0],fztime[0]],[0,fz[0]],xrange=[t1,t2],yrange=[0,max(fz)*1.1+100],xtit='time [sec]',charsize=1,xstyle=1
  ;   if (fztime[0]<t2) then oplot,[fztime[0],t2],[fz[0],fz[0]]
  ;endif else begin
  ;   plot,[fztime[0],fztime[0]],[0,fz[0]],xrange=[t1,t2],yrange=[0,max(fz)*1.1+100],xtit='time [sec]',charsize=1,xstyle=1
  ;   for in=1,nfz-2 do begin
  ;      oplot,[fztime[in-1],fztime[in]],[fz[in-1],fz[in-1]]
  ;      oplot,[fztime[in],fztime[in]],[fz[in-1],fz[in]]
  ;   endfor
  ;   oplot,[fztime[nfz-2],fztime[nfz-1]],[fz[nfz-2],fz[nfz-2]]
  ;   if (fztime[nfz-1] lt t2) then begin
  ;      oplot,[fztime[nfz-1],fztime[nfz-1]],[fz[nfz-2],fz[nfz-1]]
  ;      oplot,[fztime[nfz-1],t2],[fz[nfz-1],fz[nfz-1]]        
  ;   endif
  ;endelse

  if keyword_set(ps) then begin
     device,/close
     set_plot,'x'
  endif

end


pro wstrahl_calc_CHI,u
  ;Gindex=widget_info(u.id.CHIGmodel,/droplist_select)
  ;Hindex=widget_info(u.id.CHIHmodel,/droplist_select)
  Gindex=u.stat.coef.chiGdroplist
  Hindex=u.stat.coef.chiHdroplist

  Nr=u.stat.coef.coef_Nr
  Nt=u.stat.coef.coef_Nt
  t1=u.t1
  t2=u.t2
  x=findgen(Nr+1.)/(Nr)
  t=findgen(Nt)/float(Nt)*(t2-t1)+t1
  
  chi=u.stat.coef.chi
  gx=fltarr(Nr)
  ht=fltarr(Nt)
  ;if Gindex ge 0 then gx+=u.stat.coef.chiGparam.a0
  ;if Gindex ge 1 then gx+=u.stat.coef.chiGparam.a1*x
  ;if Gindex ge 2 then gx+=u.stat.coef.chiGparam.a2*x^2
  ;if Gindex ge 3 then gx+=u.stat.coef.chiGparam.a3*x^3
  ;if Gindex ge 4 then gx+=u.stat.coef.chiGparam.a4*x^4
  c0 = u.stat.coef.chiGparam.c0
  a0 = u.stat.coef.chiGparam.a0
  a1 = u.stat.coef.chiGparam.a1
  a2 = u.stat.coef.chiGparam.a2
  a3 = u.stat.coef.chiGparam.a3
  x0 = u.stat.coef.chiGparam.x0

  if Gindex eq 0 then begin
     tmph = fltarr(Nr) + 1.0
     if ((where(x lt x0))[0] ne -1) then tmph(where(x lt x0)) = 0.0
     gx = c0 + (a0+a1*(x-x0) + a2*(x-x0)^2 + a3*(x-x0)^3)*tmph
  endif else if Gindex eq 1 then begin
     print,'gaussian profile'
     gx = a0*exp(-(x-a1)^2/(2*a2^2))+a3+a4*x
  endif else begin
     gx = 0.0
  endelse 

  
  if Hindex ge 0 then ht+=u.stat.coef.chiHparam.b0
  index=where(t ge u.stat.coef.chiHparam.t0,count)
  if count gt 0 then ht[index]+=u.stat.coef.chiHparam.c
  if Hindex ge 1 then ht+=u.stat.coef.chiHparam.b1*(t-t(0))
  if Hindex ge 2 then ht+=u.stat.coef.chiHparam.b2*(t-t(0))^2
  if Hindex ge 3 then ht+=u.stat.coef.chiHparam.b3*(t-t(0))^3
  *u.stat.coef.chi=gx#ht
  *u.stat.coef.coefx=x
  *u.stat.coef.coeft=t
end

pro wstrahl_plot_CHI_spatial,u,ps=ps
  if keyword_set(ps) then begin
     set_plot,'ps'
     savepath = '~/strahl/widget/result/'+u.savepath
     device,/encapsulated,/color,filename=savepath+'/CHI_spatial.eps'
  endif else begin
     widget_control,u.id.draw_CHI_spatial,get_value=draw_win
     wset,draw_win
  endelse

  Nt=u.stat.coef.coef_Nt
  chi=*u.stat.coef.chi
  x=*u.stat.coef.coefx
  t=*u.stat.coef.coeft
  Np=5                          ; plot 5 profiles
  dt=(Nt-1)/(Np-1)
  index=findgen(Np)*dt
  colors=(findgen(Np)+1)*(200./Np)
  styles=findgen(Np)
  text=string(t(index))+ ' s'
  plot,[0,0],[0,0],/nodata,xrange=[0,1],yrange=[0,max(chi)*1.1],xstyle=1,ystyle=1,xtit='r',ytit='chi',charsize=1
  for it=0,Np-1 do oplot,x,chi[*,index(it)],color=colors[it],linestyle=styles[it]
  legend,text,colors=colors,linestyle=styles  

  if keyword_set(ps) then begin
     device,/close
     set_plot,'x'
  endif

end


pro wstrahl_plot_CHI_temporal,u,ps=ps
  if keyword_set(ps) then begin
     set_plot,'ps'
     savepath = '~/strahl/widget/result/'+u.savepath
     device,/encapsulated,/color,filename=savepath+'/CHI_temporal.eps'
  endif else begin
     widget_control,u.id.draw_CHI_temporal,get_value=draw_win
     wset,draw_win
  endelse
  Nr=u.stat.coef.coef_Nr
  chi=*u.stat.coef.chi
  x=*u.stat.coef.coefx
  t=*u.stat.coef.coeft
  Np=5                          ; plot 5 profiles
  dr=(Nr-1)/(Np-1)
  index=findgen(Np)*dr
  colors=(findgen(Np)+1)*(200./Np)
  styles=findgen(Np)
  text='r ='+string(x(index)) 
  plot,[0,0],[0,0],/nodata,xrange=[min(t),max(t)],yrange=[0,max(chi)*1.1],xstyle=1,ystyle=1,xtit='time [sec]',ytit='chi',charsize=1
  for ir=0,Np-1 do oplot,t,chi[index(ir),*],color=colors[ir],linestyle=styles[ir]
  legend,text,colors=colors,linestyle=styles  

  if keyword_set(ps) then begin
     device,/close
     set_plot,'x'
  endif
end


pro wstrahl_calc_V,u
  ;Gindex=widget_info(u.id.VGmodel,/droplist_select)
  ;Hindex=widget_info(u.id.VHmodel,/droplist_select)
  Gindex=u.stat.coef.VGdroplist
  Hindex=u.stat.coef.VHdroplist

  Nr=u.stat.coef.coef_Nr
  Nt=u.stat.coef.coef_Nt
  t1=u.t1
  t2=u.t2
  x=findgen(Nr+1.)/(Nr)
  t=findgen(Nt)/float(Nt)*(t2-t1)+t1
  v=u.stat.coef.v
  gx=fltarr(Nr)
  ht=fltarr(Nt)
  ;if Gindex ge 0 then gx+=u.stat.coef.vGparam.a0
  ;if Gindex ge 1 then gx+=u.stat.coef.vGparam.a1*x
  ;if Gindex ge 2 then gx+=u.stat.coef.vGparam.a2*x^2
  ;if Gindex ge 3 then gx+=u.stat.coef.vGparam.a3*x^3
  ;if Gindex ge 4 then gx+=u.stat.coef.vGparam.a4*x^4
  c0 = u.stat.coef.VGparam.c0
  a0 = u.stat.coef.VGparam.a0
  a1 = u.stat.coef.VGparam.a1
  a2 = u.stat.coef.VGparam.a2
  a3 = u.stat.coef.VGparam.a3
  x0 = u.stat.coef.VGparam.x0

  if Gindex eq 0 then begin
     tmph = fltarr(Nr) + 1.0
     if ((where(x lt x0))[0] ne -1) then tmph(where(x lt x0)) = 0.0
     gx = c0 + (a0+a1*(x-x0) + a2*(x-x0)^2 + a3*(x-x0)^3)*tmph
  endif else if Gindex eq 1 then begin
     print,'gaussian profile'
     gx = a0*exp(-(x-a1)^2/(2*a2^2))+a3+a4*x
  endif else begin
     gx = 0.0
  endelse 

  if Hindex ge 0 then ht+=u.stat.coef.vHparam.b0
  index=where(t ge u.stat.coef.vHparam.t0,count)
  if count gt 0 then ht[index]+=u.stat.coef.vHparam.c
  if Hindex ge 1 then ht+=u.stat.coef.vHparam.b1*(t-t(0))
  if Hindex ge 2 then ht+=u.stat.coef.vHparam.b2*(t-t(0))^2
  if Hindex ge 3 then ht+=u.stat.coef.vHparam.b3*(t-t(0))^3
  *u.stat.coef.v=gx#ht
  *u.stat.coef.coefx=x
  *u.stat.coef.coeft=t
end


pro wstrahl_plot_V_spatial,u,ps=ps
  if keyword_set(ps) then begin
     set_plot,'ps'
     savepath = '~/strahl/widget/result/'+u.savepath
     device,/encapsulated,/color,filename=savepath+'/V_spatial.eps'
  endif else begin
     widget_control,u.id.draw_v_spatial,get_value=draw_win
     wset,draw_win
  endelse

  Nt=u.stat.coef.coef_Nt
  v=*u.stat.coef.v
  x=*u.stat.coef.coefx
  t=*u.stat.coef.coeft
  Np=5                          ; plot 5 profiles
  dt=(Nt-1)/(Np-1)
  index=findgen(Np)*dt
  colors=(findgen(Np)+1)*(200./Np)
  styles=findgen(Np)
  text=string(t(index))+ ' s'
  plot,[0,0],[0,0],/nodata,xrange=[0,1],yrange=[min(v)-0.5,max(v)+0.5],xstyle=1,xtit='r',ytit='v',charsize=1
  for it=0,Np-1 do oplot,x,v[*,index(it)],color=colors[it],linestyle=styles[it]
  legend,text,colors=colors,linestyle=styles  

  if keyword_set(ps) then begin
     device,/close
     set_plot,'x'
  endif

end


pro wstrahl_plot_V_temporal,u,ps=ps
  if keyword_set(ps) then begin
     set_plot,'ps'
     savepath = '~/strahl/widget/result/'+u.savepath
     device,/encapsulated,/color,filename=savepath+'/V_temporal.eps'
  endif else begin
     widget_control,u.id.draw_V_temporal,get_value=draw_win
     wset,draw_win
  endelse

  Nr=u.stat.coef.coef_Nr
  v=*u.stat.coef.v
  x=*u.stat.coef.coefx
  t=*u.stat.coef.coeft
  Np=5                          ; plot 5 profiles
  dr=(Nr-1)/(Np-1)
  index=findgen(Np)*dr
  colors=(findgen(Np)+1)*(200./Np)
  styles=findgen(Np)
  text='r ='+string(x(index)) 
  plot,[0,0],[0,0],/nodata,xrange=[min(t),max(t)],yrange=[min(v),max(v)],xstyle=1,xtit='time [sec]',ytit='v',charsize=1
  for ir=0,Np-1 do oplot,t,v[index(ir),*],color=colors[ir],linestyle=styles[ir]
  legend,text,colors=colors,linestyle=styles  
  
  if keyword_set(ps) then begin
     device,/close
     set_plot,'x'
  endif
end


pro wstrahl_plotparam,u,ps=ps
  if not keyword_set(ps) then ps=0
  wstrahl_plot_binning,u,ps=ps
  wstrahl_plot_source,u,ps=ps
  wstrahl_calc_CHI,u
  wstrahl_plot_CHI_spatial,u,ps=ps
  wstrahl_plot_CHI_temporal,u,ps=ps
  wstrahl_calc_V,u
  wstrahl_plot_V_spatial,u,ps=ps
  wstrahl_plot_V_temporal,u,ps=ps
end

;==============================
pro wstrahl_genplt,array,ivec,jvec,labels=labels,cct=cct,pct=pct,ncntrs=ncntrs,set=set,$
                   dp=dp,chars=chars,thick=thick,ir=ir,jr=jr,ps=ps,prefix=prefix,savepath=savepath,norm=norm
  x=size(array)
  IF NOT keyword_set(ivec) THEN ivec=indgen(x[1])
  IF NOT keyword_set(jvec) THEN jvec=indgen(x[2])
  IF NOT keyword_set(norm) THEN norm=0

  ;set plotting variables
  ;ir=[min(ivec), max(ivec)]
  ;hard coded for now
  if not keyword_set(ir) then ir = [0.0,1.0]
  if not keyword_set(jr) then jr=[min(jvec), max(jvec)]
  IF NOT keyword_set(ncntrs) THEN  ncntrs=18.0
  extremepts=localmaxmin(array,ivec,jvec,ir,jr)
  maxplot=extremepts[1]*1.1
  minplot=0.0 < extremepts[0]*1.1 
  io=make(ir[0],ir[1],6)
  io=io[1:4]
  jo=make(jr[0],jr[1],6)
  jo=jo[1:4]
  IF NOT keyword_set(dp) THEN dp=2
  IF NOT keyword_set(chars) THEN begin
     if not keyword_set(ps) then chars = 1.0 else chars = 1.0
  endif
  IF NOT keyword_set(thick) THEN begin
     if not keyword_set(ps) then thick = 1.0 else thick = 2.0
  endif

  IF NOT keyword_set(cct) THEN cct=12
  IF NOT keyword_set(pct) THEN pct=12

  ;correnct selection of i,j ranges outside of ivec,jvec
  ;IF (ir[0] LT min(ivec)) THEN ir[0]=min(ivec)
  ;IF (ir[1] GT max(ivec)) THEN ir[1]=max(ivec)
  ;IF (jr[0] LT min(jvec)) THEN jr[0]=min(jvec)
  ;IF (jr[1] GT max(jvec)) THEN jr[1]=max(jvec)


  ;plot contour
  savepath = '~/strahl/widget/result/'+savepath
  if not keyword_set(norm) then begin
     if not keyword_set(ps) then ps = 0
     if keyword_set(ps) then begin
        if not keyword_set(prefix) then prefix = ''
        set_plot,'ps'
        device,/encapsulated,/color,filename=savepath+'/'+prefix+'contour.eps'
     endif else begin
        wset,set[0]
     endelse

  ;this is weird, but it seems to be necessary, to get good contours, 
  ;must truncate the emiss array to the plot size
     loadct,cct,/silent
     tmp1=where(jvec GE jr[0] AND jvec LE jr[1])
     tmp2=where(ivec GE ir[0] AND ivec LE ir[1])
     array_temp1=array[*,tmp1]
     array_plot=array_temp1[tmp2,*]
     
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
     levels=make(minplot,maxplot,ncntrs)
     position=[0.085,0.08,0.825,0.95]	
     ;stop
     contour,array_plot,ivec[tmp2],jvec[tmp1],min=minplot,max=maxplot,fill=1,levels=levels,$
             xtit=xtit,ytit=ytit,tit=tit,xrange=ir,yrange=jr,xstyle=1,ystyle=1,position=position,chars=chars   
     cbar=fltarr(2,n(levels)+1)
     cbar[0,*]=levels
     cbar[1,*]=levels
     position=[0.92,0.08,0.97,0.95]
     contour,cbar,[0,1],levels,fill=1,levels=levels,ytit=ztit,position=position,/noerase,chars=chars,$
             /xsty,/ysty,xticklayout=0,xticks=1,xtickname=[' ',' ']

     if keyword_set(ps) then begin
        device,/close
        set_plot,'x'
     endif
  endif


  ;plot the chosen or deafult i-coords
  if keyword_set(ps) then begin
     if not keyword_set(prefix) then prefix = ''
     set_plot,'ps'
     ;savepath = '~/strahl/widget/result/'+savepath
     if not keyword_set(norm) then begin
        device,/encapsulated,/color,filename=savepath+'/'+prefix+'iaxis.eps'
     endif else begin
        device,/encapsulated,/color,filename=savepath+'/'+prefix+'iaxis_norm.eps'
     endelse
  endif else begin
     wset,set[1]
  endelse

  loadct,pct,/silent
  iplots=fltarr(n(io)+1,n(jvec)+1)
  FOR k=0,n(io) DO iplots[k,*]=array_slice(array,ivec,jvec,i=io[k])
  colormap=colormap(io)
  maxplot=max(iplots[*,where(jvec GE jr[0] AND jvec LE jr[1])])
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
  if keyword_set(norm) then begin
     yrange=[0,1.2]
  endif else begin
     yrange=[minplot,maxplot*1.1]
  endelse 



  plot,[0],[0],xtit=xtit,ytit=ytit,tit=tit,xrange=jr,yrange=yrange,chars=chars,/xsty,/ysty
  xyouts,0.17,0.85,note, /norm
  FOR i=0,n(io) DO BEGIN
     if keyword_set(norm) then begin
        ;calculate the maxium value
        index=where(jvec ge jr[0] and jvec le jr[1])
        scale=max(iplots[i,index])*1.0
     endif else begin
        scale=1.0
     endelse
     oplot,jvec, iplots[i,*]/scale,color=colormap[i],thick=thick
     xyouts,0.2,0.80-0.05*i,num2str(io[i],dp=dp),color=colormap[i], chars=chars,/norm
  ENDFOR
  IF minplot LT 0 THEN oplot,jr,[0,0],linestyle=2
  if keyword_set(ps) then begin
     device,/close
     set_plot,'x'
  endif


  ;plot the chosen or deafult j-coords	
  if not keyword_set(norm) then begin
     if keyword_set(ps) then begin
        if not keyword_set(prefix) then prefix = ''
        set_plot,'ps'
        ;savepath = '~/strahl/widget/result/'+savepath
        device,/encapsulated,/color,filename=savepath+'/'+prefix+'jaxis.eps'
     endif else begin
        wset,set[2]
     endelse

     jplots=fltarr(n(jo)+1,n(ivec)+1)
     FOR k=0,n(jo) DO jplots[k,*]=array_slice(array,ivec,jvec,j=jo[k])
     colormap=colormap(jo)
     maxplot=max(jplots[*,where(ivec GE ir[0] AND ivec LE ir[1])])
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
     if keyword_set(ps) then begin
        device,/close
        set_plot,'x'
     endif
  endif

end


;pro wstrahl_plotdens_contour,u
;  widget_control,u.id.draw_contour_density,get_value=draw_contour
;end

;pro wstrahl_plotdens_t,u
;  widget_control,u.id.draw_time_density,get_value=draw_time
;end

;pro wstrahl_plotdens_r,u
;  widget_control,u.id.draw_radius_density,get_value=draw_radius
;end

pro wstrahl_plot_sim_dens_r_chargestate,u,ps=ps
  if not keyword_set(ps) then ps = 0


  widget_control,u.id.draw_radius_chargestate_density,get_value=set
  ;widget_control,u.id.density_time,get_value=t
  ;widget_control,u.id.zid,get_value=z
  ;widget_control,u.id.zplotid,get_value=zplot
  t = u.stat.times.density_time
  data = *u.stat.strahl.data
  csden = data.csden
  time = data.time
  rho = data.rho
  it = ipt(time,t)
  csden = csden[*,*,it]
  ncs = (size(csden,/dim))[1]
  yr=[0,max(csden)*1.2]
  
  if keyword_set(ps) then begin
     widget_control,u.id.savepath,get_value=savepath
     savepath = '~/strahl/widget/result/'+savepath
     set_plot,'ps'     
     device,/encapsulated,/color,filename=savepath+'/dens_r_chargestate.eps'
  endif else begin
     wset,set
  endelse

  plot,rho,csden[*,u.zplot],xr=u.stat.strahl.rr,yr=yr,/nodata,xtit='r/a',ytit='Charge State Density'
  for ics=0,ncs-1 do oplot,rho,csden[*,ics]
  color0 = 200 ;red
  color1 = 100 ;blue
  color2 = 50  ;green
  oplot,rho,csden[*,u.zplot],color=color0
  y1=max(csden[*,u.zplot],x1)
  xyouts,rho(x1),y1,num2elem(u.z)+'!u'+num2str(u.zplot,1)+'+!n',color=color0,charsize=2

  if (u.zplot lt ncs - 1) then begin
     oplot,rho,csden[*,u.zplot+1],color=color1
     y2 = max(csden[*,u.zplot+1],x2)
     xyouts,rho(x2),y2,num2elem(u.z)+'!u'+num2str(u.zplot+1,1)+'+!n',color=color1,charsize=2
  endif


  if (u.zplot lt ncs - 2) then begin
     oplot,rho,csden[*,u.zplot+2],color=color2
     y3 = max(csden[*,u.zplot+2],x3)
     xyouts,rho(x3),y3,num2elem(u.z)+'!u'+num2str(u.zplot+2,1)+'+!n',color=color2,charsize=2
  endif

  if keyword_set(ps) then begin
     device,/close
     set_plot,'x'
  endif

end

;==============================
;==============================
;plot experimental (strahl)


pro wstrahl_plot_sim_dens,u,ps=ps
  if not keyword_set(ps) then ps = 0
  ;wstrahl_plotdens_contour,u
  ;wstrahl_plotdens_t,u
  ;wstrahl_plotdens_r,u
  data=*u.stat.strahl.data
  labels={ilab:'r/a',jlab:'Time [sec]',klab:'Charge State Density',ctit:'Time Evolution of '$
          +num2elem(u.Z)+'!u'+num2str(u.Zplot,1)+'+!n',itit:'',jtit:''}
  help,ps
  wstrahl_genplt,reform(data.csden[*,u.Zplot,*]),data.rho,data.time,$
                 labels=labels,cct=39,ncntrs=50,set=u.stat.strahl.set_dens_abs,$
                 ir=u.stat.strahl.rr,jr=u.stat.strahl.tr,$
                 ps=ps,prefix='dens_',savepath=u.savepath
  
  wstrahl_genplt,reform(data.csden[*,u.Zplot,*]),data.rho,data.time,$
                 labels=labels,cct=39,ncntrs=50,set=u.stat.strahl.set_dens_norm,$
                 ir=u.stat.strahl.rr,jr=u.stat.strahl.tr,$
                 ps=ps,prefix='dens_',savepath=u.savepath,/norm
  wstrahl_plot_sim_dens_r_chargestate,u,ps=ps
end

;pro wstrahl_plotemiss_contour,u
  
;end

;pro wstrahl_plotemiss_t,u

;end

;pro wstrahl_plotemiss_r,u

;end



pro wstrahl_plot_sim_emiss,u,ps=ps
  if not keyword_set(ps) then ps = 0
  ;wstrahl_plotemiss_contour,u
  ;wstrahl_plotemiss_t,u
  ;wstrahl_plotemiss_r,u

  ;;plot w line sim emissivity
  emiss=(*u.stat.strahl.emiss_w).emiss
  rho=(*u.stat.strahl.emiss_w).rho
  time=(*u.stat.strahl.emiss_w).time
  labels={ilab:'r/a',jlab:'Time [sec]',klab:'W line Emissivity',ctit:'Time Evolution of '$
          +num2elem(u.Z)+'!u'+num2str(u.Zplot,1)+'+!n W line Emissivity',itit:'',jtit:''}
  ;labels={ilab:'r/a',jlab:'Time [sec]',klab:'W line Emissivity',ctit:'Time Evolution of '$
  ;        +num2elem(u.Z)+'!u'+num2str(u.Zplot,1)+'+!n',itit:'',jtit:''}
  wstrahl_genplt,emiss,rho,time,$
                 labels=labels,cct=39,ncntrs=50,set=u.stat.strahl.set_emiss_sim_w_abs,$
                 ir=u.stat.strahl.rr,jr=u.stat.strahl.tr,$
                 ps=ps,prefix='emiss_sim_w_',savepath=u.savepath
  wstrahl_genplt,emiss,rho,time,$
                 labels=labels,cct=39,ncntrs=50,set=u.stat.strahl.set_emiss_sim_w_norm,$
                 ir=u.stat.strahl.rr,jr=u.stat.strahl.tr,$
                 ps=ps,prefix='emiss_sim_w_',savepath=u.savepath,/norm

  ;;plot z line sim emissivity
  emiss=(*u.stat.strahl.emiss_z).emiss
  rho=(*u.stat.strahl.emiss_z).rho
  time=(*u.stat.strahl.emiss_z).time
  labels={ilab:'r/a',jlab:'Time [sec]',klab:'W line Emissivity',ctit:'Time Evolution of '$
          +num2elem(u.Z)+'!u'+num2str(u.Zplot,1)+'+!n Z line Emissivity',itit:'',jtit:''}
  ;labels={ilab:'r/a',jlab:'Time [sec]',klab:'Z line Emissivity',ctit:'Time Evolution of '$
  ;        +num2elem(u.Z)+'!u'+num2str(u.Zplot,1)+'+!n',itit:'',jtit:''}
  wstrahl_genplt,emiss,rho,time,$
                 labels=labels,cct=39,ncntrs=50,set=u.stat.strahl.set_emiss_sim_z_abs,$
                 ir=u.stat.strahl.rr,jr=u.stat.strahl.tr,$
                 ps=ps,prefix='emiss_sim_z_',savepath=u.savepath
  wstrahl_genplt,emiss,rho,time,$
                 labels=labels,cct=39,ncntrs=50,set=u.stat.strahl.set_emiss_sim_z_norm,$
                 ir=u.stat.strahl.rr,jr=u.stat.strahl.tr,$
                 ps=ps,prefix='emiss_sim_z_',savepath=u.savepath,/norm
end


pro wstrahl_plot_sim_bright,u,ps=ps
  if not keyword_set(ps) then ps = 0
  ;wstrahl_plotemiss_contour,u
  ;wstrahl_plotemiss_t,u
  ;wstrahl_plotemiss_r,u

  ;;plot w line sim brightness
  bright=(*u.stat.strahl.bright_w).bright
  rho=(*u.stat.strahl.bright_w).rho
  time=(*u.stat.strahl.bright_w).time
  
  ;only take half
  dum = min(rho,imax)
  bright = bright(0:imax,*)
  rho = rhO(0:imax)

  labels={ilab:'r/a',jlab:'Time [sec]',klab:'W line Brightness',ctit:'Time Evolution of '$
          +num2elem(u.Z)+'!u'+num2str(u.Zplot,1)+'+!n W line Brightness',itit:'',jtit:''}
  wstrahl_genplt,bright,rho,time,$
                 labels=labels,cct=39,ncntrs=50,set=u.stat.strahl.set_bright_sim_w_abs,$
                 ir=u.stat.strahl.rr,jr=u.stat.strahl.tr,$
                 ps=ps,prefix='bright_sim_w_',savepath=u.savepath
  wstrahl_genplt,bright,rho,time,$
                 labels=labels,cct=39,ncntrs=50,set=u.stat.strahl.set_bright_sim_w_norm,$
                 ir=u.stat.strahl.rr,jr=u.stat.strahl.tr,$
                 ps=ps,prefix='bright_sim_w_',savepath=u.savepath,/norm
  ;;plot z line sim brightness
  bright=(*u.stat.strahl.bright_z).bright
  rho=(*u.stat.strahl.bright_z).rho
  time=(*u.stat.strahl.bright_z).time
  
  ;only take half
  dum = min(rho,imax)
  bright = bright(0:imax,*)
  rho = rhO(0:imax)

  labels={ilab:'r/a',jlab:'Time [sec]',klab:'Z line Brightness',ctit:'Time Evolution of '$
          +num2elem(u.Z)+'!u'+num2str(u.Zplot,1)+'+!n Z line Brightness',itit:'',jtit:''}
  wstrahl_genplt,bright,rho,time,$
                 labels=labels,cct=39,ncntrs=50,set=u.stat.strahl.set_bright_sim_z_abs,$
                 ir=u.stat.strahl.rr,jr=u.stat.strahl.tr,$
                 ps=ps,prefix='bright_sim_z_',savepath=u.savepath
  wstrahl_genplt,bright,rho,time,$
                 labels=labels,cct=39,ncntrs=50,set=u.stat.strahl.set_bright_sim_z_norm,$
                 ir=u.stat.strahl.rr,jr=u.stat.strahl.tr,$
                 ps=ps,prefix='bright_sim_z_',savepath=u.savepath,/norm
end



pro wstrahl_plot_sim,u,ps=ps
  if not keyword_set(ps) then ps = 0
  wstrahl_plot_sim_dens,u,ps=ps
  wstrahl_plot_sim_emiss,u,ps=ps
  wstrahl_plot_sim_bright,u,ps=ps
end

;==============================
;==============================
;plot experimental (hirex)

pro wstrahl_plot_exp_emiss,u,ps=ps
  if not keyword_set(ps) then ps = 0

  ;;plot w line exp emissivity
  emiss=(*u.stat.hirex.emiss_w).emiss
  rho=(*u.stat.hirex.emiss_w).rho
  time=(*u.stat.hirex.emiss_w).time
  labels={ilab:'r/a',jlab:'Time [sec]',klab:'W line Emissivity',ctit:'Time Evolution of '$
          +num2elem(u.Z)+'!u'+num2str(u.Zplot,1)+'+!n W line Emissivity',itit:'',jtit:''}
  wstrahl_genplt,emiss,rho,time,$
                 labels=labels,cct=39,ncntrs=50,set=u.stat.hirex.set_emiss_exp_w_abs,$
                 ir=u.stat.hirex.rr,jr=u.stat.hirex.tr,$
                 ps=ps,prefix='emiss_exp_w_',savepath=u.savepath
  wstrahl_genplt,emiss,rho,time,$
                 labels=labels,cct=39,ncntrs=50,set=u.stat.hirex.set_emiss_exp_w_norm,$
                 ir=u.stat.hirex.rr,jr=u.stat.hirex.tr,$
                 ps=ps,prefix='emiss_exp_w_',savepath=u.savepath,/norm
  ;;plot z line exp emissivity
  emiss=(*u.stat.hirex.emiss_z).emiss
  rho=(*u.stat.hirex.emiss_z).rho
  time=(*u.stat.hirex.emiss_z).time
  labels={ilab:'r/a',jlab:'Time [sec]',klab:'Z line Emissivity',ctit:'Time Evolution of '$
          +num2elem(u.Z)+'!u'+num2str(u.Zplot,1)+'+!n Z line Emissivity',itit:'',jtit:''}
  wstrahl_genplt,emiss,rho,time,$
                 labels=labels,cct=39,ncntrs=50,set=u.stat.hirex.set_emiss_exp_z_abs,$
                 ir=u.stat.hirex.rr,jr=u.stat.hirex.tr,$
                 ps=ps,prefix='emiss_exp_z_',savepath=u.savepath
  wstrahl_genplt,emiss,rho,time,$
                 labels=labels,cct=39,ncntrs=50,set=u.stat.hirex.set_emiss_exp_z_norm,$
                 ir=u.stat.hirex.rr,jr=u.stat.hirex.tr,$
                 ps=ps,prefix='emiss_exp_z_',savepath=u.savepath,/norm
end


pro wstrahl_plot_exp_bright,u,ps=ps
  if not keyword_set(ps) then ps = 0
  
  ;;plot w line exp brightness
  bright=(*u.stat.hirex.bright_w).bright
  rho=(*u.stat.hirex.bright_w).rho
  time=(*u.stat.hirex.bright_w).time

  ;only take half
  dum = min(rho,imax)
  bright = bright(0:imax,*)
  rho = rhO(0:imax)

  labels={ilab:'r/a',jlab:'Time [sec]',klab:'W line Brightness',ctit:'Time Evolution of '$
          +num2elem(u.Z)+'!u'+num2str(u.Zplot,1)+'+!n W line Brightness',itit:'',jtit:''}
  wstrahl_genplt,bright,rho,time,$
                 labels=labels,cct=39,ncntrs=50,set=u.stat.hirex.set_bright_exp_w_abs,$
                 ir=u.stat.hirex.rr,jr=u.stat.hirex.tr,$
                 ps=ps,prefix='bright_exp_w_',savepath=u.savepath
  wstrahl_genplt,bright,rho,time,$
                 labels=labels,cct=39,ncntrs=50,set=u.stat.hirex.set_bright_exp_w_norm,$
                 ir=u.stat.hirex.rr,jr=u.stat.hirex.tr,$
                 ps=ps,prefix='bright_exp_w_',savepath=u.savepath,/norm
  ;;plot z line exp brightness
  bright=(*u.stat.hirex.bright_z).bright
  rho=(*u.stat.hirex.bright_z).rho
  time=(*u.stat.hirex.bright_z).time

  ;only take half
  dum = min(rho,imax)
  bright = bright(0:imax,*)
  rho = rhO(0:imax)


  labels={ilab:'r/a',jlab:'Time [sec]',klab:'Z line Brightness',ctit:'Time Evolution of '$
          +num2elem(u.Z)+'!u'+num2str(u.Zplot,1)+'+!n Z line Brightness',itit:'',jtit:''}
  wstrahl_genplt,bright,rho,time,$
                 labels=labels,cct=39,ncntrs=50,set=u.stat.hirex.set_bright_exp_z_abs,$
                 ir=u.stat.hirex.rr,jr=u.stat.hirex.tr,$
                 ps=ps,prefix='bright_exp_z_',savepath=u.savepath
  wstrahl_genplt,bright,rho,time,$
                 labels=labels,cct=39,ncntrs=50,set=u.stat.hirex.set_bright_exp_z_norm,$
                 ir=u.stat.hirex.rr,jr=u.stat.hirex.tr,$
                 ps=ps,prefix='bright_exp_z_',savepath=u.savepath,/norm
end

pro wstrahl_plot_exp,u,ps=ps
  if not keyword_set(ps) then ps = 0
  wstrahl_plot_exp_emiss,u,ps=ps
  wstrahl_plot_exp_bright,u,ps=ps
end

;==============================




pro wstrahl_cleanup,u

end

pro wstrahl_savestrahl,u
  if (not file_test('~/strahl/widget/result',/directory)) then file_mkdir,'~/strahl/widget/result'
  widget_control,u.id.savepath,get_value=savepath
  savepath = '~/strahl/widget/result/'+savepath
  if (not file_test(savepath,/directory)) then file_mkdir,savepath
  filename = savepath + '/result.sav'
  usave=u
  save,usave,filename=filename
  print,'file saved to '+filename
  wstrahl_plotparam,u,/ps
  wstrahl_plot_sim,u,/ps
  wstrahl_plot_exp,u,/ps
  ;wstrahl_plotwline,u,/ps
  ;wstrahl_plotzline,u,/ps
  print,'figures are saved to '+savepath
end


pro wstrahl_loadstrahl,u
  widget_control,u.id.loadpath,get_value=loadpath
  loadpath = 'result/'+loadpath
  filename = loadpath + '/result.sav'
  if (not file_test(filename)) then begin
     print,'no such file'
  endif else begin
     restore,filename=filename  ; usave will be restored
     wstrahl_copy,usave,u
     wstrahl_load_info,u
     wstrahl_plotparam,u
     
     catch,error_status
     if error_status NE 0 then begin
        print,'no strahl result available'
        catch,/cancel
        goto, JUMP
     endif
     print,'plotting strahl result: emissivity...'
     wstrahl_sim_plot_sim_dens,u
     if u.z eq 18 then begin
        print,'plotting strahl result: emissivity...'
        wstrahl_plot_sim_emiss,u
        wstrahl_plot_sim_bright,u
        wstrahl_plot_exp,u
     endif
     
     
     ;JUMP:
     ;catch,error_status2
     ;if error_status2 NE 0 then begin
     ;   print,'no hirex result available'
     ;   catch,/cancel
     ;   goto, JUMP2
     ;endif
     print,'loading emissivity profiles from HiReX-Sr'
     wstrahl_plotwline,u
     wstrahl_plotzline,u
     JUMP:
  endelse
end



pro w_strahl_event,event
  widget_control,event.top,get_uvalue=u
  id = u.id
  tag = tag_names(event,/st)
  button=' '
  idtags=tag_names(id)
  for i=0,n(idtags) do if id.(i) eq event.id then ename=idtags[i]
  case tag of
     "WIDGET_BASE": begin

     end
     "WIDGET_BUTTON": begin
        widget_control,event.id,get_value=button,get_uvalue=uvalue
        if not keyword_set(uvalue) then uvalue='none'
        case ename of
           "QUIT": begin
              widget_control,event.top,/destroy
              wstrahl_cleanup,u
              heap_gc
              !except=1
           end

           "LOAD":begin
              xwplot
              widget_control,/hourglass
              wstrahl_load_info,u
              wstrahl_plotparam,u
           end

           "LOAD_HIREX":begin
              print,'loading emissivity profiles from HiReX-Sr using cgao script'
              shot=u.shot
              line=0
              hirexsr_load_result,shot,profile,moment,lineint,line=line
              time=profile.tau
              rnorm=mean(profile.rnorm,dimen=2)
              rnormvar=variance(profile.rnorm,dimen=2)
              emiss=profile.prof.emiss
              index=where(time ge u.t1 and time le u.t2)
              time=time(index)
              emiss=emiss(*,index)
              *u.stat.hirex.emiss_w = {emiss:emiss,rho:rnorm,time:time}
              
              line=2
              hirexsr_load_result,shot,profile,moment,lineint,line=line
              time=profile.tau
              rnorm=mean(profile.rnorm,dimen=2)
              rnormvar=variance(profile.rnorm,dimen=2)
              emiss=profile.prof.emiss
              index=where(time ge u.t1 and time le u.t2)
              time=time(index)
              emiss=emiss(*,index)
              *u.stat.hirex.emiss_z = {emiss:emiss,rho:rnorm,time:time}
             
              wstrahl_plot_exp_emiss,u,ps=ps
           end

           "RUN":begin

              print,'STRAHL simulation is running. It will take a while. ......'
              ;if ((*u.fz)[0] eq -1) then source = -1*(*u.fztime) else source={fz:*u.fz,time:*u.fztime}
              source={fz:*u.fz,time:*u.fztime}
              ; calculate the ta, tb specified in the new /usr/local/cmod/idl/GENIE/GENTRAN/strahl/run_cmod_strahl.pro
              ; ta: start times of simulation
              ; tb: delta times of simulations
              delta = (u.t2 - u.t1)/u.nt
              ta = findgen(u.nt)*delta + u.t1
              tb = fltarr(u.nt) + delta
              RUN_CMOD_STRAHL,u.shot,u.z,ta,tb,source=source,$
                              fits=u.stat.profile.fits,$
                              nogrid=u.stat.optional.nogrid,nopp=u.stat.optional.nopp,$
                              noparam=u.stat.optional.noparam,$
                              nostrahl=u.stat.optional.nostrahl,quiet=u.stat.optional.quiet,$
                              diff={diff:*u.stat.coef.chi,rho:*u.stat.coef.coefx,time:*u.stat.coef.coeft},$
                              conv={conv:*u.stat.coef.v,rho:*u.stat.coef.coefx,time:*u.stat.coef.coeft},$
                              data=data,term=term
              *u.stat.strahl.data = data
              *u.stat.strahl.term = term
              print,'plotting strahl result: density...'
              wstrahl_plot_sim_dens,u
              print,'done'
              ;if u.z eq 18 then begin
              ;   print,'calculating emissivity and brightness profiles'         
              ;   emiss = den2emiss_cgao(data=data)
                 
              ;endif 
                                ;widget_control,event.top,set_uvalue=u
           end

           "EXP_SIM":begin
              if u.z eq 18 then begin
                 data = *u.stat.strahl.data
                 line = [0,2]   ; hard code: w & z line 
                 print,'calculating emissivity and brightness profiles'
                 result = hirexsr_genrad_argon_cgao(u.shot,line,u.t1,u.t2,csden,cserr,temp,temperr,dens,$
                                                 denserr,rhovec,tht=0,data=data)
                
                 wsim = *result.sim[0]
                 zsim = *result.sim[1]
                 wexp = *result.exp[0]
                 zexp = *result.exp[1]

                 ; for now I use mean() to conver brRho[t,r] to brRho[t].
                 ; this is not ideal. I should use interpolate in the future...
                 *u.stat.strahl.emiss_w = {emiss:wsim.emiss,rho:wsim.rho,time:wsim.time}
                 *u.stat.strahl.emiss_z = {emiss:zsim.emiss,rho:zsim.rho,time:zsim.time}
                 *u.stat.strahl.bright_w = {bright:wsim.bright,rho:mean(wsim.brRho,dim=2),time:wsim.time}
                 *u.stat.strahl.bright_z = {bright:zsim.bright,rho:mean(zsim.brRho,dim=2),time:zsim.time}
                 
                 *u.stat.hirex.emiss_w = {emiss:wexp.emiss,rho:wexp.rho,time:wexp.time}
                 *u.stat.hirex.emiss_z = {emiss:zexp.emiss,rho:zexp.rho,time:zexp.time}
                 *u.stat.hirex.bright_w = {bright:wexp.bright,rho:mean(wexp.brRho,dim=2),time:wexp.time}
                 *u.stat.hirex.bright_z = {bright:zexp.bright,rho:mean(zexp.brRho,dim=2),time:zexp.time}
                 
              
                 print,'plotting strahl result: emissivity...'
                 wstrahl_plot_sim,u
                 wstrahl_plot_exp,u
                 print,'done'           
              end
           end



           "PLOT":begin
              xwplot
              wstrahl_plot_sim,u
              wstrahl_plot_exp,u
           end
           
           
           "STOP":begin
              stop
           end
           
           "SAVEBUTTON":begin
              wstrahl_savestrahl,u
           end


           "LOADBUTTON":begin
              wstrahl_loadstrahl,u
           end


           
        endcase
     end
     
     

     "WIDGET_DROPLIST":begin
        widget_control,event.id,get_value=droplist
        case event.id of
           u.id.CHIGmodel: begin
             u.stat.coef.chiGdroplist = widget_info(u.id.CHIGmodel,/droplist_select)
           end
           u.id.CHIHmodel: begin
             u.stat.coef.chiHdroplist = widget_info(u.id.CHIHmodel,/droplist_select)
           end
           u.id.VGmodel: begin
             u.stat.coef.VGdroplist = widget_info(u.id.VGmodel,/droplist_select)
           end
           u.id.VHmodel: begin
             u.stat.coef.VHdroplist = widget_info(u.id.VHmodel,/droplist_select)
           end
        endcase
     end


     "WIDGET_TEXT_CH":begin
        widget_control,event.id,get_value=text
        case event.id of
           u.id.shotid: begin
              widget_control,u.id.shotid,get_value=shot
              u.shot=shot
           end
           u.id.t1id: begin
              widget_control,u.id.t1id,get_value=t1
              u.t1=t1
           end
           u.id.t2id: begin
              widget_control,u.id.t2id,get_value=t2
              u.t2=t2
           end
           u.id.Ntid: begin
              widget_control,u.id.Ntid,get_value=Nt
              u.Nt=Nt
           end
           u.id.Zid: begin
              widget_control,u.id.Zid,get_value=Z
              u.Z=Z
           end
           u.id.Zplotid: begin
              widget_control,u.id.Zplotid,get_value=Zplot
              u.Zplot=Zplot
           end
           u.id.fzid: begin
              widget_control,u.id.fzid,get_value=fz_text
              fz=double(strsplit(fz_text,',',/extract))
              u.fz_text=fz_text
              *u.fz=fz
              ;print,fztime
              ;print,fztime_text
           end
           u.id.fztimeid: begin
              widget_control,u.id.fztimeid,get_value=fztime_text
              fztime=double(strsplit(fztime_text,',',/extract))
              u.fztime_text=fztime_text
              *u.fztime=fztime
              ;print,fztime
              ;print,fztime_text
           end
           u.id.savepath: begin
              widget_control,u.id.savepath,get_value=savepath
              u.savepath = savepath
           end
           u.id.loadpath: begin
              widget_control,u.id.loadpath,get_value=loadpath
              u.loadpath = loadpath
           end
           
           u.id.strahltmin: begin
              widget_control,u.id.strahltmin,get_value=tmin
              u.stat.strahl.tr[0]=tmin
           end
           u.id.strahltmax: begin
              widget_control,u.id.strahltmax,get_value=tmax
              u.stat.strahl.tr[1]=tmax
           end
           u.id.strahlrmin: begin
              widget_control,u.id.strahlrmin,get_value=rmin
              u.stat.strahl.rr[0]=rmin
           end
           u.id.strahlrmax: begin
              widget_control,u.id.strahlrmax,get_value=rmax
              u.stat.strahl.rr[1]=rmax
           end
           u.id.hirextmin: begin
              widget_control,u.id.hirextmin,get_value=tmin
              u.stat.hirex.tr[0]=tmin
           end
           u.id.hirextmax: begin
              widget_control,u.id.hirextmax,get_value=tmax
              u.stat.hirex.tr[1]=tmax
           end
           u.id.hirexrmin: begin
              widget_control,u.id.hirexrmin,get_value=rmin
              u.stat.hirex.rr[0]=rmin
           end
           u.id.hirexrmax: begin
              widget_control,u.id.hirexrmax,get_value=rmax
              u.stat.hirex.rr[1]=rmax
           end           
           u.id.coef_Nr:begin
              widget_control,u.id.coef_Nr,get_value=coef_Nr
              u.stat.coef.coef_Nr=coef_Nr
           end
           u.id.coef_Nt:begin
              widget_control,u.id.coef_Nt,get_value=coef_Nt
              u.stat.coef.coef_Nt=coef_Nt              
           end
           u.id.CHIc0:begin
              widget_control,u.id.CHIc0,get_value=CHIc0
              u.stat.coef.chiGparam.c0=CHIc0
           end
           u.id.CHIa0:begin
              widget_control,u.id.CHIa0,get_value=CHIa0
              u.stat.coef.chiGparam.a0=CHIa0
           end
           u.id.CHIa1:begin
              widget_control,u.id.CHIa1,get_value=CHIa1
              u.stat.coef.chiGparam.a1=CHIa1
           end
           u.id.CHIa2:begin
              widget_control,u.id.CHIa2,get_value=CHIa2
              u.stat.coef.chiGparam.a2=CHIa2
           end
           u.id.CHIa3:begin
              widget_control,u.id.CHIa3,get_value=CHIa3
              u.stat.coef.chiGparam.a3=CHIa3
           end
           u.id.CHIx0:begin
              widget_control,u.id.CHIx0,get_value=CHIx0
              u.stat.coef.chiGparam.x0=CHIx0
           end
          u.id.CHIc:begin
             widget_control,u.id.CHIc,get_value=CHIc
             u.stat.coef.chiHparam.c=CHIc
          end
          u.id.CHIt0:begin
             widget_control,u.id.CHIt0,get_value=CHIt0
             u.stat.coef.chiHparam.t0=CHIt0
          end
          u.id.CHIb0:begin
             widget_control,u.id.CHIb0,get_value=CHIb0
             u.stat.coef.chiHparam.b0=CHIb0
          end
          u.id.CHIb1:begin
             widget_control,u.id.CHIb1,get_value=CHIb1
             u.stat.coef.chiHparam.b1=CHIb1
          end
          u.id.CHIb2:begin
             widget_control,u.id.CHIb2,get_value=CHIb2
             u.stat.coef.chiHparam.b2=CHIb2
          end
          u.id.CHIb3:begin
             widget_control,u.id.CHIb3,get_value=CHIb3
             u.stat.coef.chiHparam.b3=CHIb3
          end
          u.id.Vc0:begin
              widget_control,u.id.Vc0,get_value=Vc0
              u.stat.coef.vGparam.c0=Vc0
           end
           u.id.Va0:begin
              widget_control,u.id.Va0,get_value=Va0
              u.stat.coef.vGparam.a0=Va0
           end
           u.id.Va1:begin
              widget_control,u.id.Va1,get_value=Va1
              u.stat.coef.vGparam.a1=Va1
           end
           u.id.Va2:begin
              widget_control,u.id.Va2,get_value=Va2
              u.stat.coef.vGparam.a2=Va2
           end
           u.id.Va3:begin
              widget_control,u.id.Va3,get_value=Va3
              u.stat.coef.vGparam.a3=Va3
           end
           u.id.Vx0:begin
              widget_control,u.id.Vx0,get_value=Vx0
              u.stat.coef.vGparam.x0=Vx0
           end
          u.id.Vc:begin
             widget_control,u.id.Vc,get_value=Vc
             u.stat.coef.vHparam.c=Vc
          end
          u.id.Vt0:begin
             widget_control,u.id.Vt0,get_value=Vt0
             u.stat.coef.vHparam.t0=Vt0
          end
          u.id.Vb0:begin
             widget_control,u.id.Vb0,get_value=Vb0
             u.stat.coef.vHparam.b0=Vb0
          end
          u.id.Vb1:begin
             widget_control,u.id.Vb1,get_value=Vb1
             u.stat.coef.vHparam.b1=Vb1
          end
          u.id.Vb2:begin
             widget_control,u.id.Vb2,get_value=Vb2
             u.stat.coef.vHparam.b2=Vb2
          end
          u.id.Vb3:begin
             widget_control,u.id.Vb3,get_value=Vb3
             u.stat.coef.vHparam.b3=Vb3
          end
          u.id.density_time:begin
             widget_control,u.id.density_time,get_value=density_time
             u.stat.times.density_time = density_time
             widget_control,u.id.density_time_slider,$
                 set_value= int(2000*(u.stat.times.density_time - u.t1)/(u.t2-u.t1))
             wstrahl_plot_sim_dens_r_chargestate,u
         end
       endcase


     end


     "WIDGET_SLIDER":begin
        widget_control,event.id,get_value=text
        case event.id of
           u.id.density_time_slider:begin
              widget_control,u.id.density_time_slider,get_value = val
              time = (u.t2 - u.t1)*val/2000. + u.t1
              u.stat.times.density_time = time
              widget_control,u.id.density_time,set_value = num2str(time,dp=2)
              wstrahl_plot_sim_dens_r_chargestate,u
           end
        endcase
     end



     "WIDGET_DRAW":begin
        
     end
     
     else:
     
  endcase

  IF button NE 'QUIT' THEN begin
     widget_control,event.top,set_uvalue=u
  ENDIF

end



pro w_strahl,shot=shot,z=z,Zplot=Zplot,t1=t1,t2=t2,Nt=Nt
  user=(get_login_info()).user_name
  loadct,12,/silent
  mdim=get_screen_size()
  IF mdim[0] NE 1600 AND mdim[1] NE 1200 THEN $
     base=widget_base(title='STRAHL Analysis Widget',/row,tlb_size_events=1,/scroll,$
     x_scroll_size=mdim[0]*0.95,y_scroll_size=mdim[1]*0.90) $
  ELSE base=widget_base(title='STRAHL Analysis Widge',/row,tlb_size_events=1)  
; Parameter panel
  A=widget_base(base,/column)
; base parameters
  A1=widget_base(A,/column,xsize=500,frame=5)
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  A1a=widget_base(A1,/col,frame=1)
  dum=widget_label(A1a,value='Shot information:',/align_left)
  A1a1=widget_base(A1a,/row)
  dum=widget_label(A1a1,value='shot:')
  shotid=widget_text(A1a1,xsize=10,/edit)
  dum=widget_label(A1a1,value='t1:')
  t1id=widget_text(A1a1,xsize=5,/edit)
  dum=widget_label(A1a1,value='t2:')
  t2id=widget_text(A1a1,xsize=5,/edit)
  dum=widget_label(A1a1,value='Nt:')
  Ntid=widget_text(A1a1,xsize=5,/edit)

  A1b=widget_base(A1,/col,frame=1)
  dum=widget_label(A1b,value='Impurity information:',/align_left)
  A1b1=widget_base(A1b,/row)
  dum=widget_label(A1b1,value='Z: ')
  Zid=widget_text(A1b1,xsize=3,/edit)
  dum=widget_label(A1b1,value='Charge state to plot: ')
  Zplotid=widget_text(A1b1,xsize=3,/edit)

  A1c=widget_base(A1,/col,frame=1)
  dum=widget_label(A1c,value='Source Setup: ',/align_left)
  A1c1=widget_base(A1c,/row)
  dum=widget_label(A1c1,value='fz*1e17 (array):     ')
  fzid=widget_text(A1c1,xsize=50,/edit)
  A1c2=widget_base(A1c,/row)
  dum=widget_label(A1c2,value='fztime (array): ')
  fztimeid=widget_text(A1c2,xsize=50,/edit)

  A1d=widget_base(A1,/col,frame=1)
  dum=widget_label(A1d,value='Profile from: ',/align_left)
  A1d1=widget_base(A1d,/row,/exclusive)
  fits=widget_button(A1d1,value='FiTs file',tooltip='load T_e and N_e profile from fits result in ~/fits/fits_shotnumber.save')
  qfits=widget_button(A1d1,value='quickFit file',tooltip='use the "quick_fits" program in /usr/local/cmod/idl/ to compute the T_e and n_e profiles')
  
  A1e=widget_base(A1,/col,frame=1)
  dum=widget_label(A1e,value='Optional Parameters: ',/align_left)
  A1e1=widget_base(A1e,/row,/nonexclusive)
  nogrid=widget_button(A1e1,value='nogrid',tooltip='skips MAKE_CMOD_GRID and does not copy the output file from ~/strahl/cmod to ~/strahl/nete')
  nopp=widget_button(A1e1,value='nopp',tooltip='skips MAKE_CMOD_PP and does not copy the output file from ~/strahl/cmod to ~/strahl/nete')
  noparam=widget_button(A1e1,value='noparam',tooltip='skips MAKE_CMOD_PARAM and does not copy the output file from ~/strahl/cmod to ~/strahl/param_files')
  nostrahl=widget_button(A1e1,value='nostrahl',tooltip='skips executing STRAHL as well as the MAKE_CMOD*, useful for viewing already run data')
  quiet=widget_button(A1e1,value='quiet',tooltip='suppress terminal notifications ')

  
  A1f=widget_base(A1,/col,frame=1)
  dum=widget_label(A1f,value='Action Controls: ',/align_left)
  A1f1=widget_base(A1f,/row)
  load=widget_button(A1f1,value='LOAD',tooltip='Load and visualize the simulation parameters')
  load_hirex=widget_button(A1f1,value='LOAD HiReX',tooltip='Load the emissivity of w and z lines from HiReX (using /home/cgao/idl/hirex/hirexsr_load_result.pro)')
 
  run=widget_button(A1f1,value='RUN',tooltip='Run the STRAHL simulation with current parameters')
 exp_sim=widget_button(A1f1,value='Exp v.s. Sim',$
                            tooltip='compare experimental and simulation emiisivity and brightness profiles (using ~/strahl/hirexsr_genrad_argon)')
  plot=widget_button(A1f1,value='PLOT',tooltip='PLOT the STRAHL result')
  stop=widget_button(A1f1,value='STOP',tooltip='Stop and give control to terminal interactive mode')
  quit=widget_button(A1f1,value='QUIT',tooltip='Quit the widget')
  
  A1f2=widget_base(A1f,/row)
  dum=widget_label(A1f2,value='save to: ',/align_left)
  savepath=widget_text(A1f2,xsize=40,/edit)
  savebutton=widget_button(A1f2,value='save',tooltip='save ufile and figures to ~/strahl/widget/result')
  
  A1f3=widget_base(A1f,/row)
  dum=widget_label(A1f3,value='load from: ',/align_left)
  loadpath=widget_text(A1f3,xsize=40,/edit)
  loadbutton=widget_button(A1f3,value='load',tooltip='load ufile and figures from ~/strahl/widget/result')


  A1g=widget_base(A1,/col,frame=1)
  dum=widget_label(A1f,value='Plot Controls: ',/align_left)
  A1g1=widget_base(A1g,/row)
  dum=widget_label(A1g1,value='STRAHL ',/align_left)
  dum=widget_label(A1g1,value='tmin: ')
  strahltmin=widget_text(A1g1,xsize=5,/edit)
  dum=widget_label(A1g1,value='tmax: ')
  strahltmax=widget_text(A1g1,xsize=5,/edit)
  dum=widget_label(A1g1,value='rmin: ')
  strahlrmin=widget_text(A1g1,xsize=5,/edit)
  dum=widget_label(A1g1,value='rmax: ')
  strahlrmax=widget_text(A1g1,xsize=5,/edit)
  A1g2=widget_base(A1g,/row)
  dum=widget_label(A1g2,value='HIREX ',/align_left)
  dum=widget_label(A1g2,value='tmin: ')
  hirextmin=widget_text(A1g2,xsize=5,/edit)
  dum=widget_label(A1g2,value='tmax: ')
  hirextmax=widget_text(A1g2,xsize=5,/edit)
  dum=widget_label(A1g2,value='rmin: ')
  hirexrmin=widget_text(A1g2,xsize=5,/edit)
  dum=widget_label(A1g2,value='rmax: ')
  hirexrmax=widget_text(A1g2,xsize=5,/edit)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  A2=widget_base(A,/column,xsize=500,frame=5)
  ; profile average time binning
  A2a=widget_base(A2,/col,frame=1)
  dum=widget_label(A2a,value='Time Binning of Profile Averaging: ', /align_left)
  draw_binning=widget_draw(A2a,xsize=500,ysize=175,frame=3)
  ; source strength v.s. time
  A2b=widget_base(A2,/col,frame=1)
  dum=widget_label(A2b,value='Source v.s. Time: ', /align_left)
  draw_source=widget_draw(A2b,xsize=500,ysize=175,frame=3)

  

;
  dum=widget_base(base,/row,xsize=2)	

;===========================================================================================================
  B=widget_base(base,/column)
  Btb=widget_tab(B,value='SETUP')
  ;transport coefficient setup
  
  
  B1=widget_base(Btb,title='TRANSPORT SETUP',/column,group_leader=base,/frame)

  B1z=widget_base(B1,/row,frame=3)
  dum=widget_label(B1z,value='Nr (for Chi and V): ')
  coef_Nr=widget_text(B1z,xsize=5,/edit)
  dum=widget_label(B1z,value='Nt (for Chi and V): ')
  coef_Nt=widget_text(B1z,xsize=5,/edit)
    

  B1a=widget_base(B1,/column,xsize=500,frame=3)
  dum=widget_label(B1a,value='DIFFUSIVITY MODEL: chi(x,t) = g(x)*h(t) ',/align_left)
  B1a1=widget_base(B1a,xsize=500,/row)
  dum=widget_label(B1a1,value='Spatial Model g(x): ')
  ;values = ['0th: g(x) = a0','linear: g(x) = a0+a1*x','quadratic: g(x) = a0+a1*x+a2*x^2',$
  ;          'cubic: g(x) = a0+a1*x+a2*x^2+a3*x^3','quartic: g(x) = a0+a1*x+a2*x^2+a3*x^3+a4*x^4']
  values = [ 'poly: g(x) = c0+(a0+a1*(x-x0)+a2*(x-x0)^2+a3*(x-x0)^3)*h(x-x0)' ,$
             'gaussian: g(x) = a0*exp(-(x-a1)^2/(2*a2^2))+a3+a4*x' ,$
             'to be added',$
             'to be added'$
             ]
  
  CHIGmodel = widget_droplist(B1a1,value=values)

  B1a2 = widget_base(B1a,/row,xsize=500)
  dum = widget_label(B1a2,value='c0:')
  CHIc0 = widget_text(B1a2,xsize=5,/edit)
  dum = widget_label(B1a2,value='a0:')
  CHIa0 = widget_text(B1a2,xsize=5,/edit)
  dum = widget_label(B1a2,value='a1:')
  CHIa1 = widget_text(B1a2,xsize=5,/edit)
  dum = widget_label(B1a2,value='a2:')
  CHIa2 = widget_text(B1a2,xsize=5,/edit)
  dum = widget_label(B1a2,value='a3:')
  CHIa3 = widget_text(B1a2,xsize=5,/edit)
  dum = widget_label(B1a2,value='x0:')
  CHIx0 = widget_text(B1a2,xsize=5,/edit)

  B1a3=widget_base(B1a,xsize=500,/row)
  dum=widget_label(B1a3,value='Temporal Model h(t): ')
  values = ['0th: h(t) = c*H(t-t0)+b0','linear: h(t) = c*H(t-t0))+b0+b1*(t-t1)',$
            'quadratic: h(t) = c*H(t-t0))+b0+b1*(t-t1)+b2*(t-t1)^2',$
            'cubic: h(t) = c*H(t-t0))+b0+b1*(t-t1)+b2*(t-t1)^2+b3*(t-t1)^3']
  CHIHmodel = widget_droplist(B1a3,value=values)

  B1a4 = widget_base(B1a,/row,xsize=500)
  dum = widget_label(B1a4,value='c:')
  CHIc = widget_text(B1a4,xsize=5,/edit)
  dum = widget_label(B1a4,value='t0:')
  CHIt0 = widget_text(B1a4,xsize=5,/edit)
  dum = widget_label(B1a4,value='b0:')
  CHIb0 = widget_text(B1a4,xsize=5,/edit)
  dum = widget_label(B1a4,value='b1:')
  CHIb1 = widget_text(B1a4,xsize=5,/edit)
  dum = widget_label(B1a4,value='b2:')
  CHIb2 = widget_text(B1a4,xsize=5,/edit)
  dum = widget_label(B1a4,value='b3:')
  CHIb3 = widget_text(B1a4,xsize=5,/edit)

  B1a5=widget_base(B1a,/row,xsize=500)
  B1a5tb=widget_tab(B1a5,value='DIFFUSION')
  B1a5tbs=widget_base(B1a5tb,title='spatial')
  draw_CHI_spatial=widget_draw(B1a5tbs,xsize=500,ysize=250,frame=1)
  B1a5tbt=widget_base(B1a5tb,title='temporal')
  draw_CHI_temporal=widget_draw(B1a5tbt,xsize=500,ysize=250,frame=1)
  
  ;B1a6=widget_base(B1a,/row)
  ;dum=widget_label(B1a6,value='Time: ')
  ;CHI_time=widget_text(B1a6,xsize=5,/edit)
  ;dum=widget_label(B1a6,value='[Sec] ')
  ;CHI_time_slider=widget_slider(B1a6,xsize=350,min=0,max=2000,value=1000,/drag,/suppress)

;;;;;;;;;;;;;;;;;;;;;;;;
  dum=widget_base(base,/row,ysize=10)	

  B1b=widget_base(B1,/column,xsize=500,frame=3)
  dum=widget_label(B1b,value='CONVECTIVITY MODEL: V(x,t) = g(x)*h(t) ',/align_left)

  B1b1=widget_base(B1b,xsize=500,/row)
  dum=widget_label(B1b1,value='Spatial Model g(x): ')
  ;values = ['0th: g(x) = a0','linear: g(x) = a0+a1*x','quadratic: g(x) = a0+a1*x+a2*x^2',$
  ;          'cubic: g(x) = a0+a1*x+a2*x^2+a3*x^3','quartic: g(x) = a0+a1*x+a2*x^2+a3*x^3+a4*x^4']
  values = [ 'poly: g(x) = c0+(a0+a1*(x-x0)+a2*(x-x0)^2+a3*(x-x0)^3)*h(x-x0)' ,$
             'gaussian: g(x) = a0*exp(-(x-a1)^2/(2*a2^2))+a3+a4*x' ,$
             'to be added',$
             'to be added'$
             ]
  VGmodel = widget_droplist(B1b1,value=values)

  B1b2 = widget_base(B1b,/row,xsize=500)
  dum = widget_label(B1b2,value='c0:')
  Vc0 = widget_text(B1b2,xsize=5,/edit)
  dum = widget_label(B1b2,value='a0:')
  Va0 = widget_text(B1b2,xsize=5,/edit)
  dum = widget_label(B1b2,value='a1:')
  Va1 = widget_text(B1b2,xsize=5,/edit)
  dum = widget_label(B1b2,value='a2:')
  Va2 = widget_text(B1b2,xsize=5,/edit)
  dum = widget_label(B1b2,value='a3:')
  Va3 = widget_text(B1b2,xsize=5,/edit)
  dum = widget_label(B1b2,value='x0:')
  Vx0 = widget_text(B1b2,xsize=5,/edit)

  B1b3=widget_base(B1b,xsize=500,/row)
  dum=widget_label(B1b3,value='Temporal Model h(t): ')
  values = ['0th: h(t) = c*H(t-t0)+b0','linear: h(t) = c*H(t-t0))+b0+b1*(t-t1)',$
            'quadratic: h(t) = c*H(t-t0))+b0+b1*(t-t1)+b2*(t-t1)^2',$
            'cubic: h(t) = c*H(t-t0))+b0+b1*(t-t1)+b2*(t-t1)^2+b3*(t-t1)^3']
  VHmodel = widget_droplist(B1b3,value=values)

  B1b4 = widget_base(B1b,/row,xsize=500)
  dum = widget_label(B1b4,value='c:')
  Vc = widget_text(B1b4,xsize=5,/edit)
  dum = widget_label(B1b4,value='t0:')
  Vt0 = widget_text(B1b4,xsize=5,/edit)
  dum = widget_label(B1b4,value='b0:')
  Vb0 = widget_text(B1b4,xsize=5,/edit)
  dum = widget_label(B1b4,value='b1:')
  Vb1 = widget_text(B1b4,xsize=5,/edit)
  dum = widget_label(B1b4,value='b2:')
  Vb2 = widget_text(B1b4,xsize=5,/edit)
  dum = widget_label(B1b4,value='b3:')
  Vb3 = widget_text(B1b4,xsize=5,/edit)


  B1b5=widget_base(B1b,/row,xsize=500)
  B1b5tb=widget_tab(B1b5,value='CONVECTION')
  B1b5tbs=widget_base(B1b5tb,title='spatial')
  draw_V_spatial=widget_draw(B1b5tbs,xsize=500,ysize=250,frame=1)
  B1b5tbt=widget_base(B1b5tb,title='temporal')
  draw_V_temporal=widget_draw(B1b5tbt,xsize=500,ysize=250,frame=1)
  
  ;B1b6=widget_base(B1b,/row)
  ;dum=widget_label(B1b6,value='Time: ')
  ;V_time=widget_text(B1b6,xsize=5,/edit)
  ;dum=widget_label(B1b6,value='[Sec] ')
  ;V_time_slider=widget_slider(B1b6,xsize=350,min=0,max=2000,value=1000,/drag,/suppress)

  ;==================
  ;strahl density display
  B2=widget_base(Btb,title='STRAHL RESULT',/column,group_leader=base,/frame)

  B2a=widget_base(B2,/col,frame=3)
  dum=widget_label(B2a,value='density(r,t) contour:',/align_left)
  draw_contour_density=widget_draw(B2a,xsize=500,ysize=200)

  ;density v.s. time
  B2b=widget_tab(B2,value='density(t) at various radius')
  B2b1=widget_base(B2b,/col,title='abolute density',frame=3)
  draw_time_density_abs=widget_draw(B2b1,xsize=500,ysize=200)
  B2b2=widget_base(B2b,/col,title='normalized density',frame=3)
  draw_time_density_norm=widget_draw(B2b2,xsize=500,ysize=200)


  ;density v.s. radius
  B2c=widget_tab(B2,value='density(r) at different time')
  B2c1=widget_base(B2c,/col,title='abolute density',frame=3)
  draw_radius_density_abs=widget_draw(B2c1,xsize=500,ysize=200)
  B2c2=widget_base(B2c,/col,title='nomalized density',frame=3)
  draw_radius_density_norm=widget_draw(B2c2,xsize=500,ysize=200)

  ;charge state density
  B2d=widget_base(B2,/col,frame=3)
  draw_radius_chargestate_density=widget_draw(B2d,xsize=500,ysize=200)
  B2d1=widget_base(B2d,/row)
  dum=widget_label(B2d1,value='Time: ')
  density_time=widget_text(B2d1,xsize=5,/edit)
  dum=widget_label(B2d1,value='[Sec] ')
  density_time_slider=widget_slider(B2d1,xsize=350,min=0,max=2000,value=1000,/drag,/suppress)

  ;==================
  ;strahl emissivity (w line)
  B3=widget_base(Btb,title='W emissivity',/column,group_leader=base,/frame)

  B3a=widget_base(B3,/col,frame=3)
  dum=widget_label(B3a,value='emissivity(r,t) contour:',/align_left)
  draw_contour_emissivity_sim_w=widget_draw(B3a,xsize=500,ysize=200)

  ;emissivity v.s. time
  B3b=widget_tab(B3,value='emissivity(t) at various radius')
  B3b1=widget_base(B3b,/col,title='absolute emissivity',frame=3)
  draw_time_emissivity_sim_w_abs=widget_draw(B3b1,xsize=500,ysize=200)
  B3b2=widget_base(B3b,/col,title='normalized emissivity',frame=3)
  draw_time_emissivity_sim_w_norm=widget_draw(B3b2,xsize=500,ysize=200)

  ;emissivity v.s. radius
  B3c=widget_tab(B3,value='emissivity(r) at different time')
  B3c1=widget_base(B3c,/col,title='absolute emissivity',frame=3)
  draw_radius_emissivity_sim_w_abs=widget_draw(B3c1,xsize=500,ysize=200)
  B3c2=widget_base(B3c,/col,title='normalized emissivity',frame=3)
  draw_radius_emissivity_sim_w_norm=widget_draw(B3c2,xsize=500,ysize=200)

  ;dum=widget_label(B3c,value='emissivity(r) at slider time:',/align_left)
  ;draw_radius_chargestate_emissivity_sim_w=widget_draw(B3c,xsize=500,ysize=200)
  ;B3c1=widget_base(B3c,/row)
  ;dum=widget_label(B3c1,value='Time: ')
  ;emissivity_time_sim_w=widget_text(B3c1,xsize=5,/edit)
  ;dum=widget_label(B3c1,value='[Sec] ')
  ;emissivity_time_slider_sim_w=widget_slider(B3c,xsize=350,min=0,max=2000,value=1000,/drag,/suppress)
    ;==================
  ;strahl emissivity (z line)
  B4=widget_base(Btb,title='Z emissivity',/column,group_leader=base,/frame)

  B4a=widget_base(B4,/col,frame=3)
  dum=widget_label(B4a,value='emissivity(r,t) contour:',/align_left)
  draw_contour_emissivity_sim_z=widget_draw(B4a,xsize=500,ysize=200)

  ;emissivity v.s. time
  B4b=widget_tab(B4,value='emissivity(t) at various radius')
  B4b1=widget_base(B4b,/col,title='absolute emissivity',frame=3)
  draw_time_emissivity_sim_z_abs=widget_draw(B4b1,xsize=500,ysize=200)
  B4b2=widget_base(B4b,/col,title='normalized emissivity',frame=3)
  draw_time_emissivity_sim_z_norm=widget_draw(B4b2,xsize=500,ysize=200)

  ;emissivity v.s. radius
  B4c=widget_tab(B4,value='emissivity(r) at different time')
  B4c1=widget_base(B4c,/col,title='absolute emissivity',frame=3)
  draw_radius_emissivity_sim_z_abs=widget_draw(B4c1,xsize=500,ysize=200)
  B4c2=widget_base(B4c,/col,title='normalized emissivity',frame=3)
  draw_radius_emissivity_sim_z_norm=widget_draw(B4c2,xsize=500,ysize=200)
  ;dum=widget_label(B4c,value='emissivity(r) at slider time:',/align_left)
  ;draw_radius_chargestate_emissivity_sim_z=widget_draw(B4c,xsize=500,ysize=200)
  ;B4c1=widget_base(B4c,/row)
  ;dum=widget_label(B4c1,value='Time: ')
  ;emissivity_time_sim_z=widget_text(B4c1,xsize=5,/edit)
  ;dum=widget_label(B4c1,value='[Sec] ')
  ;emissivity_time_slider_sim_z=widget_slider(B4c,xsize=350,min=0,max=2000,value=1000,/drag,/suppress)


 ;==================
  ;strahl brightness (w line)
  B5=widget_base(Btb,title='W brightness',/column,group_leader=base,/frame)

  B5a=widget_base(B5,/col,frame=3)
  dum=widget_label(B5a,value='Brightness(r,t) contour:',/align_left)
  draw_contour_brightness_sim_w=widget_draw(B5a,xsize=500,ysize=200)

  ;brightness v.s. time
  ;B5b=widget_base(B5,/col,frame=3)
  ;dum=widget_label(B5b,value='Brightness(t) at various radius:',/align_left)
  ;draw_time_brightness_sim_w=widget_draw(B5b,xsize=500,ysize=200)
  B5b=widget_tab(B5,value='Brightness(t) at various radius')
  B5b1=widget_base(B5b,/col,title='absolute Brightness',frame=3)
  draw_time_brightness_sim_w_abs=widget_draw(B5b1,xsize=500,ysize=200)
  B5b2=widget_base(B5b,/col,title='normalized Brightness',frame=3)
  draw_time_brightness_sim_w_norm=widget_draw(B5b2,xsize=500,ysize=200)

  ;brightness v.s. radius
  ;B5c=widget_base(B5,/col,frame=3)
  ;dum=widget_label(B5c,value='brightness(r) at different time:',/align_left)
  ;draw_radius_brightness_sim_w=widget_draw(B5c,xsize=500,ysize=200)
  B5c=widget_tab(B5,value='brightness(r) at different time:')
  B5c1=widget_base(B5c,/col,title='absolute Brightness',frame=3)
  draw_radius_brightness_sim_w_abs=widget_draw(B5c1,xsize=500,ysize=200)
  B5c2=widget_base(B5c,/col,title='normalized Brightness',frame=3)
  draw_radius_brightness_sim_w_norm=widget_draw(B5c2,xsize=500,ysize=200)


  ;dum=widget_label(B5c,value='brightness(r) at slider time: ',/align_left)
  ;draw_radius_chargestate_brightness_sim_w=widget_draw(B5c,xsize=500,ysize=200)
  ;B5c1=widget_base(B5c,/row)
  ;dum=widget_label(B5c1,value='Time: ')
  ;brightness_time_sim_w=widget_text(B5c1,xsize=5,/edit)
  ;dum=widget_label(B5c1,value='[Sec] ')
  ;brightness_time_slider_sim_w=widget_slider(B5c,xsize=350,min=0,max=2000,value=1000,/drag,/suppress)

 ;==================
  ;strahl brightness (Z line)
  B6=widget_base(Btb,title='Z brightness',/column,group_leader=base,/frame)

  B6a=widget_base(B6,/col,frame=3)
  dum=widget_label(B6a,value='Brightness(r,t) contour:',/align_left)
  draw_contour_brightness_sim_z=widget_draw(B6a,xsize=500,ysize=200)

  ;brightness v.s. time
  ;B6b=widget_base(B6,/col,frame=3)
  ;dum=widget_label(B6b,value='Brightness(t) at various radius:',/align_left)
  ;draw_time_brightness_sim_z=widget_draw(B6b,xsize=500,ysize=200)
  B6b=widget_tab(B6,value='Brightness(t) at various radius')
  B6b1=widget_base(B6b,/col,title='absolute Brightness',frame=3)
  draw_time_brightness_sim_z_abs=widget_draw(B6b1,xsize=500,ysize=200)
  B6b2=widget_base(B6b,/col,title='normalized Brightness',frame=3)
  draw_time_brightness_sim_z_norm=widget_draw(B6b2,xsize=500,ysize=200)

  ;brightness v.s. radius
  ;B6c=widget_base(B6,/col,frame=3)
  ;dum=widget_label(B6c,value='brightness(r) at different time:',/align_left)
  ;draw_radius_brightness_sim_z=widget_draw(B6c,xsize=500,ysize=200)
  B6c=widget_tab(B6,value='brightness(r) at different time:')
  B6c1=widget_base(B6c,/col,title='absolute Brightness',frame=3)
  draw_radius_brightness_sim_z_abs=widget_draw(B6c1,xsize=500,ysize=200)
  B6c2=widget_base(B6c,/col,title='normalized Brightness',frame=3)
  draw_radius_brightness_sim_z_norm=widget_draw(B6c2,xsize=500,ysize=200)


  ;dum=widget_label(B6c,value='brightness(r) at slider time:',/align_left)
  ;draw_radius_chargestate_brightness_sim_z=widget_draw(B6c,xsize=500,ysize=200)
  ;B6c1=widget_base(B6c,/row)
  ;dum=widget_label(B6c1,value='Time: ')
  ;brightness_time_sim_z=widget_text(B6c1,xsize=5,/edit)
  ;dum=widget_label(B6c1,value='[Sec] ')
  ;brightness_time_slider_sim_z=widget_slider(B6c,xsize=350,min=0,max=2000,value=1000,/drag,/suppress)

;===========================================================================================================
;==================
;HiReX-St emissivity 

  C=widget_base(base,/column)
  Ctb=widget_tab(C,value='HiReX Sr')

;==============
;w line
  C1=widget_base(Ctb,title='W line emissivity',/column,group_leader=base,/frame)
  C1a=widget_base(C1,/col,frame=3)
  dum=widget_label(C1a,value='W emissivity(r,t) contour:',/align_left)
  draw_contour_emissivity_exp_w=widget_draw(C1a,xsize=500,ysize=200)

  ;emissivity v.s. time
  ;C1b=widget_base(C1,/col,frame=3)
  ;dum=widget_label(C1b,value='emissivity(t) at various radius:',/align_left)
  ;draw_time_emissivity_exp_w=widget_draw(C1b,xsize=500,ysize=200)
  C1b=widget_tab(C1,value='emissivity(t) at various radius')
  C1b1=widget_base(C1b,/col,title='absolute emissivity',frame=3)
  draw_time_emissivity_exp_w_abs=widget_draw(C1b1,xsize=500,ysize=200)
  C1b2=widget_base(C1b,/col,title='normalized emissivity',frame=3)
  draw_time_emissivity_exp_w_norm=widget_draw(C1b2,xsize=500,ysize=200)
  ;*****************


  ;emissivity v.s. radius
  ;C1c=widget_base(C1,/col,frame=3)
  ;dum=widget_label(C1c,value='emissivity(r) at slider time:',/align_left)
  ;draw_radius_emissivity_exp_w=widget_draw(C1c,xsize=500,ysize=200)
  C1c=widget_tab(C1,value='emissivity(r) at different time')
  C1c1=widget_base(C1c,/col,title='absolute emissivity',frame=3)
  draw_radius_emissivity_exp_w_abs=widget_draw(C1c1,xsize=500,ysize=200)
  C1c2=widget_base(C1c,/col,title='normalized emissivity',frame=3)
  draw_radius_emissivity_exp_w_norm=widget_draw(C1c2,xsize=500,ysize=200)

  ;dum=widget_label(C1c,value='emissivity(r) for various charge states:',/align_left)
  ;draw_radius_chargestate_emissivity_w=widget_draw(C1c,xsize=500,ysize=200)
  ;C1c1=widget_base(C1c,/row)
  ;dum=widget_label(C1c1,value='Time: ')
  ;emissivity_time_exp-w=widget_text(C1c1,xsize=5,/edit)
  ;dum=widget_label(C1c1,value='[Sec] ')
  ;emissivity_time_slider_exp_w=widget_slider(C1c,xsize=350,min=0,max=2000,value=1000,/drag,/suppress)

;==============
;z line
  C2=widget_base(Ctb,title='Z line emissivity',/column,group_leader=base,/frame)
  C2a=widget_base(C2,/col,frame=3)
  dum=widget_label(C2a,value='W emissivity(r,t) contour:',/align_left)
  draw_contour_emissivity_exp_z=widget_draw(C2a,xsize=500,ysize=200)

  ;emissivity v.s. time
  ;C2b=widget_base(C2,/col,frame=3)
  ;dum=widget_label(C2b,value='emissivity(t) at various radius:',/align_left)
  ;draw_time_emissivity_exp_z=widget_draw(C2b,xsize=500,ysize=200)
  C2b=widget_tab(C2,value='emissivity(t) at various radius')
  C2b1=widget_base(C2b,/col,title='absolute emissivity',frame=3)
  draw_time_emissivity_exp_z_abs=widget_draw(C2b1,xsize=500,ysize=200)
  C2b2=widget_base(C2b,/col,title='normalized emissivity',frame=3)
  draw_time_emissivity_exp_z_norm=widget_draw(C2b2,xsize=500,ysize=200)

  ;emissivity v.s. radius
  ;C2c=widget_base(C2,/col,frame=3)
  ;dum=widget_label(C2c,value='emissivity(r) at slider time:',/align_left)
  ;draw_radius_emissivity_exp_z=widget_draw(C2c,xsize=500,ysize=200)
  C2c=widget_tab(C2,value='emissivity(r) at different time')
  C2c1=widget_base(C2c,/col,title='absolute emissivity',frame=3)
  draw_radius_emissivity_exp_z_abs=widget_draw(C2c1,xsize=500,ysize=200)
  C2c2=widget_base(C2c,/col,title='normalized emissivity',frame=3)
  draw_radius_emissivity_exp_z_norm=widget_draw(C2c2,xsize=500,ysize=200)


  ;dum=widget_label(C2c,value='emissivity(r) for various charge states:',/align_left)
  ;draw_radius_chargestate_emissivity_z=widget_draw(C2c,xsize=500,ysize=200)
  ;C2c1=widget_base(C2c,/row)
  ;dum=widget_label(C2c1,value='Time: ')
  ;emissivity_time_exp_z=widget_text(C2c1,xsize=5,/edit)
  ;dum=widget_label(C2c1,value='[Sec] ')
  ;emissivity_time_slider_exp_z=widget_slider(C2c,xsize=350,min=0,max=2000,value=1000,/drag,/suppress)
;===========================================================================================================
 ;HiReX-St brightness 
 ;==============
;w line
  C3=widget_base(Ctb,title='W line brightness',/column,group_leader=base,/frame)
  C3a=widget_base(C3,/col,frame=3)
  dum=widget_label(C3a,value='W brightness(r,t) contour:',/align_left)
  draw_contour_brightness_exp_w=widget_draw(C3a,xsize=500,ysize=200)

  ;brightness v.s. time
  ;C3b=widget_base(C3,/col,frame=3)
  ;dum=widget_label(C3b,value='brightness(t) at various radius:',/align_left)
  ;draw_time_brightness_exp_w=widget_draw(C3b,xsize=500,ysize=200)
  C3b=widget_tab(C3,value='brightness(t) at various radius')
  C3b1=widget_base(C3b,/col,title='absolute brightness',frame=3)
  draw_time_brightness_exp_w_abs=widget_draw(C3b1,xsize=500,ysize=200)
  C3b2=widget_base(C3b,/col,title='normalized brightness',frame=3)
  draw_time_brightness_exp_w_norm=widget_draw(C3b2,xsize=500,ysize=200)

  ;brightness v.s. radius
  ;C3c=widget_base(C3,/col,frame=3)
  ;dum=widget_label(C3c,value='brightness(r) at slider time:',/align_left)
  ;draw_radius_brightness_exp_w=widget_draw(C3c,xsize=500,ysize=200)
  C3c=widget_tab(C3,value='brightness(r) at different time')
  C3c1=widget_base(C3c,/col,title='absolute brightness',frame=3)
  draw_radius_brightness_exp_w_abs=widget_draw(C3c1,xsize=500,ysize=200)
  C3c2=widget_base(C3c,/col,title='normalized brightness',frame=3)
  draw_radius_brightness_exp_w_norm=widget_draw(C3c2,xsize=500,ysize=200)
  ;dum=widget_label(C3c,value='brightness(r) for various charge states:',/align_left)
  ;draw_radius_chargestate_brightness_w=widget_draw(C3c,xsize=500,ysize=200)
  ;C3c1=widget_base(C3c,/row)
  ;dum=widget_label(C3c1,value='Time: ')
  ;brightness_time_exp_w=widget_text(C3c1,xsize=5,/edit)
  ;dum=widget_label(C3c1,value='[Sec] ')
  ;brightness_time_slider_exp_w=widget_slider(C3c,xsize=350,min=0,max=2000,value=1000,/drag,/suppress)

;==============
;z line
  C4=widget_base(Ctb,title='Z line brightness',/column,group_leader=base,/frame)
  C4a=widget_base(C4,/col,frame=3)
  dum=widget_label(C4a,value='W brightness(r,t) contour:',/align_left)
  draw_contour_brightness_exp_z=widget_draw(C4a,xsize=500,ysize=200)

  ;brightness v.s. time
  ;C4b=widget_base(C4,/col,frame=3)
  ;dum=widget_label(C4b,value='brightness(t) at various radius:',/align_left)
  ;draw_time_brightness_exp_z=widget_draw(C4b,xsize=500,ysize=200)
  C4b=widget_tab(C4,value='brightness(t) at various radius')
  C4b1=widget_base(C4b,/col,title='absolute brightness',frame=3)
  draw_time_brightness_exp_z_abs=widget_draw(C4b1,xsize=500,ysize=200)
  C4b2=widget_base(C4b,/col,title='normalized brightness',frame=3)
  draw_time_brightness_exp_z_norm=widget_draw(C4b2,xsize=500,ysize=200)

  ;brightness v.s. radius
  ;C4c=widget_base(C4,/col,frame=3)
  ;dum=widget_label(C4c,value='brightness(r) at slider time:',/align_left)
  ;draw_radius_brightness_exp_z=widget_draw(C4c,xsize=500,ysize=200)
  C4c=widget_tab(C4,value='brightness(r) at different time')
  C4c1=widget_base(C4c,/col,title='absolute brightness',frame=3)
  draw_radius_brightness_exp_z_abs=widget_draw(C4c1,xsize=500,ysize=200)
  C4c2=widget_base(C4c,/col,title='normalized brightness',frame=3)
  draw_radius_brightness_exp_z_norm=widget_draw(C4c2,xsize=500,ysize=200)
  ;dum=widget_label(C4c,value='brightness(r) for various charge states:',/align_left)
  ;draw_radius_chargestate_brightness_z=widget_draw(C4c,xsize=500,ysize=200)
  ;C4c1=widget_base(C4c,/row)
  ;dum=widget_label(C4c1,value='Time: ')
  ;brightness_time_exp_z=widget_text(C4c1,xsize=5,/edit)
  ;dum=widget_label(C4c1,value='[Sec] ')
  ;brightness_time_slider_exp_z=widget_slider(C4c,xsize=350,min=0,max=2000,value=1000,/drag,/suppress)
;===========================================================================================================
  ;
  id={base:base,$
      shotid:shotid,t1id:t1id,t2id:t2id,Ntid:Ntid,$
      Zid:Zid,Zplotid:Zplotid,$
      fzid:fzid,fztimeid:fztimeid,$
      fits:fits,qfits:qfits,$
      nogrid:nogrid,nopp:nopp,noparam:noparam,nostrahl:nostrahl,quiet:quiet,$
      load:load, load_hirex:load_hirex,run:run,exp_sim:exp_sim,plot:plot,stop:stop,quit:quit,$
      savepath:savepath,savebutton:savebutton,loadpath:loadpath,loadbutton:loadbutton,$
      strahltmin:strahltmin,strahltmax:strahltmax,strahlrmin:strahlrmin,strahlrmax:strahlrmax,$
      hirextmin:hirextmin,hirextmax:hirextmax,hirexrmin:hirexrmin,hirexrmax:hirexrmax,$
      draw_binning:draw_binning,draw_source:draw_source,$
      coef_Nr:coef_Nr,coef_Nt:coef_Nt,$
      CHIGmodel:CHIGmodel,CHIc0:CHIc0,CHIa0:CHIa0,CHIa1:CHIa1,CHIa2:CHIa2,CHIa3:CHIa3,CHIx0:CHIx0,$
      CHIHmodel:CHIHmodel,CHIc:CHIc,CHit0:CHit0,CHIb0:CHIb0,CHIb1:CHIb1,CHIb2:CHIb2,CHIb3:CHIb3,$
      draw_CHI_spatial:draw_CHI_spatial,draw_CHI_temporal:draw_CHI_temporal,$
      ;CHI_time:CHI_time,CHI_time_slider:CHI_time_slider,$
      VGmodel:VGmodel,Vc0:Vc0,Va0:Va0,Va1:Va1,Va2:Va2,Va3:Va3,Vx0:Vx0,$
      VHmodel:VHmodel,Vc:Vc,Vt0:Vt0,Vb0:Vb0,Vb1:Vb1,Vb2:Vb2,Vb3:Vb3,$
      draw_V_spatial:draw_V_spatial,draw_V_temporal:draw_V_temporal,$
      ;V_time:V_time,V_time_slider:V_time_slider,$
      draw_contour_density:draw_contour_density,$
      draw_time_density_abs:draw_time_density_abs,$
      draw_time_density_norm:draw_time_density_norm,$
      draw_radius_density_abs:draw_radius_density_abs,$
      draw_radius_density_norm:draw_radius_density_norm,$
      draw_radius_chargestate_density:draw_radius_chargestate_density,$
      density_time:density_time,density_time_slider:density_time_slider,$
      ;
      draw_contour_emissivity_sim_w:draw_contour_emissivity_sim_w,$
      draw_time_emissivity_sim_w_abs:draw_time_emissivity_sim_w_abs,$
      draw_time_emissivity_sim_w_norm:draw_time_emissivity_sim_w_norm,$
      draw_radius_emissivity_sim_w_abs:draw_radius_emissivity_sim_w_abs,$
      draw_radius_emissivity_sim_w_norm:draw_radius_emissivity_sim_w_norm,$
      ;draw_radius_chargestate_emissivity:draw_radius_chargestate_emissivity,$
      ;emissivity_time:emissivity_time,emissivity_time_slider:emissivity_time_slider,$
      ;
      draw_contour_emissivity_sim_z:draw_contour_emissivity_sim_z,$
      draw_time_emissivity_sim_z_abs:draw_time_emissivity_sim_z_abs,$
      draw_time_emissivity_sim_z_norm:draw_time_emissivity_sim_z_norm,$
      draw_radius_emissivity_sim_z_abs:draw_radius_emissivity_sim_z_abs,$
      draw_radius_emissivity_sim_z_norm:draw_radius_emissivity_sim_z_norm,$
      ;
      draw_contour_brightness_sim_w:draw_contour_brightness_sim_w,$
      draw_time_brightness_sim_w_abs:draw_time_brightness_sim_w_abs,$
      draw_time_brightness_sim_w_norm:draw_time_brightness_sim_w_norm,$
      draw_radius_brightness_sim_w_abs:draw_radius_brightness_sim_w_abs,$
      draw_radius_brightness_sim_w_norm:draw_radius_brightness_sim_w_norm,$
      ;draw_radius_chargestate_brightness:draw_radius_chargestate_brightness,$
      ;brightness_time:brightness_time,brightness_time_slider:brightness_time_slider,$
      ;
      draw_contour_brightness_sim_z:draw_contour_brightness_sim_z,$
      draw_time_brightness_sim_z_abs:draw_time_brightness_sim_z_abs,$
      draw_time_brightness_sim_z_norm:draw_time_brightness_sim_z_norm,$
      draw_radius_brightness_sim_z_abs:draw_radius_brightness_sim_z_abs,$
      draw_radius_brightness_sim_z_norm:draw_radius_brightness_sim_z_norm,$
      ;
      draw_contour_emissivity_exp_w:draw_contour_emissivity_exp_w,$
      draw_time_emissivity_exp_w_abs:draw_time_emissivity_exp_w_abs,$
      draw_time_emissivity_exp_w_norm:draw_time_emissivity_exp_w_norm,$
      draw_radius_emissivity_exp_w_abs:draw_radius_emissivity_exp_w_abs,$
      draw_radius_emissivity_exp_w_norm:draw_radius_emissivity_exp_w_norm,$
      ;draw_radius_chargestate_emissivity_w:draw_radius_chargestate_emissivity_w,$
      ;emissivity_time_w:emissivity_time_w,emissivity_time_slider_w:emissivity_time_slider_w,$
      ;
      draw_contour_emissivity_exp_z:draw_contour_emissivity_exp_z,$
      draw_time_emissivity_exp_z_abs:draw_time_emissivity_exp_z_abs,$
      draw_time_emissivity_exp_z_norm:draw_time_emissivity_exp_z_norm,$
      draw_radius_emissivity_exp_z_abs:draw_radius_emissivity_exp_z_abs,$
      draw_radius_emissivity_exp_z_norm:draw_radius_emissivity_exp_z_norm,$
      ;draw_radius_chargestate_emissivity_z:draw_radius_chargestate_emissivity_z,$
      ;emissivity_time_z:emissivity_time_z,emissivity_time_slider_z:emissivity_time_slider_z,$
      ;
      draw_contour_brightness_exp_w:draw_contour_brightness_exp_w,$
      draw_time_brightness_exp_w_abs:draw_time_brightness_exp_w_abs,$
      draw_time_brightness_exp_w_norm:draw_time_brightness_exp_w_norm,$
      draw_radius_brightness_exp_w_abs:draw_radius_brightness_exp_w_abs,$
      draw_radius_brightness_exp_w_norm:draw_radius_brightness_exp_w_norm,$
      ;draw_radius_chargestate_brightness_w:draw_radius_chargestate_brightness_w,$
      ;brightness_time_w:brightness_time_w,brightness_time_slider_w:brightness_time_slider_w,$
      ;
      draw_contour_brightness_exp_z:draw_contour_brightness_exp_z,$
      draw_time_brightness_exp_z_abs:draw_time_brightness_exp_z_abs,$
      draw_time_brightness_exp_z_norm:draw_time_brightness_exp_z_norm,$
      draw_radius_brightness_exp_z_abs:draw_radius_brightness_exp_z_abs,$
      draw_radius_brightness_exp_z_norm:draw_radius_brightness_exp_z_norm}
      ;draw_radius_chargestate_brightness_z:draw_radius_chargestate_brightness_z,$
      ;brightness_time_z:brightness_time_z,brightness_time_slider_z:brightness_time_slider_z}
  
  if not keyword_set(shot) then shot=1140221016
  if not keyword_set(z) then z=18
  if not keyword_set(Zplot) then Zplot= z gt 2? z-2:0
  if not keyword_set(t1) then t1=0.4
  if not keyword_set(t2) then t2=1.2
  if not keyword_set(Nt) then Nt=1
  time=(t1+t2)/2.
  strahltr=[0.8,1.2]
  strahlrr=[0.0,0.8]
  hirextr=[0.8,1.2]
  hirexrr=[0.0,0.8]
  
  denmax=1e18
  profile={fits:1,qfits:0}
  optional={nogrid:0,nopp:0,noparam:0,nostrahl:0,quiet:0}
  times={density_time:1.0,emissivity_time:time,brightness_time:time}
  
  binplt=[t1,t2,0,1.2]
  sourceplt=[t1,t2,0,denmax]
  chisplt=[0,1,0,5]
  chitplt=[t1,t2,0,5]
  vsplt=[0,1,-5,5]
  vcplt=[t1,t2,-5,5]
  strahl1plt=[0,1,t1,t2]
  strahl2plt=[t1,t2,0,denmax]
  strahl3plt=[0,1,0,denmax]
  strahl4plt=[0,1,0,denmax]
  chiGdroplist = 0
  chiHdroplist = 0
  VGdroplist = 0
  VHdroplist = 0
  chiGparam={c0:0.0,a0:1.0,a1:0.0,a2:0.0,a3:0.0,x0:0.0}
  chiHparam={c:0.0,t0:1.0,b0:1.0,b1:1.0,b2:1.0,b3:1.0}
  VGparam={c0:0.0,a0:15,a1:0.40,a2:0.25,a3:0.0,x0:0.0}
  VHparam={c:1.0,t0:1.0,b0:0.0,b1:0.0,b2:0.0,b3:0.0}
  

  plot={binplt:binplt,sourceplt:sourceplt,chisplt:chisplt,chitplt:chitplt,$
        vsplt:vsplt,vcplt:vcplt,strahl1plt:strahl1plt,strahl2plt:strahl2plt,$
        strahl3plt:strahl3plt,strahl4plt:strahl4plt}
  
  tmpNr=50
  tmpNt=50
  coefx=findgen(tmpNr+1.)/(tmpNr)
  coeft=findgen(tmpNt)/float(tmpNt)*(t2-t1)+t1
  chi=fltarr(tmpNr,tmpNt)
  v=fltarr(tmpNr,tmpNt)

  ;store chi and v information
  coef={coef_Nr:tmpNr,coef_Nt:tmpNt,$
        chiGdroplist:chiGdroplist,chiHdroplist:chiHdroplist,$
        chiGparam:chiGparam,chiHparam:chiHparam,$
        VGdroplist:VGdroplist,VHdroplist:VHdroplist,$
        VGparam:VGparam,VHparam:VHparam,$
        coefx:ptr_new(/allocate),coeft:ptr_new(/allocate),$
        chi:ptr_new(/allocate),v:ptr_new(/allocate)}  
  *coef.coefx=coefx
  *coef.coeft=coeft
  *coef.chi=chi
  *coef.v=v
  
  ;result of STRAHL simulation
  data = ptr_new(/allocate)
  term = ptr_new(/allocate)
  emiss_sim_w = ptr_new(/allocate)
  emiss_sim_z = ptr_new(/allocate)
  bright_sim_w = ptr_new(/allocate)
  bright_sim_z = ptr_new(/allocate)
  emiss_exp_w = ptr_new(/allocate)
  emiss_exp_z = ptr_new(/allocate)
  bright_exp_w = ptr_new(/allocate)
  bright_exp_z = ptr_new(/allocate)
  

  strahl = {data:data,term:term, $
            emiss_w:emiss_sim_w, emiss_z:emiss_sim_z, bright_w:bright_sim_w, bright_z:bright_sim_z,$
            set_dens_abs:[-1,-1,-1,-1],$
            set_dens_norm:[-1,-1,-1,-1],$
            set_emiss_sim_w_abs:[-1,-1,-1],set_emiss_sim_z_abs:[-1,-1,-1],$
            set_emiss_sim_w_norm:[-1,-1,-1],set_emiss_sim_z_norm:[-1,-1,-1],$
            set_bright_sim_w_abs:[-1,-1,-1],set_bright_sim_z_abs:[-1,-1,-1],$
            set_bright_sim_w_norm:[-1,-1,-1],set_bright_sim_z_norm:[-1,-1,-1],$
            tr:strahltr,rr:strahlrr}
  hirex = {emiss_w:emiss_exp_w, emiss_z:emiss_exp_z, bright_w:bright_exp_w, bright_z:bright_exp_z, $
           set_emiss_exp_w_abs:[-1,-1,-1],set_emiss_exp_z_abs:[-1,-1,-1],$
           set_emiss_exp_w_norm:[-1,-1,-1],set_emiss_exp_z_norm:[-1,-1,-1],$
           set_bright_exp_w_abs:[-1,-1,-1],set_bright_exp_z_abs:[-1,-1,-1],$
           set_bright_exp_w_norm:[-1,-1,-1],set_bright_exp_z_norm:[-1,-1,-1],$
           tr:hirextr,rr:hirexrr}
  stat={profile:profile,optional:optional,times:times,plot:plot,coef:coef,strahl:strahl,hirex:hirex}
  
  ;fz_text='0.0,0.0,1.0,1.0,0.0,0.0,1.0,1.0,0.0,0.0'
  ;fztime_text='0.0,0.9,0.9,0.91,0.91,1.0,1.0,1.01,1.01,2.0'
  fz_text='0.0,0.0,1.0,1.0'
  fztime_text='0.0,0.4,0.4,2.0'
  savepa = num2str(shot,1)+'_'
  loadpa = num2str(shot,1)+'_'
  u={id:id,stat:stat,shot:shot,Z:Z,Zplot:Zplot,t1:t1,t2:t2,Nt:Nt,$
     fz_text:fz_text,fz:ptr_new(double(strsplit(fz_text,',',/extract))),$
     savepath:savepa,loadpath:loadpa,$
     fztime_text:fztime_text,fztime:ptr_new(double(strsplit(fztime_text,',',/extract)))}


  wstrahl_load_info,u
  ;


  !except=0
  widget_control,base,/realize

  ;get window information
  widget_control,u.id.draw_contour_density,get_value=win1
  widget_control,u.id.draw_time_density_abs,get_value=win2
  widget_control,u.id.draw_radius_density_abs,get_value=win3
  widget_control,u.id.draw_radius_chargestate_density,get_value=win4
  u.stat.strahl.set_dens_abs=[win1,win2,win3,win4]
  ;get window information
  ;widget_control,u.id.draw_contour_density,get_value=win1
  widget_control,u.id.draw_time_density_norm,get_value=win2
  ;widget_control,u.id.draw_radius_density_norm,get_value=win3
  ;widget_control,u.id.draw_radius_chargestate_density,get_value=win4
  u.stat.strahl.set_dens_norm=[win1,win2,win3,win4]
  ;print,u.stat.strahl.set_dens



  widget_control,u.id.draw_contour_emissivity_sim_w,get_value=win1
  widget_control,u.id.draw_time_emissivity_sim_w_abs,get_value=win2
  widget_control,u.id.draw_radius_emissivity_sim_w_abs,get_value=win3
  u.stat.strahl.set_emiss_sim_w_abs=[win1,win2,win3]
  ;widget_control,u.id.draw_contour_emissivity_sim_w,get_value=win1
  widget_control,u.id.draw_time_emissivity_sim_w_norm,get_value=win2
  widget_control,u.id.draw_radius_emissivity_sim_w_norm,get_value=win3
  u.stat.strahl.set_emiss_sim_w_norm=[win1,win2,win3]
  ;print,u.stat.strahl.set_emiss

  widget_control,u.id.draw_contour_emissivity_sim_z,get_value=win1
  widget_control,u.id.draw_time_emissivity_sim_z_abs,get_value=win2
  widget_control,u.id.draw_radius_emissivity_sim_z_abs,get_value=win3
  u.stat.strahl.set_emiss_sim_z_abs=[win1,win2,win3]
  ;widget_control,u.id.draw_contour_emissivity_sim_z,get_value=win1
  widget_control,u.id.draw_time_emissivity_sim_z_norm,get_value=win2
  widget_control,u.id.draw_radius_emissivity_sim_z_norm,get_value=win3
  u.stat.strahl.set_emiss_sim_z_norm=[win1,win2,win3]
  ;print,u.stat.strahl.set_emiss
 
 
  widget_control,u.id.draw_contour_brightness_sim_w,get_value=win1
  widget_control,u.id.draw_time_brightness_sim_w_abs,get_value=win2
  widget_control,u.id.draw_radius_brightness_sim_w_abs,get_value=win3
  u.stat.strahl.set_bright_sim_w_abs=[win1,win2,win3]
  ;widget_control,u.id.draw_contour_brightness_sim_w,get_value=win1
  widget_control,u.id.draw_time_brightness_sim_w_norm,get_value=win2
  widget_control,u.id.draw_radius_brightness_sim_w_norm,get_value=win3
  u.stat.strahl.set_bright_sim_w_norm=[win1,win2,win3]
  ;print,u.stat.strahl.set_emiss

  widget_control,u.id.draw_contour_brightness_sim_z,get_value=win1
  widget_control,u.id.draw_time_brightness_sim_z_abs,get_value=win2
  widget_control,u.id.draw_radius_brightness_sim_z_abs,get_value=win3
  u.stat.strahl.set_bright_sim_z_abs=[win1,win2,win3]
  ;widget_control,u.id.draw_contour_brightness_sim_z,get_value=win1
  widget_control,u.id.draw_time_brightness_sim_z_norm,get_value=win2
  widget_control,u.id.draw_radius_brightness_sim_z_norm,get_value=win3
  u.stat.strahl.set_bright_sim_z_norm=[win1,win2,win3]
  ;print,u.stat.strahl.set_emiss


  widget_control,u.id.draw_contour_emissivity_exp_w,get_value=win1
  widget_control,u.id.draw_time_emissivity_exp_w_abs,get_value=win2
  widget_control,u.id.draw_radius_emissivity_exp_w_abs,get_value=win3
  u.stat.hirex.set_emiss_exp_w_abs=[win1,win2,win3]
  ;widget_control,u.id.draw_contour_emissivity_exp_w,get_value=win1
  widget_control,u.id.draw_time_emissivity_exp_w_norm,get_value=win2
  widget_control,u.id.draw_radius_emissivity_exp_w_norm,get_value=win3
  u.stat.hirex.set_emiss_exp_w_norm=[win1,win2,win3]


  widget_control,u.id.draw_contour_emissivity_exp_z,get_value=win1
  widget_control,u.id.draw_time_emissivity_exp_z_abs,get_value=win2
  widget_control,u.id.draw_radius_emissivity_exp_z_abs,get_value=win3
  u.stat.hirex.set_emiss_exp_z_abs=[win1,win2,win3]
  ;widget_control,u.id.draw_contour_emissivity_exp_z,get_value=win1
  widget_control,u.id.draw_time_emissivity_exp_z_norm,get_value=win2
  widget_control,u.id.draw_radius_emissivity_exp_z_norm,get_value=win3
  u.stat.hirex.set_emiss_exp_z_norm=[win1,win2,win3]

  widget_control,u.id.draw_contour_brightness_exp_w,get_value=win1
  widget_control,u.id.draw_time_brightness_exp_w_abs,get_value=win2
  widget_control,u.id.draw_radius_brightness_exp_w_abs,get_value=win3
  u.stat.hirex.set_bright_exp_w_abs=[win1,win2,win3]
  ;widget_control,u.id.draw_contour_brightness_exp_w,get_value=win1
  widget_control,u.id.draw_time_brightness_exp_w_norm,get_value=win2
  widget_control,u.id.draw_radius_brightness_exp_w_norm,get_value=win3
  u.stat.hirex.set_bright_exp_w_norm=[win1,win2,win3]

  widget_control,u.id.draw_contour_brightness_exp_z,get_value=win1
  widget_control,u.id.draw_time_brightness_exp_z_abs,get_value=win2
  widget_control,u.id.draw_radius_brightness_exp_z_abs,get_value=win3
  u.stat.hirex.set_bright_exp_z_abs=[win1,win2,win3]
  ;widget_control,u.id.draw_contour_brightness_exp_z,get_value=win1
  widget_control,u.id.draw_time_brightness_exp_z_norm,get_value=win2
  widget_control,u.id.draw_radius_brightness_exp_z_norm,get_value=win3
  u.stat.hirex.set_bright_exp_z_norm=[win1,win2,win3]


  ;set u value
  widget_control,base,set_uvalue=u

  xmanager,'w_strahl',base

end
