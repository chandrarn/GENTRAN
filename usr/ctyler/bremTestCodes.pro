;NAME:
;     integrateAllWavelengths
;
;PURPOSE:
;     Compares the the wavelength-integrated emissivity of energyEmiss
;     over a wide range of wavelengths to a known overall
;     bremsstrahlung power emission. The numbers should come out
;     close, although not exact due to the gaunt factor.
;
;CALLING SEQUENCE:
;     integrateAllWavelengths,n_e,t_e,z_eff
;
;INPUTS:
;     n_e: An arbitrary electron density to test, in m^-3
;     t_e: An arbitrary electron temperature to test, in keV
;     z_eff: An arbitrary z_eff to test
;
;OUTPUTS:
;     Prints results

PRO integrateAllWavelengths,n_e,t_e,z_eff
  energyEmissivity = fltarr(4000)
  wavelengths = fltarr(n_elements(energyEmissivity))
  FOR i=0,n_elements(energyEmissivity)-1 DO BEGIN
     energyEmissivity[i] = energyEmiss(1,n_e,t_e,i+1.)
     wavelengths[i] = i+1
  ENDFOR
 
  print,'Calculated by energyEmiss integration (W/m^3): ' + string(int_tabulated(wavelengths, energyEmissivity))


                                ;From magnetic fusion energy formulary
                                ;(does not account for guant factor,
                                ;so an exact match is not expected):
  print,'Calculated by expected bremsstrahlung total power loss (W/m^3): ' + string(5.35e-37*z_eff*n_e^2*t_e^(0.5))
END






;NAME:
;     testSeveralEmissEquations
;
;PURPOSE:
;     Plots four integrated emissivity equations for the n_e and t_e
;     data at radial points at a particular time. Two are energy
;     dependent, and two are wavelength dependent. Both integrate
;     identically showing the validity of the equation used in energyEmiss 
;
;CALLING SEQUENCE:
;     testSeveralEmissEquations
;
;INPUTS:
;     shot: provide some shot number
;     t: provide some shot time
;
;OUTPUTS:
;     provides a plot of the four integrated emissivities with respect
;     to radial position.

PRO testSeveralEmissEquations,shot,t
  loadct,5
  cmod_ts,shot,ts
  index = closest(ts.TIME,t)
  temp = reform(ts.TE[index,*])
  dens = reform(ts.N_E[index,*])
  rmaj = reform(ts.R[index,*])
  goodPlaces = WHERE(temp NE 0)
  temp = temp[goodPlaces]
  dens = dens[goodPlaces]
  rmaj = rmaj[goodPlaces]
  zeff=1
  path_filter='/usr/local/cmod/idl/atomic_physics/sxr_filters/be50.dat'
  filter=read_xray_data(path_filter)
  eph=filter.e
  lam=ev2ang(eph)
  npts=n(rmaj)+1		;number of radial points
  h=6.626068e-34                ;m^2 kg/s
  c=2.9979e8                    ;m/s
  e=1.60218e-19                 ;J/eV
  kb=8.617e-5                   ;eV/K 
  lamRev = reverse(lam)
  
  emissInt1=fltarr(npts)                             
  FOR i=0,npts-1 DO BEGIN       ;Calculation from GENIE (taken from Hutchinson book, as of 6/21/2012 has absolute magnitude issues but is qualitatively correct. Scaled in the plot)     (Orange)
     log_gff=0.355*lam^(-0.06)*alog10(lam)+0.3*lam^(-0.066)*alog10(temp[i]*1000/kb*1.0e-6/100.0)+0.0043
     emiss1=double(double(5.03e-14)/double(h)*double(e)*double(4.0)*double(!pi)*double(10^(log_gff))*double(zeff)*double(dens[i])/double(1.0e20)*double(dens[i])/double(1e20)/(double(temp[i]))^(0.5)*exp(-double(eph)/double(1000)/double(temp[i]))) ;watts/m^3/eV from Hutch
     emissInt1[i]=int_tabulated(eph,emiss1) ;W/m^3
  ENDFOR
  
  emissInt2=fltarr(npts)                          
  FOR i=0,npts-1 DO BEGIN       ;Lambda-dependent equation from http://adsabs.harvard.edu/full/1986A&AS...65..511M     (Red)
     T_MK = 11.6*temp[i]
     gff = 10^(0.355 * lamRev^(-0.06)*ALOG10(lamRev) + 0.3*lamRev^(-0.066)*ALOG10(T_MK/100) + 0.0043)
     emiss2=double((double(2.051e-22)*(double(dens[i])*double(1e-6))^2/(double(T_MK)^(0.5)*double(lamRev)^2))*exp(-double(143.9)/(double(T_MK)*double(lamRev)))*double(gff)*double(1e-1))
     emissInt2[i] = int_tabulated(lamRev, emiss2)
  ENDFOR
  
  emissInt3 = fltarr(npts)
  FOR i=0,npts-1 DO BEGIN       ;Energy-dependent equation from http://adsabs.harvard.edu/full/1986A&AS...65..511M     (Blue)
     T_MK = 11.6*temp[i]
     gff = 10^(0.355 * lam^(-0.06)*ALOG10(lam) + 0.3*lam^(-0.066)*ALOG10(T_MK/100) + 0.0043)
     emiss3 = double(double(1.032e-14)*(double(dens[i])*double(1e-6))^2/(double(11.6)*double(temp[i]))^(0.5)*exp(-double(11.6)*double(eph)/double(1000)/(double(temp[i])*double(11.6)))*double(gff)*double(1.602e-16))*1e6
     emissInt3[i] = int_tabulated(eph/1000,emiss3)
  ENDFOR

  emissInt4 = fltarr(npts)
  FOR i=0,npts-1 DO BEGIN       ;Taken from Rev. Sci. Instruments. 78, 023501 (2007)    (White)
     T_MK = 11.6*temp[i]
     gff = 10^(0.355 * lamRev^(-0.06)*ALOG10(lamRev) + 0.3*lamRev^(-0.066)*ALOG10(T_MK/100) + 0.0043)
     emiss4 = double(double(1.89e-28)*double(zeff)*(double(dens[i])*double(1e-10))^2*double(gff) / (double(temp[i])*double(1000))^(0.5) / double(lamRev)^2 * exp(-double(12400)/(double(temp[i])*double(1000)*double(lamRev))) * double(1e10) * double(1e10) * double(1e6) / double(1e12))
     emissInt4[i] = int_tabulated(lamRev,emiss4)
  ENDFOR

  plot,rmaj,emissInt4
  oplot,rmaj,emissInt3,color=50
  oplot,rmaj,emissInt2,color=100
  oplot,rmaj,emissInt1*3.3e-2,color=150
END







;NAME:
;     testMin
;
;PURPOSE:
;     shows the minimized spectrum of both XEUS and Loweus for a
;     certain shot and time
;
;CALLING SEQUENCE:
;     testMin,shot,t,width=width
;
;INPUTS:
;     shot: Provide a shot to look at
;     t: provide a time to look at
;
;OPTIONAL INPUTS:
;     width: the width of each sample for which to take the minimum
;     of. 50 is the default, and seems to work well for both XEUS and
;     Loweus
;
;OUTPUTS:
;     A plot showing the minimized spectrums

PRO testMin,shot,t,width=width
  IF NOT keyword_set(width) THEN width=50
  loadct,5
  mdsopen,'spectroscopy',shot
  lamLoweus=mdsvalue('dim_of(\SPECTROSCOPY::TOP.loweus.SPEC,0)')
  lamXeus  =mdsvalue('dim_of(\SPECTROSCOPY::TOP.XEUS.SPEC,0)')

  times = mdsvalue('dim_of(\SPECTROSCOPY::TOP.XEUS.SPEC,1)')
  
  specLoweus=mdsvalue('\SPECTROSCOPY::TOP.loweus.SPEC')
  specXeus  =mdsvalue('\SPECTROSCOPY::TOP.XEUS.SPEC')

  xeus_pos = [2.561,0.21,0.179,6.5*!pi/180]
  loweus_pos = [2.561,0.21,0.179,6.5*!pi/180]
  genpos_pos_reform,xeus_pos,[0.44,1,0.6,-0.6]
  genpos_pos_reform,loweus_pos,[0.44,1,0.6,-0.6]

  rawLoweus = specLoweus[*,closest(times,t)]
  rawXeus = specXeus[*,closest(times,t)]

  indicies=ExtractMin(rawLoweus,width)
  plot,lamLoweus,rawLoweus,xrange=[0,max([lamLoweus,lamXeus])]
  oplot,MakeSpecFromMin(lamLoweus,indicies),MakeSpecFromMin(rawLoweus,indicies),color=105

  oplot,lamXeus,rawXeus,color=255
  indicies=ExtractMin(rawXeus,width)
  oplot,MakeSpecFromMin(lamXeus,indicies),MakeSpecFromMin(rawXeus,indicies),color=105
END



;NAME:
;     testPlotEmiss
;
;PURPOSE:
;     Plots the emissivity curve from energyEmiss over a given range.
;
;CALLING SEQUENCE:
;     testPlotEmiss,n_e,t_e,range,zeff=zeff
;
;INPUTS:
;     n_e: provide an arbitrary n_e to plot for, 1/m^3
;     t_e: provide an arbitrary t_e to plot for, keV
;     range: Will plot from 0 angstroms to range angstroms
;
;OPTIONAL INPUTS:
;     zeff: if unspecified, assumed to be zero

;OUTPUTS:
;     A plot in the units of energyEmiss, which was W/m^3/ang at the
;     time of this writing

PRO testPlotEmiss,n_e,t_e,range,zeff=zeff
  IF NOT keyword_set(zeff) THEN zeff=1

  lam=fltarr(500)
  FOR i=0,499 DO BEGIN
     lam[i] = double(i)*range/500+0.1
  ENDFOR

  output = dblarr(n_elements(lam))

  FOR i=0,n_elements(lam)-1 DO BEGIN
     output[i] = energyEmiss(zeff,n_e,t_e,lam[i])
  ENDFOR
  plot,lam,output
END
