
;NAME:
;     energyEmiss
;
;PURPOSE:
;     Returns the emissivity expected for some plasma parameters
;     Function gaunt provides energyEmiss with the gaunt factor
;
;CALLING SEQUENCE
;     emissivity = energyEmiss(zeff,n_e,t_e,lambda)
;
;INPUTS:
;     zeff: zeff of the plasma
;     n_e: electron density in 1/m^3
;     t_e: electron temperature in keV
;     lambda: wavelength in angstroms
;
;OUTPUTS:
;     Returns emissivity in W/m^3/A 


FUNCTION gaunt, lambda, T_MK
  ;lambda in angstroms
  ;T in MK (valid for T>0.086KeV)

  ;Taken from http://adsabs.harvard.edu/full/1986A&AS...65..511M
  IF T_MK GT 1 THEN BEGIN
     gff = 10^(0.355 * lambda^(-0.06)*ALOG10(lambda) + 0.3*lambda^(-0.066)*ALOG10(T_MK/100) + 0.0043)
  ENDIF ELSE BEGIN
     gff = 10^(0.48*lambda^(-0.08)*ALOG10(lambda) + (0.133*ALOG10(lambda)-0.2)*ALOG10(T_MK) - 0.538)
  ENDELSE

 
  RETURN,gff
END

FUNCTION energyEmiss,zeff,n_e,t_e,lambda  
  ;lambda in angstroms
  ;t_e in KeV
  ;n_e in cm^-3
  T_MK = t_e*11.6

  ;Taken from Rev. Sci. Instruments. 78, 023501 (2007)
  ;Returns in W m^-3 A^-1
  emiss = double(1.89e-28)*double(zeff)*double(double(n_e)*double(1e-26))^2*double(gaunt(lambda,T_MK)) / (double(t_e)*double(1000))^(0.5) / double(lambda)^2 * (exp(-12400/(t_e*1000*lambda))) * double(1e20) * double(1e20)*1e6
;  if lambda LT 133 then begin
;     RETURN,emiss 
;  endif else begin
;     RETURN,emiss
;  endelse
return, emiss
END







;NAME:
;     emiss
;
;PURPOSE:
;     Returns a line integrated emissivity array for a lambda array
;     based on TS data
;
;CALLING SEQUENCE:
;     emissivity = emiss([lam],shot,t,z_eff,pos
;
;INPUTS:
;     lam: 1D array containing wavelengths to analyze (in angstroms)
;     shot: What shot to use TS data from
;     t: What time index of that shot to use TS data from
;     z_eff: Provide a z_eff to analyze
;     pos: Provide the POS vector that the line integral should be
;       performed over
;
;OUTPUTS:
;     Returns a 1D array of same dimension as lam, containing
;     emissivities in the units returned by energyEmiss, W/m^3/A at
;     the time of this writing



;A helper function that returns the index of the element in inpArr
;nearest in absolute value to valSearch:
FUNCTION closest,inpArr,valSearch
  index = -1
  IF (n_elements(inpArr) EQ 1) THEN BEGIN
     index = 0
  ENDIF ELSE BEGIN
     IF (n_elements(inpArr) GT 1) THEN BEGIN
        offset = ABS(valSearch - inpArr[0])
        index = 0
        FOR i=1,n_elements(inpArr) - 1 DO BEGIN
           IF (ABS(valSearch - inpArr[i]) LT offset) THEN BEGIN
              index = i
              offset = ABS(valSearch - inpArr[i])
           ENDIF
        ENDFOR
     ENDIF
  ENDELSE
  RETURN,index
END


FUNCTION emiss,lam,shot,t,z_eff,pos
  cmod_ts,shot,ts

  N_Eraw = reform(ts.N_E[closest(ts.TIME,t),*])
  TEraw = reform(ts.TE[closest(ts.TIME,t),*])
  goodLocations = where(N_Eraw NE 0)
  closestT = closest(ts.TIME,t)
  RadiusRaw = reform(ts.R[closestT,*])

  IF (goodLocations EQ [-1]) THEN BEGIN
     N_E = [0]
     TE = [0]
     Radius = [0]
  
  ENDIF ELSE BEGIN
  
     N_E = N_Eraw[goodLocations]
     TE = TEraw[goodLocations]
     Radius = RadiusRaw[goodLocations]
  ENDELSE
 

  emiss_arr = fltarr(n_elements(lam))


  ;precalculate efit to greatly accelerate things:
  efit={time:line_gettimes(shot),lcfs:line_getlcfs(shot),axis:line_getaxis(shot),rmid:line_getrmid(shot)}
  

  FOR i=0,n_elements(lam)-1 DO BEGIN
     EmissI=fltarr(n_elements(TE))
     
     FOR z=0,(n_elements(TE)-1) DO BEGIN
        EmissI[z]=energyEmiss(z_eff,N_E[z],TE[z],lam[i])
     ENDFOR

     emiss_arr[i] = genpos_line_br(pos,[EmissI],Radius,[t],shot,t,efit=efit)
  
     print, string(27B) + '[4D' + string(27B) + '[K', format='(A, $)'
     print, string(int(float(i)/n_elements(lam)*100)), format='(I3, "%", $)'
 
  ENDFOR
  print,""
  return,emiss_arr
END





;NAME:
;     calibrate
;
;PURPOSE:
;     Relatively calibrates a camera by comparing calculated emissivity to the
;     observed emissivity spectrum. The camera will be calibrate about
;     a provided lambda at which the calibration will be 1:1.
;
;CALLING SEQUENCE:
;     calibration = calibrate(emiss, observed_emiss, lam, known_lambda
;
;INPUTS:
;     These arrays should all be of equal dimension:
;
;     emiss: 1D float/double array containing calculated emissivities 
;       (could be calculated with emiss function above)
;     observed_emiss: 1D float/double array containing emissivities 
;       observed by a camera
;     lam: 1D array containing the wavelength that each element of
;       emiss represents
;
;     known_lambda: A wavelength to calibrate about (at this
;       wavelength, the calibration multiplier will be 1). This value
;       must be within the range of lam
;
;OUTPUTS:
;     Returns a 1D calibration array. This is of the form
;     known_emiss/observed_emiss, so camera channels should be
;     multiplied by the respective calibration value

FUNCTION calibrate,emiss,observed_emiss,lam,known_lambda
  IF ((n_elements(emiss) NE n_elements(observed_emiss)) OR (n_elements(observed_emiss) NE n_elements(lam)) OR (n_elements(emiss) NE n_elements(lam))) THEN BEGIN
     print,'Error: emiss, observed_emiss, and lam do not all have the same dimension. Terminating.'
     return,[0]
  ENDIF

  IF ((known_lambda GT MAX(lam)) OR (known_lambda LT MIN(lam))) THEN BEGIN
     print,'Error: known_lambda is not within the bounds of lam. Terminating.'
     return,[0]
  ENDIF

  lam_index = closest(lam,known_lambda)
  correction = observed_emiss[lam_index] / emiss[lam_index]
  calibration = fltarr(n_elements(lam))
  FOR i=0,n_elements(lam)-1 DO BEGIN
     calibration[i] = (emiss[i] * correction) / observed_emiss[i]
  ENDFOR
  return,calibration
END



 


;NAME:
;     MakeSpecFromMin
;
;PURPOSE:
;     Extracts the continnum from a spectrum by taking the minimum
;     value of bands of size width
;
;CALLING SEQUENCE:
;     minSpec = MakeSpecFromMin(spec,ExtractMin(spec,width))
;     minLam = MakeSpecFromMin(lam,ExtractMin(spec,width))
;
;INPUTS:
;     arrayInp: A 1D array for which to return the elements that
;       correspond to the minimum of the spectrum, typically would
;       either be the spectrum itself or wavelength/energy
;     minSpecIndicies: 1D array containing the indicies of the
;       minimums of the spectrum. Could be found by calling
;       ExtractMin(spec,width)
;
;OUTPUTS:
;     Returns a 1D array containing the values corresponding to the
;     minimums of the spectrum

FUNCTION ExtractMin,spec,width
  numberSlots=float(n_elements(spec)/width)
  minSpecIndicies=intarr(numberSlots)
  FOR i=0,numberSlots-1 DO BEGIN
     amin=min(spec[i*width:(i*width+width-1)],z)
     minSpecIndicies[i]=z+(width)*i
  ENDFOR
  return,minSpecIndicies
END

FUNCTION MakeSpecFromMin,arrayInp,minSpecIndicies
  minSpec=fltarr(n_elements(minSpecIndicies))
  FOR i=0,n_elements(minSpecIndicies)-1 DO BEGIN
     minSpec[i]=arrayInp[minSpecIndicies[i]]
  ENDFOR
  return,minSpec
END



  




;NAME:
;     Calibrate_XEUS_Loweus
;
;PURPOSE:
;     Performs a relative calibration of any spectrum from XEUS or
;     Loweus. Can return either an array, or a plot. Due to all the
;     line_br calls, this function may take some time to run. A
;     progress indicator will be shown.
;
;CALLING SEQUENCE:
;     Calibrate_XEUS_Loweus,shot,t,z_eff,calibration,lam,specIn,specOut,plotWhat=plotWhat,calibrateWhat=calibrateWhat,lambdaPoint=lambdaPoint
;
;INPUTS:
;     shot: Provide a shot to calibrate
;     t: Provide a time in that shot to calibrate
;     z_eff: Provide z_eff of the plasma. Since this is a relative
;       calibration, it should not affect the results, but some value
;       has to be provided to emiss calculation
;
;OPTIONAL INPUTS:
;     plotWhat:
;       0: Return no plots (DEFAULT)
;       1: Plot the calibration constants with respect to wavelength
;       2: Plot a calibrated spectrum in red, over the original
;         spectrum in white. The original spectrum is shifted by an
;         angstrom to make it easier to visually analyze. The returned
;         data is not shifted, only the plot.
;     calibrateWhat:
;       0: XEUS (DEFAULT)
;       1: Loweus
;     LambdaPoint:
;       Provides a wavelength to calibrate about (more information in
;       calibrate function documentation above). Since you
;       probably don't care about this, default for XEUS is 30A
;       and default for Loweus is 150A which both work well.
;
;OUTPUTS:
;     Calibration: A 1D array containing the calibration multipliers
;       from the calling of the calibrate function
;     Lam: Returns the wavelengths that were calibrated
;     specIn: Returns the input spectrum at those wavelengths
;     specOut: Returns the relatively calibrated spectrum at those
;       wavelengths
;

PRO Calibrate_XEUS_LOWEUS,shot,t,z_eff,calibration,lam,specIn,specOut,plotWhat=plotWhat,calibrateWhat=calibrateWhat,lambdaPoint=lambdaPoint
  IF NOT keyword_set(plotWhat) THEN plotWhat=0
  IF NOT keyword_set(calibrateWhat) THEN calibrateWhat=0
  IF NOT keyword_set(lambdaPoint) THEN BEGIN
     IF calibrateWhat EQ 0 THEN lambdaPoint = 30
     IF calibrateWhat EQ 1 THEN lambdaPoint = 150
  ENDIF

  mdsopen,'spectroscopy',shot
  IF (calibrateWhat EQ 0) THEN BEGIN
     lam=mdsvalue('dim_of(\SPECTROSCOPY::TOP.XEUS.SPEC,0)')
     times=mdsvalue('dim_of(\SPECTROSCOPY::TOP.XEUS.SPEC,1)')
     spec=mdsvalue('\SPECTROSCOPY::TOP.XEUS.SPEC')
  ENDIF
  IF (calibrateWhat EQ 1) THEN BEGIN
     lam=mdsvalue('dim_of(\SPECTROSCOPY::TOP.LOWEUS.SPEC,0)')
     times=mdsvalue('dim_of(\SPECTROSCOPY::TOP.LOWEUS.SPEC,1)')
     spec=mdsvalue('\SPECTROSCOPY::TOP.LOWEUS.SPEC')
  ENDIF

  ;fill pos vectors
  xeus_pos = [2.561,0.21,0.179,6.5*!pi/180]
  loweus_pos = [2.561,0.21,0.179,6.5*!pi/180]

  pos = xeus_pos
  IF calibrateWhat EQ 1 THEN BEGIN
     pos=loweus_pos
  ENDIF

  genpos_pos_reform,pos,[0.44,1,0.6,-0.6]

  ;grab spectrum at time t and find the minimum approximation of it
  t_index = closest(times, t)
  specT = spec[*,t_index]
  indicies = [0,ExtractMin(specT,50),n_elements(specT)-50]
  specT_clean = MakeSpecFromMin(specT,indicies)
  
  
  ;calibrationRaw contains only points that are a part of the minimum spectrum
  calibrationRaw = calibrate(emiss(MakeSpecFromMin(lam,indicies),shot,t,z_eff,pos) , specT_clean, MakeSpecFromMin(lam,indicies), lambdaPoint)
  

  calibration = interpol(calibrationRaw, lam[indicies], lam)
  

  ;Performs plotting operation selected in the optional plotWhat input
  IF (plotWhat EQ 1) THEN BEGIN ;requested plot of calibration curve
     plot,MakeSpecFromMin(lam,indicies),calibrationRaw,YRANGE=[0,50]
     oplot,lam,calibration
  ENDIF

  IF (plotWHAT EQ 2) THEN BEGIN ;requested a calibrated spectrum (in red over white original)
     loadct,5

     plot,lam+0.35+0.65*calibrateWhat,specT,color=255,YRANGE=[0,max([specT,specT*calibration])]
     oplot,lam,specT*calibration,color=105

     values = specT*calibration
  ENDIF

  specIn=specT
  specOut=specT*calibration
END








;NAME:
;     zeff_brem_calibrate
;
;PURPOSE:
;     Calibrates TS, zmeter, or a band of loweus to the calculated
;     emissivity for that wavelength band through the filter
;     transmission on the camera. 
;     
;     On a zeff=1 shot, set zeffset=1 and the result can be used as a
;     calibration constant in zeff_brem
;
;     Uses zeff from zave if no zeffset override
;
;     User will be prompted for low and high bounds of calibration
;     after showing a plot. Alternatively, these can be provided in
;     tlow=tlow and thigh=thigh. If these optional inputs are both set, no plots
;     will be shown
;
;CALLING SEQUENCE:
;     zeff_brem_calibrate, shot, measure, precision=precision,
;     tprecision=tprecision, time, magnitude, calibration
;
;INPUTS:
;     shot: provide a shot to analyze
;     measure: What device to operate on
;        'TS_4' = Thompson band 4
;        'zmeter' = zmeter
;        'Loweus' = A high wavelength band in Loweus
;
;OPTIONAL INPUTS:
;     precision (default 100): How many wavelengths to sum in the
;       integral
;     tprecision (default 100): How many time points to analyze
;     tlow: Lower time bound of calibration
;     thigh: Upper time bound of calibration
;     plotTrue (default 0): set to 1 to make a plot
;     zeffset (default 0): 0 reads from zave, alternate value can be
;       specified to override zave
;
;OUTPUTS:
;     Provides a plot if plotTrue=1
;        Blue = Actual emissivity data in whatever raw units TS,
;          zmeter, and Loweus have in tree
;        Red = Filtered actual emissivity data (same units as blue)
;        White = Calculated emisisivity prediction (W/m^2) within the 
;          detector's transmission function
;     time: An array of time values that match up to magnitude entries
;     magnitude: An array of actual data observed emissivities
;     timesCalc: An array of time values that match up to magnitudeCalc
;     magnitudeCalc: An array of calculated emissivities
;     calibration: Returns the calibration value in W/m^2/count;     




;### helper functions###
FUNCTION z_meter_transmission,shot,view,band,lambda
                                ;shot, view, and band are all ignored
                                ;by this function. Provided for
                                ;compatibility with the execute in
                                ;integrate_calculation_time

  center=536.
  FWHM=3.

  return, exp(-(lambda-center)^2/(2*(FWHM/2.35482)^2))
END

FUNCTION ts_transmission,shot,view,band,lambda
  mdsopen,'electrons',shot
  bandpass = mdsvalue('\CALIB:INST_FUNC')
  wavelength = mdsvalue('\CALIB:LAMBDA')
  index = closest(wavelength,lambda)
  return, bandpass[view,band,index]

End

FUNCTION loweus_transmission,shot,view,band,lambda
  return, 1
END

FUNCTION avgNoise, spec, avgWidth
  n = n_elements(spec)/avgWidth
  output = fltarr(n-2)
  for i=0, n-3 do begin
     output[i] = float(sum(spec[long(i)*long(avgWidth):long(i)*long(avgWidth) + long(avgWidth) - 1])) / avgWidth
  endfor
  return, output
end
;########################





PRO zeff_brem_calibrate,shot,measure,precision=precision,zeffset=zeffset,tprecision=tprecision,tLow=tLow,tHigh=tHigh,yplot=yplot,time,magnitude,timesCalc,magnitudeCalc,calibration
  IF (string(measure) EQ 'zmeter') THEN BEGIN
     MDSOPEN,'spectroscopy',shot
     magnitude = mdsvalue('\SPECTROSCOPY::TOP.Z_METER:Z_BRIGHT')
     time = mdsvalue('dim_of(\SPECTROSCOPY::TOP.Z_METER:Z_BRIGHT,0)')
     filter = 'z_meter_transmission'
     pos = [1,0,0.68,0]
     lamRange = [520,560]
     avgWidth = 100
  ENDIF
  IF (string(measure) EQ 'TS_4') THEN BEGIN
     MDSOPEN,'electrons',shot
     magnitude = mdsvalue('\ELECTRONS::TOP.ZEFF_YAG.SIGNALS.BAND_4:POLYC_07')
     time = mdsvalue('dim_of(\ELECTRONS::TOP.ZEFF_YAG.SIGNALS.BAND_4:POLYC_07,0)')
     filter = 'ts_transmission'
     pos = [1.69,0,0,0]
     lamRange = [10000,11000]
     avgWidth = 100
  ENDIF
  IF (string(measure) EQ 'loweus') THEN BEGIN
     filter = 'loweus_transmission'
     MDSOPEN,'cmod',shot
     specRead = mdsvalue('\CMOD::TOP.SPECTROSCOPY.LOWEUS:SPEC')
     lamRead = mdsvalue('dim_of(\CMOD::TOP.SPECTROSCOPY.LOWEUS:SPEC,0)')
     time = mdsvalue('dim_of(\CMOD::TOP.SPECTROSCOPY.LOWEUS:SPEC,1)')
     lamRange = [lamRead[70],lamRead[50]]
     magnitude = fltarr(n_elements(time))
     FOR i=0, n_elements(magnitude) - 1 DO BEGIN
        magnitude[i] = int_tabulated(reverse(lamRead[50:70]), reverse(specRead[50:70,i]))
     ENDFOR
     pos = [2.561,0.21,0.179,6.5*!pi/180]
     avgWidth = 5
  ENDIF



  IF ((string(measure) NE 'TS_4') AND (string(measure) NE 'zmeter') AND (string(measure) NE 'loweus')) THEN BEGIN
     print,'error: No valid device was selected'
  ENDIF

  cmod_ts,shot,ts
  
  IF NOT keyword_set(precision) THEN precision=100
  IF NOT keyword_Set(tprecision) THEN tprecision=100
  IF NOT keyword_set(rerun) THEN rerun=0
  IF NOT keyword_set(zeffset) THEN zeffset = 0
  
  lam = fltarr(precision)
  timeEmiss = fltarr(tprecision)
  times = fltarr(tprecision)


  filterFunc = fltarr(precision)
  FOR i=0, precision - 1 DO BEGIN
     lam[i] = lamRange[0] + float(lamRange[1] - lamRange[0]) / precision * i
     filterFunc[i] = call_function(filter, shot,7,0,lam[i])
  ENDFOR

    mdsopen,'cmod',shot
  z = mdsvalue('\CMOD::TOP.SPECTROSCOPY.Z_METER.ANALYSIS:Z_AVE')
  zt = mdsvalue('dim_of(\CMOD::TOP.SPECTROSCOPY.Z_METER.ANALYSIS:Z_AVE)')
  mdsclose



  efit={time:line_gettimes(shot),lcfs:line_getlcfs(shot),axis:line_getaxis(shot),rmid:line_getrmid(shot)}
  FOR t_loop=0,tprecision - 1 DO BEGIN
     times[t_loop] = 2*float(t_loop) / tprecision + 0.05
        
     z_index = closest(zt,times[t_loop])

     IF zeffset EQ 0 THEN z_eff = z[z_index]
     IF zeffset NE 0 THEN z_eff = zeffset
        
     emissiv = fltarr(n_elements(lam))
        
     N_Eraw = reform(ts.N_E[closest(ts.TIME,times[t_loop]),*])
     TEraw = reform(ts.TE[closest(ts.TIME,times[t_loop]),*])
     goodLocations = where(N_Eraw NE 0)
     closestT = closest(ts.TIME,times[t_loop])
     RadiusRaw = reform(ts.R[closestT,*])
        
     IF (goodLocations EQ [-1]) THEN BEGIN
        N_E = [0]
        TE = [0]
        Radius = [0]
           
     ENDIF ELSE BEGIN
        N_E = N_Eraw[goodLocations]
        TE = TEraw[goodLocations]
        Radius = RadiusRaw[goodLocations]
     ENDELSE
        
     unint = fltarr(n_elements(Radius))
     FOR q=0, n_elements(Radius)- 1 DO BEGIN
        unint[q] = int_tabulated(lam, energyEmiss(z_eff, N_E[q], TE[q], lam)*filterFunc)
     ENDFOR
     timeEmiss[t_loop] = genpos_line_br(pos, [unint], Radius,[times[t_loop]], shot, times[t_loop],efit=efit)
     
     
     print, string(27B) + '[4D' + string(27B) + '[K', format='(A, $)'
     print, string(int(float(t_loop)/tprecision*100)), format='(I3, "%", $)'
     
  ENDFOR



  
;plot:
  loadct,5

  void = min(abs(time - 1),indexTime)

  IF NOT keyword_set(yplot) THEN BEGIN
     yrangeplot=[0,max([magnitude, timeEmiss / timeEmiss[closest(times,1)]*magnitude[indexTime]])]
  ENDIF ELSE BEGIN
     yrangeplot=[0,yplot]
  ENDELSE


  if (NOT keyword_Set(tlow)) AND (NOT keyword_set(thigh)) then begin
     plot,[0],[0],xrange=[0,2],yrange=yrangeplot

     ;Plot the observed emissivity data:
     oplot,time,magnitude,color=50

     ;Plot the calculated emissivity data:
     oplot,times[WHERE(timeEmiss GT 0)],(timeEmiss[WHERE(timeEmiss GT 0)] / timeEmiss[closest(times,1)]*magnitude[indexTime])

     ;Plot the filtered observed emissivity data:
     oplot,avgNoise(time, avgWidth), avgNoise(magnitude, avgWidth), color=100
  endif

  magnitudeCalc = timeEmiss[WHERE(timeEmiss GT 0)]
  timesCalc = times[WHERE(timeEmiss GT 0)]
  
  avgTime = avgNoise(time, avgWidth)
  avgMagnitude = avgNoise(magnitude, avgWidth)
  calibration = 0
  timesReal = times[WHERE(timeEmiss GT 0)]
  timeEmissReal = timeEmiss[WHERE(timeEmiss GT 0)]

  IF NOT keyword_Set(tLow) THEN BEGIN
     read,tLowCal,prompt='beginning time: '
  ENDIF ELSE BEGIN
     tLowCal = tLow
  ENDELSE

  IF NOT keyword_Set(tHigh) THEN BEGIN
     read,tHighCal,prompt='ending time: '
  ENDIF ELSE BEGIN
     tHighCal = tHigh
  ENDELSE


  ;Perform calibration between tlow and thigh:
  timesGood = timesReal[WHERE(timesReal GT tLowCal)]
  timeEmissGood = timeEmissReal[WHERE(timesReal GT tLowCal)]
  timesGoodW = timesGood
  timesGood = timesGood[WHERE(timesGood LT tHighCal)]
  timeEmissGood = timeEmissGood[WHERE(timesGoodW LT tHighCal)]
  
  FOR i=0, n_elements(timeEmissGood)-1 DO BEGIN
     calibration += timeEmissGood[i] / avgMagnitude[closest(avgTime, timesGood[i])]
  ENDFOR
  calibration = calibration / n_elements(timeEmissGood)

  print,'calibration (W/m^2/count): ' + string(calibration)
END







;NAME: 
;     zeff_brem
;
;PURPOSE:
;     Finds the observed zeff at any particular time by comparing a band of loweus, zmeter,
;     or TS view 4 to the calculated emissivity
;
;     A calibration constant should be provided for each detector, in
;     W/m^2/count. This could be obtained by running the following on
;     a zeff=1 shot:
;       zeff_brem_calibrate,shot,'detector',zeffset=1
;     select the time region where zeff=1 when prompted
;
;CALLING SEQUENCE:
;     zeff =
;         zeff_brem(shot,measure,t,efit=efit,magnitudeInput=magnitudeInput,timeInput=timeInput,cmodts=cmodts,lamInput=lamInput,filterInput=filterInput)
;
;INPUTS:
;     shot: Which Shot?
;     measure:
;         'TS_4' - view 4 of thompson data
;         'zmeter' - zmeter
;         'loweus' - a band of loweus from approximately 271.93nm to 275.35nm
;
;OPTIONAL INPUTS:
;     All optional inputs are provided to increase speed if a script
;     loops this function many times. All of these items are
;     calculated if not provided.
;
;     efit: Provide an efit of the following form
;         efit={time:line_gettimes(shot),lcfs:line_getlcfs(shot),axis:line_getaxis(shot),rmid:line_getrmid(shot)}
;     magnitudeInput: The calculated emissivity corresponding to
;         timeInput
;     timeInput: The times corresponding to magnitudeInput
;     cmodts: The cmod_ts for this shot
;     lamInput: an array of wavelengths to look at. This range should
;         make sense for the camera's filter
;     filterInput: the corresponding filter inputs
;
;     
;     Other optional inputs not related to speed:
;     
;     n_eForce: Provide an electron density rather than using ts data
;     teForce: Provide an electron temperature rather than using ts data
;
;OUTPUTS:
;     Function returns the zeff value found in an array of 3
;     elements. [zeff, zeff_low_bound, zeff_high_bound]


FUNCTION zeff_brem,shot,measure,t,efit=efit,magnitudeInput=magnitudeInput,timeInput=timeInput,cmodts=cmodts,lamInput=lamInput,filterInput=filterInput,n_eForce=n_eForce,teForce=teForce
  IF (string(measure) EQ 'zmeter') THEN BEGIN
     ;#####
     ;Calibration constant:
     calibration=7.69286e-12   ;7.69286e-12 comes from setting zeff=1 in shot 1120522015
     ;#####


     if (NOT keyword_set(magnitudeInput)) OR (NOT keyword_set(timeInput)) then MDSOPEN,'spectroscopy',shot
     if NOT keyword_set(magnitudeInput) then magnitude = mdsvalue('\SPECTROSCOPY::TOP.Z_METER:Z_BRIGHT') else magnitude=magnitudeInput
     if NOT keyword_set(timeInput) then time = mdsvalue('dim_of(\SPECTROSCOPY::TOP.Z_METER:Z_BRIGHT,0)') else time = timeInput
     filter = 'z_meter_transmission'
     pos = [1,0,0.68,0]
    ; pos = [1.69,0,0,0]
     lamRange = [520,560]
     magnitude = avgNoise(magnitude,100)
     time = avgNoise(time,100)
  ENDIF
  IF (string(measure) EQ 'TS_4') THEN BEGIN
     ;#####
     ;Calibration constant:
     calibration=0.436225    ;0.436225 comes from setting zeff=1 in shot 1120522015
     ;#####

     if (NOT keyword_set(magnitudeInput)) OR (NOT keyword_set(timeInput)) then MDSOPEN,'electrons',shot
     if NOT keyword_set(magnitudeInput) then magnitude = mdsvalue('\ELECTRONS::TOP.ZEFF_YAG.SIGNALS.BAND_4:POLYC_07') else magnitude = magnitudeInput
     if NOT keyword_set(timeInput) then time = mdsvalue('dim_of(\ELECTRONS::TOP.ZEFF_YAG.SIGNALS.BAND_4:POLYC_07,0)') else time = timeInput
     filter = 'ts_transmission'
     pos = [1.69,0,0,0]
     lamRange = [10000,11000]
     magnitude = avgNoise(magnitude,100)
     time = avgNoise(time,100)
  ENDIF
  IF (string(measure) EQ 'loweus') THEN BEGIN
     ;#####
     ;Calibration constant:
     calibration=2886.39  ;unknown source or accuracy. Calibrating loweus currently seems to lack much meaning.
     ;#####

     filter = 'loweus_transmission'
     MDSOPEN,'cmod',shot
     specRead = mdsvalue('\CMOD::TOP.SPECTROSCOPY.LOWEUS:SPEC')
     lamRead = mdsvalue('dim_of(\CMOD::TOP.SPECTROSCOPY.LOWEUS:SPEC,0)')
     time = mdsvalue('dim_of(\CMOD::TOP.SPECTROSCOPY.LOWEUS:SPEC,1)')

  
     lamRange = [lamRead[70],lamRead[50]]
     magnitude = fltarr(n_elements(time))
     FOR i=0, n_elements(magnitude) - 1 DO BEGIN
        magnitude[i] = int_tabulated(reverse(lamRead[50:70]), reverse(specRead[50:70,i]))
     ENDFOR
     pos = [2.561,0.21,0.179,6.5*!pi/180]
     magnitude = avgNoise(magnitude,5)
     time = avgNoise(time,5)
  ENDIF
  if not keyword_set(cmodts) then cmod_ts,shot,ts else ts=cmodts


  
  ;Find calculated data:
  z_eff_dummy = 1
  precision=500


  if (NOT keyword_set(lamInput)) OR (NOT keyword_set(filterInput)) then begin
     lam = fltarr(precision)
     FOR i=0, precision - 1 DO BEGIN
        lam[i] = lamRange[0] + float(lamRange[1] - lamRange[0]) / precision * i
     ENDFOR
  endif else begin
     lam = lamInput
  endelse

  if NOT keyword_set(filterInput) then begin
     filterfunc = fltarr(precision)
     for i=0,precision-1 do begin
        filterFunc[i] = call_function(filter, shot,7,0,lam[i])
     endfor
  endif else begin
     filterFunc=filterInput
  endelse

 
  N_Eraw = reform(ts.N_E[closest(ts.TIME,t),*])
  TEraw = reform(ts.TE[closest(ts.TIME,t),*])

  DN_Eraw = reform(ts.DN_E[closest(ts.TIME,t),*])
  DTEraw = reform(ts.DTE[closest(ts.TIMe,t),*])

  goodLocations = where(N_Eraw NE 0)

  closestT = closest(ts.TIME,t)
  RadiusRaw = reform(ts.R[closestT,*])
  
  IF (goodLocations EQ [-1]) THEN BEGIN
     N_E = [0]
     TE = [0]
     Radius = [0]
     DN_E = [0]
     DTE = [0]
  ENDIF ELSE BEGIN
     N_E = N_Eraw[goodLocations]
     TE = TEraw[goodLocations]
     Radius = RadiusRaw[goodLocations]
     DN_E = DN_Eraw[goodLocations]
     DTE = DTEraw[goodLocations]
  ENDELSE 
  
  unint=fltarr(n_elements(N_E))
  uninthigh=fltarr(n_elements(N_E))
  unintlow=fltarr(n_elements(N_E))


  for i=0,n_elements(N_E)-1 do begin
     if keyword_set(n_eForce) then N_E[i] = n_eForce
     if keyword_set(teForce) then TE[i] = teForce
  endfor


  FOR q=0, n_elements(Radius)- 1 DO BEGIN
     unint[q] = int_tabulated(lam,energyEmiss(z_eff_dummy,N_E[q],TE[q],lam)*filterFunc)
     uninthigh[q] = int_tabulated(lam,energyEmiss(z_eff_dummy,N_E[q]+DN_E[q],TE[q]-DTE[q],lam)*filterFunc)
     unintlow[q] = int_tabulated(lam,energyEmiss(z_eff_dummy,N_E[q]-DN_E[q],TE[q]+DTE[q],lam)*filterFunc)
  ENDFOR
  if keyword_set(efit) then begin
     timeEmiss = genpos_line_br(pos, [unint], Radius,[t], shot, t, efit=efit)
     timeEmissHigh = genpos_line_br(pos, [uninthigh], Radius, [t], shot, t, efit=efit)
     timeEmissLow = genpos_line_br(pos, [unintlow], Radius, [t], shot, t, efit=efit)
  endif else begin
     timeEmiss = genpos_line_br(pos, [unint], Radius,[t], shot, t)
     timeEmissHigh = genpos_line_br(pos, [uninthigh], Radius, [t], shot, t)
     timeEmissLow = genpos_line_br(pos, [unintlow], Radius, [t], shot, t)
  endelse

                                ;Compare to actual data:
  magnitudeCompare = magnitude[closest(time,t)]

 ; return,magnitudeCompare*calibration/timeEmiss
  return,[magnitudeCompare*calibration/timeEmiss, magnitudeCompare*calibration/timeEmissHigh, magnitudeCompare*calibration/timeEmissLow]
END








;NAME: 
;     zeff_brem_t
;
;PURPOSE:
;     view the time evolution of zeff calculated from continuum
;     brightness
;
;CALLING SEQUENCE:
;     zeff_brem_t,shot,measure,tprecision=tprecision,plot=plot,n_eForce=n_eForce,teForce=teForce,timeRange=timeRange,zeff,t,zeffHigh,zeffLow
;
;INPUTS:
;     shot - which shot to look at
;     measure:
;         'TS_4' - view 4 of thompson data
;         'zmeter' - zmeter
;         'loweus' - a band of loweus from approximately 271.93nm to
;                    275.35nm
;     
;OPTIONAL INPUTS:
;     tprecision - sets number of time points to look at within timeRange
;         default=50
;     /plot - display plot of output
;         white trace - zave
;         red traces - zeff_brem, along with high and low error lines
;                      (from TS error data)
;     n_eForce - set an electron density rather than using TS data
;     teForce - set an electron temperature rather than using TS data
;     timeRange - [timeLow, timeHigh] for analysis and plotting
;     
;OUTPUTS:
;     t - time array of output
;     zeff - calculated zeff
;     zeffHigh,zeffLow - calculated zeff error bounds


pro zeff_brem_t,shot,measure,tprecision=tprecision,plot=plot,n_eForce=n_eForce,teForce=teForce,timeRange=timeRange,t,zeff,zeffHigh,zeffLow
  lambdaPrecision=500
  wavelengths = fltarr(lambdaPrecision)
  filterFunction = fltarr(lambdaPrecision)
  IF NOT keyword_set(timeRange) THEN timeRange = [0,2]
  IF NOT keyword_set(tprecision) THEN tprecision = 50
  IF measure EQ 'zmeter' THEN BEGIN
     ;zmeter declarations:
     lamRange = [520, 560]
     mdsopen,'spectroscopy',shot
     measuredEmiss = mdsvalue('\SPECTROSCOPY::TOP.Z_METER:Z_BRIGHT')
     measuredEmissTime = mdsvalue('dim_of(\SPECTROSCOPY::TOP.Z_METER:Z_BRIGHT,0)')
     mdsclose
     filterCall = 'z_meter_transmission'
  ENDIF
  IF measure EQ 'TS_4' THEN BEGIN
     ;TS_4 declarations:
     lamRange = [10000,11000]
     mdsopen,'electrons',shot
     measuredEmiss = mdsvalue('\ELECTRONS::TOP.ZEFF_YAG.SIGNALS.BAND_4:POLYC_07')
     measuredEmissTime = mdsvalue('dim_of(\ELECTRONS::TOP.ZEFF_YAG.SIGNALS.BAND_4:POLYC_07,0)')
     mdsclose
     filterCall = 'TS_transmission'
  ENDIF

  times = fltarr(tprecision)
  resultsHigh = fltarr(tprecision)
  resultsActual = fltarr(tprecision)
  resultsLow = fltarr(tprecision)

  cmod_ts,shot,ts

  ;pre-calculate the efit to save processing time:
  efit={time:line_gettimes(shot),lcfs:line_getlcfs(shot),axis:line_getaxis(shot),rmid:line_getrmid(shot)}

  FOR i=0, lambdaPrecision - 1 DO BEGIN
     wavelengths[i] = lamRange[0] + float(lamRange[1] - lamRange[0]) / lambdaPrecision * i
     filterFunction[i] = call_function(filterCall,shot,7,0,wavelengths[i])
  ENDFOR

  FOR i=0,n_elements(times)-1 DO BEGIN
     times[i] = float(i) / tprecision * float(timeRange[1] - timeRange[0]) + timeRange[0]
     IF NOT keyword_set(n_eForce) AND NOT keyword_set(teForce) THEN BEGIN
        zeffbrem = zeff_brem(shot,measure,times[i],efit=efit,cmodts=ts,magnitudeInput=measuredEmiss,timeInput=measuredEmissTime,lamInput=wavelengths,filterInput=filterFunction)
     ENDIF
     IF keyword_set(n_eForce) AND NOT keyword_set(teForce) THEN BEGIN
        zeffbrem = zeff_brem(shot,'zmeter',times[i],efit=efit,magnitudeInput=magnitudeZmeter,timeInput=timeZmeter,cmodts=ts,lamInput=lamzmeter,filterInput=filterfunczmeter,n_eForce=n_eForce)
     ENDIF
     IF NOT keyword_set(n_eForce) AND keyword_set(teForce) THEN BEGIN
        zeffbrem = zeff_brem(shot,'zmeter',times[i],efit=efit,magnitudeInput=magnitudeZmeter,timeInput=timeZmeter,cmodts=ts,lamInput=lamzmeter,filterInput=filterfunczmeter,teForce=teForce)
     ENDIF
     IF keyword_set(n_eForce) AND keyword_set(teForce) THEN BEGIN
        zeffbrem = zeff_brem(shot,'zmeter',times[i],efit=efit,magnitudeInput=magnitudeZmeter,timeInput=timeZmeter,cmodts=ts,lamInput=lamzmeter,filterInput=filterfunczmeter,n_eForce=n_eForce,teForce=teForce)
     ENDIF
     resultsActual[i] = zeffbrem[0]
     resultsLow[i] = zeffBrem[1]
     resultsHigh[i] = zeffBrem[2]

     print, string(27B) + '[4D' + string(27B) + '[K', format='(A, $)'
     print, string(int(float(i)/n_elements(times)*100)), format='(I3, "%", $)'
  ENDFOR

  zeff = resultsActual[WHERE( resultsActual LT 30)]
  zeffHigh = resultsLow[WHERE( resultsActual LT 30)]
  zeffLow = resultsHigh[WHERE( resultsActual LT 30)]
  t = times[WHERE( resultsActual LT 30)]

  IF keyword_set(plot) THEN BEGIN
     loadct,5
     plot,[0],[0],color=255,yrange=[0,10],xrange=timeRange
     oplot,t,zeffHigh,color=110
     oplot,t,zeffLow,color=110
     oplot,t,zeff,color=100

      mdsopen,'cmod',shot
      zave = mdsvalue('\CMOD::TOP.SPECTROSCOPY.Z_METER.ANALYSIS:Z_AVE')
      tave = mdsvalue('dim_of(\CMOD::TOP.SPECTROSCOPY.Z_METER.ANALYSIS:Z_AVE)')
      mdsclose

      oplot,tave,zave,color=255
  ENDIF

END








































;temporary script to look at spectral line ratios

PRO timeRelationship,shot,res,tlow,thigh
  time = fltarr(res)
  output = fltarr(res)
  middle = fltarr(res)
  fsubarray = fltarr(res)
  mocorrect = fltarr(res)
  for z=0,res-1 do begin
     t=float(thigh-tlow)/res*z+tlow
     zeff=1
     calibrate_xeus_loweus,shot,t,zeff,calibration,lam,specIn,specOut,calibrateWhat=1
 

     SpecOutAvg = SpecOut - Interpol(SpecOut[ExtractMin(SpecOut,50)],lam[ExtractMin(SpecOut,50)],lam)



     Elements=['Fe','Cr','Cu','Mo']
     integrals = fltarr(n_elements(elements))

     Ranges=[191.8,192.6,254.9,255.7,   222.7,223.7,279.5,280.5,   153.3,154.2,224.5,225.4,   127.2,128.3,175.9,177]
 

     for i=0,n_elements(integrals)-1 do begin
        low = Sum(SpecOutAvg[closest(lam,Ranges[(i*4)+1]):closest(lam,Ranges[i*4])])*(lam[closest(lam,Ranges[(i*4)])]-lam[closest(lam,Ranges[(i*4+1)])])
        high = Sum(SpecOutAvg[closest(lam,Ranges[(i*4)+3]):closest(lam,Ranges[(i*4)+2])])*(lam[closest(lam,Ranges[(i*4)+2])]-lam[closest(lam,Ranges[(i*4)+3])])
        mid = Sum(SpecOutAvg[closest(lam, Ranges[(i*4)+2]):closest(lam, Ranges[(i*4)+1])])*(lam[closest(lam,Ranges[(i*4)+1])]-lam[closest(lam,Ranges[(i*4+2)])])
        print,(Elements[i]) + String(low/high)
        integrals[i] = low/high
     endfor

  ;   Print,('Middle') + String(mid/low)
     time[z]=t
     
     output[z]=low/high
     middle[z]=mid/low

     Fsub = Sum(SpecOutAvg[closest(lam, 113.4):closest(lam,112.45)])*(lam[closest(lam,112.45)]-lam[closest(lam,113.4)])
     
     i=3
     molyLow = Sum(SpecOutAvg[closest(lam,Ranges[(i*4)+1]):closest(lam,Ranges[i*4])])*(lam[closest(lam,Ranges[(i*4)])]-lam[closest(lam,Ranges[(i*4+1)])])
     molyHigh = Sum(SpecOutAvg[closest(lam,Ranges[(i*4)+3]):closest(lam,Ranges[(i*4)+2])])*(lam[closest(lam,Ranges[(i*4)+2])]-lam[closest(lam,Ranges[(i*4)+3])])

     
     
     foffset = (Sum(SpecOutAvg[closest(lam,Ranges[13]):closest(lam,Ranges[12])])/(lam[closest(lam,Ranges[12])]-lam[closest(lam,Ranges[13])]))/Fsub

     print, 'F-offset:' + string(foffset)
     fsubarray[z] = string(foffset)
     print, 'MO-corrected:' + string((molyLow - 1.89*Fsub)/molyHigh)
     mocorrect[z] = (molyLow - 1.89*Fsub)/molyHigh
  endfor
  loadct,5
  plot,time,output,yrange=[-3,15]
 ; oplot,time,middle,color=50
  oplot,time,mocorrect,color=100
  print,''
;  print,sum(output[WHERE(output LT 100 AND output GT 0)])/n_elements(WHERE(output LT 100 AND output GT 0)) print,'mid average' + string(sum(middle[WHERE(middle LT 100 AND
;  middle GT 0)])/n_elements(WHERE(middle LT 100 AND middle GT 0)))

  print,'fsubavg:' + sum(fsubarray)/n_elements(fsubarray)
  
end

