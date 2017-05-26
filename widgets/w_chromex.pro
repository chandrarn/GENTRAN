; This is a modified version of the chr_view.pro to allow multiple shot viewing.
; KMarr 07-01-03
;******************************************************************************
; This program creates a widget which allows the user to analyze and print
; data from the Chromex spectrometer.  It needs a X Window display
; to run.
;
; 9/30/93: first functional version.
;
; 6/21/94: added a slot for a WIDGET_TEXT_CH event in the main event procedure.
; this makes pressing the RETURN key equivalent to pressing the APPLY button when
; entering a shot number.
;
; 4/11/95: starting to rewrite program for Chromex spectrometer system.
;
; 10 July 2000: IHH added a button to lock/unlock the yrange, using
; the parameter the_yrange added to common data and modifications to
; plot_frame, and new procedure lock_yrange.
; Added the ability to read in known wavelengths and over plot
; them. This is toggled by a button.
; Added ability to not plot an output file
; Removed the default landscape mode. Do an independent
; set_plot,'ps'; device,/landscape set_plot,'x' to set it.
; Rearranged the widget layout to allow more buttons
; Stopped the destroying and rerealizing of sliders. It is unnecessary.
;
; 7/1/03 Added the ability to see a second shot for comparison
; techniques.  Made many changes to plotting routine.  Multiple plots
; are now stored in x and y.  The size of each array is stored
; in xyDim with the respective index.  For example to plot 2 arrays of
; size 10 and 20 in color:
; num_of_plots = 2   (arrays)
; xyDim[0] = 10  (elements)
; xyDim[1] = 20  (elements)
; plotColorFlag = 1  ; plotSymbolFlag = 1 (yes color and symbol)
; plotColor[0] = 0 ; plotColor[1] = 255 (which color)
; plotSymbol[0] = -4 ; plotSymbol[1] = 6  (which symbol)
; Ways to store actual data (either works but xyDim is more robust):
; x[0:xyDim[0]-1,0] = xdata1 ; y[0:9,0] = ydata1 (first data)
; x[0:19,1] = xdata2 ; y[0:xyDim[1]-1,1] = ydata2 (second data)
; plot_anything
;
;03/11/11 Added variable pixel numbers in wavelength direction (wavepix_no) and in 
;bining direction (binpix_no) to account for change in camera type (CCD dimensions).
;Both variables were added to the common block Data.
;
;
;******************************************************************************
; 		HOW TO ADD A NEW FUNCTION TO THIS PROGRAM
;******************************************************************************
; Follow the steps below to add some functionality to this program.  Don't
; change the structure of the code unless you really know what you are doing.
; Send me mail if you have some comments, critique, or sugestions:
; blip@psfc.mit.edu
;
; 1) Insert your code (in the form of a procedure) somewhere before the Main
;    program, which is near the end of the file.
;    Plotting should be transparent to you. If you want to be able to print,
;    you have to include the common block "curplot". Examine the procedure
;    "Time_history" for an example.

; 2) Now you have to add a button to the widget so that your procedure can
;    get called. To do this, add  the following line to the code below the
;    line that reads: w2=widget_base(/lcol........
;	my_button = widget_button(w2, value='Do it')
;
; 3) Change EVERY OCCURENCE of the common block "widgets" by appending the
;    variable "my_button" to it.
;
; 4) In the procedure "OMA_VIEW_EVENT" find the clause of the first case statement,
;    the one that reads 'WIDGET_BUTTON'.
;
; 5) Below the line that reads: 'Print': Print_it, you have to add a line which
;    should look like this:
;	'Do it': my_procedure
;    Here 'Do it' is the value specified when you created the button in step 2.
;    My_procedure is the name of the procedure that you appended to the code and
;    want to be executed when the user presses the button 'Do it'.
;
; 6) You are done now.  Because of a minor bug in IDL you may have to exit IDL before
;    you can successfully recompile the entire program.  When your procedure crashes
;    it is also a good idea to exit and re-enter IDL before running it again.
;******************************************************************************
; This is a compound widget that displays a message and requests a Yes/No
; type of response from the user.  If the keyword ALERT is set, then only an
; OK-button is displayed to acknowledge the message.
;
pro msg,str
 common widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, w3c, w4, $
  w4c, draw, shot_label, shot_text, shot_button, frame_slider, $
  frameplot_button, lockyrange_button, timehist_button, movie_button, $
  from_label, from_text, from_label2, linewavelength_button, $
  to_label, to_text, to_label2, $
  xclick_label1, xclick_text, xclick_label2, yclick_text, $
  print_button, done_button, plot_radio, file_button, zoomout_button, $
  bin_slider, cal_button, savecal_button, delcal_button, spatial_button, $
  two_shot, current_button, curvie, di_box, plotlog_button

 widget_control, di_box, /append, set_value=str
 print,str
end

FUNCTION yes, message, ALERT=alert
 common return_value, return_value


 base = widget_base(/column, title='User Alert', xoffset=200, yoffset=200)
 w1 = widget_base(base, /row)
 message_text = widget_text(w1, value=message, xsize=strlen(message), ysize=1)
 w2 = widget_base(base, /row, space=25)
 if keyword_set(alert) then ok_button = widget_button(w2, value='OK') $
 else begin
   yes_button = widget_button(w2, value='Yes')
   no_button = widget_button(w2, value='No')
 endelse
 
 widget_control, base, /realize
 xmanager, 'alertbox', base, event_handler='yes_event', /modal
 
 return, return_value
END


;******************************************************************************
; This is the event handler for yes.
;
PRO yes_event, ev
 common return_value, return_value

 type = tag_names(ev, /structure)
 widget_control, ev.id, get_value=button_name

 if (type eq 'WIDGET_BUTTON') then begin
   widget_control, /destroy, ev.top
   if (button_name eq 'Yes') then return_value=-1  else return_value=0
 endif
END


;******************************************************************************
; This is a widget that emulates a Mac-style dialog box. It asks the user for
; input and returns the reply as a string.
;
FUNCTION ask, message
 common reply, reply, ok_button

 base = widget_base(/column, title='User Dialogue', xoffset=200, yoffset=200)
 w1 = widget_base(base, /row)
 message_text = widget_text(w1, value=message, xsize=strlen(message), ysize=1)
 reply_text = widget_text(w1, value='', xsize=32, ysize=1, /editable)
 w2 = widget_base(base, /column)
 ok_button = widget_button(w2, value='OK', uvalue=reply_text)

 widget_control, base, /realize
 widget_control, reply_text, /input_focus
 xmanager, 'ask', base, event_handler='ask_event', /modal

 return, reply
END


;******************************************************************************
; This is the event handler for ask.
;
PRO ask_event, ev
 common reply, reply, ok_button

 type = tag_names(ev, /structure)
 widget_control, ev.id, get_value=button_name

 if (type eq 'WIDGET_TEXT_CH') then begin
   widget_control, ok_button, get_uvalue=reply_text
   widget_control, reply_text, get_value=reply
   widget_control, /destroy, ev.top
 endif

 if (type eq 'WIDGET_BUTTON') then begin
   widget_control, ev.id, get_uvalue=reply_text
   widget_control, reply_text, get_value=reply
   widget_control, /destroy, ev.top
 endif

END


;*******************************************************************************
; Here we get the current wavelength boundaries for plotting by reading the 'To'
; and 'From' text fields.
;
PRO get_to_from, wto, wfrom, wto_i, wfrom_i
 common widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, w3c, w4, $
  w4c, draw, shot_label, shot_text, shot_button, frame_slider, $
  frameplot_button, lockyrange_button, timehist_button, movie_button, $
  from_label, from_text, from_label2, linewavelength_button, $
  to_label, to_text, to_label2, $
  xclick_label1, xclick_text, xclick_label2, yclick_text, $
  print_button, done_button, plot_radio, file_button, zoomout_button, $
  bin_slider, cal_button, savecal_button, delcal_button, spatial_button, $
  two_shot, current_button, curvie, di_box, plotlog_button
 common data, shotmarker, status_ok, data_ok, shot, frame, intensity, lambda, times, dt, $
             bin_numb, view_arr, counts, raw_data, bin, the_yrange, plottype, firstpix, $
             wavepix_no, binpix_no
 common overplot, shot_text2, frame_slider2, bin_slider2, Overplot_button, compshot,$
       twotimefrom, twotimeto, utility_button, secondScale

 widget_control, from_text, get_uvalue=wfrom
 widget_control, to_text, get_uvalue=wto

 if shotmarker eq 0 then begin
  widget_control, frame_slider, get_value=bing & frame[0]=bing
  widget_control, bin_slider, get_value=bong & bin[0] = bong
 endif else begin
  widget_control, frame_slider2, get_value=bing & frame[1]=bing
  widget_control, bin_slider2, get_value=bong & bin[1] = bong
 endelse

 widget_control, shot_label, get_uvalue = sizecount
  ;wto = strtrim(wto, 2)
  ;wfrom = strtrim(wfrom, 2)

 ;If the fields are empty then default to the full plot.
 if wfrom[0] eq '' then begin
  wfrom_i = firstpix
 endif else begin
  wfrom = float(wfrom[0])
  wfrom_i = where(wfrom le lambda[*,bin[0]-1,0] )
  if wfrom_i[0] eq -1 then wfrom_i = firstpix
  wfrom_i = min([wavepix_no-1, wfrom_i[*]])
 endelse

 if wto[0] eq '' then begin
  wto_i = firstpix+sizecount -1
 endif else begin
  wto = float(wto[0])
  wto_i = where(lambda[*,bin[0]-1,0] le wto)
  if (wto_i(0) eq -1) then wto_i = firstpix+sizecount-1
  wto_i = max([0, wto_i(*)])
 endelse

 
 if abs(wto_i - wfrom_i) le 3 then begin
  wfrom_i = firstpix
  wto_i = firstpix+sizecount - 1
 endif

 wfrom = lambda[wfrom_i,bin[0]-1,0]
 wto = lambda[wto_i,bin[0]-1,0]
 if wfrom gt wto then begin
  temp = wfrom
  wfrom = wto
  wto = temp
  temp = wfrom_i
  wfrom_i = wto_i
  wto_i = temp
 endif

 widget_control, from_text,  set_uvalue = wfrom,  set_value = string(form='(f7.2)', wfrom)
 widget_control, to_text,    set_uvalue = wto,    set_value = string(form='(f7.2)', wto)
 widget_control, from_label, set_uvalue = wfrom_i
 widget_control, to_label,   set_uvalue = wto_i
END


;******************************************************************************
; ZOOM_OUT
;
PRO zoom_out
 common widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, w3c, w4, $
  w4c, draw, shot_label, shot_text, shot_button, frame_slider, $
  frameplot_button, lockyrange_button, timehist_button, movie_button, $
  from_label, from_text, from_label2, linewavelength_button, $
  to_label, to_text, to_label2, $
  xclick_label1, xclick_text, xclick_label2, yclick_text, $
  print_button, done_button, plot_radio, file_button, zoomout_button, $
  bin_slider, cal_button, savecal_button, delcal_button, spatial_button, $
  two_shot, current_button, curvie, di_box, plotlog_button

 msg,' zoom out'
 widget_control, from_text, set_uvalue=0, set_value=''
 widget_control, to_text,   set_uvalue=0, set_value=''
 plot_frame
END

;*******************************************************************************
; Returns the wavelength calibration based on the coefficients stored in the tree.
;
PRO get_lambda
 common data, shotmarker, status_ok, data_ok, shot, frame, intensity, lambda, times, dt, $
             bin_numb, view_arr, counts, raw_data, bin, the_yrange, plottype, firstpix, $
             wavepix_no, binpix_no

 ;get the wavelength info right from the tree
 mdsopen, 'SPECTROSCOPY', shot[shotmarker]
 xlam = mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.CALIBRATION:LAMBDA_TRUE', /quiet, status = data_ok)
 cfg_bin_no = mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.SETTINGS:BIN_NO', /quiet)

 bin_numb[shotmarker]=cfg_bin_no

 if data_ok ne 1 or sum(finite(xlam) ne 1) ne 0 then begin ; estimate wavelength using simple linear dispersion model for grating
  msg,' exact wavelength calibration not available, using estimate'
  l_chr = mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.SETTINGS:LAMBDA_CTR') ;[nm]
  gr = mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.SETTINGS:GRATING') ; [grooves/mm]
  case gr of
	600  : disp = .06*600/gr ; [nm/pixel]
	1200 : disp = .065*600/gr ; [nm/pixel]
	1800 : disp = .0625*600/gr ; [nm/pixel]
	else : disp = .06*600/gr
  endcase
  msg, '  dispersion for '+num2str(gr,1)+'/mm grating is approximately '+num2str(disp,dp=4)+'nm/pixel'
  off = 0
  for k=0,bin_numb[shotmarker]-1,1 do begin
   lambda[*,k,shotmarker] = (findgen(wavepix_no) - floor(wavepix_no/2.))*disp + off + l_chr
  endfor
 endif else begin
  sizexlam = size(xlam,/dim)
  lambda[0:sizexlam[0]-1,0:sizexlam[1]-1,shotmarker] = xlam
 endelse

mdsclose
return
END

;*******************************************************************************
; Returns the time calibration based on the information stored in the tree.
;
PRO get_times
 common data, shotmarker, status_ok, data_ok, shot, frame, intensity, lambda, times, dt, $
             bin_numb, view_arr, counts, raw_data, bin, the_yrange, plottype, firstpix, $
             wavepix_no, binpix_no

 if n_elements(times) eq 0 then times = fltarr(1000,2)
 if n_elements(dt) eq 0 then dt = fltarr(2)

 ;get the time info right from the tree
 mdsopen, 'spectroscopy', shot[shotmarker]
 dummy = mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.CAMAC:DELTA_TIME', /quiet, status = data_ok)

 ;if data_ok then t_start=mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.CAMAC:TIME_START')

 ;**** DRM, 14 Feb 2012 - the (dummy eq '') was added because a
 ;                        previous version of TIMEBASE_CONSTRUCTION.RO
 ;                        entered data into the tree nodes as '' if no
 ;                        data was availble or if it was bad and this
 ;                        caused the status check in MDSVALUE to record
 ;                        it as okay (i.e., odd)
 if (not data_ok) OR (dummy eq '') then begin
  ;estimate time base
  dummy = mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.SETTINGS:NOM_EXP_TIME', /quiet, status = data_ok)
  if data_ok then dt[shotmarker] = dummy else dt[shotmarker] = 0.030
  t_start = -.09
  msg,' exact time base not available, using estimate'
 endif else begin
  dt[shotmarker] = dummy
  t_start=mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.CAMAC:TIME_START')
 endelse

 times[*,shotmarker] = t_start + dt[shotmarker] * findgen(1000)

 mdsclose
END

;*******************************************************************************
; This procedure retrieves the raw data from the tree.
; 
PRO get_raw_data
 common widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, w3c, w4, $
  w4c, draw, shot_label, shot_text, shot_button, frame_slider, $
  frameplot_button, lockyrange_button, timehist_button, movie_button, $
  from_label, from_text, from_label2, linewavelength_button, $
  to_label, to_text, to_label2, $
  xclick_label1, xclick_text, xclick_label2, yclick_text, $
  print_button, done_button, plot_radio, file_button, zoomout_button, $
  bin_slider, cal_button, savecal_button, delcal_button, spatial_button, $
  two_shot, current_button, curvie, di_box, plotlog_button
 common data, shotmarker, status_ok, data_ok, shot, frame, intensity, lambda, times, dt, $
             bin_numb, view_arr, counts, raw_data, bin, the_yrange, plottype, firstpix, $
             wavepix_no, binpix_no
 common overplot, shot_text2, frame_slider2, bin_slider2, Overplot_button, compshot,$
       twotimefrom, twotimeto, utility_button, secondScale 

 if n_elements(counts) eq 0 then counts = fltarr(1,1,1,2)
 if n_elements(view_arr) eq 0 then view_arr = strarr(25,2)
 if n_elements(frame) eq 0 then frame = long([0,0])
 if n_elements(bin) eq 0 then bin = long([0,0])
 mdsopen, 'spectroscopy', shot[shotmarker]
 xcounts = mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.ANALYSIS:RAW_SPECTRA', /quiet, status=data_ok)
 ycounts = counts

 if data_ok then begin
  cfg_periscop = mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.SETTINGS:CFG_PERISCOP', /quiet)
  cfg_bin_no = mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.SETTINGS:BIN_NO', /quiet)
  cfg_per_fibr = mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.SETTINGS:CFG_PER_FIBR', /quiet)
  slitWidth = mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.SETTINGS.SLIT_WIDTH',/quiet)
  grating = mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.SETTINGS.GRATING',/quiet)
  dummyx = n_elements(cfg_periscop)
  ;help,view_arr,cfg_periscop,cfg_per_fibr
  view_arr[0:dummyx-1,shotmarker] = strtrim(cfg_periscop,2) + ' ' + strtrim(string(cfg_per_fibr),2)

  msg,' '+strtrim(string(slitWidth),2)+'e-6m slit, '+strtrim(string(grating),2)+'/mm grating'
 endif else begin
  msg,'raw data is not available'
  return
 endelse

 mdsclose

 sizestr = size(xcounts)
 nframes = fix(sizestr[1])
 sizer   = size(ycounts)

 if not compshot then widget_control, shot_label, set_uvalue = sizestr[2]

 if sizestr[2] ne wavepix_no then begin
    dummy = ask("Not the full array.  What is the left pixel?")
    firstpix = int(dummy[0])
 endif else firstpix = 0 
 ;msg, ' firstpix='+string(firstpix)

 ;;NOTE: 02/08/12 DRM - This method of updating the counts will keep
 ;;                     the frame dimension as the largest number of
 ;;                     frames of any shot viewed during the current session
  counts = fltarr(max([sizestr[1], sizer[1]]), $
                  max([sizestr[2], sizer[2]]), $
                  max([sizestr[3], sizer[3]]), 2)

 ;set the new maximum number of frames for the frame slider widget
 if shotmarker eq 0 then begin
  widget_control, frame_slider, get_value=bing & frame[0] = min([bing, nframes])
  widget_control, bin_slider,   get_value=bong & bin[0] = bong
  widget_control, frame_slider, set_slider_max=nframes
  widget_control, bin_slider,   set_slider_max=bin_numb[0]
  counts[0:sizestr(1)-1,0:sizestr(2)-1,0:sizestr(3)-1,0] = xcounts
  counts[0:sizer(1)-1  ,0:sizer(2)-1,  0:sizer(3)-1,  1] = ycounts[*,*,*,1]
 endif else begin
  bin_numb[1]=cfg_bin_no
  widget_control, frame_slider2, get_value=bling & frame[1] = min([bling, nframes])
  widget_control, bin_slider2,   get_value=blong & bin[1] = blong
  widget_control, frame_slider2, set_slider_max=nframes
  widget_control, bin_slider2,   set_slider_max=bin_numb[1]
  counts[0:sizestr(1)-1,0:sizestr(2)-1,0:sizestr(3)-1,1] = xcounts
  counts[0:sizer(1)-1,  0:sizer(2)-1,  0:sizer(3)-1,  0] = ycounts[*,*,*,0]
 endelse 

 ;Restore curvefitting variables.  
 if slitWidth gt 26 then begin
  restore, file = '/usr/local/cmod/codes/spectroscopy/chromex/he_gaspuff/binfactors_km_50m_486nm.dat';, /verb
 endif else restore, file = '/usr/local/cmod/codes/spectroscopy/chromex/he_gaspuff/binfactors_blip_11_20_03.dat';, /verb
 restore, file = '/usr/local/cmod/codes/spectroscopy/chromex/position_periscopes.dat';, /verb

 zoom_out
END


;*******************************************************************************
; Let's retrieve a shot from the tree.  The success is reflected in the common
; block variables "status_ok" and "data_ok".
;
PRO get_shot
 common widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, w3c, w4, $
  w4c, draw, shot_label, shot_text, shot_button, frame_slider, $
  frameplot_button, lockyrange_button, timehist_button, movie_button, $
  from_label, from_text, from_label2, linewavelength_button, $
  to_label, to_text, to_label2, $
  xclick_label1, xclick_text, xclick_label2, yclick_text, $
  print_button, done_button, plot_radio, file_button, zoomout_button, $
  bin_slider, cal_button, savecal_button, delcal_button, spatial_button, $
  two_shot, current_button, curvie, di_box, plotlog_button
 common data, shotmarker, status_ok, data_ok, shot, frame, intensity, lambda, times, dt, $
             bin_numb, view_arr, counts, raw_data, bin, the_yrange, plottype, firstpix, $
             wavepix_no, binpix_no
 common overplot, shot_text2, frame_slider2, bin_slider2, Overplot_button, compshot,$
       twotimefrom, twotimeto, utility_button, secondScale 

 data_ok = 0
 status_ok = 0
 if n_elements(counts) eq 0	then counts = fltarr(1,1,1,2)
 if n_elements(raw_data) eq 0	then raw_data = [1,1]
 if n_elements(view_arr) eq 0	then view_arr = strarr(25,2)
 if n_elements(frame) eq 0	then frame = long([0,0])
 if n_elements(bin) eq 0	then bin = long([0,0])
 if n_elements(shot) eq 0	then shot = [1010101001,1010101001]

 ;If no shot number entered, then use current shot
 if shotmarker eq 0 then begin
  widget_control, shot_text, get_value=shot_str
  if shot_str eq '' then begin
    shot[0] = mdscur_shot("cmod")
    shot_str = string(form='(i10)', shot[0])
    widget_control, shot_text, set_value=shot_str
  endif else begin
    shot[0] = long(shot_str[0])
  endelse
 endif else begin
  widget_control, shot_text2, get_value=shot_str2
  shot[shotmarker] = long(shot_str2)
 endelse
 msg, 'Retrieving shot '+strtrim(string(shot[0]),2)

 widget_control, /hourglass
 mdsopen, 'spectroscopy', shot[shotmarker], /quiet, status=status_ok
 mdsclose

 if not status_ok then begin
  msg,'could not access the shot requested. enter another shot number'
  return
 endif

 get_times
 get_lambda

 ; ANJ 20120619: don't even ask for now since we haven't calibrated...
 ;raw_data[shotmarker] = yes('Do you want to look at raw data? (no -> try to load intensity data)')
 raw_data[shotmarker] = 1

 if raw_data[shotmarker] then begin
  get_raw_data
  raw_data[shotmarker] = 1
  ;if not data_ok then sens=0 else sens=1
  sens=data_ok

  widget_control, frameplot_button,	sensitive=sens
  widget_control, zoomout_button,	sensitive=sens
  widget_control, timehist_button,	sensitive=sens
  widget_control, spatial_button,	sensitive=sens
  widget_control, movie_button,	sensitive=sens
  widget_control, print_button,	sensitive=sens
  widget_control, file_button,	sensitive=sens

  return
 endif 

  ; the following part of procedure will only be run if one tries to get intensity data
  ; wavelength-calibration and time-base are already installed
  if n_elements(intensity) eq 0 then intensity = fltarr(1,1,1,2)
  mdsopen, 'spectroscopy', shot[shotmarker], /quiet, status=status_ok
  xintensity = mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.ANALYSIS:ANAL_SPECTRA', /quiet, status=data_ok)
  yintensity = intensity

  if data_ok then begin
   raw_data[shotmarker] = 0
   xcounts = mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.ANALYSIS:RAW_SPECTRA')
   sizestr = size(xcounts)
   nframes = fix(sizestr[1])
   ycounts = counts
   sizer = size(ycounts)
   intensize = size(xintensity)
   intensizer = size(yintensity)
   counts = fltarr(max([sizestr[1], sizer[1]]), $
                   max([sizestr[2], sizer[2]]), $
                   max([sizestr[3], sizer[3]]), 2)
   intensity = fltarr(max([intensize[1], intensizer[1]]), $
                      max([intensize[2], intensizer[2]]), $
                      max([intensize[3], intensizer[3]]), 2)
   cfg_periscop = mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.SETTINGS:CFG_PERISCOP', /quiet)
   cfg_bin_no = mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.SETTINGS:BIN_NO', /quiet)
   cfg_per_fibr = mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.SETTINGS:CFG_PER_FIBR', /quiet)
   dummyx = n_elements(cfg_periscop)
   view_arr[0:dummyx-1,shotmarker] = strtrim(cfg_periscop,2) + ' ' +  strtrim(string(cfg_per_fibr),2) 
   if shotmarker eq 0 then begin
    widget_control, frame_slider, get_value=bing & frame[0] = min([bing, nframes])
    widget_control, bin_slider,   get_value=bong & bin[0] = bong
    widget_control, frame_slider, set_slider_max=nframes
    widget_control, bin_slider,   set_slider_max=bin_numb[0]
    counts[0:sizestr(1)-1,0:sizestr(2)-1,0:sizestr(3)-1,0] = xcounts
    counts[0:sizer(1)-1,  0:sizer(2)-1,  0:sizer(3)-1,  1] = ycounts[*,*,*,1]
    intensity[0:intensize(1)-1, 0:intensize(2)-1, 0:intensize(3)-1, 0] = xintensity
    intensity[0:intensizer(1)-1,0:intensizer(2)-1,0:intensizer(3)-1,1] = yintensity[*,*,*,1]
   endif else begin
    bin_numb[1]=cfg_bin_no
    widget_control, frame_slider2, get_value=bling & frame[1] = min([bling, nframes])
    widget_control, bin_slider2,   get_value=blong & bin[1]=blong
    widget_control, frame_slider2, set_slider_max=nframes
    widget_control, bin_slider2,   set_slider_max=bin_numb[1]
    counts[0:sizestr(1)-1,0:sizestr(2)-1,0:sizestr(3)-1,1] = xcounts
    counts[0:sizer(1)-1,  0:sizer(2)-1,  0:sizer(3)-1,  0] = ycounts[*,*,*,0]
    intensity[0:intensize(1)-1, 0:intensize(2)-1, 0:intensize(3)-1, 1] = xintensity
    intensity[0:intensizer(1)-1,0:intensizer(2)-1,0:intensizer(3)-1,0] = yintensity[*,*,*,0]
   endelse

   zoom_out
  endif

 mdsclose

 if not data_ok then begin
  raw_data[shotmarker] = yes('Could not retrieve intensity data.  Do you want to use the raw data instead?')
  if not raw_data[shotmarker] then begin
    widget_control, shot_text, set_value=''
    widget_control, frameplot_button,	sensitive=0
    widget_control, zoomout_button,	sensitive=0
    widget_control, timehist_button,	sensitive=0
    widget_control, spatial_button,	sensitive=0
    widget_control, movie_button,	sensitive=0
    widget_control, print_button,	sensitive=0
    widget_control, file_button,	sensitive=0
  endif else begin
    get_raw_data
  endelse
 endif

 if data_ok then sens=1 else sens=0
 widget_control, frameplot_button,	sensitive=sens
 widget_control, zoomout_button,	sensitive=sens
 widget_control, timehist_button,	sensitive=sens
 widget_control, spatial_button,	sensitive=sens
 widget_control, movie_button,	sensitive=sens
 widget_control, print_button,	sensitive=sens
 widget_control, file_button,	sensitive=sens

 return
END

;*******************************************************************************
;
Pro Plot_frame
 common widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, w3c, w4, $
  w4c, draw, shot_label, shot_text, shot_button, frame_slider, $
  frameplot_button, lockyrange_button, timehist_button, movie_button, $
  from_label, from_text, from_label2, linewavelength_button, $
  to_label, to_text, to_label2, $
  xclick_label1, xclick_text, xclick_label2, yclick_text, $
  print_button, done_button, plot_radio, file_button, zoomout_button, $
  bin_slider, cal_button, savecal_button, delcal_button, spatial_button, $
  two_shot, current_button, curvie, di_box, plotlog_button
 common data, shotmarker, status_ok, data_ok, shot, frame, intensity, lambda, times, dt, $
             bin_numb, view_arr, counts, raw_data, bin, the_yrange, plottype, firstpix, $
             wavepix_no, binpix_no
 common curplot, x, y, xyDim, tit, xtit, ytit, xannot, yannot, annot, annotNum, $ 
                this_yrange, print_queue_num, view, num_of_plots, plotColor, plotColorFlag, $
                plotSymbol, plotSymbolFlag
 common overplot, shot_text2, frame_slider2, bin_slider2, Overplot_button, compshot,$
       twotimefrom, twotimeto, utility_button, secondScale
 common wchr_view_event, from_flag, format_str

 shotmark = shotmarker
 ;Run this twice for plotting asthetics.  Trust me.
 for ib=0, 1 do begin
  shotmarker = 0
  get_to_from, wto, wfrom, wto_i, wfrom_i
  from_flag = 0
 endfor
 shotmarker = shotmark

 plotSymbolFlag = 0
 plotColorFlag = 0


 tit = 'CHROMEX: Discharge ' + string(form='(i10, 2x)', shot[0]) + $
  '   Bin: ' + strtrim(string(bin[0]),2) + $
  '   [' + strtrim(string(form='(f6.3)',times[frame[0]-1,0]),2) + ' s, ' + $ 
  strtrim(string(form='(f6.3)', times[frame[0]-1,0]+dt[0]),2) +' s]'

 if not raw_data[0] then meas=intensity else meas=counts
 x1 = lambda[wfrom_i : wto_i,bin[0]-1, 0]
 y1 = meas[frame[0]-1, wfrom_i - firstpix : wto_i - firstpix, bin[0]-1, 0]
 ;the absolute intensity calibration is given PER PIXEL, but it will be plotted PER MICRON
 if not raw_data[0] then begin
  ctr = long( (wto_i + wfrom_i)/2. )
  pixel_per_micron = float(1.) / $
                     (0.001*abs( lambda[ctr+1,bin[0]-1,0] - lambda[ctr,bin[0]-1,0] ) )
  y1 = y1*pixel_per_micron[0]
 endif

 y1 = reform(y1,n_elements(y1),1)
 xyDim[0]=n_elements(y1)
 x[0:n_elements(x1)-1,0] = x1
 y[0:n_elements(x1)-1,0] = y1
 num_of_plots = 1

 if not raw_data[0] then ytit = 'mW /(cm^2 sterad 10^-6 m)' $
                    else ytit = 'Counts'
 annotNum = 0
 if bin[0]-1 GT n_elements(view_arr[*,0])-1 $
  then view = '' $
  else view = view_arr[bin[0]-1,0]

 ;Log or Linear Scales that seem to work best:
 ayrange =  [min(y1)*.8,max(y1)+0.1*(max(y1)-min(y1))]
 byrange =  [min(y1)-0.1*(max(y1)-min(y1)),max(y1)+0.1*(max(y1)-min(y1))]

 widget_control, xclick_label2, set_value='nm'

 format_str = '(f7.2)'
 if n_elements(the_yrange) ne 2 $
  then if plottype $
        then this_yrange = ayrange $
        else this_yrange = byrange $
  else this_yrange = the_yrange


 ;If the compare shots is active this loads them and stores them in to the second slot of the x y data
 if compshot then begin
  if not raw_data[1] then neas=intensity $
                     else neas=counts
  widget_control, twotimefrom, set_value=strtrim(string(times[frame[1]-1,1]),2)
  widget_control, twotimeto,   set_value=strtrim(string(times[frame[1]-1,1]+dt[1]),2)
  ;Hide opy in slot 3, but don't tell plot_anything its there.  For log/lin button
  opy = neas[frame[1]-1, wfrom_i : wto_i, bin[1]-1, 1]
  num_of_plots = 2
  xyDim[1] = n_elements(y1)
  xyDim[3] = n_elements(y1)
  y[0:xyDim[3]-1,3] = opy   
  if plottype then smallscale = max(byrange) - min(byrange) $
              else smallscale = max(this_yrange) - min(this_yrange)
  bigscale = (max(opy) - min(opy) + 0.2*(max(opy)-min(opy)))
  bigmin = (min(opy)-.1*(max(opy)-min(opy)))/bigscale*smallscale
  shift = bigmin - min(byrange)
  opy1 = opy/bigscale*smallscale - shift
  x2 = lambda[wfrom_i : wto_i,bin[1]-1,1]
  x[0:xyDim[1]-1,1] = x1
  y[0:xyDim[1]-1,1] = opy1
  plotColorFlag = 1
  plotColor[0] = !p.color
  plotColor[1] = 100
  ysty = 9
  ;plotSymbol[2] = 1
  if plottype then secondScale = [min(opy)*.8, max(opy) + .1*(max(opy)-min(opy))] $
              else secondScale = [min(opy)-.1*(max(opy)-min(opy)), max(opy) + .1*(max(opy)-min(opy))]
 endif

 xtit = 'nm'
 if not raw_data[1] then ytit1 = 'mW /(cm^2 sterad 10^-6 m)' $
                    else ytit1 = 'Counts'

 Plot_anything
end

;*******************************************************************************
;
PRO Plot_anything
 common widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, w3c, w4, $
  w4c, draw, shot_label, shot_text, shot_button, frame_slider, $
  frameplot_button, lockyrange_button, timehist_button, movie_button, $
  from_label, from_text, from_label2, linewavelength_button, $
  to_label, to_text, to_label2, $
  xclick_label1, xclick_text, xclick_label2, yclick_text, $
  print_button, done_button, plot_radio, file_button, zoomout_button, $
  bin_slider, cal_button, savecal_button, delcal_button, spatial_button, $
  two_shot, current_button, curvie, di_box, plotlog_button
 common data, shotmarker, status_ok, data_ok, shot, frame, intensity, lambda, times, dt, $
        bin_numb, view_arr, counts, raw_data, bin, the_yrange, plottype, firstpix, $
        wavepix_no, binpix_no
 common curplot, x, y, xyDim, tit, xtit, ytit, xannot, yannot, annot, annotNum, $ 
        this_yrange, print_queue_num, view, num_of_plots, plotColor, plotColorFlag, $
        plotSymbol, plotSymbolFlag
 common wchr_view_event, from_flag, format_str
 common line_data, line_wvlngths, line_idents, line_plot
 common overplot, shot_text2, frame_slider2, bin_slider2, Overplot_button, compshot,$
        twotimefrom, twotimeto, utility_button, secondScale 

 if compshot then ysty = 9 else  ysty = 1
 if plotColorFlag eq 1 then color = plotColor[0] ; 1 = yes color, Not 1 equals no color.  Isn't this fun?
 if plotSymbolFlag eq 1 then begin  ; 1 = yes symbols, Not 1 equals no symbols
    psym = plotSymbol[0] 
    syms = 2
 endif

 if plottype then begin
  if this_yrange[0] lt 1 then this_yrange[0]= 1
  plot_io, x[0:xyDim[0]-1,0], y[0:xyDim[0]-1,0], tit=tit, color=color, psym=psym, syms=syms, $
           xtit=xtit,         ytit=ytit, $
           xmargin=[10,10],   yrange=this_yrange, $
           xstyle=1,          ystyle=ysty
 endif else begin
  plot, x[0:xyDim[0]-1,0], y[0:xyDim[0]-1,0], tit=tit, color=color, psym=psym, syms=syms, $
        xtit=xtit,         ytit=ytit, $
        xmargin=[10,10],   yrange=this_yrange, $
        xstyle=1,          ystyle=ysty
 endelse

 ;Overplot the upper slots of data
 for p = 1, num_of_plots-1 do begin
  if plotSymbolFlag eq 1 then psym = plotSymbol[p]
  if plotColorFlag eq 1 then color = plotColor[p] 
  oplot, x[0:xyDim[p]-1,p], y[0:xyDim[p]-1,p], color = color, psym = psym, syms = syms
  if compshot then begin
   xyouts, .13, .85, view_arr[bin[1]-1,1],/norm, color= color
   axis, yaxis=1, yrange = secondScale, ystyle = 1, ytit = ytit1
  endif
 endfor 

 ;Annotate:
 for po = 0, annotNum-1 do xyouts, xannot[po], yannot[po], annot[po], /norm

 yfrac=.1
 if view ne '' then xyouts, .13, .87, view, /norm
 if n_elements(line_plot) ne 1 then line_plot=0
 if line_plot then plot_line_wv

END


;*******************************************************************************
; overplot known lines
PRO plot_line_wv
 common curplot, x, y, xyDim, tit, xtit, ytit, xannot, yannot, annot, annotNum, $ 
                 this_yrange, print_queue_num, view, num_of_plots, plotColor, plotColorFlag, $
                 plotSymbol, plotSymbolFlag
 common data, shotmarker, status_ok, data_ok, shot, frame, intensity, lambda, times, dt, $
              bin_numb, view_arr, counts, raw_data, bin, the_yrange, plottype, firstpix, $
              wavepix_no, binpix_no
 common wchr_view_event, from_flag, format_str
 common line_data, line_wvlngths, line_idents, line_plot

 y1 = y[0:xyDim[0]-1,0]
 x1 = x[0:xyDim[0]-1,0]

 if n_elements(the_yrange) ne 2 then begin
  if plottype then this_yrange = [min(y1)*.5, max(y1)-0.3*(max(y1)-min(y1))] $
              else this_yrange = [min(y1)-0.1*(max(y1)-min(y1)),max(y1)+0.1*(max(y1)-min(y1))]
 endif else this_yrange= the_yrange

 wline=where(line_wvlngths gt min(x1) and line_wvlngths lt max(x1))
 if wline[0] ne -1 then begin
  yfrac=.2
  for k=0,n_elements(wline)-1 do begin
   oplot,[line_wvlngths[wline[k]],line_wvlngths[wline[k]]],$
         [this_yrange[0],         this_yrange[0]*yfrac+this_yrange[1]*(1.-yfrac)],$
         lines=1
   xyouts,line_wvlngths[wline[k]],$
          this_yrange(0)*yfrac+this_yrange(1)*(1.-yfrac),$
          line_idents(wline(k)),alignment=0.4,orientation=90.
  endfor 
 endif

END

;*******************************************************************************
;
PRO lock_yrange
 common widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, w3c, w4, $
  w4c, draw, shot_label, shot_text, shot_button, frame_slider, $
  frameplot_button, lockyrange_button, timehist_button, movie_button, $
  from_label, from_text, from_label2, linewavelength_button, $
  to_label, to_text, to_label2, $
  xclick_label1, xclick_text, xclick_label2, yclick_text, $
  print_button, done_button, plot_radio, file_button, zoomout_button, $
  bin_slider, cal_button, savecal_button, delcal_button, spatial_button, $
  two_shot, current_button, curvie, di_box, plotlog_button
 common data, shotmarker, status_ok, data_ok, shot, frame, intensity, lambda, times, dt, $
              bin_numb, view_arr, counts, raw_data, bin, the_yrange, plottype, firstpix, $
              wavepix_no, binpix_no
 common curplot, x, y, xyDim, tit, xtit, ytit, xannot, yannot, annot, annotNum, $ 
                 this_yrange, print_queue_num, view, num_of_plots, plotColor, plotColorFlag, $
                 plotSymbol, plotSymbolFlag
 common wchr_view_event, from_flag, format_str

 if n_elements(the_yrange) lt 2 then begin
  ; The yrange is not locked. Lock it by setting the_yrange to a
  ; 2-element array to the last main data plotted.  Will work on
  ; timehist and spatial as well.  Be careful.
  x1 = x[0:xyDim[0]-1,0]
  y1 = y[0:xyDim[0]-1,0]

  ;  if (not raw_data[0]) then begin
  ;    ctr = long( (wto_i + wfrom_i)/2. )
  ;    pixel_per_micron = float(1.) / $
  ;      ( 0.001*abs( lambda[ctr+1,bin[0]-1,0] - lambda[ctr,bin[0]-1,0] ) )
  ;    y1 = y1*pixel_per_micron(0)
  ;  endif    
  the_yrange=[min(y1)-0.1*(max(y1)-min(y1)),max(y1)+0.1*(max(y1)-min(y1))]
  if plottype then the_yrange = [min(y1)*.8,max(y1)+0.1*(max(y1)-min(y1))]
  widget_control, lockyrange_button, set_value='Unlock Y-Range'  
 endif else begin
  ; The yrange is currently locked. Unlock it.
  the_yrange=0.
  widget_control, lockyrange_button, set_value='Lock Y-Range'  
 endelse

END


;*******************************************************************************
;
PRO line_wavelength
common widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, w3c, w4, $
  w4c, draw, shot_label, shot_text, shot_button, frame_slider, $
  frameplot_button, lockyrange_button, timehist_button, movie_button, $
  from_label, from_text, from_label2, linewavelength_button, $
  to_label, to_text, to_label2, $
  xclick_label1, xclick_text, xclick_label2, yclick_text, $
  print_button, done_button, plot_radio, file_button, zoomout_button, $
  bin_slider, cal_button, savecal_button, delcal_button, spatial_button, $
  two_shot, current_button, curvie, di_box, plotlog_button
 common data, shotmarker, status_ok, data_ok, shot, frame, intensity, lambda, times, dt, $
              bin_numb, view_arr, counts, raw_data, bin, the_yrange, plottype, firstpix, $
              wavepix_no, binpix_no
 common curplot, x, y, xyDim, tit, xtit, ytit, xannot, yannot, annot, annotNum, $ 
                 this_yrange, print_queue_num, view, num_of_plots, plotColor, plotColorFlag, $
                 plotSymbol, plotSymbolFlag
 common wchr_view_event, from_flag, format_str
 common line_data, line_wvlngths, line_idents, line_plot

 ;Toggle the wavelength annotation of known lines.
 if n_elements(line_plot) ne 1 then line_plot=0

 if line_plot then begin
  line_plot=0 
 endif else begin 
  line_plot=1 ; toggle

  GENIE_PATH=getenv('GENIE_PATH') & if GENIE_PATH eq '' then GENIE_PATH='/usr/local/cmod/idl/GENIE/'
  wfile=getenv('VIS_LINE_LIST') & if wfile eq '' then wfile=GENIE_PATH+'IMPSPEC/vis_line_list.csv'
  csv_read,wfile & print,'loaded lines from '+wfile
  
  line_wvlngths=wave
  line_idents=elem+' '+cs
 endelse

 plot_anything
END


;*******************************************************************************
;
PRO Time_history
 common widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, w3c, w4, $
  w4c, draw, shot_label, shot_text, shot_button, frame_slider, $
  frameplot_button, lockyrange_button, timehist_button, movie_button, $
  from_label, from_text, from_label2, linewavelength_button, $
  to_label, to_text, to_label2, $
  xclick_label1, xclick_text, xclick_label2, yclick_text, $
  print_button, done_button, plot_radio, file_button, zoomout_button, $
  bin_slider, cal_button, savecal_button, delcal_button, spatial_button, $
  two_shot, current_button, curvie, di_box, plotlog_button
 common data, shotmarker, status_ok, data_ok, shot, frame, intensity, lambda, times, dt, $
              bin_numb, view_arr, counts, raw_data, bin, the_yrange, plottype, firstpix, $
              wavepix_no, binpix_no
 common curplot, x, y, xyDim, tit, xtit, ytit, xannot, yannot, annot, annotNum, $ 
                 this_yrange, print_queue_num, view, num_of_plots, plotColor, plotColorFlag, $
                 plotSymbol, plotSymbolFlag
 common wchr_view_event, from_flag, format_str
 common overplot, shot_text2, frame_slider2, bin_slider2, Overplot_button, compshot,$
        twotimefrom, twotimeto, utility_button, secondScale

 widget_control, plotlog_button, set_value='Plot Logarithmic'
 get_to_from, wto, wfrom, wto_i, wfrom_i

 ;help, wfrom, wto, wfrom_i, wto_i  ;;;;;;;;;;;;;;;;;;;;;
 region = findgen(wto_i+1 - wfrom_i) + wfrom_i - firstpix
 ;help, region
 ;print, max(region)
 ans = ''
 if compshot then ans = ask('White or Blue graph? (W/B)')
 if ans eq 'b' or ans eq 'B' or ans eq 'blue' or ans eq 'Blue' then shotmarker = 1 else shotmarker = 0

 if not raw_data[shotmarker] then meas=intensity $
                             else meas=counts

 ;;02/07/12 DRM - added to retrieve frame no. in order to scale x-axis correctly
 mdsopen, 'spectroscopy', shot[shotmarker]
 rawspect = mdsvalue('\spectroscopy::top.chromex.analysis:raw_spectra')
 mdsclose
 spect_dim  = size(rawspect)
 frame_no   = spect_dim(1)
 ;help,frame_no ;;;;;;;;;;;;;;;;

 ctr = fix( (wto_i + wfrom_i)/2. )
 nm_per_pixel = abs( lambda[ctr+1,bin[shotmarker]-1,shotmarker] - lambda[ctr,bin[shotmarker]-1,shotmarker] )
 ;area = total( meas[*, region, bin[shotmarker]-1, shotmarker], 2)
 ;;02/07/12 DRM changed length of first dimension in MEAS to scale x-axis correctly
 area = total( meas[0:frame_no-1, region, bin[shotmarker]-1, shotmarker], 2)
 ;help, area ;;;;;;;;;;;;;;;;;;
 ;help, meas  ;;;;;;;;;;;;;;;;;;
 ;help, counts ;;;;;;;;;;;;;;;;;
 ;help,bin ;;;;;;;;;;;;;;;;
 left_1 = max([0, wfrom_i-1])
 left_2 = wfrom_i
 left_3 = wfrom_i+1
 right_1 = wto_i-1
 right_2 = wto_i
 right_3 = min([wto_i+1, wavepix_no])
 ;backgnd = (total( meas(*, $
 ;          [left_1, left_2, left_3, right_1, right_2, right_3], bin-1), 2) ) * $
 ;          (wto_i-wfrom_i+1.)/6.
 ;;02/07/12 DRM changed length of first dimension in MEAS to scale x-axis correctly
 backgnd = (total( meas(0:frame_no-1, $
           [left_1, left_2, left_3, right_1, right_2, right_3], bin-1), 2) ) * $
           (wto_i-wfrom_i+1.)/6.
 area(*) = area(*) - backgnd(*)

 ;check for saturation
  ;saturated = where(counts[*,region, bin[shotmarker]-1,shotmarker] gt 60000.)
  ;;02/07/12 DRM changed length of first dimension in SATURATED to scale
  ;;x-axis correctly
 saturated = where(counts[0:frame_no-1,region, bin[shotmarker]-1,shotmarker] gt 60000.)
 saturated = size(saturated)
 if (saturated(0) ne 0) then begin
  annotNum = 1
  xannot[0] = .5
  yannot[0] = .9
  annot[0] = 'Saturated pixels in region'
 endif else annotNum = 0

 dims = n_elements(area)
 x1 = times[0:dims-1,shotmarker] + 0.5* (times[1,shotmarker]-times[0,shotmarker])
 y1 = area[0:dims-1]
 x[0:n_elements(x1)-1,0] = x1
 y[0:n_elements(x1)-1,0] = y1
 num_of_plots = 1
 xyDim[0] = n_elements(x1)
 tit = 'CHROMEX Time History:  Discharge ' + string(form='(i10, 2x)', shot[shotmarker]) + $
      '  Bin: ' + strtrim(string(bin[shotmarker]),2) + $      
      '  [' + strtrim(string(form='(f7.2)',lambda[wfrom_i,bin[shotmarker]-1,shotmarker]),2) + $
      ' nm, ' + strtrim(string(form='(f7.2)', lambda[wto_i, bin[shotmarker]-1,shotmarker]),2) +' nm]'
 xtit = 'seconds'
 if (not raw_data[shotmarker]) then ytit = 'mW / (cm^2 sterad)' else ytit='Counts'

 this_yrange = [0,max(y[0:xyDim[0]-1,0])*1.05]
 format_str = '(f5.2)'
 widget_control, xclick_label2, set_value='sec'
 view = view_arr[bin[shotmarker]-1,shotmarker]
 plotSymbolFlag = 0
 plotColorFlag = 0
 plottype = 0
 compy = compshot
 compshot = 0

 plot_anything
 compshot = compy
END


;*******************************************************************************
;
PRO spatial
;wfrom_i and wto_i have to be calculated for all bins separately:
; -> modification necessary
 common widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, w3c, w4, $
  w4c, draw, shot_label, shot_text, shot_button, frame_slider, $
  frameplot_button, lockyrange_button, timehist_button, movie_button, $
  from_label, from_text, from_label2, linewavelength_button, $
  to_label, to_text, to_label2, $
  xclick_label1, xclick_text, xclick_label2, yclick_text, $
  print_button, done_button, plot_radio, file_button, zoomout_button, $
  bin_slider, cal_button, savecal_button, delcal_button, spatial_button, $
  two_shot, current_button, curvie, di_box, plotlog_button
 common data, shotmarker, status_ok, data_ok, shot, frame, intensity, lambda, times, dt, $
              bin_numb, view_arr, counts, raw_data, bin, the_yrange, plottype, firstpix, $
              wavepix_no, binpix_no
 common curplot, x, y, xyDim, tit, xtit, ytit, xannot, yannot, annot, annotNum, $ 
                 this_yrange, print_queue_num, view, num_of_plots, plotColor, plotColorFlag, $
                 plotSymbol, plotSymbolFlag
 common wchr_view_event, from_flag, format_str
 common overplot, shot_text2, frame_slider2, bin_slider2, Overplot_button, compshot,$
        twotimefrom, twotimeto, utility_button, secondScale 

 widget_control, plotlog_button, set_value='Plot Logarithmic'
 get_to_from, wto, wfrom, wto_i, wfrom_i
 region = findgen(wto_i+1 - wfrom_i) + wfrom_i - firstpix

 ans = ''
 if compshot then ans = ask('White or Blue graph? (W/B)')
 if ans eq 'b' or ans eq 'B' or ans eq 'blue' or ans eq 'Blue' then shotmarker = 1 else shotmarker = 0

 if (not raw_data[shotmarker]) then meas=intensity else meas=counts

 ctr = long( (wto_i + wfrom_i)/2. )
 nm_per_pixel = abs( lambda[ctr+1,bin[shotmarker]-1,shotmarker] - lambda[ctr,bin[shotmarker]-1,shotmarker] )
 area = total(meas(frame[shotmarker]-1, region, *), 2)

 left_1 = max([0, wfrom_i-1])
 left_2 = wfrom_i
 left_3 = wfrom_i+1
 right_1 = wto_i-1
 right_2 = wto_i
 right_3 = min([wto_i+1, wavepix_no])
 backgnd = (total( meas(frame[shotmarker]-1, $
           [left_1, left_2, left_3, right_1, right_2, right_3], *), 2) ) * $
           (wto_i-wfrom_i+1.)/6.
 area(*) = area(*) - backgnd(*)


 dims = n_elements(area)
 x1 = findgen(dims)+1
 y1 = area(0:dims-1)
 x[0:n_elements(x1)-1,0] = x1
 y[0:n_elements(x1)-1,0] = y1
 num_of_plots = 1
 xyDim[0] = n_elements(x1)

 ctr = 1
 annotNum = binpix_no
 xannot[0:binpix_no-1] = .15
 for ken = 0, binpix_no-1 do begin
  yannot[ken] = .9 - (ken)/20.
  annot[ken] = 'Bin ' + strtrim(string(long(ken+1)),2) + ': ' + view_arr[ken,shotmarker]
 endfor

 ;check for saturation
 saturated = where(counts(frame[shotmarker]-1,region, *) gt 60000.)
 saturated = size(saturated)
 if (saturated(0) ne 0) then begin
  annotNum = 17
  xannot[binpix_no] = .5
  yannot[binpix_no] = .9
  annot[binpix_no] = 'Saturated pixels in region'
 endif 
 view = ''

 tit = 'CHROMEX Spatial Scan:  Discharge ' + string(form='(i10, 2x)', shot[0]) + $
      '  [' + strtrim(string(form='(f6.3)',times[frame[0]-1,0]),2) + ' s, ' + $
      strtrim(string(form='(f6.3)', times[frame[0],0]),2) +' s]' + $      
      '  [' + strtrim(string(form='(f7.2)',lambda(wfrom_i,bin[0]-1)),2) + $
      ' nm, ' + strtrim(string(form='(f7.2)', lambda(wto_i, bin[0]-1)),2) +' nm]'
 xtit = 'Bin'
 if not raw_data[shotmarker] then ytit = 'mW / (cm^2 sterad)' $
                             else ytit='Counts'

 format_str = '(f5.2)'
 widget_control, xclick_label2, set_value=' bin'

 this_yrange = [0,max(y[0:xyDim[0]-1,0])*1.05]
 plottype = 0
 plotSymbolFlag = 0
 plotColorFlag = 0
 compy = compshot
 compshot = 0

 plot_anything
 compshot = compy
END


;*******************************************************************************
;
PRO Movie
 common widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, w3c, w4, $
  w4c, draw, shot_label, shot_text, shot_button, frame_slider, $
  frameplot_button, lockyrange_button, timehist_button, movie_button, $
  from_label, from_text, from_label2, linewavelength_button, $
  to_label, to_text, to_label2, $
  xclick_label1, xclick_text, xclick_label2, yclick_text, $
  print_button, done_button, plot_radio, file_button, zoomout_button, $
  bin_slider, cal_button, savecal_button, delcal_button, spatial_button, $
  two_shot, current_button, curvie, di_box, plotlog_button
 common data, shotmarker, status_ok, data_ok, shot, frame, intensity, lambda, times, dt, $
              bin_numb, view_arr, counts, raw_data, bin, the_yrange, plottype, firstpix, $
              wavepix_no, binpix_no

 widget_control, shot_button, sensitive=0
 widget_control, bin_slider, sensitive=0
 widget_control, frame_slider, sensitive=0
 widget_control, frameplot_button, sensitive=0
 widget_control, zoomout_button, sensitive=0
 widget_control, timehist_button, sensitive=0
 widget_control, spatial_button, sensitive=0
 widget_control, lockyrange_button, sensitive=0
 widget_control, linewavelength_button, sensitive=0
 widget_control, file_button, sensitive=0
 widget_control, print_button, sensitive=0
 widget_control, done_button, sensitive=1, get_uvalue=boop
 widget_control, plotlog_button, sensitive=0
 widget_control, curvie, sens = 0
 if  boop eq 'on' then widget_control, two_shot, sens = 0
 widget_control, current_button, sens = 0

 widget_control, movie_button, set_value='Stop movie', xsize=160
 widget_control, xclick_label1, set_value='   Time:'
 widget_control, xclick_label2, set_value='sec'

 shotmarker=0
 if not raw_data[shotmarker] then meas=intensity else meas=counts

 get_to_from, wto, wfrom, wto_i, wfrom_i
 region = findgen(wto_i - wfrom_i) + wfrom_i
 maxmeas = max(meas(*,region, bin[0]-1))
 if not raw_data[shotmarker] then ytit='mW / (cm^2 sterad 10^-6 m)' else ytit='Counts'
 plot, lambda[region,bin[0]-1,0], meas[0, region, bin[0]-1], /nodata, $
  yrange=[0, maxmeas], xtit='nm', ytit=ytit, color=!p.color

 nframes = n_elements(meas(*,0,0))
 repeat begin
 for k = 0, nframes-1 do begin
    oplot, lambda[region,bin[0]-1,0], meas[k, region, bin[0]-1], color=!p.color
    outstr = string(format='(f7.2)', times[k,0]) + ' sec'
    xyouts, .75, .9, /normal, outstr
    widget_control, xclick_text, set_value=string(format='(f5.2)', times[k,0])
    widget_control, frame_slider, set_value=k+1
    event = widget_event(movie_button, /nowait)
    if (event.id ne 0) and event.select then goto, finished
    wait, .3
    oplot, lambda[region,bin[0]-1,0], meas[k, region, bin[0]-1], color=!p.background
        xyouts, .75, .9, /normal, outstr;, color=!p.background
        axis, xaxis=0, xrange = !x.crange
 endfor
 endrep until 0

finished:
 widget_control, shot_button, sensitive=1
 widget_control, bin_slider, sensitive=1
 widget_control, frame_slider, sensitive=1
 widget_control, frameplot_button, sensitive=1
 widget_control, zoomout_button, sensitive=1
 widget_control, timehist_button, sensitive=1
 widget_control, spatial_button, sensitive=1
 widget_control, lockyrange_button, sensitive=1
 widget_control, linewavelength_button, sensitive=1
 widget_control, file_button, sensitive=1
 widget_control, print_button, sensitive=1
 widget_control, done_button, sensitive=1, get_uvalue=boop
 widget_control, plotlog_button, sensitive=1
 widget_control, movie_button, sens = 1
 widget_control, curvie, sens = 1
 if  boop eq 'on' then widget_control, two_shot, sens = 1
 widget_control, current_button, sens = 1

 widget_control, movie_button, set_value='Movie', xsize=160
 widget_control, xclick_label1, set_value='Cursor:'
 widget_control, xclick_label2, set_value='nm'
 widget_control, xclick_text, set_value=''
END

;*******************************************************************************
;
PRO Print_it, FILE=filename
 common data, shotmarker, status_ok, data_ok, shot, frame, intensity, lambda, times, dt, $
              bin_numb, view_arr, counts, raw_data, bin, the_yrange, plottype, firstpix, $
              wavepix_no, binpix_no
 common curplot, x, y, xyDim, tit, xtit, ytit, xannot, yannot, annot, annotNum, $ 
                 this_yrange, print_queue_num, view, num_of_plots, plotColor, plotColorFlag, $
                 plotSymbol, plotSymbolFlag
 common line_data, line_wvlngths, line_idents, line_plot
 common overplot, shot_text2, frame_slider2, bin_slider2, Overplot_button, compshot,$
        twotimefrom, twotimeto, utility_button, secondScale 

 widget_control, /hourglass
 if keyword_set(filename) then begin
  xdata = x[0:xydim[0]-1,0];lambda[*,bin[0]-1,0]
  if not raw_data[0] then meas=intensity else meas=counts
  ydata = y[0:xyDim[0]-1,0];reform(meas[frame[0]-1,*,bin[0]-1,0])
  name = filename[0]
  viewingArray = view_arr[*,0]
  if annotNum gt 0 $
   then begin Annotate = annot[0:annotNum-1] & AnnotateXPos = xannot[0:annotNum-1] & AnnotateYPos = yannot[0:annotNum-1] & endif $
   else begin Annotate = ''                  & AnnotateXPos = .5                   & AnnotateYPos = .5 & endelse
 
  save, xdata, ydata, tit, xtit, ytit, Annotate, AnnotateXPos, AnnotateYPos, this_yrange, view, viewingArray, $
        file=name, /verb
 endif else begin

  tmp = ask('Which print queue number? [Default 114 or previous. Zero: no print, just leave idl.ps]')
  tmp = strtrim(tmp,2)
  if tmp eq '' then begin
   if print_queue_num eq '' then print_queue_num = '114'
   tmp = print_queue_num
  endif else begin
   print_queue_num = tmp(0)
  endelse

  white = where(plotColor eq 255) 
  if white[0] ne -1 then plotColor[white] = 0.
  ;this_yrange = [0,max(y[0:xyDim[0]-1,0])*1.05]
 
  set_plot, 'ps'
  device, /landscape
  plot_anything
  device, /close
  if tmp(0) ne 0 then spawn,'lpr -PP'+strtrim(tmp[0],2)+' idl.ps'

  set_plot, 'x'

  if white[0] ne -1 then plotColor[white] = 255.
 endelse

END

;**********************************************************************
; Backcalibration procedure.
;WARNING:  If you intend to re-instate the cal programs you must
;update it for the proper dimensionalities.
;
PRO cal
 common widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, w3c, w4, $
  w4c, draw, shot_label, shot_text, shot_button, frame_slider, $
  frameplot_button, lockyrange_button, timehist_button, movie_button, $
  from_label, from_text, from_label2, linewavelength_button, $
  to_label, to_text, to_label2, $
  xclick_label1, xclick_text, xclick_label2, yclick_text, $
  print_button, done_button, plot_radio, file_button, zoomout_button, $
  bin_slider, cal_button, savecal_button, delcal_button, spatial_button, $
  two_shot, current_button, curvie, di_box, plotlog_button
 common data, shotmarker, status_ok, data_ok, shot, frame, intensity, lambda, times, dt, $
             bin_numb, view_arr, counts, raw_data, bin, the_yrange, plottype, firstpix, $
             wavepix_no, binpix_no
 common cal, c0, c1, c2, c3, lambda_text, text, lines, impurities, lambda_spec, $
            cal_num, cal_delta, cal_x, cal_y, cal_chr, cal_frame, cal_max, cal_new, $
            cal_aberr
 
getcal

meas = counts
get_to_from, wto, wfrom, wto_i, wfrom_i

;find peak data in chosen plot
peakmax = max( meas(frame-1,wfrom_i:wto_i,bin-1) )
peak_x_i = where ( meas(frame-1,wfrom_i:wto_i,bin[0]-1) eq peakmax )
peak_x_i = wfrom_i + peak_x_i( (n_elements(peak_x_i)-1) /2. )
peak_l = lambda(peak_x_i, bin[0]-1)
str = string(form='(f7.3)', peak_l)

tmp = 'Current wavelength calibration gives ' + strtrim(string(peak_l),2) + ' nm for the peak value. '


;create window
lines = [379.1, 379.8, 386.4, 387.6, 390.3, 441.6, 485.1332, 513.29, 559.24, 658.1]
impurities = ["O III", "Mo I", "Mo I", "C II", "Mo I", "O II", "H", "C II", "O III", "C II"]
text = strarr(n_elements(lines))
for k=1,n_elements(lines),1 do begin
  text(k-1) =  '  ' + string(form='(f7.3)', lines(k-1) ) + '   ' + impurities(k-1) 
endfor

c0 = widget_base(/column, title='Wavelength CAL', xoffset=200, yoffset=200)
c1 = widget_base(c0, /row)
message_text = widget_text(c1, value=tmp, xsize=strlen(tmp), ysize=1)
c2 = widget_base(c0, /row)
lambda_label = widget_label(c2, value ='Spectral Wavelength: ')
lambda_text = widget_text(c2, value=str, xsize=7, ysize=1, /editable)
lambda_label2 = widget_label(c2, value=' nm')
c3 = widget_base(c0, /row)
okcal_button = widget_button(c3, value='OK')
c4 = widget_base(c3, /column)
for k=1,n_elements(lines),1 do begin
  line_button = widget_button(c4, value=text(k-1) )
endfor 

widget_control, c0, /realize
xmanager, 'lambdabox', c0, event_handler='lambda_event', /modal


;get lambda_ctr, lambda_gen, lambda from tree
mdsopen, 'spectroscopy', shot
cfg_per = mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.SETTINGS:CFG_PERISCOP')
lambda_ctr = mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.SETTINGS:LAMBDA_CTR')

lambda = mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.CALIBRATION:LAMBDA_TRUE', /quiet, status=cal_ok)
if (not cal_ok) then begin
  msg,'backcalibration not possible since no exact wavelength data stored in tree'
  zoom_out
  return
endif

lambda_gen = mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.CALIBRATION:LAMBDA_GEN', /quiet, status=cal_ok)
if (not cal_ok) then begin
  lambda_gen = lambda
  mdsput, '\SPECTROSCOPY::TOP.CHROMEX.CALIBRATION:LAMBDA_GEN', '$', lambda_gen, /quiet, status=cal_ok
  if (not cal_ok) then begin
    msg,'The node LAMBDA_GEN does not exist for this shot. Wavelength calibration not possible'
    mdsclose
    return
  endif
endif
mdsclose

;get peak data throughout all bins
cal_num_before = cal_num
delta = 0. 
bin_pre = bin

;bin_numb=n_elements(cfg_per_fibr)
;bin_numb=n_elements(lambda_gen(0,*))
;cfg_bin_no = mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.SETTINGS:BIN_NO', /quiet)
;bin_numb=cfg_bin_no

for spat=1,bin_numb,1 do begin
  if (cfg_per(spat-1) eq cfg_per(bin_pre-1)) then begin
    widget_control, bin_slider, set_value=spat
    bin = spat
    plot_frame
    peakmax = max( meas(frame-1,wfrom_i:wto_i,spat-1) )
    peak_x_i = where ( meas(frame-1,wfrom_i:wto_i,spat-1) eq peakmax )
    peak_x_i = wfrom_i + peak_x_i( (n_elements(peak_x_i)-1) /2. )
    peak_l = lambda_gen(peak_x_i, spat-1)
    delta = lambda_spec - peak_l
    dist_flag = abs(delta)
    if (dist_flag(0) le 4.0) then begin
      cal_delta(cal_num) = delta
      cal_x(cal_num) = peak_x_i
      cal_y(cal_num) = spat-1
      cal_chr(cal_num) = lambda_ctr
      cal_frame(cal_num) = frame-1
      print, peak_x_i, spat-1, delta
      cal_num = cal_num + 1
    endif
  endif
endfor
bin = bin_pre
widget_control, bin_slider, set_value=bin_pre
plot_frame


;find average for boundary condition
ave = 0.
for k=4, cal_num-1,1 do begin
  ave = ave + cal_delta(k)
endfor
ave = ave/(cal_num-4)
cal_delta(0:3) = ave

;fit surface
TRIANGULATE, cal_x(0:cal_num-1), cal_y(0:cal_num-1), tr, b
cal_aberr = TRIGRID( cal_x(0:cal_num-1), cal_y(0:cal_num-1), $
                     cal_delta(0:cal_num-1), tr, [1,1], [0,0,wavepix_no-1,13] ) 
lambda = lambda_gen + cal_aberr

z = cal_delta(0:cal_num-1)
x= cal_x(0:cal_num-1)
y= cal_y(0:cal_num-1)
save,  x,y,z, cal_aberr, file='test.dat'


msg,strtrim(string(cal_num-cal_num_before),2) + $
      ' data points have been added to backcalibration. There exist ' + $
      strtrim(string(cal_num-4),2) + ' data points for this shot.'
zoom_out

END


;*******************************************************************************
; This is the event handling procedure for the wavelength calibration.
;
PRO lambda_event, ev
 common cal, c0, c1, c2, c3, lambda_text, text, lines, impurities, lambda_spec, $
             cal_num, cal_delta, cal_x, cal_y, cal_chr, cal_frame, cal_max, cal_new, $
             cal_aberr

 type = tag_names(ev,/structure)
 widget_control, ev.id, get_value=button_name

 if type eq 'WIDGET_BUTTON' then begin
  if button_name eq 'OK' then begin
    widget_control, lambda_text, get_value=str
    widget_control, /destroy, ev.top
    lambda_spec = float(str)
  endif

  for k=1,n_elements(lines),1 do $
    if button_name eq text(k-1) $
     then widget_control, lambda_text, set_value=string(form='(f7.3)', lines(k-1) )
 endif

END

;**********************************************************************
; Save backcalibration procedure.
;
PRO savecal

common widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, w3c, w4, $
  w4c, draw, shot_label, shot_text, shot_button, frame_slider, $
  frameplot_button, lockyrange_button, timehist_button, movie_button, $
  from_label, from_text, from_label2, linewavelength_button, $
  to_label, to_text, to_label2, $
  xclick_label1, xclick_text, xclick_label2, yclick_text, $
  print_button, done_button, plot_radio, file_button, zoomout_button, $
  bin_slider, cal_button, savecal_button, delcal_button, spatial_button, $
  two_shot, current_button, curvie, di_box, plotlog_button

common data, shotmarker, status_ok, data_ok, shot, frame, intensity, lambda, times, dt, $
             bin_numb, view_arr, counts, raw_data, bin, the_yrange, plottype, firstpix, $
             wavepix_no, binpix_no

common cal, c0, c1, c2, c3, lambda_text, text, lines, impurities, lambda_spec, $
            cal_num, cal_delta, cal_x, cal_y, cal_chr, cal_frame, cal_max, cal_new, $
            cal_aberr

if (cal_num eq 4) then cal_new=4

if (cal_num eq cal_new) then begin
  msg,'No new calibration data present.'
  return
endif

;get grating number and static tree number
mdsopen, 'spectroscopy', shot
gr = mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.SETTINGS:GRATING')
gr = strtrim(string(gr),2)
stat_tree = mdsvalue('\SPECTROSCOPY::TOP.CHROMEX:STATIC_SHOT', /quiet, status=cal_ok)
mdsclose
if (not cal_ok) then begin
  stat_tree = 1
endif


;add new data points to static tree
mdsopen, 'chr_static', stat_tree
node_name = '\CHR_STATIC::TOP:CAL' + gr 
calfile = mdsvalue(node_name, /quiet, status=cal_ok)
if (not cal_ok) then begin
  ind_start = 0
  ind_length = cal_num-cal_new
  newfile = fltarr(ind_length, 7)
endif else begin
  ind_start = n_elements( calfile(*,0) )
  ind_length = ind_start + cal_num-cal_new
  newfile = fltarr(ind_length, 7)
  newfile(0:ind_start-1,*) = calfile(*,*)
endelse

for k=ind_start,ind_length-1,1 do begin
  newfile(k,0) = cal_delta(cal_new+k-ind_start)
  newfile(k,1) = cal_x(cal_new+k-ind_start)
  newfile(k,2) = cal_y(cal_new+k-ind_start)
  newfile(k,3) = cal_chr(cal_new+k-ind_start)
  newfile(k,4) = cal_frame(cal_new+k-ind_start)
  tmp = strtrim(string(shot),2)
  tmp_len=strlen(tmp)
  day = float(strmid(tmp, 0, tmp_len-3))
  shot_num = float(strmid(tmp,tmp_len-3,3))
  newfile(k,5) = day
  newfile(k,6) = shot_num
endfor

mdsput, node_name, '$', newfile
mdsclose


;write new shot specific wavelength calibration into tree
mdsopen, 'spectroscopy', shot
mdsput, '\SPECTROSCOPY::TOP.CHROMEX.CALIBRATION:LAMBDA_TRUE', '$', lambda
mdsclose

cal_num = 4
msg,'Wavelength calibration saved.'
zoom_out

END


;**********************************************************************
; Delete backcalibration procedure.
;
PRO delcal

common widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, w3c, w4, $
  w4c, draw, shot_label, shot_text, shot_button, frame_slider, $
  frameplot_button, lockyrange_button, timehist_button, movie_button, $
  from_label, from_text, from_label2, linewavelength_button, $
  to_label, to_text, to_label2, $
  xclick_label1, xclick_text, xclick_label2, yclick_text, $
  print_button, done_button, plot_radio, file_button, zoomout_button, $
  bin_slider, cal_button, savecal_button, delcal_button, spatial_button, $
  two_shot, current_button, curvie, di_box, plotlog_button

common data, shotmarker, status_ok, data_ok, shot, frame, intensity, lambda, times, dt, $
             bin_numb, view_arr, counts, raw_data, bin, the_yrange, plottype, firstpix, $
             wavepix_no, binpix_no

common cal, c0, c1, c2, c3, lambda_text, text, lines, impurities, lambda_spec, $
            cal_num, cal_delta, cal_x, cal_y, cal_chr, cal_frame, cal_max, cal_new, $
            cal_aberr

tmp = 'Do you really want to delete this shot from the backcalibration?'
del_flag = yes(tmp)

if (not del_flag) then return

;get grating number and static tree number
mdsopen, 'spectroscopy', shot
gr = mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.SETTINGS:GRATING')
gr = strtrim(string(gr),2)
stat_tree = mdsvalue('\SPECTROSCOPY::TOP.CHROMEX:STATIC_SHOT', /quiet, status=cal_ok)
mdsclose
if (not cal_ok) then begin
  stat_tree = 1
endif


;access static tree
mdsopen, 'chr_static', stat_tree
node_name = '\CHR_STATIC::TOP:CAL' + gr 
calfile = mdsvalue(node_name, /quiet, status=cal_ok)

if (not cal_ok) then begin
  ;no calibration file saved yet
endif else begin


  tmp = strtrim(string(shot),2)
  tmp_len=strlen(tmp)
  day = float(strmid(tmp, 0, tmp_len-3))
  shot_num = float(strmid(tmp,tmp_len-3,3))

  ind_length = n_elements( calfile(*,0) )
  newfile = fltarr(ind_length,7)
  ind = 0
  for k=0,ind_length-1,1 do begin
    unequal = (sign(calfile(k,5)-day))^2 + (sign(calfile(k,6)-shot_num))^2
    if (unequal(0) ne 0) then begin
      newfile(ind,0) = calfile(k,0)
      newfile(ind,1) = calfile(k,1)
      newfile(ind,2) = calfile(k,2)
      newfile(ind,3) = calfile(k,3)
      newfile(ind,4) = calfile(k,4)
      newfile(ind,5) = calfile(k,5)
      newfile(ind,6) = calfile(k,6)
      ind = ind +1
    endif
  endfor

  if (ind ne 0) then begin
     newfile = newfile(0:ind-1,*)
     mdsput, node_name, '$', newfile
  endif else begin
     tmp = fltarr(1,7)
     tmp(*,*) = 0.
     mdsput, node_name, '$', tmp
  endelse
endelse
mdsclose

mdsopen, 'spectroscopy', shot
lambda = mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.CALIBRATION:LAMBDA_GEN')
mdsput, '\SPECTROSCOPY::TOP.CHROMEX.CALIBRATION:LAMBDA_TRUE', '$', lambda, /quiet
mdsclose

msg,'Backcalibration data points for this shot are deleted.'
cal_num = 4
zoom_out
  
END


;**********************************************************************
; Retrieve saved backcalibration data.
;
PRO getcal

common widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, w3c, w4, $
  w4c, draw, shot_label, shot_text, shot_button, frame_slider, $
  frameplot_button, lockyrange_button, timehist_button, movie_button, $
  from_label, from_text, from_label2, linewavelength_button, $
  to_label, to_text, to_label2, $
  xclick_label1, xclick_text, xclick_label2, yclick_text, $
  print_button, done_button, plot_radio, file_button, zoomout_button, $
  bin_slider, cal_button, savecal_button, delcal_button, spatial_button, $
  two_shot, current_button, curvie, di_box, plotlog_button


common data, shotmarker, status_ok, data_ok, shot, frame, intensity, lambda, times, dt, $
             bin_numb, view_arr, counts, raw_data, bin, the_yrange, plottype, firstpix, $
             wavepix_no, binpix_no

common cal, c0, c1, c2, c3, lambda_text, text, lines, impurities, lambda_spec, $
            cal_num, cal_delta, cal_x, cal_y, cal_chr, cal_frame, cal_max, cal_new, $
            cal_aberr

if cal_num gt 4 then return

;get grating number and static tree number
mdsopen, 'spectroscopy', shot
gr = mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.SETTINGS:GRATING')
gr = strtrim(string(gr),2)
stat_tree = mdsvalue('\SPECTROSCOPY::TOP.CHROMEX:STATIC_SHOT', /quiet, status=cal_ok)
mdsclose
if (not cal_ok) then begin
  stat_tree = 1
endif


;access static tree
mdsopen, 'chr_static', stat_tree
node_name = '\CHR_STATIC::TOP:CAL' + gr 
calfile = mdsvalue(node_name, /quiet, status=cal_ok)
mdsclose

if (not cal_ok) then begin
  cal_new = 4
  return
endif

tmp = strtrim(string(shot),2)
  tmp_len=strlen(tmp)
  day = float(strmid(tmp, 0, tmp_len-3))
  shot_num = float(strmid(tmp,tmp_len-3,3))

ind_length = n_elements( calfile(*,0) )
for k=0,ind_length-1,1 do begin
  unequal = (sign(calfile(k,5)-day))^2 + (sign(calfile(k,6)-shot_num))^2
  if (unequal(0) eq 0) then begin
    cal_delta(cal_num) = calfile(k,0)
    cal_x(cal_num) = calfile(k,1)
    cal_y(cal_num) = calfile(k,2)
    cal_chr(cal_num) = calfile(k,3)
    cal_frame(cal_num) = calfile(k,4)
    cal_num = cal_num +1
  endif
endfor
cal_new = cal_num

END

;******************************************************************************
;
PRO two_shot_pro, ev

common widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, w3c, w4, $
  w4c, draw, shot_label, shot_text, shot_button, frame_slider, $
  frameplot_button, lockyrange_button, timehist_button, movie_button, $
  from_label, from_text, from_label2, linewavelength_button, $
  to_label, to_text, to_label2, $
  xclick_label1, xclick_text, xclick_label2, yclick_text, $
  print_button, done_button, plot_radio, file_button, zoomout_button, $
  bin_slider, cal_button, savecal_button, delcal_button, spatial_button, $
  two_shot, current_button, curvie, di_box, plotlog_button
common overplot, shot_text2, frame_slider2, bin_slider2, Overplot_button, compshot,$
       twotimefrom, twotimeto, utility_button, secondScale 

if (n_params() eq 0) then shotstr2 = '' else shotstr2 = string(format='(i10)', shot)

widget_control, done_button, set_uvalue='off'
widget_control, /destroy, two_shot
base2 = widget_base(w4c, /row, /frame)
  base3 = widget_base(base2, /column)
    b3b1 = widget_base(base3, /column)
      one_shot = widget_button(b3b1, value='One Shot',xsize=160)
      Overplot_button = widget_button(b3b1, value='Last Shot', xsize = 160)
      Utility_button  = widget_button(b3b1, value='Utility',xsize = 160)
  base4 = widget_base(base2, /column)
     b4b1= widget_base(base4, /row)
       shot_label2 = widget_label(b4b1, value='Shot:')
       shot_text2 = widget_text(b4b1, /editable, xsize=10, ysize=1, value=shotstr2)
       shot_button2 = widget_button(b4b1, value='Apply', uvalue=shotstr2)
     b4b2 = widget_base(base4, /row)
       frametimefrom = widget_label(b4b2, value='From:')
       twotimefrom = widget_text(b4b2,  xsize=10, value='0.000')
       frametimefrom2 = widget_label(b4b2, value=' sec')
     b4b3 = widget_base(base4, /row)
       frametimeto = widget_label(b4b3, value='  To:')
       twotimeto = widget_text(b4b3,  xsize=10, value='0.000')
       frametimeto2 = widget_label(b4b3, value=' sec')
  base5 = widget_base(base2, /column)
    b5b1 = widget_base(base5, /row)
      bin_slider2 = widget_slider(b5b1, min=1, max=binpix_no, value=7, title='(Bin)', xsize=200,/drag)
    b5b2 = widget_base(base5, /row)
      frame_slider2 = widget_slider(b5b2, min=1, max=35, title='(Frame)', xsize=200,/drag)

widget_control, base2, /realize
xmanager, 'Twos', base2, event_handler='Two_shot_events'


END

;******************************************************************************
;This is the 2nd shot event handler
;
PRO Two_shot_Events, ev

common widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, w3c, w4, $
  w4c, draw, shot_label, shot_text, shot_button, frame_slider, $
  frameplot_button, lockyrange_button, timehist_button, movie_button, $
  from_label, from_text, from_label2, linewavelength_button, $
  to_label, to_text, to_label2, $
  xclick_label1, xclick_text, xclick_label2, yclick_text, $
  print_button, done_button, plot_radio, file_button, zoomout_button, $
  bin_slider, cal_button, savecal_button, delcal_button, spatial_button, $
  two_shot, current_button, curvie, di_box, plotlog_button
common overplot, shot_text2, frame_slider2, bin_slider2, Overplot_button, compshot,$
       twotimefrom, twotimeto, utility_button, secondScale 
common data, shotmarker, status_ok, data_ok, shot, frame, intensity, lambda, times, dt, $
             bin_numb, view_arr, counts, raw_data, bin, the_yrange, plottype, firstpix, $
             wavepix_no, binpix_no
common curplot, x, y, xyDim, tit, xtit, ytit, xannot, yannot, annot, annotNum, $ 
                this_yrange, print_queue_num, view, num_of_plots, plotColor, plotColorFlag, $
                plotSymbol, plotSymbolFlag

widget_control, shot_text2, get_value=shotstr2
widget_control, frame_slider2, get_value=bing
frame[1]=bing
widget_control, bin_slider2, get_value=bong
bin[1]=bong
  
type = tag_names(ev, /structure)
widget_control, ev.id, get_value=button_name
shotmarker= 1

case type of
  'WIDGET_SLIDER': plot_frame
   

  'WIDGET_BUTTON': begin
		     case button_name of
                      'Utility': Utility
                 'Overplot Off': begin
                                  num_of_plots = 1
                                  shotmarker = 0
                                  compshot = 0 
                                  Plot_frame
                                  widget_control, Overplot_button, set_value='Overplot On'
                                end
                 'Overplot On': begin
                                  compshot = 1
                                  Plot_frame
                                  widget_control, Overplot_button, set_value='Overplot Off'                                 
                              end
                   'Last Shot': begin 
                                if shot[1] ne 1010101001 then begin 
                                   shot_str = string(form='(i10)', shot[1])
                                   widget_control, shot_text2, set_value=shot_str 
                                   widget_control, Overplot_button, set_value='Overplot Off'
                                   goto, applied
                               endif else begin
                                  widget_control, Overplot_button, sensitive = 0                                  
                                  widget_control, Di_box, /append, set_value='No Last Shot'
                              endelse
                          end
                       'Apply': if (shotstr2 ne '') then begin
                             applied:
                                  compshot = 1
                                  get_shot
                                  Plot_frame
                                  if lambda[0,bin[0],0] ne lambda[0,bin[0],1] then msg,'Different Wavelength Scale.'
                                  widget_control, Overplot_button, set_value='Overplot Off', sensitive = 1
                                  endif else msg,'No Shot Entered, tough guy.'
                    'One Shot': begin
                                  num_of_plots = 1
                                  compshot = 0
                                  shotmarker= 0
                                  widget_control, /destroy, w4c
                                  w4c = widget_base(w4, /row)
                                  two_shot = widget_button(w4c, value='Compare Shots', xsize=160)
                                  widget_control, done_button, set_uvalue = 'on'
                                  widget_control, w4c, /realize
                                  plot_frame
                                end
                             else: wait, 1
                          endcase
                      end

  'WIDGET_TEXT_CH' :  goto, applied 
   else: help, ev.id, type
endcase

END

;******************************************************************************
;This can be converted to any utility necessary at any given time.
;
PRO Utility

common widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, w3c, w4, $
  w4c, draw, shot_label, shot_text, shot_button, frame_slider, $
  frameplot_button, lockyrange_button, timehist_button, movie_button, $
  from_label, from_text, from_label2, linewavelength_button, $
  to_label, to_text, to_label2, $
  xclick_label1, xclick_text, xclick_label2, yclick_text, $
  print_button, done_button, plot_radio, file_button, zoomout_button, $
  bin_slider, cal_button, savecal_button, delcal_button, spatial_button, $
  two_shot, current_button, curvie, di_box, plotlog_button
common data, shotmarker, status_ok, data_ok, shot, frame, intensity, lambda, times, dt, $
             bin_numb, view_arr, counts, raw_data, bin, the_yrange, plottype, firstpix, $
             wavepix_no, binpix_no
common cal, c0, c1, c2, c3, lambda_text, text, lines, impurities, lambda_spec, $
            cal_num, cal_delta, cal_x, cal_y, cal_chr, cal_frame, cal_max, cal_new, $
            cal_aberr
common overplot, shot_text2, frame_slider2, bin_slider2, Overplot_button, compshot,$
       twotimefrom, twotimeto, utility_button, secondScale 
common wchr_view_event, from_flag, format_str
common curplot, x, y, xyDim, tit, xtit, ytit, xannot, yannot, annot, annotNum, $ 
                this_yrange, print_queue_num, view, num_of_plots, plotColor, plotColorFlag, $
                plotSymbol, plotSymbolFlag
common curvie, pos_csp, pos_innerwall, pos_cxrs,pol_25, limits, a_data, binfactors 
common gaussy, gaussdata
common line_data, line_wvlngths, line_idents, line_plot
common flaggs, dist, instrumentalFlag


goto, otherutil


;Possible utitilies:
;restore, file='Zerofind.dat', /verb


funshot =  fltarr(40,wavepix_no,wavepix_no)
funlamb =  fltarr(wavepix_no, wavepix_no)
shiftdata = fltarr(3,wavepix_no)
for binny = 0, binpix_no-1 do begin
  

  widget_control, bin_slider, set_value = strtrim(binny+1)
  widget_control, bin_slider2, set_value = strtrim(binny+1)
  widget_control, from_text, set_value='435.70'
  widget_control, to_text, set_value='437.'
  compshot = 1
  plot_frame

 low = locate(lambda[*,binny,0],435.7)
 high = locate(lambda[*,binny,0],437.)
  n_i = high - low+1
 xyDim[0:1] = n_i
 num_of_plots = 2
 plotSymbolFlag = 1 & plotColorFlag = 1
 plotSymbol[0] = -1 & plotColor[0] = 100
 plotSymbol[1] = -4 & plotColor[1] = 200
 x1 = lambda[low:high,binny,0]
 y1 = counts[0,low:high,binny,0]
 x2 = lambda[low+5:high+5,binny,1]
 y2 = counts[0,low+5:high+5,binny,1]
 x[0:n_i-1,0] = x1
 x[0:n_i-1,1] = x2
 y[0:n_i-1,0] = y1
 y[0:n_i-1,1] = y2
 view = ''
 compshot = 0
 annotNum = 3
 xannot[0:2] = .6
 for ff = 1,3 do yannot[ff] = .9 - ff/20.
 xshift = 0.
 yshift = 0.
 ymult = 1.
 the_y = this_yrange
tweak: 
 annot[0] = 'X Shift:  ' +   strtrim(xshift,2)
 annot[1] = 'Y Shift:  ' +   strtrim(yshift,2)
 annot[2] = 'Y Mult:  ' +   strtrim(ymult,2)

  plot_anything
  caseask = ask('xShift, yIncrease, Multiply, Log/Lin, The range,  Done?')

    case caseask[0] of

               's':  begin
                        xshift = ask('Shift in x?  :')
                        x[0:n_i-1,1] = x2 - xshift[0]
                        goto, tweak
                     end
               'i':  begin
                        yshift = ask('Increase in y?  :')
                        y[0:n_i-1,1] = y2*ymult[0] + yshift[0]
                        goto, tweak
                     end
               'm':  begin
                        ymult = ask('Multiply y?  :')
                        y[0:n_i-1,1] = y2*ymult[0] + yshift[0]
                        goto, tweak
                     end
               'l':  begin
                     if plottype eq 1 then begin
                     plottype = 0
                     this_yrange = the_y
                     endif else begin 
                     plottype = 1
                     this_yrange=[100, 30000]
                     endelse
                     goto, tweak
                     end
               't':  begin
                     down = ask('Range min?')
                     up   = ask('Range max?')
                     if down[0] ne '' then this_yrange[0] = down[0]
                     if up[0] ne '' then this_yrange[1] = up[0]
                     goto, tweak
                     end
               'd':  goto, doney
               else: goto, tweak

    endcase
doney:
 shiftdata[*,binny] = [xshift,yshift,ymult]
 xarr = [x1, x2 - xshift[0]]
 yarr = [reform(y1), reform(y2*ymult[0] + yshift[0])]
 xint = xarr[sort(xarr)]
 yint = yarr[sort(xarr)]
 plot, xint, yint, psym = -2
 others = fltarr((wavepix_no - 2.* n_i)/2.)
 yinterim = [others, yint, others]
 help, yinterim
 funlamb[*,binny] = [others+430., xint, others+440.]
 for j = 0, 39 do  funshot[j,*,binny] = yinterim
 plot, funlamb[*,binny], yinterim, psym = -1
 press_return
 
endfor
save, shiftdata, funlamb, funshot, file = 'comboshot'

press_return


mdsopen, 'spectroscopy', 7204
mdsput, '\SPECTROSCOPY::TOP.CHROMEX.ANALYSIS:RAW_SPECTRA', '$', funshot,  /quiet, status=data_ok
mdsput, '\SPECTROSCOPY::TOP.CHROMEX.CALIBRATION:LAMBDA_TRUE', '$', funlamb,  /quiet, status=data_ok
mdsclose, 'spectroscopy', 7204




;otherutil:

ans = 'y'

gaussdata = fltarr(4, 30, 7, 10)

shotlist = ['1030630019','1030630020','1030630023','1030630024','1030630025','1030630026',$
            '1030630027','1030630028','1030630030','1030630031']
bins = ['10','1','2','3','4','5','6']
frames = findgen(30)+3

superframe = fltarr(wavepix_no)
widget_control, shot_text, get_value=shoty

 widget_control, shot_text, set_value = shotlist[1]
;   compshot = 1
;   shotmarker = 1
;   get_shot
;   widget_control, frame_slider2, set_value = '3'
;   widget_control, bin_slider2, set_value =  '10'
;   bin[1] = 10
;   frame[1] = 4

     cter = 0.
for q = 0, 9 do begin
   widget_control, shot_text, set_value = shotlist[q]
   compshot = 0
   shotmarker = 0
   get_shot
   
  for binny = 0, 6 do begin
;     if binny eq 0 then distanceFlag = 1 else distanceFlag = 0
     nowbin = string(bins[binny])
     widget_control, bin_slider, set_value = nowbin
  ; if binny eq 1 then ans = ask('Good?')
  ; if ans ne 'y' then goto, done 
     for framy = 0, 29 do begin
     ;   compshot = 1
         nowframe = string(frames[framy])
         widget_control, frame_slider, set_value = nowframe    
         widget_control, from_text, set_value='468.85'
         widget_control, to_text, set_value='470.15'

         plot_frame
      ;     ans = ask('Include this frame?') 
      ;   if ans[0] eq '' or ans[0] eq 'y' then begin
         ;     y[*,24] = y[*,24] + y[*,1]
         ;   cter = cter + 1.
      ;   endif 
         curvefittwo
         wait, .5
         gaussdata[*, framy, binny, q] = a_data
     endfor

  endfor 


;compshot = 0
;y[*,1] = y[*,24]/cter
;num_of_plots=1

  xtit = 'Time'
  ytit = 'Position (nm)'
  tit = 'New Position by Time, Shot:  '+ shotlist[q]
  annotNum = 7
  xannot[0:6] = .3
  ;y1 = y[0:xyDim[1]-1,1]
  ;this_yrange = [min(y1) - .1*(max(y1)-min(y1)), max(y1) + .1*(max(y1)-min(y1))]
  plotSymbolFlag = 1
  plotColorFlag = 1
  plotColor[0] = !p.color
  ;print, cter
  plotColor[3:9] = !p.color
  ;plotColor[2] = 83
  view = ''

  symbols = ['  Line', '  Plus', '  Star', '  X-Symbol', '  Diamond', '  Triangle','  Square']
  for arg = 0,6 do begin
     x[0:29, arg]=times[frames-1,0];[18,19,20,23,24,25,26,27,28,30,31,32]  
     y[0:29, arg] = [reform(gaussdata[1,*,arg,q])]
     plotSymbol[arg] =arg
     yannot[arg]= .4- arg/25.
     annot[arg] = 'Bin ' +bins[arg] + symbols[arg]
  endfor
   plotSymbol[3] = 7 
   y[*,0] = y[*,0] - .02
   this_yrange = [469.40,469.50]
 
   line_plot = 0
   plot_anything
 
;   press_return
   print_it
   save, gaussdata, shotlist, bins, file = 'newfitdata.dat', /verb
endfor

;ans = ask('Good?')
;if ans ne 'y' then goto, done 
;curvefittwo

;print_it



done:
widget_control, shot_text, set_value = shoty 

otherutil:
widget_control, di_box, set_value = 'No utility loaded.  See Pro Utility', /append

END


;******************************************************************************
;
PRO DoubleZee, x, a, f, pder

common curvie, pos_csp, pos_innerwall, pos_cxrs,pol_25, limits, a_data, binfactors 
common flaggs, dist, instrumentalFlag

fact1 = binfactors[*,0,bin[shotmarker]-1]
fact2 = binfactors[*,1,bin[shotmarker]-1]
fact3 = binfactors[*,2,bin[shotmarker]-1]

big=n_elements(x)
f=fltarr(big)
pder=fltarr(big,n_elements(a))
f[*]=0.
lowlim = limits[0]
hilim = limits[1]

;Splitting by alpha*B
alpha = 0.205206 ; From Swede paper


;if distanceFlag then dist = 90. else dist = 45.
B = 5.3 * 67./dist ; Off axis magfield


split = alpha*B/10.

if a[2] gt .3 then a[2] = 0.07
if a[2] lt 0.001 then a[2] = 0.05
if a[3] lt 1. then a[3] = 500
if a[0] lt 1. then a[0] = 500

for j = 0,2 do begin


f = f + fact1[j]*a[0]*exp(-(x-a[1]-fact2[j]+split/2)^2/(fact3[j]^2 + a[2]^2))+$
        fact1[j]*a[0]*exp(-(x-a[1]-fact2[j]-split/2)^2/(fact3[j]^2 + a[2]^2))+$
        fact1[j]*a[3]*exp(-(x-a[1]-fact2[j]        )^2/(fact3[j]^2 + a[2]^2)) 

pder[*,0] = pder[*,0] + fact1[j]*exp(-(x-a[1]-fact2[j]+split/2)^2/(fact3[j]^2 + a[2]^2)) + $
                        fact1[j]*exp(-(x-a[1]-fact2[j]-split/2)^2/(fact3[j]^2 + a[2]^2))
pder[*,3] = pder[*,3] + fact1[j]*exp(-(x-a[1]-fact2[j]        )^2/(fact3[j]^2 + a[2]^2))

pder[*,1] = pder[*,1] +$
            (x-a[1]+split/2-fact2[j])/(fact3[j]^2 + a[2]^2)*2*fact1[j]*a[0]*exp(-(x-a[1]+split/2-fact2[j])^2/(fact3[j]^2 + a[2]^2))+$
            (x-a[1]        -fact2[j])/(fact3[j]^2 + a[2]^2)*2*fact1[j]*a[3]*exp(-(x-a[1]        -fact2[j])^2/(fact3[j]^2 + a[2]^2))+$
            (x-a[1]-split/2-fact2[j])/(fact3[j]^2 + a[2]^2)*2*fact1[j]*a[0]*exp(-(x-a[1]-split/2-fact2[j])^2/(fact3[j]^2 + a[2]^2))

pder[*,2] = pder[*,2] +$ 
            (x-a[1]+split/2-fact2[j])^2/(fact3[j]^2 + a[2]^2)^2*a[2]*2*fact1[j]*a[0]*exp(-(x-a[1]+split/2-fact2[j])^2/(fact3[j]^2 + a[2]^2))+$
            (x-a[1]        -fact2[j])^2/(fact3[j]^2 + a[2]^2)^2*a[2]*2*fact1[j]*a[3]*exp(-(x-a[1]        -fact2[j])^2/(fact3[j]^2 + a[2]^2))+$
            (x-a[1]-split/2-fact2[j])^2/(fact3[j]^2 + a[2]^2)^2*a[2]*2*fact1[j]*a[0]*exp(-(x-a[1]-split/2-fact2[j])^2/(fact3[j]^2 + a[2]^2))


endfor
return
END

;******************************************************************************
;
PRO gfunct3, x, a, f, pder

common data, shotmarker, status_ok, data_ok, shot, frame, intensity, lambda, times, dt, $
             bin_numb, view_arr, counts, raw_data, bin, the_yrange, plottype, firstpix, $
             wavepix_no, binpix_no
     
common curvie, pos_csp, pos_innerwall, pos_cxrs,pol_25, limits, a_data, binfactors
common flaggs, dist, instrumentalFlag

fact1 = binfactors[*,0,bin[shotmarker]-1]
fact2 = binfactors[*,1,bin[shotmarker]-1]
if instrumentalFlag eq 0 then fact3 = [0,0,0] else  fact3 = binfactors[*,2,bin[shotmarker]-1]

big=n_elements(x)
pks=(n_elements(a))/3
f=fltarr(big)
pder=fltarr(big,n_elements(a))
f[*]=0.


lowlim = limits[0]
hilim = limits[1]
;print, fact1
for i=0,pks-1 do begin
  ;n_width = n_elements(a)-1
  ;width = a[n_width]
  ;if width lt 0. then width = 0.01
  ad=3*i
  

  if a[0+ad] lt 0 then a[0+ad] = 1000.
  if a[1+ad] lt x(lowlim) or a[1+ad] gt x(hilim) then a[1+ad] = x(round((hilim-lowlim)/2.))
  if a[2+ad] gt 3.0 then a[2+ad] = 1.05
  if a[2+ad] lt 0.0 then a[2+ad] = 0.1

  for j = 0, instrumentalFlag do begin
   f = f + fact1[j]*a[0+ad]*exp(-(x-a[1+ad]-fact2[j])^2/(fact3[j]^2 + a[2+ad]^2))
   pder[*,0+ad]= pder[*,0+ad] + fact1[j]*exp(-(x-a[1+ad]-fact2[j])^2/(fact3[j]^2 + a[2+ad]^2))
   pder[*,1+ad]= pder[*,1+ad] + fact1[j]*a[0+ad]*exp(-(x-a[1+ad]-fact2[j])^2/(fact3[j]^2 + a[2+ad]^2))*2*(x-a[1+ad]-fact2[j])/(fact3[j]^2 + a[2+ad]^2)
   pder[*,2+ad]= pder[*,2+ad] + fact1[j]*a[0+ad]*exp(-(x-a[1+ad]-fact2[j])^2/(fact3[j]^2 + a[2+ad]^2))*2*(x-a[1+ad]-fact2[j])^2/(fact3[j]^2 + a[2+ad]^2)^2*a[2+ad]
  endfor
endfor


return

END

;******************************************************************************
;This will fit either of two curves present in the window
;
PRO curvefittwo, ev

common widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, w3c, w4, $
  w4c, draw, shot_label, shot_text, shot_button, frame_slider, $
  frameplot_button, lockyrange_button, timehist_button, movie_button, $
  from_label, from_text, from_label2, linewavelength_button, $
  to_label, to_text, to_label2, $
  xclick_label1, xclick_text, xclick_label2, yclick_text, $
  print_button, done_button, plot_radio, file_button, zoomout_button, $
  bin_slider, cal_button, savecal_button, delcal_button, spatial_button, $
  two_shot, current_button, curvie, di_box, plotlog_button
common data, shotmarker, status_ok, data_ok, shot, frame, intensity, lambda, times, dt, $
             bin_numb, view_arr, counts, raw_data, bin, the_yrange, plottype, firstpix, $
             wavepix_no, binpix_no
common overplot, shot_text2, frame_slider2, bin_slider2, Overplot_button, compshot,$
       twotimefrom, twotimeto, utility_button, secondScale 
common curvie, pos_csp, pos_innerwall, pos_cxrs,pol_25, limits, a_data, binfactors 
common curplot, x, y, xyDim, tit, xtit, ytit, xannot, yannot, annot, annotNum, $ 
                this_yrange, print_queue_num, view, num_of_plots, plotColor, plotColorFlag, $
                plotSymbol, plotSymbolFlag
common flaggs, dist, instrumentalFlag


restore, file = '/usr/local/cmod/codes/spectroscopy/chromex/he_gaspuff/binfactors_blip_11_20_03.dat', /verb


dist = 90.
case strmid(view_arr[bin[0]-1,0],0,3) of
                      'A_S': dist = 90.
                      'CSP': dist = pos_csp[long(strmid(view_arr[bin[shotmarker]-1],8))-1]
                      'INN': dist = pos_innerwall[long(strmid(view_arr[bin[shotmarker]-1],15))-1]
                      'CXR': dist = pos_cxrs_pol_25[long(strmid(view_arr[bin[shotmarker]-1],12))-1]
else: msg,'No distance information.'
endcase




plot_frame

if n_elements(the_yrange) then lock_yrange 

ans = ''
if compshot then ans = ask('White or Blue graph? (W/B)')
if ans eq 'b' or ans eq 'B' or ans eq 'blue' or ans eq 'Blue' then shotmarker = 1 else shotmarker = 0

compy = compshot
compshot = 0

get_to_from, wto, wfrom, wto_i, wfrom_i
g = fltarr(75)
n_data = [0,0]

ans = 2
ans = ask('Would you like 1) Regular or 2) Zeeman fit? (Just the number)')
if ans eq '' then ans = 1

qans = ''
qans = ask('Would you like with instrumental triplets? [Y/N]')
if qans eq 'n' or qans eq 'N' then instrumentalFlag = 0 else instrumentalFlag = 2


  if (not raw_data[shotmarker]) then meas=intensity else meas=counts
    x2 = lambda[wfrom_i : wto_i,bin[shotmarker]-1,shotmarker]
    y2 = reform(meas[frame[shotmarker]-1, wfrom_i : wto_i, bin[shotmarker]-1, shotmarker])

xyDim[0] = n_elements(x2)
x[0:xyDim[0]-1,0] = x2
y[0:xyDim[0]-1,0] = y2

tit = 'CHROMEX: Discharge ' + string(form='(i10, 2x)', shot[shotmarker]) + $
  '   Bin: ' + strtrim(string(bin[shotmarker]),2) + $
  '   [' + strtrim(string(form='(f6.3)',times[frame[shotmarker]-1,shotmarker]),2) + ' s, ' + $ 
  strtrim(string(form='(f6.3)', times[frame[shotmarker]-1,shotmarker]+dt[shotmarker]),2) +' s]'

view = view_arr[bin[shotmarker]-1,shotmarker]

if shotmarker eq 1 then this_yrange = secondScale else this_yrange = the_yrange
if shotmarker eq 1 then plotColorFlag = 1
if shotmarker eq 1 then plotColor[0] = 100
num_of_plots = 1
plot_anything

fact1 = binfactors[*,0,bin[shotmarker]-1]
fact2 = binfactors[*,1,bin[shotmarker]-1]
if instrumentalFlag eq 0 then fact3 = [0,0,0] else  fact3 = binfactors[*,2,bin[shotmarker]-1]


 widget_control, di_box, /append, set_value='Select Peak(s).' 
   
  ctr = -1
  label1:
    wait,.25
    

  if ans ne 2 then begin
    cursor,aaa,bbb & msg,string([aaa,bbb])
    wait,.25
    xpk = aaa
    if xpk gt x2[xyDim[0]-1] then goto, label2
    location = string(xpk)
    ctr=ctr+1
    add=3*ctr
    g[0+add]= y2(locate(x2,xpk))*.6
    g[1+add]= xpk
    g[2+add]= .04
   widget_control, di_box, /append, set_value='Peak at '+location+' nm (click to right when DONE)'
    goto,label1
  label2:
    if n_elements(add) eq 0 then goto, label1
    a = g[0:(ctr+1)*3 - 1]
    gfunky = 'gfunct3'
  endif else begin
    cursor,aaa,bbb & msg,string([aaa,bbb])
    wait, .25
    xpk = aaa
    a = fltarr(4)
    a[0] = y2(locate(x2,xpk))/6
    a[3] = y2(locate(x2,xpk))/10
    a[1] = xpk
    a[2] = .03
    gfunky = 'DoubleZee'
 endelse

  max_i = where(y2 eq max(y2))
  lowlim = locate(y2[0:max_i[0]], max(y2)*.02)
  if lowlim eq -1 then lowlim = 0
  dummy = where(y2[max_i[0]:n_elements(y2)-1] lt max(y2)*.02)
  if dummy[0] eq -1 then dummy[0] = (n_elements(y2) - 1 - max_i[0])
  hilim = max_i[0] + dummy[0]
  limits = [lowlim, hilim]

;Non-Auto Background
  widget_control, di_box, /append, set_value= 'Choose Left Background Point.'
  cursor, aaa, bbb
  wait, .25

  widget_control, di_box, /append, set_value= 'Choose Right Background Point.'
  cursor, ccc, ddd
  wait, .25  


;Automatic backgrounding...  
;puh = n_elements(x2)
;aaa = x2[0] & bbb = y2[0] & ccc = x2[puh-1] & ddd = y2[puh-1]
;for i=0, 2 do begin
; if y2[i+1] - y2[i] gt 100 then goto, skippy
; if y2[puh-4+i] - y2[puh-3+i] gt 100 then goto, skippy
;endfor
;aaa = avg(x2[0:3])
;bbb = avg(y2[0:3])
;ccc = avg(x2[puh-4:puh-1])
;ddd = avg(y2[puh-4:puh-1])
skippy:

  mmm = (ddd-bbb)/(ccc-aaa)
  bbbb = ddd - mmm*ccc

  back = mmm*x2 + bbbb
  oplot, x2, back
  wait, 1
  ybacked = y2 - back

  w = ybacked
; w = fltarr(n_elements(x2))+1.
  w[*]=0.
  w[lowlim:hilim]=10.

  ;The FIT**********

  r = curvefit(x2, ybacked, w, a, function_name=gfunky, itmax = 60, iter = iter, tol = .000001)
  print, iter,min(x2)

  xyDim[0:1] = n_elements(x2)
;  x[0:n_elements(x2)-1,0] = x2
 ; y[0:n_elements(x2)-1,0]=  y2
  x[0:n_elements(x2)-1,1] = x2
  y[0:n_elements(x2)-1,1]= r + back
;Store symbol and color in the dummy y[*,0] and x[*,0] slots
  plotSymbolFlag = 1
  plotColorFlag = 1
  plotColor[0] = !p.color
  plotSymbol[0] = 2
  plotColor[1] = 76
  plotSymbol[1] = 5

  x1 = findgen(1000)/999.*(max(x2)-min(x2))+min(x2)
  backspline =  mmm*x1 + bbbb

  Print, a
  a_data = a
  plotctr = 2
  if ans ne 2 then begin
  ;  width = a[n_elements(a)-1]
    add = -1
    for k=0, ((n_elements(a))/3)-1 do begin
        add = add+1
        addd = add*3
        y[*,plotctr] = 0.
     for j = 0, instrumentalFlag do begin
        y[0:n_elements(x1)-1, plotctr] =  y[0:n_elements(x1)-1, plotctr] + fact1[j]*a[0+addd]*exp(-(x1-a[1+addd]-fact2[j])^2/(fact3[j]^2 + a[2+addd]^2))
        plotColor[plotctr] = k*30. + 50.
     endfor
        annot[addd+0] = 'Line Height ' + strtrim(add+1,2) + ': ' + strtrim(string(max(y[0:n_elements(x1)-1, plotctr])),2)  
        annot[addd+1] = 'Line Position ' + strtrim(add+1,2) + ': ' + strtrim(string(a[1+addd]),2)  
        annot[addd+2] = 'Width Factor ' + strtrim(add+1,2) + ': ' + strtrim(string(a[2+addd]),2)     ;sqrt(a[2+addd]^2*4*alog(2))),2)  
        plotctr = plotctr + 1 
  endfor
  annotNum = n_elements(a)
  xannot[*] = .65
  for hi = 0, annotNum-1  do yannot[hi] = (.9 - (hi-1.)/25.)
  endif else begin
    ;For Zeeman plot
    ;Splitting by alpha*B
    alpha = 0.205206 ; From Swede paper
   ; if distanceFlag then dist = 90. else dist = 45.
    B = 5.3 * 67./dist ; Off axis magfield
    split = alpha*B/10.

    y[*,2:4] = 0.
for j = 0, 2 do begin
    y[0:n_elements(x1)-1,2] =  y[0:n_elements(x1)-1,2] + fact1[j]*a[0]*exp(-(x1-a[1]+split/2-fact2[j])^2/(fact3[j]^2 + a[2]^2))
    plotColor[2] = 42
    y[0:n_elements(x1)-1,3] =  y[0:n_elements(x1)-1,3] + fact1[j]*a[3]*exp(-(x1-a[1]-fact2[j])^2/(fact3[j]^2 + a[2]^2))
    plotColor[3] = 60
    y[0:n_elements(x1)-1,4] =  y[0:n_elements(x1)-1,4] +  fact1[j]*a[0]*exp(-(x1-a[1]-fact2[j]-split/2)^2/(fact3[j]^2 + a[2]^2))
    plotColor[4] = 42
    plotctr = plotctr + 3
endfor
    annotNum = 4
    xannot[0:3] = .65
    for hi = 0, 3 do yannot[hi] = (.9 - (hi-1.)/20.)
    annot[0] = 'Sigma Height: '+strtrim(string(a[0]),2)
    annot[1] = 'Pi Height: '+strtrim(string(a[3]),2)
    annot[2] = 'Position: '+strtrim(string(a[1]),2)
    annot[3] = 'Full Width: '+ strtrim(string(sqrt(a[2]^2*4*alog(2))),2)
  ;  print, sqrt(a[2]^2*4*alog(2))
  endelse


pixpeak = locate(x2, a[1])
pixslope = 1/(x2[pixpeak+1]-x2[pixpeak])
pixb = pixpeak - pixslope*x2[pixpeak]
pixfrac = pixslope*a[1]+pixb
print, wfrom_i
print, pixfrac + wfrom_i

annot[annotNum] = 'Pixel of First Peak: ' + strtrim(pixfrac + wfrom_i,2)
xannot[annotNum] = .65
yannot[annotNum] = .9 - (hi-1.)/25.
annotNum = annotNum+1


if plotctr gt 3 then  y[*, plotctr] = total(y[*,2:plotctr-1],2) else  y[*, plotctr] = y[*,2]
plotColor[plotctr] = 180
x[0:n_elements(x1)-1, plotctr] = x1 
plotctr = plotctr + 1
num_of_plots = plotctr
for o = 2, plotctr-1 do begin
  y[0:999,o] = y[0:999,o] + backspline
  x[0:999,o] = x1
endfor
xyDim[2:num_of_plots-1] = n_elements(x1)
plotSymbol[2:num_of_plots-1] = 0.

for n = 0, n_elements(a) -1 do begin 
       widget_control, di_box, /append, set_value= strtrim(string(a[n]),2)+'  ',/no_newline
endfor
       widget_control, di_box, /append, set_value= ''  

plot_anything

compshot = compy

the_yrange = 0
widget_control, lockyrange_button, set_value= 'Lock Y-Range'

end

;******************************************************************************
; This is the event handling procedure for the entire widget.
;
PRO chr_view_event, ev
 common widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, w3c, w4, $
  w4c, draw, shot_label, shot_text, shot_button, frame_slider, $
  frameplot_button, lockyrange_button, timehist_button, movie_button, $
  from_label, from_text, from_label2, linewavelength_button, $
  to_label, to_text, to_label2, $
  xclick_label1, xclick_text, xclick_label2, yclick_text, $
  print_button, done_button, plot_radio, file_button, zoomout_button, $
  bin_slider, cal_button, savecal_button, delcal_button, spatial_button, $
  two_shot, current_button, curvie, di_box, plotlog_button
 common data, shotmarker, status_ok, data_ok, shot, frame, intensity, lambda, times, dt, $
             bin_numb, view_arr, counts, raw_data, bin, the_yrange, plottype, firstpix, $
             wavepix_no, binpix_no
 common cal, c0, c1, c2, c3, lambda_text, text, lines, impurities, lambda_spec, $
            cal_num, cal_delta, cal_x, cal_y, cal_chr, cal_frame, cal_max, cal_new, $
            cal_aberr
 common nudge, nudge_to, nudge_from
 common wchr_view_event, from_flag, format_str
 common curplot, x, y, xyDim, tit, xtit, ytit, xannot, yannot, annot, annotNum, $ 
                 this_yrange, print_queue_num, view, num_of_plots, plotColor, plotColorFlag, $
                 plotSymbol, plotSymbolFlag
 common overplot, shot_text2, frame_slider2, bin_slider2, Overplot_button, compshot,$
        twotimefrom, twotimeto, utility_button, secondScale 

 type = tag_names(ev, /structure)
 widget_control, ev.id, get_value=button_name

 case type of
  'WIDGET_SLIDER': plot_frame
                               
  'WIDGET_BUTTON': begin
	case button_name of
	 'Print': print_it
	 'Quit': begin
          msg,'quitting...'
          widget_control, /destroy, ev.top 
         end
	 'Curvefit': curvefittwo
	 'Plot Frame / Zoom': Plot_frame
	 'Lock Y-Range': lock_yrange
	 'Unlock Y-Range': lock_yrange
	 'Plot Linear': begin
		widget_control, plotlog_button, set_value='Plot Logarithmic'
		this_yrange = [min(y[0:xyDim[0]-1,0])-.1*(max(y[0:xyDim[0]-1,0])-min(y[0:xyDim[0]-1,0])),$
		               max(y[0:xyDim[0]-1,0])+.1*(max(y[0:xyDim[0]-1,0])-min(y[0:xyDim[0]-1,0]))]
		if compshot then secondScale = [min(y[0:xyDim[3]-1,3])-.1*(max(y[0:xyDim[3]-1,3])-min(y[0:xyDim[3]-1,3])),$
		                                max(y[0:xyDim[3]-1,3])+.1*(max(y[0:xyDim[3]-1,3])-min(y[0:xyDim[3]-1,3]))]
		plottype = 0
		plot_anything
	  end
	 'Plot Logarithmic': begin
		widget_control, plotlog_button, set_value='Plot Linear' 
		this_yrange = [min(y[0:xyDim[0]-1,0])*.8,$
		               max(y[0:xyDim[0]-1,0])+.1*(max(y[0:xyDim[0]-1,0])-min(y[0:xyDim[0]-1,0]))]
		if compshot then secondScale =  [min(y[0:xyDim[3]-1,3])*.8,$
		                                 max(y[0:xyDim[3]-1,3])+.1*(max(y[0:xyDim[3]-1,3])-min(y[0:xyDim[3]-1,3]))]
		plottype = 1
		plot_anything
	  end
	 'Identify Lines': line_wavelength
	 'Compare Shots': two_shot_pro
	 'Time history': Time_history
	 'Spatial Scan': Spatial
	 'Movie': Movie
	 'Apply': begin
		applied:
		 plottype = 0
		 widget_control, plotlog_button, set_value='Plot Logarithmic'
		 shotmarker=0
		 widget_control, bin_slider, sensitive=1
		 widget_control, frame_slider, sensitive=1
		 widget_control, frameplot_button, sensitive=1
		 widget_control, zoomout_button, sensitive=1
		 widget_control, timehist_button, sensitive=1
		 widget_control, spatial_button, sensitive=1
		 widget_control, lockyrange_button, sensitive=1
		 widget_control, linewavelength_button, sensitive=1
		 widget_control, file_button, sensitive=1
		 widget_control, print_button, sensitive=1
		 widget_control, done_button, sensitive=1, get_uvalue=boop
		 widget_control, plotlog_button, sensitive=1
		 widget_control, movie_button, sens = 1
		 widget_control, curvie, sens = 1
		 if  boop eq 'on' then widget_control, two_shot, sens = 1
		 widget_control, current_button, sens = 1

		 get_shot
		 if data_ok then Plot_frame
		 cal_num = 4
	  end
	 'Set to Current': begin
		shot[0] =  mdscur_shot("cmod")
		shot_str = string(form='(i10)', shot[0])
		widget_control, shot_text, set_value=shot_str 
		goto, applied
	  end
	 'CAL': cal
	 'Save CAL': savecal
	 'Delete CAL': delcal
	 'Save as...' :	begin
		filename = ask('File Name?  (Include .dat): ')
		Print_it, FILE=filename[0]
	  end
	 '<-- Zoom Out -->' : zoom_out
	endcase
  end

  'WIDGET_DRAW' : begin ;process a button event in the draw window
	if data_ok then begin	;only deal with the draw widget if there is data in it

	 if (ev.press ne 0) or (ev.release ne 0) and status_ok then begin
	  click_loc = convert_coord(ev.x, ev.y, /device, /to_data)
	  if ev.press eq 2 then nudge_from = click_loc(0)	;2=middle mouse button shifts the spectrum
	  if ev.release eq 2 then begin
	   nudge_to = click_loc(0)
	   msg,' shift wavelength by '+string(nudge_to-nudge_from)+'nm'
	   lambda[0:wavepix_no-1,*,0] = lambda[0:wavepix_no-1,*,0] + (nudge_to-nudge_from)
	   Plot_frame
	  endif
	  if ev.press eq 4 then begin		;right mouse button undoes all shifts
	   msg,' remove wavelength shift'
	   get_lambda
	   cal_num = 4
	   Plot_frame
	  endif
       	  if ev.press eq 1 then $		;left mouse button for zooming
	   if not from_flag then begin
	    widget_control, from_text, set_uvalue=click_loc[0], set_value=string(form='(f7.2)', click_loc(0))
	    msg,' zoom from '+strtrim(string(click_loc(0)),2)
            from_flag = 1
	   endif else begin
	    widget_control, to_text, set_uvalue=click_loc[0], set_value=string(form='(f7.2)', click_loc(0))
	    msg,'        to '+strtrim(string(click_loc(0)),2)
            Plot_frame
            from_flag = 0
	   endelse
         endif ;"if ((ev.press ne 0) or (....."

	 if (not ev.press) and (not ev.release) and status_ok then begin
	  click_loc = convert_coord(ev.x, ev.y, /device, /to_data)
	  widget_control, xclick_text, set_value=string(form=format_str, click_loc(0))
	  widget_control, yclick_text, set_value=string(click_loc(1))
	 endif

	endif ;"if data_ok"
  end
  'WIDGET_TEXT_CH' :  goto, applied
  else : help, ev, /struct
 endcase
END

;**********************************************************************
;			M A I N   P R O G R A M
;**********************************************************************
PRO chr_view, shotx
 !p.multi = 0
 window,color=2
 wdelete
 common widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, w3c, w4, $
  w4c, draw, shot_label, shot_text, shot_button, frame_slider, $
  frameplot_button, lockyrange_button, timehist_button, movie_button, $
  from_label, from_text, from_label2, linewavelength_button, $
  to_label, to_text, to_label2, $
  xclick_label1, xclick_text, xclick_label2, yclick_text, $
  print_button, done_button, plot_radio, file_button, zoomout_button, $
  bin_slider, cal_button, savecal_button, delcal_button, spatial_button, $
  two_shot, current_button, curvie, di_box, plotlog_button
 common wchr_view_event, from_flag, format_str
 common data, shotmarker, status_ok, data_ok, shot, frame, intensity, lambda, times, dt, $
              bin_numb, view_arr, counts, raw_data, bin, the_yrange, plottype, firstpix, $
              wavepix_no, binpix_no
 common curplot, x, y, xyDim, tit, xtit, ytit, xannot, yannot, annot, annotNum, $
                 this_yrange, print_queue_num, view, num_of_plots, plotColor, plotColorFlag, $
                 plotSymbol, plotSymbolFlag
 common overplot, shot_text2, frame_slider2, bin_slider2, Overplot_button, compshot,$
        twotimefrom, twotimeto, utility_button, secondScale 
 common curvie, pos_csp, pos_innerwall, pos_cxrs_pol_25

 shot = [1010101001,1010101001]
 if n_params() eq 0 $
  then shot[0] = mdscur_shot("cmod") $
  else shot[0] = shotx
 shotstr = string(format='(i10)', shot[0])
 shotmarker= 0

 ;Determine frame sizes
 print, 'Opening Shot '+strtrim(string(shot[0]),2)+' for dimensioning'
 mdsopen, 'spectroscopy', shot[0]
 rawspect = mdsvalue('\spectroscopy::top.chromex.analysis:raw_spectra')
 mdsclose
 spect_dim  = size(rawspect)
 frame_no   = spect_dim(1)
 wavepix_no = spect_dim(2)
 binpix_no  = spect_dim(3)  
 print, 'Wavelength: ', strtrim(string(wavepix_no),2), ' pixel & LOS: ', strtrim(string(binpix_no),2), ' pixel'
 print, 'Frames: ', strtrim(string(frame_no),2), ' frames'


 the_yrange = 0
 compshot = 0
 print_queue_num = ''
 lambda = fltarr(wavepix_no,binpix_no,2)
 counts = fltarr(1,wavepix_no,binpix_no,2)
 y = fltarr(wavepix_no,25)
 x = fltarr(wavepix_no,25)
 annot = strarr(28)
 xannot = fltarr(28)
 yannot = fltarr(28)
 xyDim = fltarr(25)
 plotColor = fltarr(25)
 plotSymbol = fltarr(25)
 plotColorFlag = 0
 plotSymbolFlag = 0
 annotNum = 0
 bin_numb = [binpix_no,0]
 firstpix = 0

 base = widget_base(title='CHROMEX view (v12.GENIE)', /row)
 lcol = widget_base(base, /column, space=30)
 rcol = widget_base(base, /frame, /column)

 ; Left Column.
 w1 = widget_base(lcol, /column, /frame)
 w1a = widget_base(w1, /row)
 shot_label = widget_label(w1a, value='Shot:', uvalue = wavepix_no)
 shot_text = widget_text(w1a, /editable, xsize=10, ysize=1, value=shotstr)
 shot_button = widget_button(w1a, value='Apply')

 current_button = widget_button(w1, value='Set to Current')
 w1b = widget_base(w1, /row)
 bin_slider = widget_slider(w1b, min=1, max=binpix_no, value=1, title='(Bin)', xsize=200,/drag)
 w1c = widget_base(w1, /row)
 frame_slider = widget_slider(w1c, min=1, max=frame_no, title='(Frame)', xsize=200,/drag)
 
 w2 = widget_base(lcol, /column, /frame)
 frameplot_button = widget_button(w2, value='Plot Frame / Zoom', xsize=200)
 lockyrange_button = widget_button(w2, value='Lock Y-Range', xsize=200)
 linewavelength_button = widget_button(w2, value='Identify Lines', xsize=200)
 zoomout_button = widget_button(w2, value='<-- Zoom Out -->', xsize=200)
 Curvie = widget_button(w2, value='Curvefit', xsize=200)
 plotlog_button = widget_button(w2, value='Plot Logarithmic', xsize=200)

 w5 = widget_base(lcol, /column, /frame)
 timehist_button = widget_button(w5, value='Time history', xsize=200)
 spatial_button = widget_button(w5, value='Spatial Scan', xsize=200)
 movie_button = widget_button(w5, value='Movie', xsize=200)

 w3 = widget_base(lcol, /column, /frame)
 w3a = widget_base(w3, /row)
 from_label = widget_label(w3a, uvalue = 1, value='    From:')
 from_text = widget_text(w3a, value='', uvalue = 0, xsize=7, ysize=1, /editable)
 from_label2 = widget_label(w3a, value='nm  ')
 w3b = widget_base(w3, /row)
 to_label = widget_label(w3b, uvalue = wavepix_no-1, value='        To:')
 to_text = widget_text(w3b, value='', uvalue = 0, xsize=7, ysize=1, /editable)
 to_label2 = widget_label(w3b, value='nm  ')
 w3c = widget_base(w3, /row)
 xclick_label1 = widget_label(w3c, value=' Cursor:')
 xclick_text = widget_text(w3c, value='', xsize=7, ysize=1)
 xclick_label2 = widget_label(w3c, value='nm  ')
 w3d = widget_base(w3, /row)
 yclick_label1 = widget_label(w3d, value='            ')
 yclick_text = widget_text(w3d, value='', xsize=11, ysize=1)

 w7 = widget_base(lcol, /column)
 print_button = widget_button(w7, value='Print', xsize=200)
 file_button = widget_button(w7, value='Save as...', xsize=200)
 done_button = widget_button(w7, value='Quit', xsize=200, uvalue= 'on')

 ;Right Column

 ;Drawing pane.
 draw = widget_draw(rcol, xsize= 700, ysize=550, retain=2, /motion_events, /button_events)
 w4 = widget_base(rcol, /column)
 ;w4a = widget_base(w4, /column)
 w4b = widget_base(w4, /column)
 w4c = widget_base(w4, /row)

 ;cal_button = widget_button(w4a, value='CAL',xsize=160)
 ;savecal_button = widget_button(w4a, value='Save CAL',xsize=160)
 ;delcal_button = widget_button(w4a, value='Delete CAL',xsize=160)

 di_box = widget_text(w4b, xsize = 110, ysize = 5, /wrap, /scroll,  value= '')
 two_shot = widget_button(w4c, value='Compare Shots', xsize=160)

 ;before anything is drawn, deactivate calibration buttons
 ;  widget_control, cal_button, sensitive=0
 ;  widget_control, savecal_button, sensitive=0
 ;  widget_control, delcal_button, sensitive=0
; widget_control, bin_slider, sensitive=0
; widget_control, frame_slider, sensitive=0
; widget_control, frameplot_button, sensitive=0
; widget_control, zoomout_button, sensitive=0
; widget_control, timehist_button, sensitive=0
; widget_control, spatial_button, sensitive=0
; widget_control, lockyrange_button, sensitive=0
; widget_control, linewavelength_button, sensitive=0
; widget_control, file_button, sensitive=0
; widget_control, print_button, sensitive=0
; widget_control, done_button, sensitive=1
; widget_control, movie_button, sens = 0
; widget_control, curvie, sens = 0
; widget_control, two_shot, sens = 0
; widget_control, plotlog_button, sensitive=0
 ;widget_control, current_button, sens = 0
 widget_control, base, /realize
 widget_control, get_value=draw_window, draw
 wset, draw_window
 ;print, draw_window

 plottype = 0
 data_ok = 0
 get_shot
 if data_ok then Plot_frame

 ;configure backcalibration
 common cal, c0, c1, c2, c3, lambda_text, text, lines, impurities, lambda_spec, $
             cal_num, cal_delta, cal_x, cal_y, cal_chr, cal_frame, cal_max, cal_new, $
             cal_aberr

 cal_num = 4
 cal_max=100
 cal_delta = fltarr(cal_max)
 cal_x = fltarr(cal_max)
 cal_y = fltarr(cal_max)
 cal_chr = fltarr(cal_max)
 cal_frame = fltarr(cal_max)

 ;boundary condition for fit
 cal_delta(0:3) = 0.
 cal_x(0:3) = [0,wavepix_no-1,0,wavepix_no-1]
 cal_y(0:3) = [0,0,13,13]

 from_flag = 0
 format_str = '(f7.2)'	;format for displaying the cursor location in Angstrom
 xmanager,'chr_view',base, event_handler='chr_view_event'

END


;Main Program
chr_view
END

