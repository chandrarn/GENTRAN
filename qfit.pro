;translates the output from quick_fits.pro
;written by A. Dominguez (9/10)
;generalized to use Te - M.L. Reinke (9/21/10)

pro qfit,shot,dens,temp,rmaj,time,ne_deg=ne_deg,te_deg=te_deg,data=data

if not(keyword_set(ne_deg)) then ne_deg=2
if not(keyword_set(te_deg)) then te_deg=2

quick_fit,shot,data,ne_deg=ne_deg,te_deg=te_deg,/nosave	;run quick fit

;reorganize output data [space,time] arrays
time        = data.t_ts			;time in [sec]
r0_fit      = data.r_out
r_fit       = data.ax
tvec  = time*0+1.
rvec  = r_fit*0+1.
x     = r0_fit#rvec
y     = tvec#r_fit
rmaj  = rotate((x+y),4)			;outboard midplane major radius [m]
f_fit       = data.afp
dens = rotate(f_fit(*,*,0),4)		;density in 10^20 [m-3]
temp = rotate(f_fit(*,*,1),4)/1.0e3	;temperature in [keV]


end
