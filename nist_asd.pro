pro nist_asd,spec,lowwl,uppwl,asd_struct,order=order,verb=verb,msg=msg,result=result
; nist_lines.pro - Retrieves data from the online NIST Atomic Spectra Database
; eg: nistln=nist_asd('He',150,1500)

; ANJ NOTES TO SELF:
; 1. document this damn code
; 2. occasionally NIST changes the column order for whatever reason. parse out data lines based on what was parsed from the header instead of hard coding it as presently done.

; manage order
if keyword_set(order) eq 0 then order=1
lowwl=lowwl/order
uppwl=uppwl/order

; replace spaces with '+'
spec=strtrim(spec,2)
spec_plus=spec
tmp=where(spec eq ' ')
if tmp ne -1 then for i=0,n(tmp) do strput,spec_plus,'+',tmp[i]

; build the web request
nisturl='http://physics.nist.gov/cgi-bin/ASD/lines1.pl'
postdata= $
 'encodedlist=XXT1XXR0q0qVqVIII'+'&' $ ; some key to make it work?
+'spectra='+spec_plus+'&' $ ; eg 'He' or 'He+I' or 'He+II', no spaces
+'low_wl='+num2str(lowwl,dp=3)+'&' $
+'upp_wl='+num2str(uppwl,dp=3)+'&' $
+'unit=0'+'&' $ ; wl unit 0=Angstroms, 1=nm, 2=um
+'en_unit=1'+'&' $ ; energy unit 0 cm^-1, 1 eV, 2 Rydberg
+'low_wn='+'&' $
+'upp_wn='+'&' $
+'submit=Retrieve+Data'+'&' $
+'temp='+'&' $
+'doppler='+'&' $
+'eden='+'&' $
+'iontemp='+'&' $
+'java_window=3'+'&' $
+'java_mult='+'&' $
+'tsb_value=0'+'&' $
+'format=1'+'&' $ ; 0 HTML output, 1 ascii output
+'remove_js=on'+'&' $ ; cleans up output for easier parsing
+'output=0'+'&' $ ; 0 return all output, 1 return output in pages
+'page_size=15'+'&' $
+'line_out=0'+'&' $ ; 0 return all lines, 1 only w/trans probs, 2 only w/egy levl, 3 only w/obs wls
+'order_out=0'+'&' $ ; output ordering: 0 wavelength, 1 multiplet
+'show_av=2'+'&' $ ; show wl in Vacuum (<2000A) Air (2000-20000A) Vacuum (>20,000A)
+'max_low_enrg='+'&' $ ; maximum lower level energy
+'max_upp_enrg='+'&' $ ; maximum upper level energy
+'min_str='+'&' $ ; minimum transition strength
+'max_str='+'&' $ ; maximum transition strength
+'min_accur='+'&' $ ; minimum line accuracy, eg AAA AA A B C
+'min_intens='+'&' $ ; minimum relative intensity to return
+'show_obs_wl=1'+'&' $ ; show observed wavelength
+'show_calc_wl=1'+'&' $ ; show calculated wavelength
+'A_out=0'+'&' $ ; show $
+'intens_out=on'+'&' $ ; show relative intensity
+'allowed_out=1'+'&' $ ; show allowed transitions
+'forbid_out=1'+'&' $ ; show forbidden transitions
+'conf_out=on'+'&' $ ; show electron configuration
+'term_out=on'+'&' $ ; show terms
+'enrg_out=on'+'&' $ ; show transition energies
+'J_out=on'+'&' $ ; show J
+'g_out=on' ; show g

; issue wget to pull the data from nist and use sed to split off the desired info
;  -q 'quiet' suppresses wget messages
;  -O - directs results to standard output
cmd='wget -q -O - ''' $
	+nisturl+'?'+postdata+''' ' $ ; This issues as a GET instead of POST, but it works ok anyway
	+'| sed -n ''/<pre*/,/<\/pre>/p'' ' $ ; select lines between <pre> tags
	+'| sed ''/<*pre>/d'' ' $ ; remove <pre> lines
	+'| iconv -f ISO-8859-1 -t ASCII' ; convert the web encoding to something IDL can understand...
	;'| sed ''/----*/d''' $ ; remove ---- lines
if keyword_set(verb) then print,cmd
spawn,cmd,result,err,exit_status=stat


;if keyword_set(verb) then print,result,err,stat

if stat then begin
	msg=result+string(10B)+'Error retrieving NIST data.'+string(10B)+err
	print,msg
 	nln=0
endif else if result[0] eq '' then begin
	msg='No lines in range.'
	print,msg
 	nln=0
endif else begin


 ; Parse header info
 ; this header parsing doesn't actually work, so the below fields are parsed by hard-coding. maybe not the best approach, but it works for now.
 hd1=strsplit(result[1],'|',/extract)
  i1=strsplit(result[1],'|',length=l1)
  n1=size(hd1,/n_el)
 hd2=strarr(n1)
  for i=0,n1-1 do hd2[i]=strmid(result[2],i1[i],l1[i])
 hd3=strarr(n1)
  for i=0,n1-1 do hd3[i]=strmid(result[3],i1[i],l1[i])

 hdr=strarr(n1)
  for i=0,n(hd1)-1 do hdr[i]=strtrim(strcompress(hd1[i]+' '+hd2[i]+' '+hd3[i]),2)
  if keyword_set(verb) then print,hdr

 if hdr[0] eq 'Spectrum' $
  then isp=1 $
  else isp=0

 ; Parse line data from remaining result
 res=result[4:n(result)]
 nln=0
 for i=0,n(res) do begin
  tok=res[i]
  if keyword_set(verb) then print,tok

  if strmid(tok,0,4) eq '----' then continue ; skip divider lines

  tok=strtrim(strsplit(tok,'|',/extract),2)
  if tok[0+isp] eq '' and tok[1+isp] eq '' then continue ; skip lines with no wavelength defined

  if tok[0+isp] eq '' then ichk=1+isp else ichk=0+isp
  byt=byte(strmid(tok[ichk],0,1)) ; skip if header lines are repeated
  if byt lt 48 or byt gt 57 then continue

  if isp $
   then spec=strtrim(tok[0],2) $
   else tok=[spec,tok]
  if tok[1] eq '' then obs=-1.  else obs=float(tok[1])*order
  if tok[2] eq '' then ritz=-1. else ritz=float(tok[2])*order
  if tok[3] eq '' then begin
  	rint=''
	int=-1.
  endif else begin
  	rint=tok[3]
	iint=0
	for j=0,strlen(rint)-1 do begin
	 bint=byte(strmid(rint,j,1))
	 if bint ge 48 and bint le 57 then begin
		if iint eq 0 then iint0=j
		iint+=1 ; string(byte(48)) = '0', string(byte(57)) = '9'
	 endif
	endfor
	if iint gt 0 $
	 then int=float(strmid(rint,iint0,iint)) $
	 else int=0.0
  endelse
  if tok[4] eq '' then Aki=-1.  else Aki=float(tok[4])
  if tok[5] eq '' then Acc=''   else Acc=tok[5]
  if tok[6] eq '' then begin
   Ei=-1.
   Ek=-1.
  endif else begin
   tmp=strtrim(strsplit(tok[6],'-',/extract),2)
   if strmid(tmp[0],0,1) eq '(' then tmp[0]=strmid(tmp[0],1,strlen(tmp[0])-2)
   if strmid(tmp[1],0,1) eq '(' then tmp[1]=strmid(tmp[1],1,strlen(tmp[1])-2)
   if strmid(tmp[0],0,1) eq '[' then tmp[0]=strmid(tmp[0],1,strlen(tmp[0])-2)
   if strmid(tmp[1],0,1) eq '[' then tmp[1]=strmid(tmp[1],1,strlen(tmp[1])-2)
   Ei=float(tmp[0])
   Ek=float(tmp[1])
  endelse
  if tok[7] eq '' then conf0='' else conf0=tok[7]
  if tok[8] eq '' then term0='' else term0=tok[8]
  if tok[9] eq '' then J0=''  else J0= tok[9]
  if tok[10] eq '' then conf1='' else conf1=tok[10]
  if tok[11] eq '' then term1='' else term1=tok[11]
  if tok[12] eq '' then J1=''  else J1= tok[12]
  if tok[13] eq '' then begin
	gi=-1.
	gk=-1.
  endif else begin
   	tmp=strtrim(strsplit(tok[13],'-',/extract),2)
	if size(tmp,/n_el) gt 1 then begin
		gi=float(tmp[0])
		gk=float(tmp[1])
	endif else begin
		gi=-1.
		gk=-1.
	endelse
  endelse
  if tok[14] eq '' then Type='' else Type= tok[14]
  if obs ne -1 $
   then wl=obs $
   else wl=ritz

  ln0={spec:spec, obs:obs,ritz:ritz,wl:wl, rint:rint,int:int, Aki:Aki,Acc:Acc, Ei:Ei,Ek:Ek, $
       conf0:conf0,conf1:conf1, term0:term0,term1:term1, J0:J0,J1:J1,gi:gi,gk:gk, Type:Type, order:order }

  if order gt 1 then ln0.spec=ln0.spec+' (order='+num2str(order,1)+')'
  if keyword_set(verb) then begin print,'nln=',nln & help,/str,ln0 & endif

  if nln eq 0 $
   then ln=ln0 $
   else ln=[ln,ln0]
  nln=nln+1;
 endfor ; 0,n(hd)
endelse ; if stat

if nln eq 0 then begin
 ln={spec:'none', obs:-1.,ritz:-1.,wl:-1., rint:'',int:-1., Aki:-1.,Acc:'', Ei:-1.,Ek:-1., $
     conf0:'',conf1:'', term0:'',term1:'', J0:'',J1:'', gi:0.,gk:0., Type:'', order:-1. }
 hdr=''
endif

msg=num2str(nln,1)+' lines loaded from NIST ASD.'
print,msg

asd_struct={nisturl:nisturl,postdata:postdata,stat:stat, $
   hdr:hdr,result:result,lines:ln,nln:nln,msg:msg}

end

