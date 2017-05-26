pro csv_convert,file1,file2,mlr=mlr
; written by C S Pitcher 9 June 1999
; modified by C S Pitcher 19 April 2000

; converts a comma-delimited spread-sheet produced by
; the edge data base to an idl save set

; file1 = input, comma delimited file
; file2 = output, idl save set with data, program adds *.dat
;
; Note that the '.dat' file is formed from the '.csv' file by the program
; 'csv_convert.pro' with the following caveat - blank entries are handled as
;       -99999999 for long variables
;       '' for strins
;       3e38 for floats
;       -9999 for integers
;
line=''

get_lun,u1
openr,u1,file1

;get info
nd = 0
while not eof(u1) do begin
  readf,u1,line
  temp = str_sep(line,':')
  if temp(0) eq 'Data' then begin
    nd = nd + 1
    endif
  if temp(0) eq 'Variables' then begin
    temp = str_sep(temp(1),',')
    ii = where( temp ne '' and temp ne ' ')
    var = temp(ii)
    IF keyword_set(mlr) THEN FOR i=1,n(var) DO var[i]=var[i]+num2str(i-1,1)
    nv = n_elements(var)
    endif
  if temp(0) eq 'Type' then begin
    temp = str_sep(temp(1),',')
    ii = where( temp ne '' and temp ne ' ')
    type = temp(ii)
    endif
  endwhile

; correct some of the names
for i = 0,nv-1 do begin
  i1 = strpos(var(i),'(')
  if i1 ne -1 then begin
    i2 = strpos(var(i),')')
    var(i) = strmid( var(i),0,i1) + '_' + strmid(var(i),i1+1,i2-i1-1)
    endif
  endfor

; make arrays
for i = 0,nv-1 do begin
  if type(i) eq 'Long' then begin
    temp = execute( var(i) + '= lonarr(' + string(nd) + ')' )
    endif
  if type(i) eq 'String' then begin
    temp = execute( var(i) + '= strarr(' + string(nd) + ')' )
    endif
  if type(i) eq 'Float' then begin
    temp = execute( var(i) + '= fltarr(' + string(nd) + ')' )
    endif
  if type(i) eq 'Integer' then begin
    temp = execute( var(i) + '= intarr(' + string(nd) + ')' )
    endif
  endfor

close,u1
free_lun,u1


get_lun,u1
openr,u1,file1
 
; read data into arrays

i = -1
while not eof(u1) do begin
  readf,u1,line
  temp = str_sep(line,':')
  if temp(0) eq 'Data' then begin
    i = i + 1
j1:
    temp = str_sep(line,',')
; make correction if the correct number of commas not there
    if n_elements(temp) lt nv + 1 then begin
      line = line + ','
      goto,j1
      endif
    temp = temp(1:nv)

; put in blanks or data, where appropriate
    for j = 0,nv-1 do begin
      if temp(j) eq '' then begin
        if type(j) eq 'Long' then begin
          tt = execute( var(j) + '(' + string(i) + ') = -99999999' )
          endif
        if type(j) eq 'String' then begin
          tt = execute( var(j) + '(' + string(i) + ') = temp(j)' )
          endif
        if type(j) eq 'Float' then begin
          tt = execute( var(j) + '(' + string(i) + ') = 3e38' )
          endif
        if type(j) eq 'Integer' then begin
          tt = execute( var(j) + '(' + string(i) + ') = -9999' )
          endif
      endif else begin
        if type(j) eq 'Long' then begin
          tt = execute( var(j) + '(' + string(i) + ') = long(temp(j))' )
          endif
        if type(j) eq 'String' then begin
          tt = execute( var(j) + '(' + string(i) + ') = temp(j)' )
          endif
        if type(j) eq 'Float' then begin
          tt = execute( var(j) + '(' + string(i) + ') = float(temp(j))' )
          endif
        if type(j) eq 'Integer' then begin
          tt = execute( var(j) + '(' + string(i) + ') = fix(temp(j))' )
          endif
        endelse  
      endfor

    endif
  endwhile

temp = var(0)
for j = 1,nv-1 do begin
  temp = temp + ',' + var(j)
  endfor

tt = execute( 'save,' + temp + ',file=file2,/verbose' )

;close,u1
free_lun,u1 ; automatically closes the file as well

return

end


