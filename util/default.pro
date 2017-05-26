Function Default, x, y, type = typ

;+
; NAME:
;	DEFAULT
; PURPOSE:
;	Provides an automatic default value for nondefined parameters.
; CATEGORY:
;	Programming.
; CALLING SEQUENCE:
;	Result = DEFAULT( X, Y [, TYPE = TYP])
; INPUTS:
;    X, Y
;	Arbitrary, at least one needs to be defined.
; OPTIONAL INPUT PARAMETERS:
;	None.
; KEYWORD PARAMETERS:
;    TYPE
;	Number between 1 to 8.  If provided, it is compared to the type-value 
;	of X and X is considered defined only if the two values match.  For
;	example, if TYPE = 7 then X has to be a character value.
; OUTPUTS:
;	X if it is defined, otherwise Y.
; OPTIONAL OUTPUT PARAMETERS:
;	None.
; COMMON BLOCKS:
;	None.
; SIDE EFFECTS:
;	None.
; RESTRICTIONS:
;	None.
; PROCEDURE:
;	Uses the function TYPE from MIDL.
; MODIFICATION HISTORY:
;	Created 15-JUL-1991 by Mati Meron.
;-

    if n_elements(typ) ne 0 then begin
	if typ gt 0 and Type(x) eq typ then return, x  else return, y
    endif else if n_elements(x) ne 0 then return, x  else return, y
end
