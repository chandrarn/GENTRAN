Function Int, x

;+
; NAME:
;	INT
; PURPOSE:
;	Gives the integer part of X, i.e. the largest integer(s) not 
;	exceeding X.  The result is the same as FIX(X) for positive numbers, 
;	but for negative non-integer X it is FIX(X) - 1.
; CATEGORY:
;	Mathematical Function (General) / Type Conversion.
; CALLING SEQUENCE:
;	Result = INT(X)
; INPUTS:
;    X
;	Numerical, otherwise arbitrary.
; OPTIONAL INPUT PARAMETERS:
;	None.
; KEYWORD PARAMETERS:
;	None.
; OUTPUTS:
;	Returns the value of INT(X), see above, as a long integer.
; OPTIONAL OUTPUT PARAMETERS:
;	None.
; COMMON BLOCKS:
;	None.
; SIDE EFFECTS:
;	None.
; RESTRICTIONS:
;	For complex X the result is INT(REAL(X)), the imaginary part is ignored
; PROCEDURE:
;	Using the system function LONG and the function SIGN from MIDL.
; MODIFICATION HISTORY:
;	Created 15-JUL-1991 by Mati Meron.
;-

    return, long(2*x - 1 + Sign(x))/2
end
