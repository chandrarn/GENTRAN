pro ssvdc,a,s,u,v
; Single precision SVDecomposition. Uses a LINPACK  routine.
; No tolerance, please. See USER10:[HORNE]DSVDC.TXT for documentation. 
; GT 11/19/93. IHH 30 Aug 2001.
;c     dsvdc is a subroutine to reduce a double precision nxp matrix x
;c     by orthogonal transformations u and v to diagonal form.  the
;c     diagonal elements s(i) are the singular values of x.  the
;c     columns of u are the corresponding left singular vectors,
;c     and the columns of v the right singular vectors.
; That is, u^T a v = diag(s), u diag(s) v^T = a.

ss=size(a)
if (ss(0) ne 2) then begin
        print,'Array must be 2 dimensional, not',ss(0),'dimensional.'
        return
endif
; Correct dimensions, let's go.

x=float(a)    ; array must be float precision
ldx=long(ss(1))
n=long(ss(1))
p=long(ss(2))
ldu=n
ldv=p
work=fltarr(n)
job=long(11)
s=fltarr(min([n+1,p]))
e=fltarr(p)
u=fltarr(ldu,n)
v=fltarr(ldv,p)
info=long(111)
;err=call_external('/home/hutch/idld/linpack/liblinpackidl.so','ssvdc',x,ldx,n,p,s,e,u,ldu,v,ldv,work,job,info)
err=call_external('/usr/local/cmod/lib/liblinpackidl.so','ssvdc',x,ldx,n,p,s,e,u,ldu,v,ldv,work,job,info)
;print,'INFO = ',info
return
end

