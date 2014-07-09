function BProc_obj,n
; 
; generates a test triangular, centered "object" of length 1/4*n in a vector
; with the total length of n pix
; 
; use: obj = BProc_obj(256)

   y = fltarr(n)
   y[0:63] = findgen(n/4)
   y = shift(y,-5*n/8)
   y = y/total(y)

return, y
end

   