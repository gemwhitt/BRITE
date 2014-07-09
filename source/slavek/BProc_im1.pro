function BProc_im1,n,m
;
; test "object" defined here over n pix
; shifted from centre by m 
;
; use: im = BProc_im1(32,1)   ; for 32 pix long with one-pix shift to right


  y = fltarr(n)
  y[0:7] = findgen(n/4)
  y = shift(y,-20+m)        ; m=0: center, m=-1: one pix to left
  y = y/total(y)

return,y
end

   