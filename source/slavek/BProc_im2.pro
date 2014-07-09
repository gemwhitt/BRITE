function BProc_im2,obj,s1,s2, sh
;
; simulation of an observation, assume the test "object" sampled 8x better
;    than the "observed" one
; 
; takes the theoretical test image obj (256 pix) and 
;     1. shifts by sh big (32) pixels (sh>0, right), Gaussian distributed (disp=s1)
;     2. adds Gauss noise with disp=s2 std error per small pixel (256)
;            now noise only where the image is present, not in the sky 
;
; use: obs = BProc_im2(obj,4.0,0.01, sh)  
;            for disp~=size of image (8 big, 64 small pix)
  
   n = n_elements(obj)         ; i.e. 256

R:
   g = randomn(seed,n+1)       ; errors, 0-th used for the shift, sig=1
   g0 = g[0]
   sh = s1*g0               ; shift in big (32) 
   if (abs(sh) gt 2.5*s1) then goto,R

   m = round(sh*8.)        ; shift in small pixels
   y = shift(obj,m)        ; m=0: center, m=+1: 1-pix right, m=-1: 1-pix left
   
; now noise only within the image, this can be removed
   w = where(y gt 0.)   ; comment out this line if noise in the sky as well
;   w = indgen(n)        ; then un-comment this line

   nw = n_elements(w)
   y[w] = y[w] + s2*g[1:nw]    ; noise added to small (256) pixels

   y1 = 8.*rebin(y,32)          ; rebinned to 32 pix

return,y1
end

   