pro BProc_im3, n, obj, obs,sh
;    generation of n observations
; run:
;    BProc_im3, 100, obj, obs,sh

   sig_shift = 4.0    ; disp shifts, pixels
   sig_ph = 0.01      ; disp phot per small pix (1/8 big)

   sh = fltarr(n)
   obs = fltarr(n,32) ; 32 big pix long
   sh1 = 0.
   obs1 = fltarr(32)

   for i=0,n-1 do begin
      obs1 = BProc_im2(obj,4.0,0.01, sh1)  
      obs[i,*] = obs1
      sh[i] = sh1
   end

end

   