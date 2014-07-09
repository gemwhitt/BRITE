pro BProc_im4, obj,obs, sh1
;    determination of shifts by FFT
; run:
;    BProc_im4, obj,obs, sh1


   n = n_elements(obs[*,0])   ; number of observations
   sh1 = fltarr(n)

   f0=fft(obj,-1)             ; FFT of template for now 

   for i=0,n-1 do begin
      x = reform(obs[i,*])
      x = rebin(x,256)/8.     ; same size as template, normalized area
      f1 = fft(x,-1)          ; FFT 
      y = float(fft(f1*conj(f0),+1))      ; CCF
      mx = max(y,mm)                      ; where max
      if mm gt 128 then mm = mm - 256
      sh1[i] = mm/8.          ; in big pix
   end

   sh1 = sh1 + 0.5            ; systematic shift by 1/2 big pix 
                              ; why this bias?
return
end

   