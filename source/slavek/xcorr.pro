PRO xcorr, x1, y1, x2, y2, ccfx, ccfy

; Compute cross-correlation function for functions y1(x1) and y2(x2)
; The result is returned as ccfy(ccfx)

; NB.  Currently assumes that x1 and x2 are identical.

     nccf = n_elements(y1)

;    construct fast fourier transforms
     ff1 = fft(y1)
     ff2 = fft(y2)

;    convolve transforms
     ffxff = ff1 * conj(ff2) * nccf/8

;    ccf is real part of inverse transform
     ccfy = float (fft ( ffxff, +1 ))

;    swap order around
     work = ccfy (nccf/2:nccf-1)
     ccfy(nccf/2:nccf-1) = ccfy(0:nccf/2-1)
     ccfy(0:nccf/2-1)=work

;    compute ordinate
     dw = (x1(n_elements(x1)-1) - x1(0))  / (nccf-1)
     mw =  x1(nccf/2)
     ccfx = (findgen(nccf) - nccf/2) * dw/mw

;    scaling factors to produce identity with DIPSO
     ccfx = ccfx * 4 
     ccfy = ccfy / 2

END
