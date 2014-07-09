function B2_im7b, psf,dat,shx,shy,dud, e, yfit2
  ; SMR Oct.2013
  ;
  ; determination of the flux and background using an experimental PSF
  ;
  ; input:
  ;    psf = PSF used for fits
  ;    dat = input images for averaging after shifts, 256x256
  ;    shx,shy = shifts in small pixels
  ;    dud = idx array with images which should NOT be used marked by 1
  ; output:
  ;    result of function: [flux, background]
  ;          coefficients of fit as 2-column array
  ;    e = errors to coeff.
  ; run:
  ;   fit = B2_im7(psf3,dat1,shx2,shy2,dud, fite)
  
  n = n_elements(dat[0,0,*])     ; how many images in sequence
  z = fltarr(2,n)                ; 2-col array of coefficients
  e = z                          ; 2-col array of errors
  e1 = [0.,0.]
  yfit1=1
  
  x = psf
  m = n_elements(x)   ; linear pixels in 2-D PSF image
  x = reform(x,m)     ; PSF converted from 2-D to 1-D for SVDfit.pro
  
  ;print,'Execution may be slow...'
  
  for i=0,n-1 do begin
  
    ; print,f='(a,i4)','Image: ',i
    if dud[i] eq 1 then goto,E     ; both coefficients set to 0.
    y = reform(dat[*,*,i])         ; take one image
    y = shift(reform(dat[*,*,i]),-shx[i],-shy[i])
    ; i-th 2-D image shifted to origin
    y1 = reform(y,m)     ; converted to a vector
    z1 = svdfit(x,y1,2,f='svd_psf',sig=e1, /double, yfit=yfit1)
    ; linear function fit: "svd_psf.pro"
    z[*,i] = z1          ; 2-el vector: [intensity,background]
    e[*,i] = e1          ; 2-el vector: [err int, err bckg]
    E:
  end
  
  yfit2=reform(yfit1, nbin*32, nbin*32)
  
  return,z
end

