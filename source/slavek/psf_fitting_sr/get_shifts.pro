pro get_shifts, dat,psf,nbin,xsize,ysize, shx,shy
  ; SMR Oct.2013
  ;
  ; determination of shifts by FFT for series of images
  ;     no rotation of the field!
  ;
  ; input:
  ;   dat = input image
  ;   psf = assumed PSF
  ; output:
  ;   shx,shy = shifts in pixels
  ; uses:
  ;   where2d.pro
  ; run:
  ;   for a single image i.e. 0th:
  ;      B2_im5, reform(dat1[*,*,0]),psf0, shx,shy
  ;   for a series of images:
  ;      B2_im5, dat1,psf0, shx,shy
  
  
  n = n_elements(dat[0,0,*])          ; how many images
  shx = fltarr(n) & shy = shx         ; shifts

  f0 = fft(psf,-1)                    ; FFT of the PSF

  for i=0,n-1 do begin
    im = reform(dat[*,*,i])
    f1 = fft(im,-1)                     ; FFT of the image
    y = float(fft(f1*conj(f0),+1))      ; CCF
    mx = max(y,mm)                      ; where max, mm in 2-D
    where2d, im, mm, wx,wy              ;    wx,wy idx in 1-D
    if wx gt xsize/2. then wx = wx - xsize
    if wy gt ysize/2. then wy = wy - ysize
    shx[i] = wx
    shy[i] = wy
   stop
  end
 
  return
end

