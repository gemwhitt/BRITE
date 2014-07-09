function rebin_data, dat, nbins
  ; SMR Oct.2013
  ;
  ; re-bins into series of 256x256 pix clean images - MODIFICATION OF B2_im5.pro (Slavek)
  ;
  ; input:
  ;   dat = input data
  ;   nbins = number of bins to use
  ;
  ; output:
  ;   dat1 = rebinned data, note: the total changed 64x
  ; run:
  ;   dat1 = B2_im4(dat,rbin[hh])
  
  n = n_elements(dat[0,0,*])          ; how many images
  nrows=(size(dat, /dim))[1]
  ncols=(size(dat, /dim))[0]
  x = fltarr(ncols*nbins,nrows*nbins,n)
  
  for i=0,n-1 do begin
    y = reform(dat[*,*,i])
    y = rebin(y,ncols*nbins,nrows*nbins)       ;      and re-binned
    x[*,*,i] = y
  end
  
  return,x
end

