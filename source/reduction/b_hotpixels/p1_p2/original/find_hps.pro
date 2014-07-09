function find_hps, dat,cr1,cr2, ww,wx,wy
  ; SMR Oct.2013
  ; 
  ; Modified from b2_im1b
  ;
  ; creates an average image in pixel coordinates and finds bad pixels
  ;    produces cleaned image and idx values where bad pixels
  ;
  ; criteria can be changed, here just an example
  ;
  ; input:
  ;    dat = 3-D stack of images, last index is the image number
  ;    cr1 = first criterion for rejection, just high pixel
  ;    cr2 = above local median
  ;    cr3 = how many bad pixels
  ; output:
  ;    ima = result of function, suggested variable name
  ;              this average image is not cleaned
  ;    imc = cleaned image
  ;    ww = indices of bad pixels in 2-D (linear count)
  ; uses:
  ;    where2d.pro     for conversion: linear count -> 2-D count
  ; run:
  ;    ima = B2_im1(dat,150.,100.,2, imc,ww)
  
  ; get dimensions of frame
  xdim=(size(dat, /dim))[0]
  ydim=(size(dat, /dim))[1]

  x = fltarr(xdim,ydim)               ; average image
  
  m = n_elements(dat[0,0,*])          ; how many images
  for i=0,m-1 do x = x + dat[*,*,i]   ; adding images
  x = x/float(m)                             ; average image
  
  img_med=median(x)
  img_sig=robust_sigma(x)

  x2=(lonarr(xdim+4,ydim+4)*0.)+img_med
  x2[2:xdim+1,2:ydim+1]=x
    

  ;---------------------------------------------------------------------
 ; ww = where(x2 gt img_med+(cr1*img_sig), nww)      ; cr1*img_med = criterion for high pixels
  ww = where(x2 gt cr1, nww)      ; cr1*img_med = criterion for high pixels
  where2d, x2, ww, wx,wy     ;    all hot
  n = n_elements(ww)        ;    1-st cut
  
  if nww eq 0 then stop
  
  ; analyse area around bad pixels, find out how many elevated
  nx = intarr(n)
  
  for i=0,n-1 do begin
    k1 = wx[i] & k2 = wy[i]
        
    xx = x2[k1-1:k1+1,k2-1:k2+1]     ; 3x3 neighbourhood
    ord1=reverse(sort(xx))
    mx = median(xx[ord1[0:7]])
    sigxx=robust_sigma(xx[ord1[0:7]])
    
    j = where((xx-mx) gt 5*sigxx, nj)       ; by cr2 above local median
        
    cj=where(j eq 4)  ; check that central pixel in square is the highest
            
    if nj eq 0 OR cj[0] eq -1 OR nj gt 5 then begin
      nx[i]=1
      goto, next_hp
    endif
  
  ; check j[4] against median of j's    
    if nj gt cr2 then begin
      ; define a larger area to test....
      xx = x2[k1-2:k1+2,k2-2:k2+2]     ; 3x3 neighbourhood
      mx = median(xx)
      sigxx=robust_sigma(xx)
      j = where((xx-mx) gt 3*sigxx, nj)       ; by cr2 above local median
        
      cj=where(j eq 12)  ; check that central pixel in square is the highest
            
      if nj eq 0 OR cj[0] eq -1 then begin
        nx[i]=1
        goto, next_hp
      endif 
    endif
    
  next_hp:
      
  endfor
    
  xx=where(nx eq 0, nhp)
  
  if nhp eq 0 then begin
    undefine, ww
    goto, nohp
  endif

  ww = ww[where(nx eq 0)]     ; bad pixels, 2-nd cut
  n = n_elements(ww)
  where2d, x2, ww, wx,wy
  
  ; trim of ww, wx, wy
  wx=wx-2
  wy=wy-2
  ww=(wy*ydim)+wx
  
  nohp:
  return,x          ; average image
end

