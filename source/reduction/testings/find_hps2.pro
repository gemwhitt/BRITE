function find_hps2, dat,cr1,cr2,ii, ww,wx,wy
 
  ; get dimensions of frame
  xdim=(size(dat, /dim))[0]
  ydim=(size(dat, /dim))[1]
  
  x=dat
  
  x2=(lonarr(xdim+2,ydim+2)*0.)+median(x)
  x2[1:xdim,1:ydim]=x
  
  ;---------------------------------------------------------------------
  ww = where(x2 gt cr1, nww)      ; cr1*img_med = criterion for high pixels
  where2d, x2, ww, wx,wy     ;    all hot
  n = n_elements(ww)        ;    1-st cut
  
;  plot_image, bytscl(x2, 20, 500)
;  oplot, wx, wy, color=cgcolor('orange'), psym=2

  if nww eq 0 then begin
    undefine, ww
    goto, nohp
  endif
  
  ; analyse area around bad pixels, find out how many elevated
  nx = intarr(n)
  
  for i=0, n-1 do begin
    k1 = wx[i] & k2 = wy[i]
   
  ; print, k1-1, k2-1 
;    if k2 ge 14 then stop
    
    xx = x2[k1-1:k1+1,k2-1:k2+1]     ; 3x3 neighbourhood
    mx = median(xx)                  ; local median - should be < cr1
    
    totxx=total(xx)
    avgdn=float(totxx)/9.
   
    j = where((xx-mx) ge cr1, nj)       ; by cr2 above local median
    
    cj=where(j eq 4)  ; check that central pixel in square is the highest
;  stop  
    if nj eq 0 OR cj[0] eq -1 OR nj gt cr2 then begin
      nx[i]=1
      goto, next_hp
    endif

    if nj gt 1 AND nj le cr2 then begin ; there may be 3 HPs
    
      if mx le cr1 then goto, next_hp else nx[i]=1
            
    endif
    
    next_hp:
;    if k2 ge 7 then stop
  endfor
  
  xx=where(nx eq 0, nhp)
  
  if nhp eq 0 then begin
    undefine, ww
    undefine, wx
    undefine, wy
    goto, nohp
  endif
  
  ww = ww[where(nx eq 0)]     ; bad pixels, 2-nd cut
  n = n_elements(ww)
  where2d, x2, ww, wx,wy
  
  ; trim of ww, wx, wy
  wx=wx-1
  wy=wy-1
  ww=(wy*ydim)+wx
  
  nohp:
  return,x          ; average image
end

