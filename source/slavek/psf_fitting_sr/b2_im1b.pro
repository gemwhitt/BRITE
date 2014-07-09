function B2_im1b, dat,cr1,cr2,cr3, imc,ww
  ; SMR Oct.2013
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
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; added by gemma
  ; FOR PLOTTING
  ;devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
  ;imagelib    ; adds system variable !IMAGE needed for plotting
  ;!p.background=cgcolor('white')
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  x = fltarr(32,32)                   ; average image
  
  m = n_elements(dat[0,0,*])          ; how many images
  for i=0,m-1 do x = x + dat[*,*,i]   ; adding images
  x = x/m                             ; average image
  
  med_x=median(x)
  x2=(lonarr(34,34)*0.)+med_x
  x2[1:32,1:32]=x
  
  ;---------------------------------------------------------------------
  ww = where(x2 gt cr1)      ; cr1 = criterion for high pixels
  where2d, x2, ww, wx,wy     ;    all hot
  n = n_elements(ww)        ;    1-st cut
  
  ; analyse area around bad pixels, find out how many elevated
  k1 = 0 & k2 = 0 & xx = fltarr(3,3) & nx = intarr(n)
  
  for i=0,n-1 do begin
    k1 = wx[i] & k2 = wy[i]
   ; if k1 eq 0 then  k1=1
   ; if k1 eq 31 then k1=30
   ; if k2 eq 0 then  k2=1
   ; if k2 eq 31 then k2=30
    xx = x2[k1-1:k1+1,k2-1:k2+1]     ; 3x3 neighbourhood
    mx = median(xx)
    j = where((xx-mx) gt mx, nj)       ; by cr2 above local median
    
    if nj eq 0 OR nj gt cr3 then nx[i]=1 else begin
      
      cj=where(j eq 4)  ; check that central pixel in square is the highest
      
      if cj[0] eq -1 then nx[i]=1
      
    endelse
    
  endfor
;  stop
  ww = ww[where(nx eq 0)]     ; bad pixels, 2-nd cut
  n = n_elements(ww)
  where2d, x2, ww, wx,wy
  
  ; trim of ww, wx, wy 
  wx=wx-1
  wy=wy-1
  ww=(wy*32)+wx
  
  
  ;---------------------------------------------------------------------------
  ; clean bad pixels in the average image and produce cleaned image imc
  
  imc = x                             ; cleaned image
  for i=0,n-1 do begin                ; 2-nd cut bad pix replaced
    k1 = wx[i] & k2 = wy[i]
    if k1 eq 0 then  k1=1
    if k1 eq 31 then k1=30
    if k2 eq 0 then  k2=1
    if k2 eq 31 then k2=30
    xx = x[k1-1:k1+1,k2-1:k2+1]     ; 3x3 neighbourhood
    mx = median(xx)
    imc[wx[i],wy[i]] = mx           ; median cleaning
  end
  
  return,x          ; average image
  ;stop
end

