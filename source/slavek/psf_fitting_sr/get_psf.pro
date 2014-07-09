function get_psf, dat,shx,shy,nbin,xsize,ysize, dud,bck,w
  ; SMR Oct.2013
  ;
  ; determination of PSF from the average image, after the shifts are applied
  ;
  ; input:
  ;    dat = input images for averaging after shifts, 256x256
  ;    shx,shy = shifts in small pixels
  ; output:
  ;    psf1 = result of function
  ;    dud = idx array with images NOT used marked by 1
  ;          criterion of rejection: shifts larger than 1/4 window = 64 small pix
  ;    bck = background value which was subtracted
  ; run:
  ;   psf1 = B2_im6(dat1,shx,shy, dud)
  
;  stop
  n = n_elements(dat[0,0,*])          ; how many images
  x = fltarr(xsize,ysize) & y = x
  dud = intarr(n)
  
  for i=0,n-1 do $
    if (sqrt((shx[i])^2 + (shy[i])^2) gt nbin*8.) then dud[i]=1   ; arbitrary
    
  w = where(dud eq 0, m)                 ; good images
  
;  if (m ne n) then print,'Rejected ',n-m,' images'
;  stop
  for i=0,m-1 do begin
    y = reform(dat[*,*,w[i]])
    x = x + shift(y,-fix(shx[w[i]]),-fix(shy[w[i]]))
  end
  x = x/m                             ; average image
  
  ;r1 = nbin*32./3.    ; inner and outer ring for background
  ;r2 = r1*1.5
  ;t = findgen(nbin*32.)-nbin*32./2.
  ;rx = t#replicate(1.,nbin*32.)
  ;ry = replicate(1.,nbin*32.)#t
  ;r = sqrt(rx^2+ry^2)       ; distance in pixels from centre
  ;w = where((r ge r1) and (r le r2))
  
  ;xw=(array_indices(findgen(256,256), w))[0,*]
  ;yw=(array_indices(findgen(256,256), w))[1,*]
  
  ;temp_x=x
  ;temp_x[0*nbin:(2*nbin)-1,*]=-9999 & temp_x[(xsize-2)*nbin:(xsize*nbin)-1,*]=-9999 & temp_x[*,0*nbin:(2*nbin)-1]=-9999 & temp_x[*,(ysize-2.)*nbin:(ysize*nbin)-1]=-9999
  
  ;bk_pix=array_indices(temp_x, where(temp_x eq -9999))
  
  ;plot_image, x
  ;oplot, bk_pix[0,*], bk_pix[1,*], color=cgcolor('orange'), psym=1

  ;print,'Background from: ',n_elements(w),' points'  ; can be commented out
  
  ;bck = median(x[bk_pix])
  ;x = x - bck
  
;  plot_image, x
;  oplot, bk_pix[0,*], bk_pix[1,*], color=cgcolor('orange'), psym=1
  
  return,x
end

