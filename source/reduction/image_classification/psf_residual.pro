pro psf_residual, psf,sdat2,shx,shy, res

nimg=n_elements(shx)


stop
xdim=(size(psf, /dim))[0]
ydim=(size(psf, /dim))[1]

for im=0, nimg-1 do begin
  
  img=sdat2[*,*,im]
  
  sx=abs(shx[im])
  sy=abs(shy[im])

  if shx[im] lt 0 then begin
  
    img1=img[0:xdim-1-sx,*]
    model=psf[sx:xdim-1,*]
    
  endif else begin
    
    img1=img[sx:xdim-1,*]
    model=psf[0:xdim-1-sx,*]
    
  endelse
    
  if shy[im] lt 0 then begin
  
    img1=img1[*,0:ydim-1-sy]
    model=model[*,sy:ydim-1]
    
  endif else begin
  
    img1=img1[*,sy:ydim-1]
    model=model[*,0:ydim-1-sy]
    
  endelse
  
  res=img1-model
  
stop

endfor


end