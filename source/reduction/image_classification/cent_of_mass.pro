pro cent_of_mass, data2,iloc, xy_cen

; Called by modelpsf_ba1
; Calculate the centroid coordinates for each "good" image 

nimg=n_elements(iloc)

for im=0, nimg-1 do begin
  
  im0=data2[*,*,iloc[im]]
  
  s=size(im0, /dim)
  
  imtot=total(im0)
  
  x_cen=total(total(im0, 2)*indgen(s[0]))/imtot
  y_cen=total(total(im0, 1)*indgen(s[1]))/imtot
  
  xy_cen[0,iloc[im]]=x_cen
  xy_cen[1,iloc[im]]=y_cen
  
endfor

end