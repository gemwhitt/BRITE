pro clean_hps, data1, hpmap2, flag

s=size(data1, /dim)

nimg=s[2]

for im=0, nimg-1 do begin
  if flag[im] eq 0 then continue
  
  hpmap=hpmap2[*,*,im]
  
  hp=where(hpmap eq 1, nhp)
  
  if nhp eq 0 then continue
  
  im0=data1[*,*,im]
  im1=intarr(s[0]+2,s[1]+2)
  im1[1:s[0],1:s[1]]=im0
  
  for h=0, nhp-1 do begin
    
    xi=(array_indices(hpmap, hp))[0,h]+1
    yi=(array_indices(hpmap, hp))[1,h]+1
    
    x1=xi-1
    x2=xi+1
    y1=yi-1
    y2=yi+1
    
    cutout=im1[x1:x2,y1:y2]
    
    im1[xi,yi]=median(cutout)
    if xi lt s[0] then im1[xi+1,yi]=median(im1[(x1+1):(x2+1),y1:y2])
    
  endfor
  
  im3=im1[1:s[0],1:s[1]]
  
  
;  plot_image, bytscl(im0, 0, 500)
;  stop
;  data1[*,*,im]=im3
  ;stop
endfor

end