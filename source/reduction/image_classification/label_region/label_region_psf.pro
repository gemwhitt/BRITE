pro label_region_psf

; Program to determine frames with an adequate number of PSf pixels - use label_region.pro 
; Use p2 images - which have been mostly - but not entirely cleaned of HPs
;
Compile_opt idl2

sat='BA'
field='ORION'

indir='~/BRITE/'+sat+'/'+field+'/data/p2/'
filesin=file_search(indir+'*.sav', count=nf)

for f=0, nf-1 do begin
  
  print, file_basename(filesin[f], '_p2.sav')

  obj=obj_new('IDL_Savefile', filesin[f])
  obj->restore, 'jd'
  obj->restore, 'data1'
  obj->restore, 'flag'
  obj->restore, 'ccd_temp'
  
  jd1=jd-jd[0]
  
  nimg=n_elements(jd)
    
  ; set threshold for finding illuminated pixels
  thr=50
  
  nhp=intarr(nimg)
  npix_psf=intarr(nimg)
  
  for im=65, nimg-1 do begin
    ; skip over flag=0 images
    if flag[im] eq 0 then continue
    
    dat=data1[*,*,im]
    
    ; make all pixels in dat=1 equal to 0 instead
    xx=where(dat eq 1, nxx)
    if nxx gt 0 then begin
      twod=array_indices(dat, xx)
      dat[twod[0,*],twod[1,*]]=0
    endif
    
    s=size(dat, /dim)
    
    ; add a border to the image
    dat2=lonarr(s[0]+2,s[1]+2)
    dat2[1,1]=dat
    
    ; use label region to determine number of illuminated regions and number of PSF pixels    
    r2=label_region(dat2 ge thr, /all_neighbors)
    r2=r2[1:s[0],1:s[1]]                         ; trim off border     
                         
    hist2=histogram(r2, binsize=1, locations=loc2, reverse_indices=ri2)
    
    if n_elements(hist2) eq 1 then begin
      flag[im]=1
      continue
    endif
    
    ;xx=where(hist2 eq 1, nxx)
    ;plot_image, bytscl(dat, 20, 200)
    ;for ii=0, nxx-1 do begin
    ;  ind=ri2[ri2[xx[ii]]:ri2[xx[ii]+1]-1]
    ;  ind2=array_indices(dat, ind)
    ;  oplot, [ind2[0,*]], [ind2[1,*]], color=cgcolor('green'), psym=2
    ;endfor
        
    ; the number of results in hist2 are the number of regions (blobs) > thr, plus the background pixels
    ; so if n_elements(hist2) > 2 then - HP!
    ; record n_elements(hist2)-1 + number of pixels in second largest "blob" - i.e. PSF pixels
    nhp[im]=n_elements(hist2)-2
    
  ;  wset, 0
  ;  plot_image, bytscl(dat, 20, 400)
    ; now overplot 2nd biggest group of pixels
    ; sort hist1 in decreasing order
    sort2=reverse(sort(hist2))
    ind2=ri2[ri2[sort2[1]]:ri2[sort2[1]+1]-1]
    i2=array_indices(dat, ind2)
  ;  oplot, i2[0,*], i2[1,*], color=cgcolor('green'), psym=2
    
    npix_psf[im]=hist2[sort2[1]]  
  
  endfor  ; end loop over image
  
  ; remove zeros
  n2=npix_psf[where(npix_psf gt 5)]
  
  plot, n2, color=cgcolor('black'), psym=2
  plot, n2-smooth(n2, 45, /edge_truncate), color=cgcolor('purple'), psym=2
stop
  
endfor  ; end loop over file


print, 'end of program'
end