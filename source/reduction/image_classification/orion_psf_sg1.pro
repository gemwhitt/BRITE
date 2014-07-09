pro orion_psf_sg1

; progam to start to characterise the Orion PSFs using the sg1 dataset
; 
Compile_opt idl2
!p.background=cgcolor('white')
plotsym, 0, /fill, 1.5
devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
imagelib    ; adds system variable !IMAGE needed for plotting

indir='~/BRITE/data/UB/p2/ORION/sg1/sav/'

filesin=file_search(indir+'*.sav', count=nfl)

for ff=0, nfl-1 do begin
  
  restore, filesin[ff]  ;vmag, hdname, jd, data1, ccd_temp
 
  ; select files < 10 days
  jd1=jd-jd[0]
  
  xx=where(jd1 lt 10., nxx)
  
  jd=jd[xx]
  jd1=jd1[xx]
  data1=data1[*,*,xx]
  ccd_temp=ccd_temp[xx]
  
  nfrm=nxx
  
  xdim=(size(data1, /dim))[0]
  ydim=(size(data1, /dim))[1]
  
  good_bad=lonarr(nfrm)
  
  jd2=jd1[1:nfrm-1]
  jdiff=jd2-jd1
  gap=where(jdiff gt 0.015, ngap) ; ngap+1 orbits in this sequence
  gap=[-1,gap,nfrm-1]
  
  ; variables to determine in sr_fitting_sg1
  rbin=10.
  xy_psf=fltarr(2,nfrm) ; for all frames 
  mod1=fltarr(xdim*rbin,ydim*rbin,ngap+1)
  mod2=fltarr(xdim*rbin,ydim*rbin)  
  fl1=dblarr(1)    
  fl2=dblarr(1)   
  npix1=intarr(ngap+1)
  npix2=intarr(nfrm)
  avg_temp=fltarr(ngap+1)
  
  for orbit=0, ngap do begin
    
    iloc=indgen(gap[orbit+1]-gap[orbit])+gap[orbit]+1
        
    if n_elements(iloc) lt 2 then good_bad[iloc]=1
    if n_elements(iloc) lt 2 then continue
    
    avg_temp[orbit]=median(ccd_temp[iloc])
    
    data2=data1[*,*,iloc]
      
    ; get the PSF centers
    sr_fitting_sg1, data2,iloc,orbit,outdir,hdname, xy_psf,mod1,mod2,fl1,fl2,npix1,npix2
                
  endfor
  
xx=where(avg_temp ge 18 AND avg_temp le 22, nxx, complement=yy)

npix=npix1[xx]
avg_temp1=avg_temp[xx]

model=reform(mod1[*,*,0])*0.0

for jj=0, nxx-1 do model=model+mod1[*,*,xx[jj]]

model=model/float(nxx)

model2=model/max(model)

;wset, 0
;plot_image, bytscl(model2, 0, 0.5)

sigres=fltarr(nfrm)

for im=0, nfrm-1 do begin
  
  data2=reform(data1[*,*,im])
  data2=rebin_data(data2, rbin)
  data2=data2/max(data2)
  
  ;wset, 1
  ;plot_image, bytscl(data2, 0, 0.5)
  
  xc=xy_psf[0,im]
  yc=xy_psf[1,im]
  
  ; shift the model to match the center of the image PSF
  ccd_xc=(size(data2, /dim))[0]/2.
  ccd_yc=(size(data2, /dim))[1]/2.
  
  mod2=model2*0.0 ; new model - shifted
  ; get index locations of model PSF in model2 - ixpsf and iypsf
  ixpsf=(array_indices(model2, where(model2 gt 0)))[0,*]
  iypsf=(array_indices(model2, where(model2 gt 0)))[1,*]
  shx=xc-ccd_xc
  shy=yc-ccd_yc
  
  mod2[ixpsf+shx,iypsf+shy]=model2[ixpsf,iypsf]
  
  ;wset, 0
  ;plot_image, bytscl(mod2, 0, 0.5)
  
  res=data2-mod2
  
  sigres[im]=robust_sigma(res)

endfor

wset, 0
plot, ccd_temp, sigres, color=cgcolor('black'), psym=2

wset, 1
plot, jd1, sigres, color=cgcolor('black'), psym=2

stop
  
endfor

print, 'end of program'
end
