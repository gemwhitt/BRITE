pro orion_psf_all

  ; progam to start to characterise the Orion PSFs using the sg1 dataset
  ;
  Compile_opt idl2
  !p.background=cgcolor('white')
  plotsym, 0, /fill, 1.5
  devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
  imagelib    ; adds system variable !IMAGE needed for plotting
  
  indir='~/BRITE/data/UB/p2/ORION/sg1/sav/'
  
  outdir='~/BRITE/data/UB/p2/ORION/sg1/models/'
  
  filesin=file_search(indir+'*.sav', count=nfl)
  
  for ff=0, nfl-1 do begin
  
    restore, filesin[ff]  ;vmag, hdname, jd, data1, ccd_temp
    
    jd0=jd
    data0=data1
    ccd_temp0=ccd_temp
    
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
    rbin=1.
    xy_psf=fltarr(2,nfrm) ; for all frames
    mod1=fltarr(xdim*rbin,ydim*rbin,ngap+1)
    mod2=fltarr(xdim*rbin,ydim*rbin)
    fl1=dblarr(nfrm)
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
      sr_fitting_unbin, data2,iloc,orbit,outdir,hdname, xy_psf,mod1,mod2,fl1,fl2,npix1,npix2
      
    endfor
    
    xx=where(avg_temp ge 18 AND avg_temp le 22, nxx, complement=yy)
    
    npix=npix1[xx]
    avg_temp1=avg_temp[xx]
    
    model=reform(mod1[*,*,0])*0.0
    
    for jj=0, nxx-1 do model=model+mod1[*,*,xx[jj]]
    
    model=model/float(nxx)
    
    model2=model/max(model)
    
   ; wset, 0
   ; plot_image, bytscl(model2, 0, 0.5)
    
    ; do data2-model - shift the model using xy_psf
   ; mod2=model2*0.0
   ; xymod=where(model2 gt 0)
   ; xmod=(array_indices(model2, xymod))[0,*]
   ; ymod=(array_indices(model2, xymod))[1,*]
   ; 
   ; xc=xdim/2.
   ; yc=ydim/2.
   ; xdif=xc-xy_psf[0,0]
   ; ydif=yc-xy_psf[1,0]
   ; 
   ; mod2[xmod-xdif,ymod-ydif]=model2[xmod,ymod]
   ; 
   ; wset, 1
   ; plot_image, bytscl(mod2, 0, 0.5)
   ; stop
    
    ; calculate 95% flux and number of pixels containing it
    totflx=round(total(model))
    flx=round(0.95*totflx)
    flx50=round(0.5*totflx)
    
    mod1d=reform(model, xdim*ydim)
    mod1d=mod1d[reverse(sort(mod1d))]
    tempflx=fltarr(xdim*ydim)
    for jj=0, xdim*ydim-1 do tempflx[jj]=total(mod1d[0:jj])
    xx=where(tempflx ge flx)  ; xx[0]=number of pixels containing flx
    
    xx50=where(tempflx ge flx50)  ; xx[0]=number of pixels containing flx
    
    yy=where(model2 gt 0.0, nyy)
    
    ; save output
    fileout=outdir+hdname+'_mod1.sav'
    save, filename=fileout, model, model2
    
    ; write out stats
    ;statfile=outdir+'models.txt'
    ;openw, lun, statfile, /get_lun, /append
    ;printf, lun, hdname, vmag, xx50[0], flx50, xx[0], flx, nyy, totflx, format='(a8,x,f7.3,x,i3,x,i8,x,i3,x,i8,x,i3,x,i8)'
    ;free_lun, lun 
    
    ; save x and y positions 
    jd=jd0
    data1=data0
    ccd_temp=ccd_temp0
    save, filename=filesin[ff], vmag, hdname, jd, data1, ccd_temp, xy_psf
    
   
  endfor
  
  print, 'end of program'
end
