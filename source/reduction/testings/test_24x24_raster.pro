pro test_24x24_raster

  ; created on 15/1/13 - to calculate percentage of points which would be lost with 24x24 image rasters
  
  ;
  Compile_opt idl2
  
  ;
  ; FOR PLOTTING
  devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
  imagelib    ; adds system variable !IMAGE needed for plotting
  !p.background=cgcolor('white')
  ;
  ;
  indir='~/BRITE/data/UB/p3/'
  filesin=file_search(indir+'HD33111*.sav', count=nfiles)
  fname=file_basename(filesin, '_p3.sav')
  
  outdir='~/BRITE/data/UB/p4_24x24/'
  
  for bb=0, nfiles-1 do begin
    ;stop
    print, fname[bb]
    print, 'Start time ', systime()
    
    missed=0
    good=0
    
    restore, filesin[bb] ;roi_name, exp_num, ra_dec1, jd, data1, roi_dim, xc, yc, ccd_temp, simbad_radec, $
    ;simbad_mags, parlax, otype, sptype, p2_trend, med_trend
    
    ; print the magnitude
    vmag=strmid(simbad_mags, strpos(simbad_mags, 'V=')+2)
    print, 'vmag='+vmag
    ;  stop
    
    dark_current=0. ; change this when dark current is known  ; from Jake's MSC thesis
    ;- tests show that dark current is negligible at 1s exposure
    
    rn=0  ; 15 photoelectrons per pixel
    
    dat=data1 ; no effect to data1
    
    ; get total number of frames in this file
    nfrm=(size(dat, /dim))[2]
    
    flux=fltarr(nfrm)
    backgd=fltarr(nfrm)
    backer=fltarr(nfrm)
    mean_resid=fltarr(nfrm)
    sig_resid=fltarr(nfrm)
    dresid=fltarr(nfrm)
    nsubpix=intarr(nfrm)
    xy_psf=fltarr(2,nfrm)
    
    ;window, 0, xsize=500, ysize=450, xpos=100, ypos=200
    ;window, 1, xsize=500, ysize=450, xpos=700, ypos=200
    ;window, 2, xsize=500, ysize=450, xpos=1300, ypos=200
    ;    stop
    for img=0, nfrm-1 do begin
    
      ;print, img
      
      dat1=dat[*,*,img]
      
      jd1=jd[img]-jd[0]
      ;if jd1 gt 5.5202 and jd1 lt 5.57 then stop
      
      ; find approx center of PSF
      mthr=2.5*median(dat1)
      spr=2
      pks=brite_findpeaks(dat1, minthresh=mthr, spread=spr)
      
      if pks[0] eq -999.9 then CONTINUE
      
      npks=float(n_elements(pks))/3.
      
      if npks gt 1 then begin
        ;stop  ; continue
        brightest=(where(pks[2,*] eq max(pks[2,*])))[0]
        pks=pks[*,brightest]
      endif
      
      pks1=round(pks)
      
      ; determine the background - away from the PSF center - subtract this from the image
      rbin=1
      get_bkgd, dat1,pks1,rbin, bkgd,bk_err
      
      ; subtract the background measurement from the image
      dat1=dat1-bkgd
      
      ; record background measurements
      backgd[img]=bkgd
      backer[img]=bk_err
      
      xdim=(size(dat1, /dim))[0]
      ydim=(size(dat1, /dim))[1]
      
      ; make cutout - centered on PSF center
      if pks1[0]-(8) lt 0 then x1=0 else x1=pks1[0]-(8)
      if pks1[0]+(8) gt xdim-1 then x2=xdim-1 else x2=pks1[0]+(8)
      if pks1[1]-(8) lt 0 then y1=0 else y1=pks1[1]-(8)
      if pks1[1]+(8) gt ydim-1 then y2=ydim-1 else y2=pks1[1]+(8)
      cutout0=dat1[x1:x2,y1:y2]
      ; stop
      ;wset, 0
      ;plot_image, bytscl(cutout0, 50, 500)
      
      ; rebin the cutout - 4x4, 8x8, 16x16...
      rbin=8. ;8.
      cutout=rebin_data(cutout0,rbin)
      
      cutout2=cutout
      
      for i=0, (size(cutout0, /dim))[0]-1 do begin      ; column
        for j=0, (size(cutout0, /dim))[1]-1 do begin    ; row
        
          cutout2[i*rbin:((i+1)*rbin)-1,j*rbin:((j+1)*rbin)-1]=cutout0[i,j] ; this array is equivalent in dimension to cutout
          
        endfor
      endfor
      
      cutout2=cutout2/(rbin^2)
      
      ; find center of PSF again...
      ;mthr=2.5*median(cutout)
      ;spr=2*rbin
      ;pks2=brite_findpeaks(cutout, minthresh=mthr, spread=spr)
      ;if pks2[0] eq -999.9 then stop ;CONTINUE
      
      ; npks=float(n_elements(pks2))/3.
      
      ; if npks gt 1 then begin
      ;stop  ; continue
      ;   brightest=(where(pks2[2,*] eq max(pks2[2,*])))[0]
      ;   pks2=pks2[*,brightest]
      ; endif
      
      ; get location of PSF center
      xy_psf[0,img]=pks[0]  ;(pks2[0]/rbin)+x1
      xy_psf[1,img]=pks[1]  ;(pks2[1]/rbin)+y1
      
      ;wset, 1
      ;plot_image, bytscl(cutout, 50, 5000)
      ;oplot, ([pks[0]]-x1)*rbin, ([pks[1]]-y1)*rbin, psym=2, color=cgcolor('purple')
      
      
      ncut=n_elements(cutout)
      
      ;reorder pixels by data number
      order1=reverse(sort(cutout))
      temp_arr=order1
      signal1=cutout[order1]
      
      signal2=cutout2[order1]
      
      count=2
      snr1=0
      temp_arr[0:count]=-9999
      snr_cum=[]
      ;stop
      repeat begin
      
        count=count+1
        
        if count ge float(ncut)/2.  then goto, endofloop
        
        snr=snr1
        
        sig=total(signal1[0:count])
        
        ;snr1=sig/sqrt(sig + (count+1)*((bkgd + dark_current)^2))   ; original
        ;snr1=sig/sqrt(signal1[count] + (count+1)*((bkgd*bk_err)^2))
        ;gn=3.25   ; photoelectron per ADU
        snr1=sig/sqrt(signal1[count] + bkgd*count)
        snr_cum=[snr_cum,snr1]
        
        temp_arr[count]=-9999
        
      endrep until snr1 lt snr      ; temp removed ;count eq ncut-1
      
      endofloop:
      if count ge float(ncut)/2. or count lt 600 then CONTINUE
      
      nsubpix[img]=count
      
      ; get mean and scatter of residuals in the ROI....
      temp_dat2=dat1
      temp_dat2[x1:x2,y1:y2]=-9999
      
      respix=[signal1[count+1:ncut-1]]
      
      mean_resid[img]=robust_mean(respix,2)
      sig_resid[img]=robust_sigma(respix)
      dresid[img]=max(respix)-min(respix)
      
      ;flux[img]=sig  ; old - wrong
      flux[img]=total(signal2[0:count])
      
      flux2=total(signal2[0:count])
      
      xx=where(temp_arr eq -9999)
      
      loc_2d=array_indices(cutout, order1)
            
      loc_rb=float(loc_2d[*,xx])/8.
      loc_img=loc_rb
      loc_img[0,*]=loc_img[0,*]+x1
      loc_img[1,*]=loc_img[1,*]+y1
      
      psf_dim=[min(loc_img[0,*]), max(loc_img[0,*]), min(loc_img[1,*]), max(loc_img[1,*])]
      
      
      ; check psf_dim against hypothetical roi margins
      outside=where(psf_dim lt 4. OR psf_dim gt 28, nout)
      
      if nout gt 0 then stop  ;missed=missed+1 else good=good+1
      
      
      ;wset, 0
      ;plot_image, bytscl(cutout, 50, 500), title=hdname, charsize=0.7
      ;oplot, loc_2d[0,xx], loc_2d[1,xx], psym=4, color=cgcolor('purple'), symsize=0.8
      ;stop
      
    endfor ; END LOOP OVER THIS SUB-SECTION OF ORBITS - IS USING
 ; stop  
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    print, 'missed=',missed
    print, 'good=',good
    
    ;save, filename=outdir+fname[bb]+'_8bin.sav', flux, jd, bkgd, $
     ; exp_num, roi_dim, ccd_temp, xy_psf, dresid, counts, vmag
      
  endfor  ; end loop over file
  
  print, 'end of program'
  
end