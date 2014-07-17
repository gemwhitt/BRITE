pro p3_p4_gem1

; modified from snr_investigate_target
; produce p4 light curves
; this version treats saturated stars and gives 4 versions of light curve using a ratio of snr/xx
  
; 31 Jan 2014
;
Compile_opt idl2
!p.background=cgcolor('white')
plotsym, 0, 0.8, /fill
devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
imagelib    ; adds system variable !IMAGE needed for plotting
 
indir='~/BRITE/data/UB/p3/'
filesin=file_search(indir+'HD*.sav', count=nfiles)
fname=file_basename(filesin, '_p3.sav')
  
outdir='/Users/gemmawhittaker/BRITE/data/UB/p4_gemphot1/'
  
for bb=0, nfiles-1 do begin
  
  ;stop
  print, fname[bb]
    
  restore, filesin[bb] ;roi_name, exp_num, ra_dec1, jd, data1, roi_dim, xc, yc, ccd_temp, simbad_radec, $
  ;simbad_mags, parlax, otype, sptype, p2_trend, med_trend
    
  jd1=jd-jd[0]
  keep=where(jd1 ge 0.0)
  
  jd1=jd1[keep]
  data1=data1[*,*,keep]
  roi_name=roi_name[keep]
  exp_num=exp_num[keep]
  ra_dec1=ra_dec1[*,keep]
  ccd_temp=ccd_temp[keep]
  p2_trend=p2_trend[keep]
  
  dat=data1 ; no effect to data1
    
  jd2=jd1[1:n_elements(jd1)-1]
  jdiff=jd2-jd1
  cadence=robust_mean(jdiff,2)  ; in days
  cadence=cadence*24.*60.       ; in minutes
  n_img_per_orbit=fix(15./cadence)
    
  ; print the magnitude
  vmag=strmid(simbad_mags, strpos(simbad_mags, 'V=')+2)
  print, 'vmag='+vmag
    
  ;get total number of frames in this file
  nfrm=(size(dat, /dim))[2]
  
  ; get x and y pixel dimensions
  xdim=(size(dat, /dim))[0]
  ydim=(size(dat, /dim))[1]
  
  ; set up arrays for collecting results
  nsubpix=intarr(4,nfrm)      ; number of subpixels in PSF
  max_snr=fltarr(4,nfrm)      ; SNR at nsubpix
  max_dn=fltarr(nfrm)         ; maximum data number on image        
  xy_psf=fltarr(2,nfrm)       ; x and y coords of PSF center
  flux1=fltarr(4,nfrm)         ; flux array - using signal1/64.
  flux2=fltarr(4,nfrm)         ; flux array - using signal2
  nsat=intarr(nfrm)           ; number of saturated pixels - left out of flux determination
    
  for img=0, nfrm-1 do begin  ;nfrm-1 do begin
    
    dat1=dat[*,*,img]
      
    ; find approx center of PSF
    mthr=2.5*median(dat1)
    spr=2
    pks=brite_findpeaks(dat1, minthresh=mthr, spread=spr)
      
    if pks[0] eq -999.9 then CONTINUE
      
    npks=float(n_elements(pks))/3.
      
    if npks gt 1 then begin
      brightest=(where(pks[2,*] eq max(pks[2,*])))[0]
      pks=pks[*,brightest]
    endif
    
    ; use centroid function to get new x and y location?            ??????
    ; centroid, dat1, XCEN, YCEN, XY_PEAK=[pks[0],pks[1]], FWHM=1.5, /PEAK_LOC
      
    xy_psf[0,img]=pks[0]  ;xcen
    xy_psf[1,img]=pks[1]  ;ycen
      
    pks1=round(pks) ; round of pks to nearest whole pixel
      
    ; determine the background - away from the PSF center - subtract this from the image
    rbin=1
    get_bkgd, dat1,pks1,rbin, bkgd,bk_err
    ; subtract the background measurement from the image
    dat1=dat1-bkgd
   
    ; make cutout - centered on PSF center
    if pks1[0]-(8) lt 0 then x1=0 else x1=pks1[0]-(8)
    if pks1[0]+(8) gt xdim-1 then x2=xdim-1 else x2=pks1[0]+(8)
    if pks1[1]-(8) lt 0 then y1=0 else y1=pks1[1]-(8)
    if pks1[1]+(8) gt ydim-1 then y2=ydim-1 else y2=pks1[1]+(8)
    cutout0=dat1[x1:x2,y1:y2]
      
    ; rebin the cutout - 4x4, 8x8, 16x16...
    rbin=8. ;8.
    cutout1=rebin_data(cutout0,rbin)
      
    cutout2=cutout1
      
    for i=0, (size(cutout0, /dim))[0]-1 do begin      ; column
      for j=0, (size(cutout0, /dim))[1]-1 do begin    ; row
        
        cutout2[i*rbin:((i+1)*rbin)-1,j*rbin:((j+1)*rbin)-1]=cutout0[i,j] ; this array is equivalent in dimension to cutout
          
      endfor
    endfor
    
    ; deal with saturated pixels
    sat=where(cutout2 gt 9000, numsat, complement=lin)
      
    if numsat gt 0 then begin
      nsat[img]=numsat/(rbin^2)
      ;  print, nsat
      cutout1[sat]=-9999
    endif else nsat[img]=0
    
    ; divide pixel counts by e.g. 64 to add together later for flux values  
    cutout2=cutout2/(rbin^2)
      
    ncut=n_elements(cutout1)
      
    ;reorder pixels by data number
    order1=reverse(sort(cutout1))
    temp_arr=order1
    signal1=cutout1[order1]
    signal2=cutout2[order1]
      
    count=2 ; start on count=2, not count=0
    snr1=0
    snr_cum=[]
    sig_cum=[]
      
    repeat begin
      
      count=count+1
        
      snr=snr1
        
      sig=total(signal1[0:count])
      sig_cum=[sig_cum,sig]
        
      snr1=sig/sqrt(signal1[count] + (bkgd*count))  ; testing   ;2  - less noise
      snr_cum=[snr_cum,snr1]
        
    endrep until count eq long(8000) ;snr1 lt snr
      
    ;endofloop:   ; only need this if using endrep until snr1 lt snr
      
    xx=where(snr_cum eq max(snr_cum), nxx) ; location of maximum SNR - i.e. number of subpixels
    
    if nxx gt 1 then xx=xx[nxx-1] ; choose last snr - highest number of pixels
    
    if xx lt 500 OR xx gt 7000. then continue ; reject this point
      
    max_snr[0,img]=snr_cum[xx]
    
    ; for plotting subpixels       
    loc_2d=array_indices(cutout1, order1[0:xx])
    resid=array_indices(cutout1, order1[xx+1:ncut-1])
      
    xpix=loc_2d[0,*]
    ypix=loc_2d[1,*]
      
    ;window, 0, xsize=500, ysize=500, xpos=300, ypos=200
    ;plot_image, bytscl(cutout1, 50, 300), color=cgcolor('black')
    ;oplot, xpix, ypix, psym=8, color=cgcolor('purple')
    ;oplot, resid[0,*], resid[1,*], psym=8, color=cgcolor('green')  
    ;stop
      
    ; save other variables  
    nsubpix[0,img]=xx
      
    flux1[0,img]=total(signal1[0:xx])
    flux2[0,img]=total(signal2[0:xx])
      
    max_dn[img]=signal1[0]
      
    ; GET RATIO - max_snr/nsubpix
    ratio=max_snr[img]/float(xx)
    scale1=1+(ratio/100.)
    scale2=1+((2*ratio)/100.)
    scale3=1+((4*ratio)/100.)
    
    ; get 3 more locations for cutoffs of nsubpix
    yy1=round(scale1*xx)
    yy2=round(scale2*xx)
    yy3=round(scale3*xx)
      
    nsubpix[1,img]=yy1
    nsubpix[2,img]=yy2
    nsubpix[3,img]=yy3
      
    max_snr[1,img]=snr_cum[yy1]
    max_snr[2,img]=snr_cum[yy2]
    max_snr[3,img]=snr_cum[yy3]
      
    flux1[1,img]=total(signal1[0:yy1])
    flux1[2,img]=total(signal1[0:yy2])
    flux1[3,img]=total(signal1[0:yy3])
    
    flux2[1,img]=total(signal2[0:yy1])
    flux2[2,img]=total(signal2[0:yy2])
    flux2[3,img]=total(signal2[0:yy3])
    
    ; for plotting.....  
    ;for gem=1, 3 do begin
    ;
    ;loc_2d=array_indices(cutout, order1[0:max_count[gem,img]])
      
    ;xpix=loc_2d[0,*]
    ;ypix=loc_2d[1,*]
    
    ;plot_image, bytscl(cutout1, 50, 300)
    ;oplot, xpix, ypix, psym=8, color=cgcolor('purple')
      
    ;stop
     
    ;endfor
            
  endfor  ; end loop over images
    
  ; reject points with really low or really high counts
  mean_count=robust_mean(nsubpix[0,*],2)
  rej=where(reform(nsubpix[0,*]) le mean_count/2. OR reform(nsubpix[0,*]) ge mean_count*2., nrej, complement=keep1)
    
  ;plot, jd1, max_count[0,*], color=cgcolor('black'), psym=8
  ;oplot, jd1[rej], max_count[0,rej], color=cgcolor('red'), psym=8
  ;stop
  nsubpix=nsubpix[*,keep1]
    
  ; fit the curve - MODIFY THIS WHEN OBSERVATIONS ARE SHORTER/LONGER THAN 1S
  width=n_img_per_orbit
  sm1=smooth(reform(nsubpix[0,*]), width, /edge_truncate)
  norm=reform(nsubpix[0,*])-sm1
  scat_norm=robust_sigma(norm)
    
  xx=where(norm ge 5*scat_norm OR norm le (-5)*scat_norm OR nsubpix[0,*] lt 500, nxx, complement=keep2)
    
  ;plot, jd1[keep1], max_count[0,*], color=cgcolor('black'), psym=8
  ;if nxx gt 0 then  oplot, jd1[keep1[xx]], max_count[0,xx], color=cgcolor('red'), psym=8
  ; stop

  ; modify results - rejecting bad points...
  nsubpix=nsubpix[*,keep2]
  max_snr=max_snr[*,keep1[keep2]]
  max_dn=max_dn[keep1[keep2]]
  xy_psf=xy_psf[*,keep1[keep2]]
  flux1=flux1[*,keep1[keep2]]
  flux2=flux2[*,keep1[keep2]]
  jd1=jd1[keep1[keep2]]
  ccd_temp=ccd_temp[keep1[keep2]]
  roi_name=roi_name[keep1[keep2]]
  exp_num=exp_num[keep1[keep2]]
  ra_dec1=ra_dec1[*,keep1[keep2]]
  p2_trend=p2_trend[keep1[keep2]]
  nsat=nsat[keep1[keep2]]
    
  ; save results
  fileout=outdir+fname[bb]+'_p4_gem.sav'
  save, filename=fileout, vmag, nsubpix, max_snr, max_dn, xy_psf, flux1, flux2, jd1, nsat, $
  roi_name, exp_num, ra_dec1, jd, data1, roi_dim, xc, yc, ccd_temp, simbad_radec, $
  simbad_mags, parlax, otype, sptype, p2_trend, med_trend
    
endfor
  
  
print, 'end of program'
end


