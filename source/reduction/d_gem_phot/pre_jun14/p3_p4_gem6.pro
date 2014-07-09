pro p3_p4_gem6

; modified from p3_p4_gem5 - 20 March 2014
; 
; Add the following routines:
; 1. get average flux - for each orbit
; 2. get average number of PSF pixels and minimum number of PSF pixels - i.e. which contain about 60-75% target flux
; 3. get accurate PSF centers using Slavek programs - fourier space convolution with idealised PSF
; 
version=strtrim('2',2) ; use interpolated signal 1 with no filter 

gain=3.25
qe=0.3

; number of apertures
nap=4
;
Compile_opt idl2
!p.background=cgcolor('white')
plotsym, 0, 0.8, /fill
devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
imagelib    ; adds system variable !IMAGE needed for plotting

nstk=1
  
indir='~/BRITE/data/UB/p2/CENTAURUS/'+strtrim(nstk,2)+'stk/'

filesin=file_search(indir+'*.sav', count=nfiles)

fname=file_basename(filesin, '_p2_'+strtrim(nstk,2)+'.sav')
  
outdir='/Users/gemmawhittaker/BRITE/data/UB/p4/CENTAURUS/'+strtrim(nstk,2)+'stk/'
  
plotdir='/Users/gemmawhittaker/BRITE/data/UB/P4/CENTAURUS/'+strtrim(nstk,2)+'stk/'+'plots/'
 
for b=0, nfiles-1 do begin
  
  plotpsf=0
    
  print, fname[b]
        
  restore, filesin[b] ;roi_name, exp_num, ra_dec, jd, data1, roi_dim, xc, yc, ccd_temp, exp_time, exp_ttl, $
                      ;medcols, medimg, ndead, rightcol, simbad_radec, vmag, bmag, parlax, otype, sptype
                       
  hdname=roi_name[0]               
                     
  jd1=jd-jd[0]
  keep=where(jd1 ge 0.0)
    
  bad=where(jd1 lt 0.0, nbad)
  if nbad gt 0 then stop
  
  dat=data1 ; no effect to data1 OR times 10. for 0.1s exp
  
  ;get total number of frames in this file
  nfrm=(size(dat, /dim))[2]
    
  jd2=jd1[1:nfrm-1]
  jdiff=jd2-jd1
  gap=where(jdiff gt 0.015, ngap)
  ; calculate cadence of observations
  cadence=robust_mean(jdiff,2)  ; in days
  cadence_m=cadence*24.*60.       ; in minutes
  cadence_s=cadence*24.*3600.       ; in seconds
  n_img_per_orbit=fix(15./cadence)
  
  xy_psf=fltarr(2,nfrm)
  minpix=fltarr(2)  ; value plux sigma
  minflux=dblarr(2) ; value plux sigma
  maxpix=fltarr(2)  ; value plux sigma
  maxflux=dblarr(2) ; value plux sigma
  fwhm=fltarr(2)
  
  ; get x and y pixel dimensions
  xdim=(size(dat, /dim))[0]
  ydim=(size(dat, /dim))[1]
  
  ; call external program - get PSF centers, average flux per orbit and number of pixels in PSF
  get_psf_info, dat,gap,ngap,nfrm, xy_psf,minpix,minflux,maxpix,maxflux,fwhm
       
  
stop    
  ; set up arrays for collecting results
  nsubpix=intarr(nap,nfrm)      ; number of subpixels in PSF
  max_snr=fltarr(nfrm)        ; SNR at nsubpix
  max_dn=fltarr(nfrm)         ; maximum data number on image
  flux=fltarr(nap,nfrm)         ; flux array - using signal2
  nsat=intarr(nfrm)           ; number of saturated pixels - left out of flux determination
    
  for img=0, nfrm-1 do begin
    
    dat1=dat[*,*,img]
       
    pks=[xy_psf[0,img], xy_psf[1,img], max(dat1)]
    
    if pks[0] lt 0 OR pks[0] gt xdim OR pks[1] lt 0 OR pks[1] gt ydim then CONTINUE
    
    pks1=round(pks) ; round of pks to nearest whole pixel
      
    ; determine the background - away from the PSF center
    rbin=1
    get_bkgd, dat1,pks1,rbin, bkgd,bk_err
          
    ; make cutout - centered on PSF center
    if pks1[0]-(10) lt 0 then x1=0 else x1=pks1[0]-(10)
    if pks1[0]+(10) gt xdim-1 then x2=xdim-1 else x2=pks1[0]+(10)
    if pks1[1]-(10) lt 0 then y1=0 else y1=pks1[1]-(10)
    if pks1[1]+(10) gt ydim-1 then y2=ydim-1 else y2=pks1[1]+(10)
    cutout0=dat1[x1:x2,y1:y2]
           
    ; rebin the cutout - 4x4, 8x8, 16x16...
    rbin=8. ;8.
    cutout1=rebin_data(cutout0,rbin)
    
    max_thr=(maxpix[0]*rbin^2)
    min_thr=(minpix[0]*rbin^2)-(minpix[1]*minpix[0]*rbin^2)
    
    if (size(cutout1, /dim))[0]*(size(cutout1, /dim))[1] lt max_thr then CONTINUE
       
    cutout2=cutout1

    for i=0, (size(cutout0, /dim))[0]-1 do begin      ; column
      for j=0, (size(cutout0, /dim))[1]-1 do begin    ; row
        
        cutout2[i*rbin:((i+1)*rbin)-1,j*rbin:((j+1)*rbin)-1]=cutout0[i,j] ; this array is equivalent in dimension to cutout
          
      endfor
    endfor
    
      ; deal with saturated pixels
      sat=where(cutout2 ge 9000, numsat, complement=lin)
      
     ; if numsat gt 0 then begin
        ;stop
        nsat[img]=numsat/(rbin^2)
        ;  print, nsat
     ;   cutout1[sat]=-9999
     ;   cutout2[sat]=-9999
     ; endif else nsat[img]=0
      
      ; divide pixel counts by e.g. 64 to add together later for flux values
      cutout2=cutout2/(rbin^2)
      
      ncut=n_elements(cutout1)
      
      ;reorder pixels by data number
      order1=reverse(sort(cutout1))
      temp_arr=order1
      signal1=(cutout1[order1])
      signal2=cutout2[order1]
      
      count=1 ; start on count=2, not count=0
      snr1=0
      sig_cum=[signal1[0]*gain,total(signal1[0:1]*gain)]
      snr_cum=[sig_cum[0]/sqrt(signal1[0]*gain), sig_cum[1]/sqrt((signal1[1]*gain)+(p2_trend[img]*gain))]    
      
      repeat begin
      
        count=count+1
        
        if count gt max_thr then goto, skipover
        
        snr=snr1
        
        sig=total(signal1[0:count]*gain)
        sig_cum=[sig_cum,sig]
        
        snr1=sig/sqrt(signal1[count]*gain + (p2_trend[img]*count*gain))  ; original
        snr_cum=[snr_cum,snr1]
        
      endrep until snr1 lt snr  ;ORIGINAL
            
     xx=count-1 ; location of maximum SNR - i.e. number of subpixels  ;ORIGINAL
     
     if xx lt min_thr then goto, skipover
            
     max_snr[img]=snr_cum[xx]
            
     ; get the ratio: snr/xx
     ratio1=max_snr[img]/float(xx)
     ; get 3 scale factors
     sc1=ratio1
     sc2=ratio1*2.
     sc3=ratio1*4.
    stop 
     ; get 3 more apertures
     x1=((1.+sc1)*xx) < max_thr
     x2=((1.+sc2)*xx) < max_thr
     x3=((1.+sc3)*xx) < max_thr
     ;x4=2632.  ; based on graph (snr/xx vs snr)  ; indpendent of Vmag
     
     if nap gt 1 then nsubpix[*,img]=[xx,x1,x2,x3] else nsubpix[*,img]=[xx]
     
     for gem=0, nap-1 do flux[gem,img]=(total(signal1[0:nsubpix[gem,img]]))/64.
     
     max_dn[img]=signal1[0]
     
     if plotpsf eq 0 then begin
      
       plotout=plotdir+fname[b]+'_psf.ps'
       if nap gt 1 then !p.multi=[0,2,2,0,0] else !p.multi=0
       ps_on, plotout, xsize=18, ysize=18
       
       loc2d=array_indices(cutout1, order1[0:xx])   
       xpix=loc2d[0,*]
       ypix=loc2d[1,*]
       
       plot_image, bytscl(cutout1, 0, 500), color=cgcolor('black'), title=hdname+',V='+strtrim(vmag,2)+', scale= 0-500', charsize=0.5
       oplot, xpix, ypix, psym=8, color=cgcolor('purple')
       
       ;if x4 gt xx then begin
       ; 
       ; loc2d=array_indices(cutout1, order1[xx+1:x4])
        
       ; xpix=loc2d[0,*]
       ; ypix=loc2d[1,*]
        
       ; oplot, xpix, ypix, psym=8, color=cgcolor('orange')
       ;endif
       
       if nap gt 1 then begin
        
       for gem=1, 3 do begin
       
         loc2d=array_indices(cutout1, order1[0:nsubpix[gem,img]])
         
         xpix=loc2d[0,*]
         ypix=loc2d[1,*]
         
         plot_image, bytscl(cutout1, 0, 500), color=cgcolor('black'), title=hdname+',V='+strtrim(vmag,2)+', scale= 0-500', charsize=0.5
         oplot, xpix, ypix, psym=8, color=cgcolor('purple')
         
       endfor 
       endif
             
       ps_off
      
       plotpsf=1
       !p.background=cgcolor('white')
     endif

      skipover:
    endfor  ; end loop over images
    
    ; remove 0s
    zero=where(max_dn eq 0, nzero, complement=nonzero)
    
    if nzero gt 0 then begin
      
      ; modify results - rejecting bad points...
      nsubpix=nsubpix[*,nonzero]
      max_snr=max_snr[nonzero]
      max_dn=max_dn[nonzero]
      xy_psf=xy_psf[*,nonzero]
      flux=flux[*,nonzero]
      jd=jd[nonzero]
      jd1=jd1[nonzero]
      ccd_temp=ccd_temp[nonzero]
      roi_name=roi_name[nonzero]
      exp_num=exp_num[nonzero]
      ra_dec1=ra_dec1[*,nonzero]
      p2_trend=p2_trend[nonzero]
      nsat=nsat[nonzero]     
      
    endif
     
    ; reject points with really low or really high: nsubpix, max_dn, max_snr....
    mean_count=robust_mean(nsubpix[0,*],3)
    mean_dn=robust_mean(max_dn,3)
    mean_snr=robust_mean(max_snr,3)
    
    rej=where(reform(nsubpix[0,*]) le mean_count*0.5 OR reform(nsubpix[0,*]) ge mean_count*1.5 OR $
      max_dn le mean_dn*0.5 OR max_dn ge mean_dn*1.5 OR $
      max_snr lt mean_snr*0.5 OR max_snr ge mean_snr*1.5, nrej, complement=keep)
    
    
    !p.multi=[0,1,2,0,0]
    !p.background=cgcolor('white')
    plotsym, 0, /fill, 0.3
    plotout=plotdir+fname[b]+'_clips.ps'
    ps_on, plotout, xsize=16, ysize=24
    
    ;window, 0, xsize=700, ysize=400
    plot, jd1, flux[0,*], color=cgcolor('black'), psym=8, title=hdname+', Vmag='+strtrim(vmag,2), xtitle='Time (days)', $
      ytitle='Flux (DN)', charsize=0.8, xrange=[min(jd1)-0.1, max(jd1)+0.1], xstyle=1
    if nrej gt 0 then  oplot, jd1[rej], flux[0,rej], color=cgcolor('red'), psym=8

    ; modify results - rejecting bad points...
    nsubpix=nsubpix[*,keep]
    max_snr=max_snr[keep]
    max_dn=max_dn[keep]
    xy_psf=xy_psf[*,keep]
    flux=flux[*,keep]
    jd=jd[keep]
    jd1=jd1[keep]
    ccd_temp=ccd_temp[keep]
    roi_name=roi_name[keep]
    exp_num=exp_num[keep]
    ra_dec1=ra_dec1[*,keep]
    p2_trend=p2_trend[keep]
    nsat=nsat[keep]
    
    
    
    ; save results
    txtout=outdir+'txt/'+fname[b]+'_p4gem'+version+'.txt'
    savout=outdir+'sav/'+fname[b]+'_p4gem'+version+'.sav'
    
    ; save .sav file
    save, filename=savout, vmag, nsubpix, max_snr, max_dn, xy_psf, flux, jd, nsat, $
    roi_name, exp_num, ra_dec1, roi_dim, xc, yc, ccd_temp, simbad_radec, $
    simbad_mags, parlax, otype, sptype, p2_trend, med_trend
    
    ; find best LC and save data
    sigflux=fltarr(nap)
    for i=0, nap-1 do begin
      
      n=n_elements(jd)
      jd2=jd1[1:n_elements(jd1)-1]
      jdiff=jd2-jd1
      gap=where(jdiff gt 0.015, ngap)
      gap=[0,gap,n-1]
      
      sigt=fltarr(ngap+1)
      
      for j=0, ngap do begin
        
        if j eq 0 then iloc=indgen(gap[j+1]+1) else iloc=indgen(gap[j+1]-gap[j])+gap[j]+1
        
        fl=flux[i,iloc]
        sigt[j]=robust_sigma(fl/robust_mean(fl,2))
      
      endfor
        
      sigflux[i]=(robust_mean(sigt,2))*100.
        
    endfor
    
    xx=(where(sigflux eq min(sigflux)))[0]
    
    flux=reform(flux[xx,*])
    mags=-2.5*alog10(flux)
    
    npix_in_ap=robust_mean(nsubpix[xx,*],2)/(rbin^2)
    err_npix=round(robust_sigma(nsubpix[xx,*]/(rbin^2)))
    
      
    plot, jd1, flux, color=cgcolor('black'), psym=8, title='After clips', xtitle='Time (days)', $
      ytitle='Flux (DN)', charsize=0.8, xrange=[min(jd1)-0.1, max(jd1)+0.1], xstyle=1, $
      ystyle=1
      
    ps_off

    
    ; save .txt file
    openw, lun, txtout, /get_lun
    printf, lun, hdname, format='(a)'
    printf, lun, 'UB', format='(a)'
    printf, lun, fname[b], format='(a)'
    printf, lun, 'Processed by GNW', format='(a)'
    printf, lun, 'Version '+version, format='(a)'
    printf, lun, 'Notes:', format='(a)'
    printf, lun, 'Vmag='+strtrim(vmag,2), ', Otype: '+strtrim(otype,2), ', SPtype: '+strtrim(sptype,2), format='(a,x,a,x,a)'
    printf, lun, 'Best aperture: '+strtrim(xx,2), format='(a)'
    printf, lun, 'Number of pixels in aperture: '+strtrim(npix_in_ap,2)+' +/- '+strtrim(err_npix,2)
    printf, lun, 'Average FWHM: '+strtrim(fwhm[0],2)+', '+strtrim(fwhm[1],2)+' pixels'
    printf, lun, 'Average sigma across each orbit is: '+strtrim(min(sigflux),2)+'%', format='(a)'
    printf, lun, 'Data below: JD, Magnitudes, DN, ExpTime, CCDtemp, x_psf, y_psf', format='(a)'
    for i=0, n-1 do printf, lun, jd[i], mags[i], flux[i], exp_time, ccd_temp[i], xy_psf[0,i], xy_psf[1,i], $
      format='(d20.8, x, f, x, f, x, i, x, f7.3, x, f, x, f)'
    free_lun, lun
    

  endfor
  
  
  print, 'end of program'
  stop
end


