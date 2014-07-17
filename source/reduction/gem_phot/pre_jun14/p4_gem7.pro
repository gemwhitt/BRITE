pro p4_gem7

; modified from p3_p4_gem6 on 28th April 2014 - for CENTAURUS TEST DATA
; 
; use "info.dat" file which is already made by get_psf_info.pro and sr_fitting.pro
; 
; 
;
Compile_opt idl2
!p.background=cgcolor('white')
plotsym, 0, 0.8, /fill
devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
imagelib    ; adds system variable !IMAGE needed for plotting

nstk=1  ; number of stacked frames
nap=4 ; number of apertures to test
version='1'
  
indir='~/BRITE/data/UB/p2/CENTAURUS/'+strtrim(nstk,2)+'stk/'

filesin=file_search(indir+'*120307*.sav', count=nfiles)

outdir='~/BRITE/data/UB/p4/CENTAURUS/'+strtrim(nstk,2)+'stk/'
  
plotdir='~/BRITE/data/UB/p4/CENTAURUS/'+strtrim(nstk,2)+'stk/plots/'

psfdir='~/BRITE/data/UB/p4/CENTAURUS/'+strtrim(nstk,2)+'stk/psf/'
  
for bb=0, nfiles-1 do begin
  
  plotpsf=0
  
  fname=file_basename(filesin[bb], '_p2_'+strtrim(nstk,2)+'.sav')
    
  print, fname
    
  restore, filesin[bb] ;roi_name, exp_num, ra_dec, jd, data1, roi_dim, xc, yc, ccd_temp, exp_time, exp_ttl, $
                       ;medcols, medimg, ndead, rightcol, simbad_radec, vmag, bmag, parlax, otype, sptype 
    
  hdname=roi_name[0]
    
  jd1=jd-jd[0]
  keep=where(jd1 ge 0.0)
    
  jd1=jd1[keep]
  jd=jd[keep]
  data1=data1[*,*,keep]
  ccd_temp=ccd_temp[keep]
  medimg=medimg[keep]
  
  dat=data1 ; no effect to data1 OR times 10. for 0.1s exp
    
  ;get total number of frames in this file
  nfrm=(size(dat, /dim))[2]
    
  jd2=jd1[1:n_elements(jd1)-1]
  jdiff=jd2-jd1
  ; get number of orbits = ngap+1
  gap=where(jdiff gt 0.015, ngap)
  ; calculate cadence of observations
  cadence=robust_mean(jdiff,2)  ; in days
  cadence_m=cadence*24.*60.       ; in minutes
  cadence_s=cadence*24.*3600.       ; in seconds
  n_img_per_orbit=fix(15./cadence)
    
  ; get x and y pixel dimensions
  xdim=(size(dat, /dim))[0]
  ydim=(size(dat, /dim))[1]
    
  ; set up arrays for collecting results
  nsubpix=lonarr(nap,nfrm)      ; number of subpixels in PSF
  max_snr=fltarr(nfrm)        ; SNR at nsubpix
  flux=dblarr(nap,nfrm)         ; flux array - using signal2
  nsat=intarr(nfrm)           ; number of saturated pixels - left out of flux determination
  max_dn=lonarr(nfrm)
  xy_psf=fltarr(2,nfrm)
    
  for img=0, nfrm-1 do begin
    
    dat1=dat[*,*,img]
    
    ; get xy_psf from info file
    readcol, psfdir+fname+'_psf.txt', skipline=8+img, numline=1, maxdn, xpsf, ypsf, format='l,f,f', /silent
   
    ; get number of pixels in PSF from info file
    readcol, psfdir+fname+'_psf.txt', skipline=1, numline=1, npix0, npix1, npix2, format='x,i,i,i', /silent
    npix1=2*npix0
    npix2=3*npix0
    
    ; get fwhm from info file
    readcol, psfdir+fname+'_psf.txt', numline=1, fwhmx, fwhmy, format='x,f,f', /silent

    pks=[xpsf, ypsf, maxdn]
    
    xy_psf[0,img]=xpsf
    xy_psf[1,img]=ypsf
      
    pks1=floor(pks) ; round of pks to nearest whole pixel
      
    ; determine the background - away from the PSF center
    rbin=1
    get_bkgd, dat1,pks1,rbin, bkgd,bk_err
    
    dat1=dat1-bkgd
    
      
    ; make cutout - centered on PSF center
    if pks1[0]-(11) lt 0 then x1=0 else x1=pks1[0]-(11)
    if pks1[0]+(11) gt xdim-1 then x2=xdim-1 else x2=pks1[0]+(11)
    if pks1[1]-(11) lt 0 then y1=0 else y1=pks1[1]-(11)
    if pks1[1]+(11) gt ydim-1 then y2=ydim-1 else y2=pks1[1]+(11)
    ;x1=0
    ;x2=xdim-1
    ;y1=0
    ;y2=ydim-1
    
    cutout0=dat1[x1:x2,y1:y2]
      
    ; rebin the cutout - 4x4, 8x8, 16x16...
    rbin=10. ;8.
    np0=npix0*(rbin^2)
    np1=npix1*(rbin^2)
    np2=npix2*(rbin^2)
    
    cutout1=rebin_data(cutout0,rbin)
          
    cutout2=cutout1
      
    for i=0, (size(cutout0, /dim))[0]-1 do begin      ; column
      for j=0, (size(cutout0, /dim))[1]-1 do begin    ; row
        
        cutout2[i*rbin:((i+1)*rbin)-1,j*rbin:((j+1)*rbin)-1]=cutout0[i,j] ; this array is equivalent in dimension to cutout
          
      endfor
    endfor
      
    ; deal with saturated/nonlinear pixels
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
    
    if np2 gt ncut then np2=np1+10.
      
    ;reorder pixels by data number
    order1=reverse(sort(cutout1))
    temp_arr=order1
    signal1=(cutout1[order1])
    signal2=cutout2[order1]
      
    count=1 ; start on count=2, not count=0
    snr1=0
    sig_cum=[signal1[0],total(signal1[0:1])]
    snr_cum=[sig_cum[0]/sqrt(signal1[0]), sig_cum[1]/sqrt((signal1[1])+medimg[img])]
      
      repeat begin
      
        count=count+1
        
        if count gt np2*2 then goto, skipover
        
        snr=snr1
        
        sig=total(signal1[0:count])
        sig_cum=[sig_cum,sig]
        
        snr1=sig/sqrt(signal1[count] + (medimg[img]*count))  ; original
        snr_cum=[snr_cum,snr1]
        
      endrep until snr1 lt snr  ;ORIGINAL
      
      xx=count-1 ; location of maximum SNR - i.e. number of subpixels  ;ORIGINAL
      
      if xx lt np0 then goto, skipover
      
      max_snr[img]=snr_cum[xx]
      
      if nap gt 1 then nsubpix[*,img]=[xx,np0,np1,np2] else nsubpix[*,img]=[xx]
      
      for gem=0, nap-1 do flux[gem,img]=(total(signal1[0:nsubpix[gem,img]]))/(rbin^2)
      
      max_dn[img]=signal1[0]
      
      if plotpsf eq 0 then begin
      
        plotout=plotdir+fname+'_aper.ps'
        if nap gt 1 then !p.multi=[0,2,2,0,0] else !p.multi=0
        ps_on, plotout, xsize=18, ysize=18
        
        loc2d=array_indices(cutout1, order1[0:xx])
        xpix=loc2d[0,*]
        ypix=loc2d[1,*]
        
        plot_image, bytscl(cutout1, 0, 500), color=cgcolor('black'), $
          title=hdname+',V='+strtrim(vmag,2)+', NPIX='+strtrim(float(nsubpix[0,img])/(rbin^2),2), charsize=0.5
        oplot, xpix, ypix, psym=8, color=cgcolor('purple')
             
        if nap gt 1 then begin
        
          for gem=1, 3 do begin
          
            loc2d=array_indices(cutout1, order1[0:nsubpix[gem,img]])
            
            xpix=loc2d[0,*]
            ypix=loc2d[1,*]
            
            plot_image, bytscl(cutout1, 0, 500), color=cgcolor('black'), $
              title=hdname+',V='+strtrim(vmag,2)+', NPIX='+strtrim(float(nsubpix[gem,img])/(rbin^2),2), charsize=0.5
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
      ccd_temp=ccd_temp[nonzero]
      medimg=medimg[nonzero]
      nsat=nsat[nonzero]
      jd1=jd1[nonzero]
      
    endif
    
    ; reject points with really low or really high: nsubpix, max_dn, max_snr....
    mean_count=median(nsubpix[0,*])
    mean_dn=median(max_dn)
    mean_snr=median(max_snr)
    
    rej=where(reform(nsubpix[0,*]) lt mean_count*0.5 OR reform(nsubpix[0,*]) gt mean_count*2. OR $
      max_dn lt mean_dn*0.5 OR max_dn gt mean_dn*2. OR $
      max_snr lt mean_snr*0.5 OR max_snr gt mean_snr*2., nrej, complement=keep)
      
    !p.multi=[0,1,2,0,0]
    !p.background=cgcolor('white')
    plotsym, 0, /fill, 0.3
    plotout=plotdir+fname+'_clips.ps'
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
    ccd_temp=ccd_temp[keep]
    medimg=medimg[keep]
    nsat=nsat[keep]
    jd1=jd1[keep]

    ; save results
    txtout=outdir+'lc_txt/'+fname+'_p4gem.txt'
    savout=outdir+'lc_sav/'+fname+'_p4gem.sav'
    
    ; save .sav file
        save, filename=savout, roi_name, exp_num, ra_dec, jd, data1, roi_dim, xc, yc, ccd_temp, exp_time, exp_ttl, $
      medcols, medimg, ndead, rightcol, simbad_radec, vmag, bmag, parlax, otype, sptype, $
      nsubpix, max_snr, max_dn, xy_psf, flux, nsat
    
    ; find best LC and save data
    sigflux=fltarr(nap)
    for i=0, nap-1 do begin
    
      n=n_elements(jd)
      jd1=jd-jd[0]
      jd2=jd1[1:n-1]
      jdiff=jd2-jd1
      gap=where(jdiff gt 0.015, ngap)
      gap=[0,gap,n-1]
      
      sigt=fltarr(ngap+1)
      
      for j=0, ngap do begin
      
        if j eq 0 then iloc=indgen(gap[j+1]+1) else iloc=indgen(gap[j+1]-gap[j])+gap[j]+1
        
        fl=reform(flux[i,iloc])
        sigt[j]=robust_sigma(fl/robust_mean(fl,2))
        
      endfor
      
      sigflux[i]=median(sigt)*100.
      
    endfor
    
    xx=(where(sigflux eq min(sigflux)))[0]
    
    flux=reform(flux[xx,*])
    mags=-2.5*alog10(flux)
    
    npix_in_ap=median(nsubpix[xx,*])/(rbin^2)
    err_npix=round(robust_sigma(nsubpix[xx,*]/(rbin^2)))
    
    plot, jd1, flux, color=cgcolor('black'), psym=8, title='After clips', xtitle='Time (days)', $
      ytitle='Flux (DN)', charsize=0.8, xrange=[min(jd1)-0.1, max(jd1)+0.1], xstyle=1, $
      ystyle=1
      
    ps_off
    
    nsat=robust_mean(nsat,2)
    
    exp_time=float(1)
    
    bminv=bmag-vmag
    fwhm=[fwhmx,fwhmy]
  
    ; save .txt file
    openw, lun, txtout, /get_lun
    printf, lun, hdname, format='(a)'
    printf, lun, 'UB', format='(a)'
    printf, lun, fname, format='(a)'
    printf, lun, 'Processed by GNW', format='(a)'
    printf, lun, 'Version '+version, format='(a)'
    printf, lun, 'Notes:', format='(a)'
    printf, lun, 'Vmag='+strtrim(vmag,2), format='(a)'
    printf, lun, 'BminV='+strtrim(bminv,2), format='(a)'
    printf, lun, 'Num saturated pixels'+strtrim(nsat,2), format='(a)'
    printf, lun, 'Best aperture: '+strtrim(xx,2), format='(a)'
    printf, lun, 'Number of pixels in aperture: '+strtrim(npix_in_ap,2)+' +/- '+strtrim(err_npix,2)
    printf, lun, 'Average FWHM: '+strtrim(fwhm[0],2)+', '+strtrim(fwhm[1],2)+' pixels'
    printf, lun, 'Average sigma across each orbit is: '+strtrim(min(sigflux),2)+'%', format='(a)'
    printf, lun, 'Data below: JD, Magnitudes, DN, ExpTime, CCDtemp, x_psf, y_psf, max_dn, medimg', format='(a)'
    for i=0, n-1 do printf, lun, jd[i], mags[i], flux[i], exp_time, ccd_temp[i], xy_psf[0,i], xy_psf[1,i], max_dn[i], medimg[i], $
      format='(d20.8, x, f, x, f, x, f, x, f7.3, x, f, x, f, x, f, x, f)'
    free_lun, lun
    
    
  endfor
  
  
  print, 'end of program'
  stop
end


