pro total_energy

  ; 15 March 2014
  ; Calculate the total energy coming from each target - calculate the number of pixels enclosing 95% of the total energy
  ; get the flux, snr and scatter of the light curve using a set number of pixels in the PSF
  ; 
  ;
  ;
  Compile_opt idl2
  !p.background=cgcolor('white')
  plotsym, 0, 0.8, /fill
  devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
  imagelib    ; adds system variable !IMAGE needed for plotting
  
  indir='~/BRITE/data/UB/p3/'
  filesin=file_search(indir+'Orion-CF1-2*.sav', count=nfiles)
  fname=file_basename(filesin, '_p3.sav')
  
  outdir='/Users/gemmawhittaker/BRITE/data/UB/total_energy/'
    
  for bb=0, nfiles-1 do begin
  
    print, fname[bb]
    
    restore, filesin[bb] ;roi_name, exp_num, ra_dec1, jd, data1, roi_dim, xc, yc, ccd_temp, simbad_radec, $
    ;simbad_mags, parlax, otype, sptype, exp_time, exp_ttl, p2_trend, med_trend, sig_trend    
    
    jd1=jd-jd[0]
    keep=where(jd1 ge 0.0)
    
    jd1=jd1[keep]
    jd=jd[0]
    data1=data1[*,*,keep]
    roi_name=roi_name[keep]
    exp_num=exp_num[keep]
    ra_dec1=ra_dec1[*,keep]
    ccd_temp=ccd_temp[keep]
    p2_trend=p2_trend[keep]
    
    jd=jd+2450000D  ; print, jd[0], format='(d20.7)'
    
    dat=data1 ; no effect to data1 OR times 10. for 0.1s exp
    
    jd2=jd1[1:n_elements(jd1)-1]
    jdiff=jd2-jd
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
    max_snr=fltarr(nfrm)        ; SNR at nsubpix
    max_dn=fltarr(nfrm)         ; maximum data number on image
    xy_psf=fltarr(2,nfrm)       ; x and y coords of PSF center
    flux=fltarr(4,nfrm)         ; flux array - using signal2
    nsat=intarr(nfrm)           ; number of saturated pixels - left out of flux determination
    
    for img=0, nfrm-1 do begin
    
      dat1=dat[*,*,img]
      
      ; reform dat1 to a 1d vector
      dat2=reform(dat1, 32*32)
      
      ; order in order of decreasing flux
      order1=reverse(sort(dat2))
      
      loc2d=array_indices(dat1, order1)
      xl=reform(loc2d[0,*])
      yl=reform(loc2d[1,*])
      
      npix=n_elements(dat2)
      
      signal2=dat2[order1]
      
      thresh=20;50.
      
      count=0
      sig_cum=[]
      
      repeat begin
               
        sig_cum=[sig_cum,total(signal2[0:count])]
        
        count=count+1
        
      endrep until signal2[count] le thresh
      
      count=count-1
      
      ; calculate 90% of the total flux
      sig90=sig_cum[count]*0.9  ; this is total target flux
      
      xx=(where(sig_cum ge sig90))[0]
      
      plot_image, bytscl(dat1, 20, 500)
      oplot, xl[0:xx], yl[0:xx], psym=2, color=cgcolor('orange')
      
      sig95=0.95*sig90  ; this is 95% of target flux
      
      xx=(where(sig_cum ge sig95))[0]
      
      oplot, xl[0:xx], yl[0:xx], psym=2, color=cgcolor('green')
       
      stop
   
    endfor  ; end loop over images
    
  
  endfor
  
  
  print, 'end of program'
  stop
end


