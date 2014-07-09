pro plot_temp_trend

  ; input directory
  indir='~/BRITE/data/UB/testing/nostack/p2/'
  
  outdir='/Users/gemmawhittaker/BRITE/mini_meeting/plots/'
  
  filesin=file_search(indir+'*.sav', count=nsav)
  
  dark_rate=fltarr(nsav)
  mags=fltarr(nsav)
  
  for kk=0, nsav-1 do begin
  
    restore, filesin[kk]  ; roi_name, exp_num, ra_dec1, jd, data1, roi_dim, xc, yc, ccd_temp, simbad_radec, $
   ; simbad_mags, parlax, otype, sptype
   
   vpos=strpos(simbad_mags, 'V=')
   vmag=float(strmid(simbad_mags, vpos+2))
   mags[kk]=vmag
    
    
    hdname=file_basename(filesin[kk], '_p2_test.sav')
    print, hdname
    
    ; determine number of orbits - 1 orbit=15 mins (0.01 days), 1 gap = 85 mins (0.059 days)
    npts=n_elements(jd)
    
    jd1=jd[1:npts-1]
    jdiff=jd1-jd
    
    gap=where(jdiff gt 0.02, ngap)  ; there are ngap+1 observations
    gap=[gap,npts-1]
    
    ; define two floating-point arrays in which to record values for the change in temp and data number across the frame
    chg_temp=fltarr(ngap+1)
    chg_dn=fltarr(ngap+1)
    
    expn=intarr((ngap+1)*43)
    flx=fltarr((ngap+1)*43)
    temp1=fltarr((ngap+1)*43)
    
    ;fileout=outdir+'stack_orbit.ps'
    ;ps_on, fileout, xsize=15, ysize=11
    ;plotsym, 0, /fill, 0.5
     
    ; loop over each observation window
    for i=-1, ngap-1 do begin
      print, i+1
      if i eq -1 then begin
        obs_data=data1[*,*,0:gap[i+1]]
        obs_jd=jd[0:gap[i+1]]
        obs_temp=ccd_temp[0:gap[i+1]]
        obs_exp=exp_num[0:gap[i+1]]
      endif else begin
        obs_data=data1[*,*,gap[i]+1:gap[i+1]]
        obs_jd=jd[gap[i]+1:gap[i+1]]
        obs_temp=ccd_temp[gap[i]+1:gap[i+1]]
        obs_exp=exp_num[gap[i]+1:gap[i+1]]
      endelse
      
      obs_exp=obs_exp-obs_exp[0]
      expn[(i+1)*43]=obs_exp
      
      opts=n_elements(obs_jd)
      
      ; calculate median of all pixels in this frame
      med_frm=fltarr(opts)
      for j=0, opts-1 do med_frm[j]=median(obs_data[*,*,j])
      
      flx[(i+1)*43]=med_frm
      
      temp1[(i+1)*43]=obs_temp
      
      ; record the range in temperature and the range in median values
      chg_temp[i+1]=max(obs_temp)-min(obs_temp)
      
     ; if i eq -1 then plot, obs_exp, med_frm, psym=8, color=cgcolor('black'), $
     ;   xtitle='Exposure number in sequence' , ytitle='Image Median - stacked', $
     ;   charsize=0.8, yrange=[80,110] else $
     ;   oplot, obs_exp, med_frm, psym=8, color=cgcolor('black')
      
    endfor
    
    
    
    
    ; calculate average slope across each orbit
    nbin1=43
    width1=21
    orbit_trend, expn,flx,temp1,nbins1,width1,  avg_flx,sm_flx,avg_temp
    
    smtemp=smooth(avg_temp,21)
    
    ;oplot, expn, sm_flx, thick=5, color=cgcolor('orchid')
    ;ps_off
    ;stop
    
    ; calculate average dark current rate
    delta_temp=robust_mean(chg_temp,2)
    delta_dn=avg_flx[42]-avg_flx[0]
    
   ; dark_rate=delta_dn/delta_temp ; number of data counts per degree
   
   fileout0=outdir+'trend_removed.ps'
   flux2=(avg_flx-sm_flx)+robust_mean(sm_flx,2)
   ps_on, fileout0, xsize=15, ysize=11
   plotsym, 0, /fill, 0.5
   plot, expn, flux2, psym=8, color=cgcolor('black'), xtitle='Exposure Number', ytitle='Average Image Median - Corrected', $
     charsize=0.8, yrange=[80,110]
   ;oplot, expn, sm_flx, psym=0, thick=4, color=cgcolor('orchid')
   ps_off
      
    fileout1=outdir+'avg_flx_orbit.ps'
    ps_on, fileout1, xsize=15, ysize=11
    plotsym, 0, /fill, 0.5
    plot, expn, avg_flx, psym=8, color=cgcolor('black'), xtitle='Exposure Number', ytitle='Average Image Median', $
      charsize=0.8, yrange=[80,110]
    ;oplot, expn, sm_flx, psym=0, thick=4, color=cgcolor('orchid')
    ps_off
   stop 
    fileout1=outdir+'avg_temp_orbit.ps'
    ps_on, fileout1, xsize=15, ysize=11
    plotsym, 0, /fill, 0.5
    plot, expn, avg_temp, psym=8, color=cgcolor('black'), xtitle='Exposure Number', ytitle='Average CCD Temperature', $
      charsize=0.8, /ynozero
  ;  oplot, expn, smtemp, psym=0, thick=4, color=cgcolor('orchid')
    ps_off
    
    
    delta_flx=avg_flx[42]-avg_flx[0]
    delta_temp=avg_temp[42]-avg_temp[0]
    
    dark_rate[kk]=(delta_flx/delta_temp)*3.2  ; in flx per degree
    stop
    endfor
    
    print,  mags
    print, dark_rate
    
stop

end