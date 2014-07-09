pro remove_orbit_trend

; run a test on a segment of observations to show the effect of removing the mean trend in light level
; test section is from 4 to 6 days into the CF1-2 observations

Compile_opt idl2
!p.background=cgcolor('white')
  
; input directory
indir='~/BRITE/data/UB/testing/nostack/p2/'

outdir='~/BRITE/data/UB/testing/nostack/p3/'
    
filesin=file_search(indir+'*.sav', count=nsav)

plot0='n'
plot1='n'
plot2='n'
plot3='n'
plot4='n'

for kk=0, nsav-1 do begin
    
  restore, filesin[kk]  ; roi_name, exp_num, ra_dec1, jd, data1, roi_dim, xc, yc, ccd_temp, simbad_radec, $
  ;simbad_mags, parlax, otype, sptype

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
  chg_sm2=fltarr(ngap+1)
  
  expn=intarr((ngap+1)*43)
  flx=fltarr((ngap+1)*43)
  temp1=fltarr((ngap+1)*43)
  
   
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
 
    
    ; compute short and long smooth functions over the orbit
    sm1=smooth(med_frm, 3, /edge_truncate)
    sm2=smooth(sm1, opts, /edge_truncate)
    
    ; plot median of each frame versus time with a smoothed function overplotted
    if plot0 eq 'y' then begin
    if i eq -1 then begin 
      window, 0, xsize=700, ysize=500, xpos=100, ypos=200
      plotsym, 0, 1., /fill    
    plot, obs_jd-obs_jd[0], med_frm, psym=8, color=cgcolor('black'), xtitle='Time (days)', ytitle='Median Image', charsize=0.7, $
    title=file_basename(filesin[kk], '.sav'), yrange=[90,120] 
    endif else begin
      oplot, obs_jd-obs_jd[0], med_frm, psym=8, color=cgcolor('black')
      endelse
    ;oplot, obs_jd, sm2, color=cgcolor('purple'), thick=2
    endif 
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    ; plot ccd_temp versus time ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    if plot1 eq 'y' then begin
    window, 1, xsize=700, ysize=500, xpos=500, ypos=300
    plotsym, 0, 1., /fill
    plot, obs_jd, obs_temp, psym=8, color=cgcolor('black'), xtitle='Time (days)', ytitle='CCD Temp', charsize=0.7, $
    title=file_basename(filesin[kk], '.sav'), yrange=[18,28]
    endif
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     
    ; plot flux with increasing trend removed ;;;;;;;;;;;;;;;;;;;;;;;;;;;
    if plot2 eq 'y' then begin
    flux2=(med_frm-sm2)+robust_mean(med_frm,2)  
    window, 2, xsize=700, ysize=500, xpos=1000, ypos=400
    plotsym, 0, 1., /fill
    plot, obs_jd, flux2, color=cgcolor('black'), psym=8, yrange=[80,120]
    endif
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    ; plot auto-correlation versus time-lag  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    if plot3 eq 'y' then begin
    lag=indgen(opts)
    result=a_correlate(sm1, lag);, /covariance, /double)
    window, 3, xsize=700, ysize=500, xpos=150, ypos=400
    plotsym, 0, 1., /fill
    plot, lag, result, color=cgcolor('black'), psym=8, yrange=[-1,1]
    oplot, [0,50], [0,0], color=cgcolor('green'), thick=2
    endif
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ;stop   
  endfor
; stop 
  ; calculate average slope across each orbit
  nbin1=43
  width1=21
  orbit_trend, expn,flx,temp1,nbins1,width1,  avg_flx,sm_flx,avg_temp
  
  ; calculate average dark current rate 
  delta_temp=robust_mean(chg_temp,2)
  delta_dn=avg_flx[42]-avg_flx[0]
  
  dark_rate=delta_dn/delta_temp ; number of data counts per degree
  
 ; window, 0
 ; plotsym, 0, /fill, 0.4
 ; plot, expn, avg_flx, psym=8, color=cgcolor('black'), yrange=[90,120]
 ; oplot, expn, sm_flx, psym=8, color=cgcolor('orange')
 ;stop 
  ; plot results for each orbit.....
  if plot4 eq 'y' then begin
  
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
    
 
     ;plot sm_flx over each individual orbit - show the residual
    wset, 0
  plot, obs_exp, med_frm, color=cgcolor('black'), psym=8, title='Exposure Number', charsize=0.7, $
    ytitle='Avg Count', /ynozero
    oplot, obs_exp, sm_flx, color=cgcolor('orange'), thick=2
    
    res1=(med_frm-sm_flx)+robust_mean(med_frm,2)
    window, 1
    plot, obs_exp, res1, color=cgcolor('black'), psym=8, title='Exposure Number', charsize=0.7, $
      ytitle='Trend removed', /ynozero
 stop
      endfor
      
      endif
      
      ; save p3 files with trend data included in new file - for info
      for i=-1, ngap-1 do begin
        print, i+1
        if i eq -1 then begin
          obs_data=data1[*,*,0:gap[i+1]]
          for j=0, 42 do obs_data[*,*,j]=(obs_data[*,*,j]-sm_flx[j])+robust_mean(sm_flx,2)
          data1[*,*,0:gap[i+1]]=obs_data
        endif else begin
          obs_data=data1[*,*,gap[i]+1:gap[i+1]]
          for j=0, 42 do obs_data[*,*,j]=(obs_data[*,*,j]-sm_flx[j])+robust_mean(sm_flx,2)
          data1[*,*,gap[i]+1:gap[i+1]]=obs_data
        endelse

        endfor
        
        fileout=outdir+hdname+'_p3_test.sav'
        save, filename=fileout, roi_name, exp_num, ra_dec1, jd, data1, roi_dim, xc, yc, ccd_temp, sm_flx, $
          avg_flx, expn, simbad_radec, $
          simbad_mags, parlax, otype, sptype




 ; stop
  
endfor


stop
end

