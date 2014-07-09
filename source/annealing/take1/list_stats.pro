pro list_stats

  ; program to analyse the residuals in the columns and pixels after annealing
  
  Compile_opt idl2
  
  plotout='y'
  
  ; FOR PLOTTING
  devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
  imagelib    ; adds system variable !IMAGE needed for plotting
  !p.background=cgcolor('white')
  ;
  ; input directory
  indir1='~/BRITE/data/anneal_test_data/2013Sept24_ThermalTest_Triumf_Proton_test_CCD/'
  indir2='~/BRITE/data/anneal_test_data/2013Oct11_AnnealingTest_Triumf_Proton_test_CCD/'
  
  ; save results to....
  outdir='~/BRITE/results/annealing_151013/'
  
  ;times=['60ms1', '60ms2', '1s1', '1s2', '10s']
  times=['1s1']
  ntimes=n_elements(times)
  
  ;read_temp=['0','10','20','30']
  read_temp=['20','30']
  ntemp=n_elements(read_temp)
  
  for i=0, ntemp-1 do begin
  
    pre_fits=file_search(indir1+read_temp[i]+'_'+times+'.fits')
    pre_data=mrdfits(pre_fits, 0, header1)
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; get temp data
    pre_temps=indir1+'Temperatures_'+read_temp[i]+'.txt'
    readcol, pre_temps, skipline=2, numline=5, read_time, t1, t2, t3, t4, format='a,f,f,f,f'
    
    ; modify read_time array for comparison to 'times'
    dot_pos=strpos(read_time, '.')
    for k=0, n_elements(dot_pos)-1 do read_time[k]=strmid(read_time[k], 0, dot_pos[k])
    
    ; get average of all temps at this read time
    loc1=where(read_time eq times[0], nloc1)
    
    if nloc1 ne 1 then stop
    
    ;pre_at= average pre temp
    all_temps=[t1[loc1],t2[loc1],t3[loc1],t4[loc1]]
    pre_at=strtrim(robust_mean(all_temps,2),2)
    pre_at=strmid(pre_at, 0, 4)
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    ;window, 0, xsize=600, ysize=500, xpos=1500, ypos=150
    ;plot_image, bytscl(pre_data, 50, 150), color=cgcolor('black'), title='Pre Anneal '+strtrim(read_temp[i],2)
    
    post_fits=file_search(indir2+'*_'+times+'.fits')
    post_data=mrdfits(post_fits, 0, header2)
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; get temp data
    pre_temps=indir2+'Temperatures_20.txt'
    readcol, pre_temps, skipline=2, numline=5, read_time, t1, t2, t3, t4, format='a,f,f,f,f'
    
    ; modify read_time array for comparison to 'times'
    dot_pos=strpos(read_time, '.')
    for k=0, n_elements(dot_pos)-1 do read_time[k]=strmid(read_time[k], 0, dot_pos[k])
    
    ; get average of all temps at this read time
    loc1=where(read_time eq times[0], nloc1)
    
    if nloc1 ne 1 then stop
    
    ;pre_at= average pre temp
    all_temps=[t1[loc1],t2[loc1],t3[loc1],t4[loc1]]
    post_at=strtrim(robust_mean(all_temps,2),2)
    post_at=strmid(post_at, 0, 4)
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    fileout=outdir+'stats_'+times+'.txt'
    openw, lun, fileout, /get_lun
    printf, lun, 'Statistic', 'Pre-Ann', 'Post_Ann', format='(a20,x,a20,x,a20)'
    printf, lun, 'Avg Temp', pre_at, post_at, format='(a20,x,a20,x,a20)'
    
    
    ; compute data medians and sigmas
    pre_med=strmid(strtrim(robust_mean(pre_data,2),2), 0, 5)
    post_med=strmid(strtrim(robust_mean(post_data,2),2), 0, 5)
    
    pre_sig=strmid(strtrim(robust_sigma(pre_data), 2), 0, 5)
    post_sig=strmid(strtrim(robust_sigma(post_data),2), 0, 5)
    
    printf, lun, 'Median', pre_med, post_med, format='(a20,x,a20,x,a20)'
    printf, lun, 'Scatter', pre_sig, post_sig, format='(a20,x,a20,x,a20)'
    
    thr1=pre_med+5*pre_sig
    thr2=pre_med+20*pre_sig
    thr3=pre_med+80*pre_sig
    
    totp=float((size(pre_data, /dim))[0]*(size(pre_data, /dim))[1])
    
    pre_norm=where(pre_data lt thr1, nprenorm)
    pcprenorm=strmid(strtrim(float(nprenorm)/totp*100, 2), 0, 5)
    
    post_norm=where(post_data lt thr1, npostnorm)
    pcpostnorm=strmid(strtrim(float(npostnorm)/totp*100, 2), 0, 5)
    
    pre_int=where(pre_data gt thr1 and pre_data lt thr3, npreint)
    pcpreint=strmid(strtrim(float(npreint)/totp*100, 2), 0, 5)
    
    post_int=where(post_data gt thr1 and post_data lt thr3, npostint)
    pcpostint=strmid(strtrim(float(npostint)/totp*100, 2), 0, 5)
    
    pre_hp=where(pre_data gt thr3, nprehp)
    pcprehp=strmid(strtrim(float(nprehp)/totp*100, 2), 0, 5)
    
    post_hp=where(post_data gt thr3, nposthp)
    pcposthp=strmid(strtrim(float(nposthp)/totp*100, 2), 0, 5)
    
    
    
    printf, lun, 'Percent lt '+strmid(strtrim(thr1,2), 0, 5), pcprenorm, pcpostnorm, format='(a20,x,a20,x,a20)'
    printf, lun, 'Percent: '+strmid(strtrim(thr1,2), 0, 5)+'-'+strmid(strtrim(thr3,2), 0, 5), pcpreint, pcpostint, format='(a20,x,a20,x,a20)'
    printf, lun, 'Percent gt '+strmid(strtrim(thr3,2), 0, 5), pcprehp, pcposthp, format='(a20,x,a20,x,a20)'
    
    free_lun, lun
    
    
    
    
stop

endfor
  
  print, 'End of Program'
end