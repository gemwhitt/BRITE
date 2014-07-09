pro therm_test_plot_histo
  
  Compile_opt idl2
  
  ; FOR PLOTTING
  devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
  imagelib    ; adds system variable !IMAGE needed for plotting
  ;!p.background=cgcolor('white')
  
  ; first program to analyse thermal test results
  ;
  ; input directory
  indir='~/BRITE/data/Radiation_Test_Data/2013Sept24_ThermalTest_Triumf_Proton_test_CCD/'
  
  ; save results to....
  outdir='~/BRITE/results/ThermalTest_240913/plots/'
  
  toplot='y'
  
  temps=['0','10','20','30']
  ntemp=n_elements(temps)
  
  times=['10s','1s1', '1s2', '60ms1', '60ms2']
  ntimes=n_elements(times)
  
  for i=0, ntimes-1 do begin
  
    fitsfiles=file_search(indir+'*_'+times[i]+'.fits', count=nfits)
    
    fileout2=outdir+'data_hist/data_hist_alltemp_'+times[i]+'.ps'

    
    if toplot eq 'n' then window, 0, xsize=600, ysize=600, xpos=800, ypos=100 else $
      ps_on, fileout2, xsize=17, ysize=17
      
    for j=0, nfits-1 do begin
    
      test_temp=file_basename(fitsfiles[j], '_'+times[i]+'.fits')
      test_time=times[i]
           
      print, 'Temp='+strtrim(test_temp,2)+'  Time='+test_time
      
      ; get fits info
      fits_info, fitsfiles[j], n_ext=next, extname=extnames, /silent
      
      ; read header info
      data=mrdfits(fitsfiles[j], 0, header)  ; data is 2672 x 4048
      
      total_pix=2672*4048
      fifty_pc=0.5*total_pix
      
      ; do IQR analysis
      rstat, data, med, hinge1, hinge2, ifence1, ifence2, ofence1, ofence2, mind, maxd, /noprint, /descrip
      
      if i eq 0 and j eq 0 then begin
        scale1=hinge1
        scale2=hinge2
      endif else begin
        scale1=scale1
        scale2=scale2
      endelse
      
      cghistoplot, data, binsize=5, maxinput=150, max_value=fifty_pc, backcolorname=cgcolor('white'), $
        datacolorname=cgcolor('blue'), /line_fill, orientation=45, xtitle='Data Count', ytitle='Number of Pixels', $
        title=test_time+', '+strtrim(test_temp,2)+' deg, median='+strtrim(fix(med),2), axiscolorname=cgcolor('black'), $
        charsize=0.7, polycolor=cgcolor('blue')
    endfor
    
    if toplot eq 'y' then ps_off
    
  endfor
  
  stop
end