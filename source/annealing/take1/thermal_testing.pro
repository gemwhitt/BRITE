pro thermal_testing

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
outdir='~/BRITE/results/ThermalTest_240913/'

temps=['0','10','20','30']
ntemp=n_elements(temps)

times=['10s','1s1', '1s2', '60ms1', '60ms2']
ntimes=n_elements(times)

for i=0, ntimes-1 do begin
  
  fitsfiles=file_search(indir+'*_'+times[i]+'.fits', count=nfits)
  
  ; define a 1x4 array for the columns/rows medians for each temperature
  med_cols=fltarr(4, 4048)
  med_rows=fltarr(4,2672)
  test_temp=fltarr(4)
  
  for j=0, nfits-1 do begin
    
    test_temp[j]=file_basename(fitsfiles[j], '_'+times[i]+'.fits')
    test_time=times[i]
    
    print, 'Temp='+strtrim(test_temp[j],2)+'  Time='+test_time
    
    ; get fits info
    fits_info, fitsfiles[j], n_ext=next, extname=extnames, /silent
    
    ; read header info
    data=mrdfits(fitsfiles[j], 0, header)  ; data is 2672 x 4048
   
    xx=where(data eq 0, nxx)
    
    print, nxx
    
    print, robust_mean(data[*,0:100],2)
    
    print, robust_mean(data[*,2500:2600],2)
    
    stop
    
    ; do IQR analysis
    rstat, data, med, hinge1, hinge2, ifence1, ifence2, ofence1, ofence2, mind, maxd, /noprint, /descrip
    
    ;window, 0, xsize=600, ysize=600, xpos=0, ypos=100
    ;plot_image, data, color=cgcolor('white')
    
    window, 1, xsize=600, ysize=600, xpos=800, ypos=100
    plot_image, bytscl(data, ifence1, ifence2), color=cgcolor('white')
    plot_image, bytscl(data, 5.5, 81.5), color=cgcolor('white')
    stop
    ; get median of columns and median of rows
    ncols=(size(data, /dim))[0]
    nrows=(size(data, /dim))[1]
    
    for jj=0, ncols-1 do med_cols[j,jj]=robust_mean(data[jj,*],2)
    for jj=0, nrows-1 do med_rows[j,jj]=robust_mean(data[*,jj],2)
    
  endfor
    
    file1=outdir+'cols_'+test_time+'.txt'
    file2=outdir+'rows_'+test_time+'.txt'
    
   ; openw, lun, file1, /get_lun
   ; printf, lun, test_temp
   ; for jj=0, ncols-1 do printf, lun, med_cols[*,jj]
   ; free_lun, lun
    
   ; openw, lun, file2, /get_lun
   ; printf, lun, test_temp
   ; for jj=0, nrows-1 do printf, lun, med_rows[*,jj]
   ; free_lun, lun
    
      
endfor

stop
end