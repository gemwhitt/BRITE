pro thermal_testing3

  Compile_opt idl2
  
  ; FOR PLOTTING
  devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
  imagelib    ; adds system variable !IMAGE needed for plotting
  ;!p.background=cgcolor('white')
  
  ; 2nd program to analyse thermal test results - look at individual hot pixels -
  ;                 - what affect does the temp have on these at certain integration times
  ;
  ;
  ; input directory
  indir='~/BRITE/data/Radiation_Test_Data/2013Sept24_ThermalTest_Triumf_Proton_test_CCD/'
  
  ; save results to....
  outdir='~/BRITE/results/ThermalTest_240913/hot_pixels/'
  
  times=['10s','1s1'] ;, '1s2']  ;, '60ms1', '60ms2']
  ntimes=n_elements(times)
 ; times=['1s1']
  
  sat_temp='20'
  
  for i=0, ntimes-1 do begin
  
    ; first get the locations of the 10 hottest pixels at room temp and a threshold for evaluating a hot pixel
    fitsfile=file_search(indir+sat_temp+'_'+times[i]+'.fits')
    
    data=mrdfits(fitsfile, 0, header)
    
    ; find hottest pixels
    ; what are the hottest pixels and where - need to rebin the data and then get multi-dimensional subscripts
    reordered_data=reverse(sort(data))
    hot_loc=reordered_data[0:9]
    
    ; get thresh
    mean_at_20=robust_mean(data,2)
    sig_at_20=robust_sigma(data)
    ;thresh_at_20=mean_at_20+(5*sig_at_20)
    thresh_at_20=180
    
    ; get 6 random hot pixels at 20 degrees - just hot, moderately hot and very hot
    hot1=200
    hot2=2000
    hot3=10000
    
    if i eq 0 then begin
    x1=(where(abs(data-hot1) lt 1, nx1))[0:1]
    x2=(where(abs(data-hot2) lt 10, nx2))[0:1]
    x3=(where(data gt hot3, nx3))[0:1]
    endif else begin
      x1=x1
      x2=x2
      x3=x3
    endelse
      
    ; now get stats for all files
    fitsfiles=file_search(indir+'*_'+times[i]+'.fits', count=nfits)
    
    ; define arrays for robust mean and sigma for the CCD
    mean_data=fltarr(4)
    sigma_data=fltarr(4)
    test_temp=fltarr(4)
    min_data=fltarr(4)
    max_data=fltarr(4)
    pchot_pix=fltarr(4)
    hot_int=fltarr(4,6)
    ;thresh=fltarr(4)
    
    for j=0, nfits-1 do begin
    
      test_temp[j]=file_basename(fitsfiles[j], '_'+times[i]+'.fits')
      test_time=times[i]
      
      print, 'Temp='+strtrim(test_temp[j],2)+'  Time='+test_time
      
      ; get fits info
      fits_info, fitsfiles[j], n_ext=next, extname=extnames, /silent
      
      ; read header info
      data=mrdfits(fitsfiles[j], 0, header)  ; data is 2672 x 4048
      
      ; do IQR analysis
      rstat, data, med, hinge1, hinge2, ifence1, ifence2, ofence1, ofence2, mind, maxd, /noprint
      
      ;window, 0, xsize=600, ysize=600, xpos=1500, ypos=60
      ;plot_image, data, color=cgcolor('white')
      
      ;window, 1, xsize=600, ysize=600, xpos=2300, ypos=60
      ;plot_image, bytscl(data, ofence1, ofence2), color=cgcolor('white')
      
      ; calculate robust mean and sigma of the data
      mean_data[j]=robust_mean(data,2)
      sigma_data[j]=robust_sigma(data)
      
      ; find hot pixels
      ; threshold is 5sigma above the mean
      ;thresh[j]=sigma_data[j]*5+mean_data[j]
      
      ; number of pixels above the threshold?
      hotpix=where(data gt thresh_at_20, nhot)
      
      total=(size(data, /dim))[0]*(size(data, /dim))[1]
      
      pcnt_hot=float(nhot)/float(total)*100.
      
      pchot_pix[j]=float(nhot)/float(total)*100.
      
      ; get min and max data values
      min_data[j]=min(data)
      max_data[j]=max(data)
      
      ; what are the hottest pixels and where - need to rebin the data and then get multi-dimensional subscripts
      reordered_data=reverse(sort(data))
      
      hot_loc=[x1,x2,x3]
      
      hot_int[j,*]=data[hot_loc]
      
      
      
    endfor
    
    fileout=outdir+'hotpix2_'+times[i]+'.txt'
    
    openw, lun, fileout, /get_lun
    printf, lun, 'Temp', 'Mean', 'Sigma', 'Thresh', 'Min',  'Max', 'PC_HP', $
      format='(a10,x,a10,x,a10,x,a10,x,a10,x,a10,x,a10)'
    for jj=0, 3 do printf, lun,  test_temp[jj], mean_data[jj], sigma_data[jj], thresh_at_20, min_data[jj], $
      max_data[jj], pchot_pix[jj], $
      format='(i10,x,d10.2,x,d10.3,x,d10.2,x,d10.2,x,d10.2,x,d10.2)'
    printf, lun, 'HP_Locations'
    for jj=0, 5 do printf, lun, hot_loc[jj]
    printf, lun, 'HP_Intensities'
    for jj=0, 3 do printf, lun, hot_int[jj,*], format='(d10.2,x,d10.2,x,d10.2,x,d10.2,x,d10.2,x,d10.2)'
    free_lun, lun
    
    print, 'see output file '+fileout
    
    ;stop
    
  endfor
  
  stop
end