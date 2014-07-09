pro ccd_anneal_hps

; Investigate very hot, medium hot and just hot pixels
; 
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
 
  post_fits=file_search(indir2+'*_'+times+'.fits')
  post_data=mrdfits(post_fits, 0, header2)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  ; identify Hot, med hot and just hot pixels in orignal data
  
  ; record mean and sigma of data
  mean_pre=robust_mean(pre_data,2)
  sig_pre=robust_sigma(pre_data)
  
  thr1=mean_pre+5*sig_pre
  thr2=mean_pre+20*sig_pre
  thr3=mean_pre+80*sig_pre
  
  hp0=(where(abs(pre_data-mean_pre) lt 1, nhp0))[0:7200]
  hp1=where(abs(pre_data-thr1) lt 1, nhp1)
  hp2=where(abs(pre_data-thr2) lt 16, nhp2)
  hp3=where(abs(pre_data-thr3) lt 72, nhp3)
  
  ; compute residuals for these 3 groups
  res1=pre_data[hp1]-post_data[hp1]
  res2=pre_data[hp2]-post_data[hp2]
  res3=pre_data[hp3]-post_data[hp3]
  res0=pre_data[hp0]-post_data[hp0]
  
  if plotout eq 'n' then window, 0, xsize=600, ysize=500, xpos=1500, ypos=150 else $
    ps_on, outdir+'hps_thrsh_hist.ps', xsize=15, ysize=15
  cghistoplot, res1, datacolorname='blue', thick=1, backcolorname='white', axiscolorname='black', maxinput=200, mininput=-200, $
    /fill, polycolor='sky blue', orientation=45, histdata=result, $
    title='Three types of hot pixel, '+read_temp[i]+'deg', xtitle='', ytitle='', binsize=10, charsize=0.7
  ;cghistoplot, res1, datacolorname='orange', /oplot, binsize=10, thick=2
    cghistoplot, res2, datacolorname='green', maxinput=200, mininput=-200, /oplot, binsize=10, /line_fill, $
      orientation=45, polycolor='green', thick=2
    
    cghistoplot, res3, datacolorname='purple', maxinput=200, mininput=-200, /oplot, binsize=10, /line_fill, orientation=135, $
      polycolor='purple', thick=2
      
     al_legend, ['just hot (180)','med hot (500)','very hot (1800)'], color=[cgcolor('sky blue'), cgcolor('green'), cgcolor('purple')], psym=[2,2,2], textcolor=cgcolor('black'), /right, charsize=0.7
    if plotout eq 'y' then ps_off
      if plotout eq 'y' then ps_off
      
  ;oplot, [robust_mean(resid,2), robust_mean(resid,2)], [0, max(result)+100], color=cgcolor('purple'), thick=2
  ;al_legend, 'mean residuals is '+strtrim(robust_mean(resid,2),2), textcolor=cgcolor('black'), /right
  
  stop
endfor

end