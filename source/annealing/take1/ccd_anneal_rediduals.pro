pro ccd_anneal_rediduals

; program to analyse the residuals in the columns and pixels after annealing

Compile_opt idl2

plotout='y'

; FOR PLOTTING
devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
imagelib    ; adds system variable !IMAGE needed for plotting
!p.background=cgcolor('white')
;
; input directory
indir='/Users/gemmawhittaker/BRITE/data/anneal_test_data/sav_files/'
prefile=file_search(indir+'pre.sav')
postfile=file_search(indir+'post.sav')

; save results to....
outdir='~/BRITE/results/annealing/date_8nov13/'

readtimes=[0.06,1.,10.]
ntime=n_elements(readtimes)

toplot1='y'

if toplot1 eq 'y' then begin
fileout=outdir+'plots/residuals.ps'
ps_on, fileout, xsize=27, ysize=9, /landscape
!p.multi=[0,3,1,0,0]

for i=0, ntime-1 do begin
  
  restore, prefile  ;fdate, temp, data, readtime  ;DATA=LONG_Array[4048, 2672, 3]
  
  pre_dat=data[*,*,i]
  pre_temp=temp[i]
  
  ; identify 3 sub-groups of pixels
  med_ccd=median(pre_dat)
  sig_ccd=robust_sigma(pre_dat)
  
  grp1=(where(abs(pre_dat-med_ccd) gt sig_ccd-(sig_ccd*0.1) AND abs(pre_dat-med_ccd) lt sig_ccd+(sig_ccd*0.1), n1))[0:1999]
  thr1=med_ccd+sig_ccd
  
  grp2=(where(abs(pre_dat-med_ccd) gt 10*sig_ccd-(10*sig_ccd*0.1) AND abs(pre_dat-med_ccd) lt 10*sig_ccd+(10*sig_ccd*0.1), n2))[0:1999]
  thr2=med_ccd+sig_ccd*10.
  
  grp3=(where(abs(pre_dat-med_ccd) gt 20*sig_ccd-(20*sig_ccd*0.1) AND abs(pre_dat-med_ccd) lt 20*sig_ccd+(20*sig_ccd*0.1), n3))[0:1999]
  thr3=med_ccd+sig_ccd*20.
 ; stop
  ; now restore post-data
  restore, postfile  ;fdate, temp, data, readtime  ;DATA=LONG_Array[4048, 2672, 3]
  
  post_dat=data[*,*,i]
  post_temp=temp[i]
 
  res1=pre_dat[grp1]-post_dat[grp1]
  res2=pre_dat[grp2]-post_dat[grp2]
  res3=pre_dat[grp3]-post_dat[grp3]
  
  med1=median(res1)
  med2=median(res2)
  med3=median(res3)
  
  x1=where(res1 gt 0, nx1)
  x2=where(res2 gt 0, nx2)
  x3=where(res3 gt 0, nx3)
  
  ; calculate improvements
  imp1=median(res1[x1])/thr1*100.
  imp2=median(res2[x2])/thr2*100.
  imp3=median(res3[x3])/thr3*100.
  
  pci1=strmid(strtrim(float(nx1)/float(n_elements(res1))*100.,2),0,3)
  pci2=strmid(strtrim(float(nx2)/float(n_elements(res2))*100.,2),0,3)
  pci3=strmid(strtrim(float(nx3)/float(n_elements(res3))*100.,2),0,3)
 
  
  openw, lun, '~/Desktop/anneal_stats.txt', /get_lun, /append
  printf, lun, strtrim(readtimes[i],2)+'seconds', format='(a)'
  printf, lun, 'group 1, average reduction is '+strtrim(median(res1[x1])/thr1*100.,2), format='(a)'
  printf, lun, 'group 2, average reduction is '+strtrim(median(res2[x2])/thr2*100.,2), format='(a)'
  printf, lun, 'group 3, average reduction is '+strtrim(median(res3[x3])/thr3*100.,2), format='(a)'
  
  x1=where(res1 lt 0, nx1)
  x2=where(res2 lt 0, nx2)
  x3=where(res3 lt 0, nx3)
  
  ; calculate damages
  dm1=(-1)*median(res1[x1])/thr1*100.
  dm2=(-1)*median(res2[x2])/thr2*100.
  dm3=(-1)*median(res3[x3])/thr3*100.
  
  pcd1=strmid(strtrim(float(nx1)/float(n_elements(res1))*100.,2),0,3)
  pcd2=strmid(strtrim(float(nx2)/float(n_elements(res2))*100.,2),0,3)
  pcd3=strmid(strtrim(float(nx3)/float(n_elements(res3))*100.,2),0,3)
  
  printf, lun, strtrim(readtimes[i],2)+'seconds', format='(a)'
  printf, lun, 'group 1, average increase is '+strtrim(median(res1[x1])/thr1*100.,2), format='(a)'
  printf, lun, 'group 2, average increase is '+strtrim(median(res2[x2])/thr2*100.,2), format='(a)'
  printf, lun, 'group 3, average increase is '+strtrim(median(res3[x3])/thr3*100.,2), format='(a)'
  
  free_lun, lun

 ; window, 0, xsize=700, ysize=700, xpos=1500, ypos=200
  cghistoplot, res1, backcolorname='white', axiscolorname='black', datacolorname='blue', polycolor='blue', /line_fill, $
    orientation=45, xtitle='Residuals. Pre-anneal DN='+strmid(strtrim(thr1,2),0,4), ytitle='Number of Pixels', charsize=1., binsize=10, $
    title='Group 1: "Not hot", Readtime = '+strmid(strtrim(readtimes[i],2),0,4)+' s', $
    xrange=[(-1)*thr1,thr1], yrange=[0,1000], xstyle=1
    oplot, [0,0], [0, 2000], thick=3, color=cgcolor('purple')
    xyouts, 10, 900, 'In '+pci1+'% of pixels', charsize=0.65, color=cgcolor('purple')
    xyouts, ((-1)*thr1+10), 900, 'In '+pcd1+'% of pixels', charsize=0.65, color=cgcolor('purple')
    xyouts, 10, 850, 'DC reduces by '+strmid(strtrim(imp1,2),0,3)+'%', charsize=0.65, color=cgcolor('purple')
    xyouts, ((-1)*thr1+10), 850, 'DC increases by '+strmid(strtrim(dm1,2),0,3)+'%', charsize=0.65, color=cgcolor('purple')
    
    cghistoplot, res2, backcolorname='white', axiscolorname='black', datacolorname='blue', polycolor='blue', /line_fill, $
    orientation=45, xtitle='Residuals. Pre-anneal DN='+strmid(strtrim(thr2,2),0,4), ytitle='Number of Pixels', charsize=1., binsize=10, $
    title='Group 2: "Hot", Readtime = '+strmid(strtrim(readtimes[i],2),0,4)+' s', $
    xrange=[(-1)*thr2,thr2], yrange=[0,1000], xstyle=1
     oplot, [0,0], [0, 2000], thick=3, color=cgcolor('purple')
     xyouts,  10, 900, 'In '+pci2+'% of pixels', charsize=0.65, color=cgcolor('purple')
     xyouts, ((-1)*thr2+10), 900, 'In '+pcd2+'% of pixels', charsize=0.65, color=cgcolor('purple')
     xyouts,  10, 850, 'DC reduces by '+strmid(strtrim(imp2,2),0,3)+'%', charsize=0.65, color=cgcolor('purple')
     xyouts, ((-1)*thr2+10), 850, 'DC increases by '+strmid(strtrim(dm2,2),0,3)+'%', charsize=0.65, color=cgcolor('purple')
    
    cghistoplot, res3, backcolorname='white', axiscolorname='black', datacolorname='blue', polycolor='blue', /line_fill, $
    orientation=45, xtitle='Residuals. Pre-anneal DN='+strmid(strtrim(thr3,2),0,4), ytitle='Number of Pixels', charsize=1., binsize=10, $
    title='Group 3: "Very hot", Readtime = '+strmid(strtrim(readtimes[i],2),0,4)+' s', $
    xrange=[(-1)*thr3,thr3], yrange=[0,1000], xstyle=1
     oplot, [0,0], [0, 2000], thick=3, color=cgcolor('purple')
     xyouts, 10, 900, 'In '+pci3+'% of pixels', charsize=0.65, color=cgcolor('purple')
     xyouts, 10, 850, 'DC reduces by '+strmid(strtrim(imp3,2),0,3)+'%', charsize=0.65, color=cgcolor('purple')
     xyouts, ((-1)*thr3+10), 900, 'In '+pcd3+'% of pixels', charsize=0.65, color=cgcolor('purple')
     xyouts, ((-1)*thr3+10), 850, 'DC increases by '+strmid(strtrim(dm3,2),0,3)+'%', charsize=0.65, color=cgcolor('purple')
  
  ;stop
endfor

ps_off

endif

; now make plot with all data
fileout2=outdir+'plots/residuals_all.ps'
!p.multi=[0,3,1,0,0]
ps_on, fileout2, xsize=27, ysize=9, /landscape

for i=0, ntime-1 do begin
  
  restore, prefile
  pre_dat=data[*,*,i]
  
  restore, postfile
  post_dat=data[*,*,i]
  
  res0=pre_dat-post_dat
  nx=n_elements(res0)
  
  med_res=median(res0)
  
  x1=where(res0 lt 0, nx1)
  val1=float(nx1)/float(nx)*100.
  
  x2=where(res0 gt 0, nx2)
  val2=float(nx2)/float(nx)*100.
  
  
  
  cghistoplot, res0, backcolorname='white', axiscolorname='black', datacolorname='blue', polycolor='blue', /line_fill, $
    orientation=45, xtitle='Residuals', ytitle='Number of Pixels', charsize=0.8, binsize=10, $
    title='Readtime = '+strmid(strtrim(readtimes[i],2),0,4)+' s', $
    xrange=[-200,200], yrange=[0,4000000]
    oplot, [med_res,med_res], [0,4000000], thick=2, color=cgcolor('purple')
    xyouts, 80, 3000000, 'median at '+strmid(strtrim(med_res,2),0,3), charsize=0.6, color=cgcolor('blue')
    xyouts, -180, 3000000, strmid(strtrim(val1,2),0,3)+'% worse', charsize=0.6, color=cgcolor('red')
    xyouts, -180, 2500000, strmid(strtrim(val2,2),0,3)+'% better', charsize=0.6, color=cgcolor('red')
    
    
endfor

ps_off
spawn, 'open '+fileout+' &' 


spawn, 'open '+fileout2+' &' 

stop
print, 'End of Program'
end