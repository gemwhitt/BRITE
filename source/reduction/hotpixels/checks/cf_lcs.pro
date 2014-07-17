pro cf_lcs

  ; progam to compare light curves from Gemma and Slavek photometry....
  ;
  Compile_opt idl2
  
  !p.background=cgcolor('white')
  plotsym, 0, /fill, 0.8
  
  ingem='~/BRITE/data/UB/p4_gemphot1/'
  inslr='~/BRITE/data/UB/p4_sr1/'
  
  fgem=file_search(ingem+'*.sav', count=nfiles)
  fslv=file_search(inslr+'*.sav', count=nfiles1)
  
  cols1=[cgcolor('purple'), cgcolor('sea green')]
  cols2=['purple', 'sea green']
  
  outdir='/Users/gemmawhittaker/BRITE/reports/myreports/'
  
  mags=fltarr(nfiles)
  
  mf1=fltarr(nfiles)
  mf2=fltarr(nfiles)
  scat1=fltarr(nfiles)
  scat2=fltarr(nfiles)
  
  fileout=outdir+'cf_lcs.ps'
  ps_on, fileout, xsize=18, ysize=25
  !p.multi=[0,2,3,0,0]
  
  for i=0, nfiles-1 do begin
    
    hdname=file_basename(fgem[i],'_p4_gem.sav')
  
    restore, fgem[i]  ;flux, jd, bkgd, exp_num, roi_dim, ccd_temp, xy_psf, sig_res, counts, vmag
    t1=jd
    
    fscat=fltarr(4)
    for j=0, 3 do begin
      meanf=robust_mean(flux[j,*],2)
      normf=flux[j,*]/meanf
      fscat[j]=robust_sigma(normf)
    endfor
    
    xx=(where(fscat eq min(fscat)))[0]
    f1=flux[xx,*]
    
    n1=n_elements(t1)
    
    
    restore, fslv[i]  ;flux, jd, bkgd, exp_num, roi_dim, ccd_temp, xy_psf, sig_res, counts, vmag
    f2=flux
    t2=jd
    n2=n_elements(t2)
    
    plotsym, 0, /fill, 0.4
;    window, 0, xsize=550, ysize=400, xpos=1500, ypos=200
    plot, t1, f1/robust_mean(f1,2), color=cgcolor('black'), xtitle='Time (days)', $
      ytitle='Normalized DN', charsize=0.8, yrange=[0.8,1.2], psym=8, title='Gem-phot: '+hdname
      
    
;    window, 1, xsize=550, ysize=400, xpos=2300, ypos=200
    plot, t2, f2/robust_mean(f2,2), color=cgcolor('black'), xtitle='Time (days)', $
      ytitle='Normalized DN', charsize=0.8, yrange=[0.8,1.2], psym=8, title='S.R.-phot: '+hdname
       

    
;stop    
  endfor
  ps_off
  
  spawn, 'open '+fileout+' &'
  

  stop
  print, 'End of Program'
end