pro xy_shifts

; plot xy psf shifts - compare gemphot to sr-phot
; 
; Compile_opt idl2

devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
imagelib    ; adds system variable !IMAGE needed for plotting

!p.background=cgcolor('white')
plotsym, 0, /fill, 0.8

ingem='~/BRITE/data/UB/p4_gemphot1/'
inslr='~/BRITE/data/UB/p4_sr1/'

fgem=file_search(ingem+'*.sav', count=nfiles)
fslv=file_search(inslr+'*.sav', count=nfiles1)

outdir='/Users/gemmawhittaker/BRITE/reports/myreports/'

fileout=outdir+'xy_shifts.ps'

;ps_on, fileout, xsize=18, ysize=25
;!p.multi=[0,2,3,0,0]

for i=0, nfiles-1 do begin
  
  restore, fgem[i]
  print, fgem[i]
  
  xgem=xy_psf[0,*]
  ygem=xy_psf[1,*]
  
  restore, fslv[i]  ; also has data1
  print, fslv[i]
  
  xslv=xy_psf[0,*]
  yslv=xy_psf[1,*]
  
  plotsym, 0, /fill, 0.6
  
  meanx=robust_mean(xgem,2)
  meany=robust_mean(ygem,2)
  
  window, 0, xsize=550, ysize=450, xpos=1500, ypos=200
  plot, xgem, ygem, psym=8, color=cgcolor('black'), xrange=[meanx-1.,meanx+1.], yrange=[meany-1.,meany+1.], $
    xtitle='X-pixel', ytitle='Y-pixel', charsize=0.7, title='PSF shifts - gem'

  window, 1, xsize=550, ysize=450, xpos=2300, ypos=200
  plot, xslv, yslv, psym=8, color=cgcolor('black'), xrange=[meanx-1.,meanx+1.], yrange=[meany-1.,meany+1.], $
    xtitle='X-pixel', ytitle='Y-pixel', charsize=0.7, title='PSF shifts - S.R.'
    
  window, 2, xsize=550, ysize=550, xpos=2000, ypos=-400
  plot_image, data1[*,*,0]
  oplot, [xgem[0]], [ygem[0]], psym=2, color=cgcolor('orange')
  oplot, [xslv[0]], [yslv[0]], psym=2, color=cgcolor('green')
  

stop
endfor

;ps_off

;spawn, 'open '+fileout+' &'

stop
print, 'End of program'
end