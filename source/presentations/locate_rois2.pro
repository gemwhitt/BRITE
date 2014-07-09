pro locate_rois2

  ; Program to locate the ROI's on the whole CCD, to determine central, lower left, top left, lower right, top right....
  ;
  Compile_opt idl2
  !p.background=cgcolor('white')
  devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
  imagelib    ; adds system variable !IMAGE needed for plotting
  
  targets=['HD36861', 'HD31237', 'HD34085', 'HD36486', 'HD39801']
  
  indir='~/BRITE/data/UB/roi_proc1_sav/'
  
  savfiles=file_search(indir+'Orion-CF1-2_'+targets+'_p1.sav', count=nsav)
  
;  stop
  
  for i=0, nsav-1 do begin
    
    print, i
  
    restore, savfiles[i]  ;roi_name, exp_num, ra_dec1, helio_jd, data1, roi_dim, xc, yc
    
    xpos=xc
    ypos=yc
    
    x1=roi_dim[0]
    y1=roi_dim[2]
    
    roi=roi_name[0]
     
  outdir='/Users/gemmawhittaker/BRITE/mini_meeting/plots/'
  fileout=outdir+roi+'_on_ccd.ps'
  ps_on, fileout, xsize=17, ysize=15
  plot, [0,4008,4008,0,0], [0,0,2672,2672,0], color=cgcolor('black'), thick=3, xrange=[0,4008], yrange=[0,2672], $
    xstyle=1, ystyle=1
  plotsym, 0, /fill, 2.
  ;oplot, xpos, ypos, psym=8, color=cgcolor('purple')
  oplot, [0,4008], [2672./2.,2672./2.], linestyle=2, color=cgcolor('blue'), thick=2
  oplot, [4008./2.,4008./2.], [0,2672], linestyle=2, color=cgcolor('blue'), thick=2
  tvbox, 32, xpos, ypos, color=cgcolor('purple'), fill=1
  xyouts, x1, ypos+25, roi, charsize=0.7, color=cgcolor('black'), charthick=4
  ps_off
  
;  stop
  
  endfor
  
 
  print, 'End of program'
  
  
end