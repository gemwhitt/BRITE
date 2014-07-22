pro plot_centroid_results

  Compile_opt idl2
  
  !p.background=cgcolor('white')
  
  indir='/Users/gemmawhittaker/BRITE/results/UB/pointing'
  
  obsdate='20130717'
  
  output='ps' ; 'x' or 'ps' depending on result
  
  filesin=file_search(indir+'/datfiles/*'+obsdate+'.dat', count=nfiles)
  
  for i=0, nfiles-1 do begin
    
    fname=file_basename(filesin[i], '.dat')
     
    readcol, filesin[i], sat, exp_num, jd_obs, ra, dec, xc, yc, xcentroid, ycentroid, x_diff, y_diff, xy_diff, $
      format='(a,i,d,d,d,f,f,f,f,f,f,f)'
      
    nelem1=n_elements(sat)
    
    keep=where(xy_diff lt xc, nelem2)
    
    xy_diff2=xy_diff[keep]
    x_diff2=x_diff[keep]
    y_diff2=y_diff[keep]
    exp_num2=exp_num[keep]
    
    ; calculate average of xy_diff and plot this
    mean_absdiff=strtrim(mean(xy_diff2), 2)
    dotpos=strpos(mean_absdiff, '.')
    mean_xydiff=strmid(mean_absdiff, 0, dotpos+4)
    
    maxdiff=max([x_diff2,y_diff2])
    mindiff=min([x_diff2,y_diff2])
    
    cols=[cgcolor('blue'),cgcolor('green')]
    syms=[-2, -4]
        
    fileout=indir+'/plots/'+fname+'.ps'
    !p.multi=[0,1,2,0,0]
    
    if output eq 'ps' then ps_on, fileout, xsize=18, ysize=18 else window, 0, xsize=600, ysize=600
    
    plot, exp_num2, xy_diff2, color=cgcolor('black'), yrange=[min(xy_diff2)-2,max(xy_diff2)+2], ystyle=1, xrange=[-1, 60], xstyle=1, $
      title=fname+' - Difference between peak and image center', xtitle='Exposure Number', ytitle='Absolute differece (pixels)', charsize=0.7
    oplot, exp_num2, xy_diff2, psym=syms[0], symsize=0.5, color=cgcolor('black'), linestyle=3
    oplot, exp_num, replicate(mean_absdiff, nelem1), linestyle=1, color=cgcolor('blue')
    xyouts, 0.5, max(xy_diff2)+0.5, 'Mean difference is '+mean_xydiff, charsize=0.6, color=cgcolor('blue')
    
    plot, exp_num2, xy_diff2, color=cgcolor('black'), yrange=[min(xy_diff2)-2,max(xy_diff2)+2], ystyle=1, xrange=[59, nelem1+1], xstyle=1, $
      title=fname, xtitle='Exposure Number', ytitle='Absolute differece (pixels)', charsize=0.7
    oplot, exp_num2, xy_diff2, psym=syms[0], symsize=0.5, color=cgcolor('black'), linestyle=3
    oplot, exp_num, replicate(mean_absdiff, nelem1), linestyle=1, color=cgcolor('blue')
    xyouts, 0.5, max(xy_diff2)+0.5, 'Mean difference is '+mean_xydiff+' pixels', charsize=0.75, color=cgcolor('blue')
    
    
    ; plot x and y difference separately
    ;plot, exp_num2, x_diff2, color=cgcolor('black'), yrange=[mindiff-3,maxdiff+3], ystyle=1, xrange=[-1, nelem1+1], xstyle=1, $
    ;  title=fname, xtitle='Exposure Number', ytitle='Abs diff from center PSF to center raster', charsize=0.7
    
    ;al_legend, ['X diff','Y diff'], colors=cols, psym=syms, /box
    ;oplot, exp_num2, x_diff2, psym=syms[0], symsize=0.5, color=cols[0]
    ;oplot, exp_num2, y_diff2, psym=syms[1], symsize=0.5, color=cols[1]
    ;oplot, exp_num, replicate(mean_absdiff, nelem1), linestyle=1, color=cols[0]
    ;xyouts, 0.5, max(xy_diff)+0.5, 'Mean difference is '+mean_xydiff, charsize=0.6, color=cgcolor('black')
    
    if output eq 'ps' then ps_off else stop
    
  ;  if i eq 4 then stop
  endfor
  
  print, 'End of Program'
end