pro plot_tracking_1507

Compile_opt idl2

!p.background=cgcolor('white')

indir='/Users/gemmawhittaker/BRITE/results/UB/pointing'

obsdate='20130705'

filesin=file_search(indir+'/datfiles/*'+obsdate+'.dat', count=nfiles)

; get date for output file
fnames=file_basename(filesin, '.dat')
strlen1=strlen(fnames)
dates=strarr(nfiles)
for i=0, nfiles-1 do dates[i]=strmid(fnames[i], strlen1[i]-8)

fileout=indir+'/plots/'+dates[0]+'xy_diff.ps'
ps_on, fileout, xsize=17, ysize=13

for i=0, nfiles-1 do begin
  
  readcol, filesin[i], sat, exp_num, jd_obs, ra, dec, xc, yc, xcentroid, ycentroid, x_diff, y_diff, xy_diff, $
    format='(a,i,d,d,d,f,f,f,f,f,f,f)'
    
    nelem=n_elements(sat)
    
    time1=float(jd_obs-jd_obs[0])
    ; convert time to minutes
    time=time1*24.*60.
    
    ; calculate average of xy_diff and plot this
    mean_absdiff=strtrim(robust_mean(xy_diff,2), 2)
    dotpos=strpos(mean_absdiff, '.')
    mean_xydiff=strmid(mean_absdiff, 0, dotpos+4)
    
    allcolors=[cgcolor('red'), cgcolor('blue'), cgcolor('green')]
    
    if i eq 0 then begin
      plot, exp_num, xy_diff, color=cgcolor('black'), /nodata, yrange=[0,max(xy_diff)+4], ystyle=1
      legend, ['ext=0', 'ext=1', 'ext=2'], PSym=[2,2,2], Color=allcolors, $
        Position=[min(exp_num)+0.5,max(xy_diff)+3], /box, textcolors=[cgcolor('black'), cgcolor('black'), cgcolor('black')], $
        outline_color=cgcolor('black'), charsize=0.75
      oplot, exp_num, xy_diff, psym=-2, color=allcolors[0]
      ;oplot, exp_num, x_diff, psym=-2, color=cgcolor('blue')
      ;oplot, exp_num, y_diff, psym=-2, color=cgcolor('green')
      oplot, exp_num, replicate(mean_absdiff,nelem), linestyle=1, color=allcolors[0]
      ;oplot, [exp_num[0],exp_num[max(exp_num)]+1], [0,0], color=cgcolor('black'), linestyle=2, thick=2
      xyouts, 0.3, mean_xydiff, mean_xydiff, color=cgcolor('black'), charsize=0.5
    endif else begin
      oplot, exp_num, xy_diff, psym=-2, color=allcolors[i]
      ;oplot, exp_num, x_diff, psym=-2, color=cgcolor('blue')
      ;oplot, exp_num, y_diff, psym=-2, color=cgcolor('green')
      oplot, exp_num, replicate(mean_absdiff,nelem), linestyle=1, color=allcolors[i]
      ;oplot, [exp_num[0],exp_num[max(exp_num)]+1], [0,0], color=cgcolor('black'), linestyle=2, thick=2
      xyouts, 0.3, mean_xydiff, mean_xydiff, color=cgcolor('black'), charsize=0.5
    endelse
endfor

ps_off
stop

print, 'End of Program'
end