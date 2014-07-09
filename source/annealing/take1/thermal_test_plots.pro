pro thermal_test_plots

; first program to plot results from thermal_testing.pro
; 
Compile_opt idl2
;
!p.background=cgcolor('white')
; 
indir='~/BRITE/results/ThermalTest_240913/rows_cols/'

; for output
outdir='~/BRITE/results/ThermalTest_240913/plots/'
toplot='y'

opt='rows'  ; cols or rows

filesin=file_search(indir+opt+'*.txt', count=nfiles)

fileout=outdir+'histoplot_'+opt+'.ps'

if toplot eq 'n' then window, 1, xsize=600, ysize=600, xpos=2300, ypos=60 else $
  ps_on, fileout, xsize=17, ysize=17

for i=0, nfiles-1 do begin
  
  readcol, filesin[i], med_t1, med_t2, med_t3, med_t4
  
  print, file_basename(filesin[i],'.txt')
  
  test_temp=[med_t1[0],med_t2[0],med_t3[0],med_t4[0]]
  test_temp=strtrim(fix(test_temp),2)
  
  cols1=['blue','green','orange','purple']
  cols2=['cornflower blue','light sea green','pale goldenrod','thistle']
  
  if opt eq 'cols' then begin
    med_cols=[[med_t1],[med_t2],[med_t3],[med_t4]]
    
    ; find the maximum histogram y-value to plot from all 4 temps
    maxresult=fltarr(4)
    for j=0, 3 do begin
      result=histogram(med_cols[*,j], binsize=2)
      maxresult[j]=max(result)
    endfor
    
    ;
    xx=reform(med_cols, 4*4049)
    xy=xx[sort(xx)]
    maj_data=round(0.99*4*4049)
    maxdatavalue=xy[maj_data]
    
    
    cghistoplot, med_cols[*,0], datacolorname=cols1[0], axiscolorname='black', backcolorname='white', $
      maxinput=maxdatavalue, max_value=max(maxresult)+20, binsize=2, /fill, polycolor=cols2[0], $
      ytitle='Number of columns', xtitle='Robust mean of data column', title=file_basename(filesin[i],'.txt'), charsize=0.8
    for j=1,3 do cghistoplot, med_cols[*,j], datacolorname=cols1[j], /oplot, maxinput=200, binsize=2, $
      /fill, polycolor=cols2[j]
      for j=0, 3 do oplot, [robust_mean(med_cols[*,j],2),robust_mean(med_cols[*,j],2)], [0,max(maxresult)+20], $
      linestyle=0, color=cgcolor(cols1[j]), thick=2
    al_legend, test_temp+cgsymbol('deg'), colors=cols1, psym=[replicate(6,4)], textcolors=cgcolor('black'), /right, box=1, $
      outline_color=cgcolor('black'), charsize=0.7
     

  endif else begin
    med_rows=[[med_t1],[med_t2],[med_t3],[med_t4]]
    
    ; find the maximum histogram y-value to plot from all 4 temps
    maxresult=fltarr(4)
    for j=0, 3 do begin
      result=histogram(med_rows[*,j], binsize=2)
      maxresult[j]=max(result)
    endfor
    
    ;
    xx=reform(med_rows, 4*2673)
    xy=xx[sort(xx)]
    maj_data=round(0.99*4*2673)
    maxdatavalue=xy[maj_data]
    
   
    cghistoplot, med_rows[*,0], datacolorname=cols1[0], axiscolorname='black', backcolorname='white', $
      maxinput=maxdatavalue, max_value=max(maxresult)+20, binsize=2, /fill, polycolor=cols2[0], $
      ytitle='Number of rows', xtitle='Robust mean of data row', title=file_basename(filesin[i],'.txt'), charsize=0.8
    for j=1,3 do cghistoplot, med_rows[*,j], datacolorname=cols1[j], /oplot, maxinput=200, binsize=2, $
      /fill, polycolor=cols2[j]
    for j=0, 3 do oplot, [robust_mean(med_rows[*,j],2),robust_mean(med_rows[*,j],2)], [0,max(maxresult)+20], $
      linestyle=0, color=cgcolor(cols1[j]), thick=2
    al_legend, test_temp+cgsymbol('deg'), colors=cols1, psym=[replicate(6,4)], textcolors=cgcolor('black'), /right, box=1, $
      outline_color=cgcolor('black'), charsize=0.7
  endelse
  
endfor
 if toplot eq 'y' then ps_off
stop
end