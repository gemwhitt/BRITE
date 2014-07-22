pro locate_rois

; Program to locate the ROI's on the whole CCD, to determine central, lower left, top left, lower right, top right....
; 
Compile_opt idl2

sat='BA'

field='CENTAURUS'

indir='~/BRITE/'+sat+'/'+field+'/data/raw_sav/'

outdir='~/BRITE/'+sat+'/'+field+'/data/plots/'

fileout=outdir+sat+'-'+field+'_ccdroi.ps'

; each target is stored in its own directory - just get dirname + restore roi_loc from one file in the directory

tardirs=file_search(indir+'HD*/', count=ntar)

targets=file_basename(tardirs)

; define array variables to record xc and yc for each roi
xpos=lonarr(ntar)
ypos=lonarr(ntar)
roi=strarr(ntar)
x1=lonarr(ntar)
y1=lonarr(ntar)

for i=0, ntar-1 do begin

  ; no need to restore whole file
  savfile=(file_search(tardirs[i]+'/*.sav', count=nsav))[0]
  
  obj=obj_new('IDL_Savefile', savfile)
  obj->restore, 'roi_loc'
  
  ; calculate xc and yc from roi_loc
  xc=(roi_loc[1]-roi_loc[0])/2.+roi_loc[0]
  yc=(roi_loc[3]-roi_loc[2])/2.+roi_loc[2]
  
  xpos[i]=xc
  ypos[i]=yc
  
  x1[i]=roi_loc[0]
  y1[i]=roi_loc[2]
  
  roi[i]=targets[i]
  
  print, roi[i], xpos[i], ypos[i]
  
endfor



; plot diagram showing locations of rois on the ccd
;window, 0, xsize=800, ysize=600, xpos=1500, ypos=200
ps_on, fileout, xsize=19, ysize=17
plot, [0,4008,4008,0,0], [0,0,2672,2672,0], color=cgcolor('black'), thick=3, xrange=[0,4008], yrange=[0,2672], $
  xstyle=1, ystyle=1, title=sat+' '+field
plotsym, 0, /fill, 2.
;oplot, xpos, ypos, psym=8, color=cgcolor('purple')
oplot, [0,4008], [2672./2.,2672./2.], linestyle=2, color=cgcolor('blue'), thick=2
oplot, [4008./2.,4008./2.], [0,2672], linestyle=2, color=cgcolor('blue'), thick=2
;for i=0, nsav-1 do rectangle, x1[i], y1[i], 32, 32, color=cgcolor('purple'), thick=2, fill=1, fcolor=cgcolor('purple')
for i=0, ntar-1 do tvbox, 24, xpos[i], ypos[i], color=cgcolor('purple'), fill=1;, fcolor=cgcolor('purple')
; check disances in x&y - if too close then shift up ypos bt 2*25
for kk=0, ntar-1 do begin
  
  xdiff=abs(x1-x1[kk])
  ydiff=abs(ypos-ypos[kk])
  
  xx=where(xdiff lt 400 AND ydiff lt 25, nxx)
  if nxx gt 1 then begin
  
    xx=xx[where(xx ne kk)]
    
    if n_elements(xx) gt 1 then stop
    
    ; change the ypos
    ypos[kk]=ypos[kk]+50
  endif
endfor
for i=0, ntar-1 do xyouts, x1[i], ypos[i]+25, roi[i], charsize=0.5, color=cgcolor('black'), charthick=3
ps_off


spawn, 'open '+fileout+' &'

print, 'End of program'
end