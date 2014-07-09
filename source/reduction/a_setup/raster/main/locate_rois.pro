pro locate_rois

; Program to locate the ROI's on the whole CCD, to determine central, lower left, top left, lower right, top right....
; 
Compile_opt idl2

sat='BA'

field='ORION'

indir='~/BRITE/data/'+sat+'/roi_raw_sav/'+field+'/'

outdir='~/BRITE/data/'+sat+'/plots/'
fileout=outdir+sat+'-'+field+'_ccdloc.ps'

savfiles=file_search(indir+'*.sav', count=nsav)
targets=file_basename(savfiles, '_p0.sav')

; define array variables to record xc and yc for each roi
xpos=lonarr(nsav)
ypos=lonarr(nsav)
roi=strarr(nsav)
x1=lonarr(nsav)
y1=lonarr(nsav)

for i=0, nsav-1 do begin

  restore, savfiles[i]  ;roi_name, exp_num, ra_dec1, helio_jd, data1, roi_loc
  
  ; calculate xc and yc from roi_loc
  xc=(roi_loc[1]-roi_loc[0])/2.+roi_loc[0]
  yc=(roi_loc[3]-roi_loc[2])/2.+roi_loc[2]
  
  xpos[i]=xc
  ypos[i]=yc
  
  x1[i]=roi_loc[0]
  y1[i]=roi_loc[2]
  
  roi[i]=roi_name[0]
  
  print, roi[i], xpos[i], ypos[i]
  
endfor

; find ROI closest to the center of the CCD
ccd_xc=4008./2.
ccd_yc=2672./2.

dist=sqrt((ccd_xc-xpos)^2+(ccd_yc-ypos)^2)
reorder=sort(dist)
; center is at reorder[0]
cnt_roi=reorder[0]

ex1=where(xpos eq 1040 AND ypos eq 2067)
ex1=where(xpos eq 2065 AND ypos eq 2067)


; plot diagram showing locations of rois on the ccd
;window, 0, xsize=800, ysize=600, xpos=1500, ypos=200
ps_on, fileout, xsize=17, ysize=15
plot, [0,4008,4008,0,0], [0,0,2672,2672,0], color=cgcolor('black'), thick=3, xrange=[0,4008], yrange=[0,2672], $
  xstyle=1, ystyle=1
plotsym, 0, /fill, 2.
;oplot, xpos, ypos, psym=8, color=cgcolor('purple')
oplot, [0,4008], [2672./2.,2672./2.], linestyle=2, color=cgcolor('blue'), thick=2
oplot, [4008./2.,4008./2.], [0,2672], linestyle=2, color=cgcolor('blue'), thick=2
;for i=0, nsav-1 do rectangle, x1[i], y1[i], 32, 32, color=cgcolor('purple'), thick=2, fill=1, fcolor=cgcolor('purple')
for i=0, nsav-1 do tvbox, 24, xpos[i], ypos[i], color=cgcolor('purple'), fill=1;, fcolor=cgcolor('purple')
for i=0, nsav-1 do xyouts, x1[i], ypos[i]+25, roi[i], charsize=0.6, color=cgcolor('black'), charthick=3
ps_off

print, 'ROI closest to center of CCD is '+(roi[cnt_roi])[0]

spawn, 'open '+fileout+' &'

stop
print, 'End of program'


end