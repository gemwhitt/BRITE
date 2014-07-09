pro check_roi_shifts

Compile_opt idl2

; program to check if, when and where any target ROIs shift location - use roi_dim
; 
indir='/Users/gemmawhittaker/BRITE/data/UB/p1/ORION/'

outdir='/Users/gemmawhittaker/BRITE/data/UB/reduction/p1_info/'

fileout=outdir+'Orion_rasterShifts.txt'

filesin=file_search(indir+'HD33111*.sav', count=nf)

if nf eq 0 then begin
  print, 'No files at this location'
  stop
endif

xshifts=[]
yshifts=[]

for f=0, nf-1 do begin
  
  fname=file_basename(filesin[f], '_p1.sav')
  
  obj=obj_new('IDL_savefile', filesin[f])
  obj->restore, 'roi_dim'
  obj->restore, 'roi_name'
  obj->restore, 'data1'
  obj->restore, 'jd'
  destroy, obj
  
 ; for i=1970, 2070 do begin
 ;   plot_image, bytscl(data1[*,*,i],20,200)
    
 ;   stop
 ; endfor
 ; stop
  
  ; check ROI names match with filename
  uroiname=roi_name[uniq(roi_name, sort(roi_name))]
  
  if n_elements(uroiname) gt 1 then stop
  
  match, fname, uroiname, suba, subb, count=match1
  
  if match1 ne 1 then stop
    
  nfrm=(size(roi_dim, /dim))[1]
  
  ; get corner locations
  x1=reform(roi_dim[0,*])
  y1=reform(roi_dim[2,*])
  
  jd1=jd-jd[0]
  jd2=jd1[1:nfrm-1]
  jdiff=jd2-jd1
  gap=where(jdiff gt 0.015, ngap)
  if ngap eq 0 then stop
  gap=[-1,gap,nfrm-1]
  stop
  for orbit=0, ngap do begin
    
    iloc=indgen(gap[orbit+1]-gap[orbit])+gap[orbit]+1
    
    x1_y1=strtrim(x1[iloc],2)+'_'+strtrim(y1[iloc],2)
    
    uloc=uniq(x1_y1, sort(x1_y1))
    nuniq=n_elements(uloc)
    
    if nuniq gt 1 then stop ; shift within orbit
    
    ; check location with previous orbit
    if orbit ne 0 then begin
      
      match, x1y1, x1_y1[uloc], suba, subb, count=count1
      
      if count1 ne 1 then begin
        
        xdiff=x1[iloc[0]]-x1[iloc[0]-1]
        ydiff=y1[iloc[0]]-y1[iloc[0]-1]
        
      openw, lun, fileout, /get_lun, /append
      printf, lun, fname, roi_name[0], orbit, iloc[0], x1y1, x1_y1[uloc], xdiff, ydiff, format='(a10,x,a10,x,i4,x,i7,x,a10,x,a10,x,i7,x,i7)'
      free_lun, lun
        
      endif
      
    endif 
    
    x1y1=x1_y1[uloc]
      
  endfor
  

endfor


spawn, 'open '+fileout+' &'

print, 'End of Program'

stop
end