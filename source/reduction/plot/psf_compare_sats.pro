pro psf_compare_sats

; PURPOSE: Compare PSF shapes for the same targets at the same temperatures - but different satellites

; Output 1: plot an example ROI from each sat at a similar temperature - give temp as a parameter
; Output 2: Statistics: number of pixels holding 50% flux, 90% flux AND...
;   fit a Gaussian model and measure residual

Compile_opt idl2
  
sat='BA'
field='CENTAURUS'
temp=19

; outdir
outdir='~/BRITE/'+sat+'/'+field+'/plots/p1/PSFs/'
  
; input directory
savdir='~/BRITE/'+sat+'/'+field+'/data/p1/medcol2/'
  
filesin=file_search(savdir+'*.sav', count=nf)
  
fname=file_basename(filesin, '.sav')
  
for ff=0, 0 do begin  ;nf-1 do begin
  obj=obj_new('IDL_Savefile', filesin[ff])
  obj->restore, 'data1'
  obj->restore, 'vmag'
  obj->restore, 'bmag'
  obj->restore, 'roi_name'
  obj->restore, 'ccd_temp'
  
  ; find images in temperature range
  tdiff=abs(ccd_temp[0,*]-temp)
  ; sort in increasing order of difference
  sort1=sort(tdiff)
  tdiffs=tdiff[sort1]
  
  if tdiffs[0] lt 1. then im=sort1[30] else stop
  
  target=roi_name[0]
    
  bminv=bmag-vmag
    
  vmag=strmid(strtrim(vmag,2),0, 4)
  bminv=strmid(strtrim(bminv,2),0, 5)
  
 ; data2=data1[*,*,im]
 ; plot_image, bytscl(data2, 20, 500)
  
 ; im=0
  
 ; stop
    
  data2=data1[*,*,im]*(-1.)
  
 ; stop
    
  fileout=outdir+fname[ff]+'_psf_t'+strtrim(temp,2)+'.ps'
  ps_on, fileout, xsize=19, ysize=17
  plot_image, bytscl(data2, -200, 0), xtitle='x-pixel', ytitle='y-pixel', $
    title=sat+', '+field+', '+target+', V = '+vmag+', B-V = '+bminv+', Temp='+strtrim(temp,2), charsize=0.7, background=cgcolor('white')
  ps_off
  endfor
  
  print, 'end of program'
  print, 'See output files in '+outdir
end

