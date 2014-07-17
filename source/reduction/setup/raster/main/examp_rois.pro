pro examp_rois

; PURPOSE: Produce ROI images to show different PSF shapes....
  
Compile_opt idl2

sat='UB'
field='CENTAURUS'

; input directory
indir='~/BRITE/data/'+sat+'/p2/'+field+'/'
  
savfiles=file_search(indir+'*.sav', count=nsav)
  
fnames=file_basename(savfiles, '_p2.sav')
targets=fnames
  
outdir='~/BRITE/data/'+sat+'/plots/'+field+'/PSFs/'
    
for kk=0, nsav-1 do begin
  
  obj=obj_new('IDL_Savefile', savfiles[kk])
  obj->restore, 'data1'
  obj->restore, 'vmag'
  obj->restore, 'bmag'
  obj->restore, 'roi_name'
  
  bminv=bmag-vmag
  
  vmag=strmid(strtrim(vmag,2),0, 4)
  bminv=strmid(strtrim(bminv,2),0, 5)
  
    
  hdname=roi_name[0]
    
  ; determine number of images
  npts=n_elements(jd)
    
  data2=data1[*,*,0]*(-1.)
    
  fileout=outdir+hdname+'_psf.ps'
  ps_on, fileout, xsize=19, ysize=17
  plot_image, bytscl(data2, -200, 0), xtitle='x-pixel', ytitle='y-pixel', $
    title=hdname+', V = '+vmag+', B-V = '+bminv, charsize=0.7, background=cgcolor('white')
  ps_off
    
endfor
  
print, 'end of program'
end

