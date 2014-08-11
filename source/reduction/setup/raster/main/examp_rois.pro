pro examp_rois

; PURPOSE: Produce ROI images to show different PSF shapes....
  
Compile_opt idl2

sat='BA'
field='CENTAURUS'

im=0

; outdir
outdir='~/BRITE/'+sat+'/'+field+'/plots/p1/PSFs/'

; input directory
savdir='~/BRITE/'+sat+'/'+field+'/data/p1/medcol2/'

filesin=file_search(savdir+'*.sav', count=nf)

fname=file_basename(filesin, '.sav')
  
for ff=0, nf-1 do begin
    obj=obj_new('IDL_Savefile', filesin[ff])
    obj->restore, 'data1'
    obj->restore, 'vmag'
    obj->restore, 'bmag'
    obj->restore, 'roi_name'
    obj->restore, 'ccd_temp'
    
    stop
    
    target=roi_name[0]
    
    bminv=bmag-vmag
    
    vmag=strmid(strtrim(vmag,2),0, 4)
    bminv=strmid(strtrim(bminv,2),0, 5)
       
    ; determine number of images
    nimg=n_elements(jd)
    
    data2=data1[*,*,im]*(-1.)
    
    fileout=outdir+fname[ff]+'_psfL.ps'
    ps_on, fileout, xsize=19, ysize=17
    plot_image, bytscl(data2, -200, 0), xtitle='x-pixel', ytitle='y-pixel', $
      title=sat+', '+field+', '+target+', V = '+vmag+', B-V = '+bminv, charsize=0.7, background=cgcolor('white')
    ps_off
  endfor
  
print, 'end of program'
print, 'See output files in '+outdir
end

