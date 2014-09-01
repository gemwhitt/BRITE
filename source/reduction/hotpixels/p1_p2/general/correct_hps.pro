pro correct_hps

; Correct images for HPs and CPs using the hpmaps and a simple median replacement from the surrounding pixels
; 
; Also append the output files with flag - to denote good/bad images
; 
; Do not "clean" bad images (with flag=0)
; 
Compile_opt idl2

sat='TOR'
field='CENTAURUS'

mapdir='~/BRITE/'+sat+'/'+field+'/reduction/hp_maps/'
maps=file_search(mapdir+'*hpmaps2.sav', count=nmaps)

if nmaps eq 0 then print, 'No hp_maps found - run map_hps2.pro first!'

savdir='~/BRITE/'+sat+'/'+field+'/data/p1/medcol2/'

outdir='~/BRITE/'+sat+'/'+field+'/data/p2/'

for f=0, nmaps-1 do begin

  obj=obj_new('IDL_Savefile', maps[f])
  obj->restore, 'hps01'
  obj->restore, 'flag'
  
  target_file=file_basename(maps[f],'_hpmaps2.sav')
  
  restore, savdir+target_file+'.sav'  ;roi_name, exp_num, ra_dec, jd, roi_loc, ccd_temp, exp_time, exp_ttl, $
                                      ;simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, ndead, nsat, 
                                      ;medcol1, medcol2, data1
  stop     
  nfrm=n_elements(jd)
  
  data2=data1*0
  
  sdata=size(data2, /dim)
  
  for im=0, nfrm-1 do begin
    
    if flag[im] eq 0 then continue
    
    im2=data1[*,*,im]
    hps=hps01[*,*,im]
    
    ind=where(hps gt 0, nxy)
    if nxy eq 0 then continue
    
    xy=array_indices(hps, ind)
    
    xi=xy[0,*]
    yi=xy[1,*]
    
    for ii=0, nxy-1 do begin
      x1=xi[ii]-1 > 0
      x2=xi[ii]+1 < (sdata[0]-1)
      y1=yi[ii]-1 > 0
      y2=yi[ii]+1 < (sdata[1]-1)
      
      im2[xi[ii],yi[ii]]=median(im2[x1:x2,y1:y2])
      
      ; and do CPs
      if xi[ii] le (sdata[0]-3) then im2[xi[ii]+1,yi[ii]]=median(im2[x1+1:x2+1,y1:y2])  
    endfor
  
  data2[*,*,im]=im2 
  endfor
  
  ; save cleaned images!
  data1=data2 
  
  fileout=outdir+target_file+'_p2.sav'
  
  save, filename=fileout, roi_name, exp_num, ra_dec, jd, roi_loc, ccd_temp, exp_time, exp_ttl, $
    simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, ndead, nsat, medcol1, medcol2, data1, flag
  
endfor

print, 'end of program'
end
