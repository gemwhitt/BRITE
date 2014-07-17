pro fix_hps_orion

; Program to take level p1-images (median subtracted), and produce p2-images, which are cleaned of HPs and CP  
; Use maps 
;
Compile_opt idl2

sat='BA'
field='ORION'

; input directory
indir='~/BRITE/data/'+SAT+'/p1/'+field+'/'
  
outdir='~/BRITE/data/'+sat+'/p2/'+field+'/'
  
hpdir='~/BRITE/data/'+sat+'/reduction/hot_pixel_maps/'+field+'/'
  
filesin=file_search(indir+'*.sav', count=nsav)
fname=file_basename(filesin, '.sav')
  
for kk=0, nsav-1 do begin
  
    ; restore observation file
    restore, filesin[kk]  ;roi_name, exp_num, ra_dec, jd, roi_loc, ccd_temp, exp_time, exp_ttl, $
                          ;simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, ndead, nsat, medcol1, data1
   
    hdname=file_basename(filesin[kk],'_p1.sav')
    
    ;restore hp map
    hpmap=hpdir+hdname+'_hpmap.sav'
    
    restore, hpmap  ; hps[32,32,nfrm]
    
    nfrm=n_elements(jd) ; number of frames
    
    ; cycle through nfrm and replace -'ve pixels with 0
    for im=0, nfrm-1 do begin
      data2=data1[*,*,im]
      zero=where(data2 lt 0, nzero)
      if nzero gt 0 then begin
        xl=reform((array_indices(data2, zero))[0,*])
        yl=reform((array_indices(data2, zero))[1,*])
        data2[xl, yl]=0
        data1[*,*,im]=data2
      endif
    endfor
    
    ; create iflag=0, 1, 2 for bad, blank or target - make everything 2 to begin with
    iflag=lonarr(nfrm)+2
    
    ; miss out "bad" frames - flag 
    totdn=lonarr(nfrm)
    for img=0, nfrm-1 do totdn[img]=total(data1[*,*,img])
    rej1=where(medimg0 gt 5000 OR totdn lt 5000, nrej1, complement=keep1) ; bad images - discard
    iflag[rej1]=0
    
    nfrm2=n_elements(keep1) ; number of frames
;    stop
    ; clean the image of HPs
    for j=0, nfrm2-1 do begin
      
      img0 = reform(data1[*,*,keep1[j]])
      
      hp1 = reform(hps[*,*,keep1[j]])
      
      img1 = clean_img(img0,hp1,j)           ; image cleaned of HPs and CPs
      
      data1[*,*,keep1[j]]=img1

    
    endfor  ; end loop over images
    
    ; save new file
    
    fileout=outdir+hdname+'_p2.sav'

    save, filename=fileout, roi_name, exp_num, ra_dec, jd, roi_loc, ccd_temp, exp_time, exp_ttl, $
    simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, ndead, nsat, medcol1, data1, iflag
      
  endfor  ; end loop over file
  
  print, 'end of program'
end
