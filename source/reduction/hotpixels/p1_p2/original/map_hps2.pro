pro map_hps2

; PURPOSE - Use the LABEL_REGION method to locate and record the majority of HPs in each frame
; Also assign a "flag" to indicate bad frames (based on medimg0 and totdn), which subsequently DO NOT have HP maps 
; flag=0,1,2 to indicate: bad, blank or target
; 
; Save map_01 and map_int - to represent a 0/1 map and the recording of the actual data values for all HPs and their CPs
; Save flag - add this to .sav files in HP removal program
; 
; Also save the medcol1 OR medcol2 array
; 
; Output: map_01, map_int, flag, medcol...1/2
  
Compile_opt idl2

sat='LEM'
field='CENTAURUS'
thr=50  ; threshold above which HPs are identified
  
indir='~/BRITE/'+sat+'/'+field+'/data/p1/medcol2/'
filesin=file_search(indir+'*.sav', count=nf)
  
outdir='~/BRITE/'+sat+'/'+field+'/reduction/hp_maps/'
  
for f=0, nf-1 do begin
  
  obj=obj_new('IDL_Savefile', filesin[f])
  obj->restore, 'jd'
  obj->restore, 'data1'
  obj->restore, 'medimg0'
  obj->restore, 'medcol2'
    
  jd1=jd-jd[0]
    
  nfrm=n_elements(jd)
    
  fname=file_basename(filesin[f],'.sav')
    
  totdn=lonarr(nfrm)
  for img=0, nfrm-1 do totdn[img]=total(data1[*,*,img])
    
  s=size(data1, /dim)
  xdim=s[0]
  ydim=s[1]
    
  hps01=lonarr(xdim, ydim, nfrm)
  hps_int=lonarr(xdim, ydim, nfrm)
  
  flag=intarr(nfrm)+2 ; 0=bad, 1=blank, 2=target - to begin assume everything is target
  
  ; find out where orbits begin and end  
  jd2=jd1[1:nfrm-1]
  jdiff=jd2-jd1
  gap=where(jdiff gt 0.015, ngap) ; ngap+1 orbits
  gap=[-1,gap,nfrm-1]
        
  for orbit=0, ngap do begin                                  ; LOOP OVER EACH ORBIT
    
    iloc=indgen(gap[orbit+1]-gap[orbit])+gap[orbit]+1         ; INDICES OF IMAGES IN THIS ORBIT
      
    imgs=data1[*,*,iloc]    
    nimg=n_elements(iloc)
      
    rej1=where(medimg0[iloc] gt 5000 OR totdn[iloc] lt 5000, nrej1, complement=keep1) ; bad images - discard
      
    ngood=n_elements(keep1)
      
    if nrej1 gt 0 then begin                                  ;stop - analyse discarded images
      
      ; change flag for this image to 0=bad
      flag[iloc[rej1]]=0
        
      if ngood eq 0 then goto, next_orbit else iloc=iloc[keep1]
      
    endif
      
    ; Start
    hpo=lonarr(xdim, ydim, ngood)                         ; array to collect HPs from each frame in this orbit
      
    for img=0, ngood-1 do begin
      
      im=imgs[*,*,keep1[img]]
             
      im1=lonarr(xdim+2,ydim+2)   ; append array with extra row/column (for label_region.pro)
      im1[1,1]=im      
        
      r1=label_region(im1 gt thr, /all_neighbors)  ; find groups of pixels above thr
      r1=r1[1:s[0],1:s[1]]                         ; trim off border
        
      ; do a histogram of r1 and save indices in ri
      result=histogram(r1, binsize=1, reverse_indices=ri, locations=loc)
        
      i=0
      repeat begin
        
        ind=ri[ri[i]:ri[i+1]-1] ; 1d index of HPs - convert to 2d
          
        if n_elements(ind) le 3 then begin
          ind2=array_indices(r1, ind)
            
          hpo[ind2[0,*],ind2[1,*],img]=1
            
        endif
        i=i+1
      endrep until i eq n_elements(result)
        
    endfor
      
    ; check frequency of HP occurance for this orbit
    for ii=0, xdim-1 do begin  ; check every pixel in x and y
      for jj=0, ydim-1 do begin
        pix=hpo[ii,jj,*]      
        xx=where(pix gt 0, nxx)
        if nxx ge 0.2*ngood then begin
          hps01[ii,jj,iloc]=1                   ; save 01 array
          
          ; save hps_int array - for later analysis - record HP & CP data values 
          for kk=0, ngood-1 do begin
            hps_int[ii,jj,iloc[kk]]=data1[ii,jj,iloc[kk]]
            if ii le xdim-2 then hps_int[ii+1,jj,iloc[kk]]=data1[ii+1,jj,iloc[kk]]
          endfor
        endif
      endfor
    endfor
            
    next_orbit:
  endfor ; end loop over this orbit
    
  fileout=outdir+fname+'_hpmaps2.sav'
    
  save, filename=fileout, hps01, hps_int, flag, medcol2
    
endfor  ;end loop over this file
  
print, 'end of program'
end