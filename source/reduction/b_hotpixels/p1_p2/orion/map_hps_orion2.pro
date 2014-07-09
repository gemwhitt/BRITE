pro map_hps_orion2

; PURPOSE - first attempt at correcting HPs in teh Orion field - in the absense of information  

Compile_opt idl2
!p.background=cgcolor('white')
plotsym, 0, /fill, 1.3

sat='BA'
field='ORION'

indir='~/BRITE/data/'+sat+'/p1/'+field+'/'
  
outdir='~/BRITE/data/'+sat+'/reduction/hot_pixel_maps/'+field+'/'
  
filesin=file_search(indir+'*.sav', count=nf)
  
for f=10, 10 do begin ; nf-1 do begin
  
  obj=obj_new('IDL_Savefile', filesin[f])
  obj->restore, 'jd'
  obj->restore, 'data1'
  obj->restore, 'ccd_temp'
  obj->restore, 'vmag'
  obj->restore, 'medimg0'
    
  jd1=jd-jd[0]
  
  nfrm=n_elements(jd)
  
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
    
  fname=file_basename(filesin[f],'_p1.sav')
    
  totdn=lonarr(nfrm)
  for img=0, nfrm-1 do totdn[img]=total(data1[*,*,img])
    
  xdim=(size(data1, /dim))[0]
  ydim=(size(data1, /dim))[1]
    
  hps=lonarr(xdim, ydim, nfrm)
  hps2=lonarr(xdim, ydim, nfrm)
    
  jd2=jd1[1:nfrm-1]
  jdiff=jd2-jd1
  gap=where(jdiff gt 0.015, ngap) ; ngap+1 orbits
  gap=[-1,gap,nfrm-1]
  
  num_hp=intarr(ngap+1
    
  for orbit=0, ngap do begin
            
    iloc=indgen(gap[orbit+1]-gap[orbit])+gap[orbit]+1
      
    dat=data1[*,*,iloc]
    
    nimg=n_elements(iloc)
           
    rej1=where(medimg0[iloc] gt 5000 OR totdn[iloc] lt 5000, nrej1, complement=keep1) ; bad images - discard
    
   ; if nrej1 gt 0 then stop ;- analyse discarded images
      
    ngood=n_elements(keep1)
    
    if ngood lt 2 then begin    ; save previous HP values and go to next orbit
        
      for p=0, nimg-1 do hps[*,*,iloc[p]]=hps[*,*,iloc[0]-1]
      
      goto, next_orbit
        
    endif
      
    hpo=lonarr(xdim, ydim, ngood)
         
    for img=0, ngood-1 do begin
      
      data2=dat[*,*,keep1[img]]
      
      ; make any negative pixels 0.0
      zero=where(data2 lt 0, nzero)
      if nzero gt 0 then begin
        xl=reform((array_indices(data2, zero))[0,*])
        yl=reform((array_indices(data2, zero))[1,*])
        data2[xl,yl]=0
      endif
         
      cr1=50.   
      cr2=3
                
      ima=find_hps2(data2,cr1,cr2,iloc[keep1[img]], ww,wx,wy)
        
      ;plot_image, bytscl(data2, 0, 200)
      ;oplot, [wx], [wy], color=cgcolor('orange'), psym=2
      ;stop
         
      nnew=n_elements(ww)
      
      print, nnew
        
      if nnew eq 0 then goto, next_image  ; skip this image
        
      ; record HP locations and values in each image - hps
      ;hps[wx,wy,iloc[keep1[img]]]=data1[wx,wy,iloc[keep1[img]]]
        
      temp_hps=lonarr(xdim, ydim)
      temp_hps[wx,wy]=data2[wx,wy]
        
      hpo[*,*,img]=temp_hps
             
      next_image:
    endfor  ; end loop over this image
      
    ; check frequency of HP occurance for this orbit 
    for ii=0, xdim-1 do begin  ; check every pixel
      for jj=0, ydim-1 do begin
        
        pix=hpo[ii,jj,*]
          
        xx=where(pix gt cr1, nxx)
          
        if nxx gt 2 then $ ; HP!
        hps[ii,jj,iloc]=data1[ii,jj,iloc]
        
      endfor       
    endfor
  
  next_orbit:
  
endfor ; end loop over this orbit
    
; check frequency of HPs within 1-day
result=histogram(jd1, binsize=1, locations=loc, reverse_indices=ri)
    
;nonz=where(result gt 0, nresult)
nresult=n_elements(result)
temp=-1
    
for freq=0, nresult-1 do begin
      
      if result[freq] eq 0 then continue  
      
      if temp[0] ne -1 then ind=[temp,ri[ri[freq]:ri[freq+1]-1]] else ind=[ri[ri[freq]:ri[freq+1]-1]]
      
      nind=n_elements(ind)
      
      if nind lt 100 then begin
        temp=ind
        
        goto, next_freq
      endif else temp=-1
      
      xloc=[]
      yloc=[]
      for ii=0, xdim-1 do begin  ; check every pixel
        for jj=0, ydim-1 do begin
        
          pix=hps[ii,jj,ind]
          
          xx=where(pix gt cr1, nxx)
          
          if float(nxx)/float(nind) gt 0.2 then begin $ ; HP!
            hps2[ii,jj,ind]=data1[ii,jj,ind]
            
            xloc=[xloc,ii]
            yloc=[yloc,jj]
            
          endif
            
        endfor
      endfor
      
     ; for kk=0, nind-1 do begin
     ;   plot_image, bytscl(data1[*,*,ind[kk]], 20, 200)
     ;   oplot, xloc, yloc, color=cgcolor('blue'), psym=8
     ;   
     ;   wait, 0.6
        
     ; endfor

    next_freq:
    endfor
    
    hps=hps2
 
    fileout=outdir+fname+'_hpmap.sav'
    print, fileout
    
    save, filename=fileout, hps
    
  endfor  ;end loop over this file

  print, 'end of program'
end