pro analyse_medcols

; Purpose: Analyse medcol variations and establish which is best for the current set of data 
; Use statistics - make plots and print to file which medcol to remove from data1

Compile_opt idl2
  
sat='UB'
field='CENTAURUS'
  
indir='~/BRITE/'+sat+'/'+field+'/data/p1/2014_0607/' ; CONTAINS ALL FOLDERS FOR THE TARGETS
  
outdir='~/BRITE/'+sat+'/'+field+'/reduction/medcols/2014_0607/'

; check dir exists
chk=file_search(outdir, count=nc)
if nc eq 0 then spawn, 'mkdir -p '+outdir

fileout=outdir+'medcol_results.txt'

filein=file_search(indir+'HD*/*.sav', count=nf)
fname=file_basename(filein, '.sav')
  
for ff=0, nf-1 do begin  ; begin loop over each file
    
    ; restore the whole file - because we are re-saving
    restore, filein[ff] ;roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl, $
    ;simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, medcol1, medcol2, ndead, nnlin, medcol3, flag, pix50
              
    nimg=n_elements(jd)
    
    medres3=[]
    sigres3=[]
    
    medres1=[]
    sigres1=[]
      
    ; loop over each image and each column
    for im=0, nimg-1 do begin
        
        if flag[im] eq 0 then continue  ; established from update_medcol
          
        im0=data1[*,*,im]
        
        s=size(im0, /dim)
        
        ;      plot_image, bytscl(im0, 0, 500)
        ;      stop
        
        ; now remove each medcol in turn, make pix50 pixels NaN's and calculate the mean and sigma as well as the extra number of negative pixels
        
        ; medcol1 - whole frame
        im1=im0*0.
        if s[0] ne (size(medcol1, /dim))[0] then stop
        
        for c=0, s[0]-1 do im1[c,*]=im0[c,*]-medcol1[c,im]
        
        ; medcol1 - raster
        im3=im0*0.
        if s[0] ne (size(medcol3, /dim))[0] then stop
        
        for c=0, s[0]-1 do im3[c,*]=im0[c,*]-medcol3[c,im]
        
        ; location 2d indices of pix50 (pixels to ignore in the calculation)
        bright=where(finite(pix50[*,*,im]) eq 0, nbright)
        if nbright gt 0 then begin
          bright2d=array_indices(im0, bright)
          for pix=0, nbright-1 do im1[bright2d[0,pix], bright2d[1,pix]]=!Values.F_NAN
          for pix=0, nbright-1 do im3[bright2d[0,pix], bright2d[1,pix]]=!Values.F_NAN
        endif
        
        medres1=[medres1,median(im1)]
        sig1=sqrt((total((im1-median(im1))^2, /nan))/float(nimg)-1.)
        sigres1=[sigres1,sig1]
        
        medres3=[medres3,median(im3)]
        sig3=sqrt((total((im3-median(im3))^2, /nan))/float(nimg)-1.)
        sigres3=[sigres3,sig3]
    
    endfor ; end loop over image
    
    medres1=median(medres1)
    medres3=median(medres3)
    
    sigres1=median(sigres1)
    sigres3=median(sigres3)
    
    choice=['medcol1','medcol3']
    
    best=where([abs(medres1),abs(medres3)] eq min([abs(medres1),abs(medres3)]))
    
    if n_elements(best) gt 1 then $ ; choose the one with the lowest sigma
      best=where([sigres1,sigres3] eq min([sigres1,sigres3]))      
   
    result=choice[best]
      
    ; print the results to file
    ; check file exists
    chk=file_search(fileout, count=nchk)
    if nchk eq 0 then begin
      openw, lun, /get_lun, fileout
      printf, lun, 'filename', 'medres1', 'medres3', 'sigres1', 'sigres3', 'best', $
      format='(a20,x,a10,x,a10,x,a10,x,a10,x,a10)'
    free_lun, lun
   endif
      
      openw, lun, /get_lun, /append, fileout
      printf, lun, fname[ff], medres1, medres3, sigres1, sigres3, result, $
        format='(a20, x, d10.5, x, d10.5, x, d10.5, x, d10.5, a10)'
      free_lun, lun
    
  
  endfor  ; end loop over file
        

print, 'End of Progrmam'
print, 'Run remove_medcol_1s with choice of 1,2, or 3
end