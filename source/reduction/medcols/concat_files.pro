pro concat_files

; Purpose: concatenate files together for the same target, same roi - after med columns removal

; destroy the directories in the process

sat='BA'
field='ORION'
  
indir='~/BRITE/'+sat+'/'+field+'/data/p1/medcol2/'
  
tardir=file_search(indir+'HD*/', count=ntar)

outdir=indir

for tar=0, ntar-1 do begin
  
  filesin=file_search(tardir[tar]+'/*.sav', count=nf)
  
  fname=file_basename(filesin, '.sav')
  
  dash=strsplit(fname, '_')
  
  hd_roi=strarr(nf)
  ind=strarr(nf)
  
  for ff=0, nf-1 do hd_roi[ff]=strmid(fname[ff], 0, dash[ff,3]-1)
  for ff=0, nf-1 do ind[ff]=strmid(fname[ff], dash[ff,4])
  
  ; get number of unique output files to make
  uf=hd_roi[uniq(hd_roi, sort(hd_roi))]
  nuf=n_elements(uf)
  
  for ii=0, nuf-1 do begin
    
    files=filesin[where(hd_roi eq uf[ii], nf)]  ; files to be restored and concatenated
    
    ; set up temporary arrays
    exp_num1=[]
    ra_dec1=[]
    jd1=[]
    roi_loc1=[]
    ccd_temp1=[]
    medimg01=[]
    medcol11=[]
    medcol21=[]
    data11=[]
    
    for jj=0, nf-1 do begin
      
    restore, files[jj]  ; roi_name, exp_num, ra_dec, jd, roi_loc, ccd_temp, exp_time, exp_ttl, $
                        ; simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, ndead, nsat, medcol1, medcol2, data1
                          
    exp_num1=[exp_num1,exp_num]
    ra_dec1=[[ra_dec1],[ra_dec]]
    jd1=[jd1,jd]
    roi_loc1=[[roi_loc1],[roi_loc]]
    ccd_temp1=[[ccd_temp1],[ccd_temp]]
    medimg01=[medimg01,medimg0]
    medcol11=[[medcol11],[medcol1]]
    medcol21=[[medcol21],[medcol2]]
    data11=[[[data11]],[[data1]]]     
               
    endfor
    
    exp_num=exp_num1
    ra_dec=ra_dec1
    jd=jd1
    roi_loc=roi_loc1
    ccd_temp=ccd_temp1
    medimg0=medimg01
    medcol1=medcol11
    medcol2=medcol21
    data1=data11
    
    roi_name=roi_name[0]
    exp_time=exp_time[0]
    exp_ttl=exp_ttl[0]
    ;simbad_radec, vmag, bmag, parlax, otype, sptype
    
    ; save output
    fileout=outdir+uf[ii]+'.sav'
    
    save, filename=fileout, roi_name, exp_num, ra_dec, jd, roi_loc, ccd_temp, exp_time, exp_ttl, $
      simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, ndead, nsat, medcol1, medcol2, data1
          
  endfor  ; end loop over this unique output file

endfor  ; end loop over this target


for ii=0, ntar-1 do spawn, 'rm -r '+tardir[ii]
 
print, 'End of program'
print, 'deal with HPs and CPs'

end


