pro concat_files3

  ; Purpose: concatenate files together for the same target, same roi ....
  
  ; when files have been added for the same field.
  
  ; Use concat_files1 if there are no previous files for that field
  
  sat='UB'
  field='ORION'
  
  indir='~/BRITE/'+sat+'/'+field+'/data/p2/'
  
  storage=indir+'storage/'
  
  ; check output directory exists and make accordingly
  chk=file_search(storage, count=nchk)
  if nchk eq 0 then spawn, 'mkdir -p '+storage
  
  outdir=indir
  
  filesin=file_search(indir+'*.sav', count=nf)
  
  fname=file_basename(filesin, '.sav')
  
  ; find dash positions
  dashpos=strsplit(fname, '_')
  
  fname2=strarr(nf)
  for ff=0, nf-1 do fname2[ff]=strmid(fname[ff], 0, dashpos[ff,3]-1)
  
  ; get number of unique output files to make
  ufile=fname2[uniq(fname2, sort(fname2))]
  nuf=n_elements(ufile)
  
  for ii=0, nuf-1 do begin
  
    files=filesin[where(fname2 eq ufile[ii], nf)]  ; files to be restored and concatenated
    
    ; set up temporary arrays
    exp_num1=[]
    ra_dec1=[]
    jd1=[]
    roi_loc1=[]
    ccd_temp1=[]
    medimg1=[]
    medcol1=[]
    data2=[]
    flag1=[]
    ndead1=[]
    nnlin1=[]
    
    for jj=0, nf-1 do begin
    
      restore, files[jj]  ; roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl, $
      ;simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, medcol, ndead, nnlin, flag
      
      exp_num1=[exp_num1,exp_num]
      ra_dec1=[[ra_dec1],[ra_dec]]
      jd1=[jd1,jd]
      roi_loc1=[[roi_loc1],[roi_loc]]
      ccd_temp1=[[ccd_temp1],[ccd_temp]]
      medimg1=[medimg1,medimg0]
      medcol1=[[medcol1],[medcol]]
      data2=[[[data2]],[[data1]]]
      flag1=[flag1,flag]
      ndead1=[ndead1,ndead]
      nnlin1=[nnlin1,nnlin]
      
      spawn, 'mv '+files[jj]+' '+storage
      
    endfor
    
    exp_num=exp_num1
    ra_dec=ra_dec1
    jd=jd1
    roi_loc=roi_loc1
    ccd_temp=ccd_temp1
    medimg0=medimg1
    medcol=medcol1
    data1=data2
    flag=flag1
    ndead=ndead1
    nnlin=nnlin1
    
    roi_name=roi_name[0]
    exp_time=exp_time[0]
    exp_ttl=exp_ttl[0]
    ;simbad_radec, vmag, bmag, parlax, otype, sptype
    
    ; save output
    fileout=outdir+ufile[ii]+'.sav'
    
    save, filename=fileout, roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl, $
      simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, medcol, ndead, nnlin, flag
      
      
      
  endfor  ; end loop over this unique output file
  
  
  
  print, 'End of program'
  print, 'Check output files'
  print, 'Run map_hps2.pro, then correct_hps.pro'
  
end


