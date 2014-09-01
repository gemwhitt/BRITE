pro concat_files2

  ; Purpose: concatenate files together for the same target, same roi ....
  
  ; when files have been added for the same field. 
  
  ; Use concat_files1 if there are no previous files for that field
  
  sat='BA'
  field='CENTAURUS'
  
  indir=['~/BRITE/'+sat+'/'+field+'/data/p1/medcol2/','~/BRITE/'+sat+'/'+field+'/data/p1/medcol2/2014_0607/']
    
  outdir='~/BRITE/'+sat+'/'+field+'/data/p1/medcol2/concat/'
  ; check outdir exisits and if not - make it
  chk=file_search(outdir, count=nchk)
  if nchk eq 0 then spawn, 'mkdir -p '+outdir
  
  filesin=file_search(indir+'*.sav', count=nf)
  
  fname=file_basename(filesin, '.sav')
  
  ; get number of unique output files to make
  uf=fname[uniq(fname, sort(fname))]
  nuf=n_elements(uf)
    
    for ii=0, nuf-1 do begin
    
      files=filesin[where(fname eq uf[ii], nf)]  ; files to be restored and concatenated
      
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
    
  
   
  print, 'End of program'
  print, 'Check output files'
  print, 'Run map_hps2.pro, then correct_hps.pro'
  
end


