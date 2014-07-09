pro stack_files

; program to stack up N-number of images from single exposures and save separate .sav file .....
;  .... process accordingly
;  
Compile_opt idl2

indir='~/BRITE/data/UB/roi_raw_sav/CENTAURUS/1stk/'

nstk=10

outdir='~/BRITE/data/UB/roi_raw_sav/CENTAURUS/'+strtrim(nstk,2)+'stk/'

filesin=file_search(indir+'*_p0.sav', count=nf)

for i=0, nf-1 do begin
  
  restore, filesin[i]  ; roi_name, exp_num, ra_dec, jd, data1, roi_dim, xc, yc, ccd_temp, exp_time, exp_ttl
  
  ; get x and y dimensions of data1
  xdim=(size(data1, /dim))[0]
  ydim=(size(data1, /dim))[1]
  
  ; find number of orbits - only stack within an orbit
  nfrm=n_elements(jd)
  
  jd1=jd-jd[0]
  jd2=jd1[1:nfrm-1]
  
  jdiff=jd2-jd1
  gap=where(jdiff gt 0.015, ngap)
  gap=[0,gap, nfrm-1]
    
  for j=0, ngap do begin
    
    if j eq 0 then iloc=indgen(gap[j+1]+1) else iloc=indgen(gap[j+1]-gap[j])+gap[j]+1
    
    n=n_elements(iloc)
    
    data2=data1[*,*,iloc]
    jd2=jd[iloc]
    head_temp=ccd_temp[iloc]
    
    ; find number of stacks in orbit
    ns=floor(float(n)/float(nstk))
    
    count=0
    for k=0, ns-1 do begin
      
      x=lonarr(xdim,ydim)
      for l=0, nstk-1 do begin
        x=x+data2[*,*,count]
        count=count+1
      endfor
      
      x=x/float(nstk)
      
      if j eq 0 AND k eq 0 then begin
        data3=x 
        jd3=jd2[count-5]
        ccdt=head_temp[count-5]
      endif else begin
        data3=[[[data3]],[[x]]]
        jd3=[jd3,jd2[count-5]]
        ccdt=[ccdt,head_temp[count-5]]
      endelse
      
    endfor  ; end loop over stacks
    
  endfor  ; end loop over orbit
  
  data1=data3
  exp_ttl=exp_ttl*nstk
  jd=jd3
  ccd_temp=ccdt
  
  nfrm=n_elements(jd)
  
  ; save file
  fname=file_basename(filesin[i], '.sav')
  fileout=outdir+fname+'_'+strtrim(nstk,2)+'.sav'
  
  save, filename=fileout, roi_name, exp_num, ra_dec, jd, data1, roi_dim, xc, yc, ccd_temp, exp_time, exp_ttl
  
endfor  ; end loop over file

print, 'end of program'
end