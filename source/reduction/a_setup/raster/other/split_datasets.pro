pro split_datasets

; program to separate any stacked data from non-stacked and move to a different (stacked) file
; 
Compile_opt idl2

indir='/Users/gemmawhittaker/BRITE/data/UB/roi_raw_sav/ORION/'

outdir='/Users/gemmawhittaker/BRITE/data/UB/roi_raw_sav/ORION/all_frames/'

filesin=file_search(indir+'*.sav', count=nf)

for i=1, nf-1 do begin
  
  restore, filesin[i] ;roi_name, exp_num, ra_dec, jd, data1, roi_dim, xc, yc, ccd_temp, exp_time, exp_ttl
  
  fname=file_basename(filesin[i], '_p0.sav')
  
  ; select good data
  non_stk=where(exp_ttl eq 1., ngood)
  
  if ngood eq n_elements(jd) then CONTINUE
  
  ; establish new arrays
  roi_name=roi_name[0]
  exp_num=exp_num[non_stk]
  ra_dec=ra_dec[*,non_stk]
  jd=jd[non_stk]
  data1=data1[*,*,non_stk]
  roi_dim=roi_dim[*,non_stk]
  ccd_temp=ccd_temp[non_stk]
  exp_time=exp_time[non_stk]
  exp_ttl=exp_ttl[non_stk]
  
  ; save as new file and move old file
  spawn, 'mv '+filesin[i]+' '+outdir+fname+'_allfrm_p0.sav'
  
  fileout=filesin[i]
  
  save, filename=fileout, roi_name, exp_num, ra_dec, jd, data1, roi_dim, xc, yc, ccd_temp, exp_time, exp_ttl
  
endfor

print, 'end of program'
end