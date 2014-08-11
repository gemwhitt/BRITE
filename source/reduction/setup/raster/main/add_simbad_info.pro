pro add_simbad_info

; PURPOSE: Append ALL *_p0*.sav files with SIMBAD info
; OUTPUT: _p0.sav ; roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl, $
                  ; simbad_radec, vmag, bmag, parlax, otype, sptype
                  
; REQUIREMENTS: This requires that target info has been obtained from the SIMBAD website - use the script for this

Compile_opt idl2

sat='BA'
field='CENTAURUS'

if field eq 'ORION' then $
  template_file='~/IDLWorkspace82/BRITE/resource/simbad_templates/simbad_brite_orion.sav' 
  
if field eq 'CENTAURUS' then $
  template_file='~/IDLWorkspace82/BRITE/resource/simbad_templates/simbad_brite_centaurus.sav'
  
restore, template_file  ; template1=['ID', 'COO', 'FLX', OTYPE', 'PLX', 'SPTYPE']

simb_results='~/BRITE/simbad/returns/'+field+'.txt'  ; file with simbad results

savdir='~/BRITE/'+sat+'/'+field+'/data/raw_sav/2014_0607/'

result=read_ascii(simb_results, template=template1) ; read the results file using the template

; remove spaces from ID name:
id1=STRCOMPRESS(result.id, /REMOVE_ALL)

; manually modify 1 entry:
for kk=0, n_elements(id1)-1 do if strpos(id1[kk], ',') gt 0 then id1[kk]=strmid(id1[kk], strpos(id1[kk], ',')+1)

; check results are all unique
hd_num=strmid(id1,2)
uniq_id=uniq(hd_num, sort(hd_num))

if n_elements(uniq_id) gt n_elements(id1) then stop ; address accordingly - if multiple entries.. 

; add info to .save files
for i=0, n_elements(id1)-1 do begin  
  
  filesin=file_search(savdir+id1[i]+'/*.sav', count=nsav)
  
  if nsav eq 0 then print, 'No saved directory with the ID: '+id1[i]
  
  for j=0, nsav-1 do begin
     
  restore, filesin[j]  ;roi_name, exp_num, ra_dec, jd, data1, roi_dim, ccd_temp, exp_time, exp_ttl
 
  simbad_radec=STRCOMPRESS(result.coo[i])
  simbad_mags=STRCOMPRESS(result.flx[i])
  parlax=STRCOMPRESS(result.plx[i])
  otype=STRCOMPRESS(result.otype[i])
  sptype=STRCOMPRESS(result.sptype[i])
  
  vmag=float(strmid(simbad_mags, strpos(simbad_mags, 'V=')+2))
  bmag=float(strmid(simbad_mags, strpos(simbad_mags, 'B=')+2, strpos(simbad_mags, 'V=')-strpos(simbad_mags, 'B=')-2))
  
  fileout=filesin[j]
  
  save, filename=fileout, roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl, $
  simbad_radec, vmag, bmag, parlax, otype, sptype
    
  endfor
 

endfor

print, 'End of program'
print, 'Run save_medcol and do analysis on medimg0, medcol1, medcol2 etc'

end