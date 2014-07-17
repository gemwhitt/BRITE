pro simbad_info

; PURPOSE: Append _p0.sav files with SIMBAD info
; OUTPUT: _p0.sav ; roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl, $
                  ; simbad_radec, vmag, bmag, parlax, otype, sptype

Compile_opt idl2

template_file='~/IDLWorkspace82/BRITE/resource/simbad_templates/simbad_brite_orion1.sav'  ; change template for each field 
restore, template_file  ; template1=['ID', 'COO', 'FLX', OTYPE', 'PLX', 'SPTYPE']

simb_results='~/BRITE/data/UB/simbad/returns/orion1.txt'  ; file with simbad results

sat='BA'
field='ORION'

savdir='~/BRITE/data/'+sat+'/roi_raw/'+field+'/sav/'

result=read_ascii(simb_results, template=template1) ; read the results file using the template

; remove spaces from ID name:
id1=STRCOMPRESS(result.id, /REMOVE_ALL)

; check results are all unique
hd_num=strmid(id1,2)
uniq_id=uniq(hd_num, sort(hd_num))

if n_elements(uniq_id) gt n_elements(id1) then stop ; address accordingly - if multiple entries.. 

; add info to .save files
for i=0, n_elements(id1)-1 do begin  
  
  savfile=file_search(savdir+'*'+id1[i]+'*.sav', count=nsav)
  
  if nsav eq 0 then stop 
  
  for j=0, nsav-1 do begin
     
  restore, savfile[j]  ;roi_name, exp_num, ra_dec, jd, data1, roi_dim, ccd_temp, exp_time, exp_ttl
 
  simbad_radec=STRCOMPRESS(result.coo[i])
  simbad_mags=STRCOMPRESS(result.flx[i])
  parlax=STRCOMPRESS(result.plx[i])
  otype=STRCOMPRESS(result.otype[i])
  sptype=STRCOMPRESS(result.sptype[i])
  
  vmag=float(strmid(simbad_mags, strpos(simbad_mags, 'V=')+2))
  bmag=float(strmid(simbad_mags, strpos(simbad_mags, 'B=')+2, strpos(simbad_mags, 'V=')-strpos(simbad_mags, 'B=')-2))
  
  fileout=savfile[j]
  
  save, filename=fileout, roi_name, exp_num, ra_dec, jd, data1, roi_loc, ccd_temp, exp_time, exp_ttl, $
  simbad_radec, vmag, bmag, parlax, otype, sptype
    
  endfor
 

endfor

print, 'End of program'
print, 'Do median column removal - with median_column2.pro'
end