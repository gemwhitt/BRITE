pro newtargets,indir, targets

; Created 11/04/2013 - program to get the identifiers of all new targets in the new field.
; need this for fits_to_sav

filesin=file_search(indir+'/*.fits', count=nf)

;filesin=filesin[0:nf-1]
;nf=nf-387

uroi=[]

for i=0, nf-1, 10 do begin
  
  fits_info, filesin[i], extname=extnames, n_ext=next, /silent
  
  hdnames=extnames[1:next-2]
  
  data0=mrdfits(filesin[i], 0, header, /status, /silent)
  num_roi=sxpar(header, 'NUM_ROI')
  
  ; check num_roi matches with extnames[1:next-2]
  if num_roi ne n_elements(hdnames) then stop
  
  uroi=[uroi,hdnames[uniq(hdnames, sort(hdnames))]]
  
  uroi=uroi[uniq(uroi, sort(uroi))]
  
endfor

targets=uroi
;targets=strcompress(uroi, /remove_all) ; old

end

