; This is code to grow or shrink data1 or new data to accomodate changing raster sizes.....
; 
if n_elements(data1) eq 0 then data1=data else begin

; check size of new data
sn=size(data, /dim)

; check size of previous data
sp=size(data1, /dim)

;if dimensions agree then just save as normal
if (sn[0] eq sp[0]) AND (sn[1] eq sp[1]) then data1=[[[data1]],[[data]]] else begin

  ; if new data SMALLER than old data then expand new data to fit
  if sn[0] lt sp[0] then begin
  
    ; grow new data...
    temp_data=lonarr(sp[0],sp[1])
    temp_data[0:sn[0]-1,0:sn[1]-1]=data
    
    data1=[[[data1]],[[temp_data]]]
    
  endif else begin
  
    ; grow old data
    temp_data=lonarr(sn[0],sn[1],so[2])
    
    for ss=0, so[2]-1 do temp_data[0:so[0]-1,0:so[1]-1,ss]=data1[*,*,ss]
    
    data1=[[[temp_data]],[[data]]]
    
  endelse
  
endelse

endelse