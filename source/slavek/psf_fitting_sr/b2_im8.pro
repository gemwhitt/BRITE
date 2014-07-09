function B2_im8, psf,dat,wgt,shx,shy,dud
; SMR Jan.2014
;
; determination of the flux and background using an experimental PSF
; 
; direct 2D LSQ fits using  TwoDfit1.pro
;
; input:
;    psf = PSF used for fits, 256x256 small pixels
;    dat = input images for averaging after shifts, 256x256
;    wgt = array of weights, 256x256
;    shx,shy = star shifts in small pixels
;    dud = idx array: images which should NOT be used are marked by 1
; output:
;    result of function: coefficients of fit as 5-column array
;          [flux, background, error flux, err-backgr, err-per-pixel]
; run: 
;   fit = B2_im8(psf3,dat1,wgt,shx2,shy2,dud2)
     
    n = n_elements(dat[0,0,*])     ; how many images in sequence
    z = fltarr(5,n)                ; 5-col array of results
    x = psf               
;stop
    for i=0,n-1 do begin
       ;print,f='(a,i4)','Image: ',i
       if dud[i] eq 1 then goto,E      ; results set to 0. (all 5 numbers)
       y = reform(dat[*,*,i])          ; take one image
;       stop
       y = shift(y,-shx[i],-shy[i])   
                            ; i-th 2-D image shifted to stellar centre
;stop
                            
       z1 = TwoDfit1(x,y,wgt)
       z[*,i] = z1          ; 5-el vector of results 
    E:
    end
;stop

return,z
end

   