function b2_im2b, im,hp_xy
; SMR Oct.2013
;
; cleans the image by medianing in small regions 3x3 pix at locations ww
;      takes care of the edges
; input:
;   im = image 32x32
;   ww = indices (2-D) where bad pixels found
; output:
;   im1 = cleaned image
; uses: 
;   where2d.pro
; run: 
;   im1 = B2_im2(im,ww)
  
    x = im

    wx=hp_xy[0,*]
    wy=hp_xy[1,*]
    
    n = (size(hp_xy, /dim))[1]
    
    ; add on border columns and rows
    med_frm=median(x)
    x2=(lonarr(34,34)*0.)+med_frm
    x2[1:32,1:32]=x

    for i=0,n-1 do begin
          k1 = wx[i]+1 & k2 = wy[i]+1         
          ;if k1 eq 0 then  k1=1
          ;if k1 eq 31 then k1=30 
          ;if k2 eq 0 then  k2=1
          ;if k2 eq 31 then k2=30 
    x[wx[i],wy[i]] = median(x2[k1-1:k1+1,k2-1:k2+1]) ; median in 3x3
    end

return,x
end

   