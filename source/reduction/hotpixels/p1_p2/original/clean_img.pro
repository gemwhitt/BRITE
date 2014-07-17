function clean_img, im,hp1,img
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
  
  ; called by p1_p2
  
  x = im
  
  xdim=(size(x, /dim))[0]
  ydim=(size(x, /dim))[1]
  
  ;get locations of HPs in this frame from hp1
  ;hp=where(hp1[*,*,img] eq 1, nhp)  ; or
  hp=where(hp1[*,*] ne 0, nhp)  ;new
  
  if nhp eq 0 then goto, endofclean
  
  hp_xy=lonarr(2,nhp)
  hp_xy[0,*]=(array_indices(x, hp))[0,*]
  hp_xy[1,*]=(array_indices(x, hp))[1,*]
  
  wx=hp_xy[0,*]
  wy=hp_xy[1,*]
  
  n = n_elements(wx)
  
  ; add on border columns and rows
  med_frm=median(x)
  x2=(lonarr(xdim+2,ydim+2)*0.)+med_frm
  x2[1:xdim,1:ydim]=x
  
  for i=0,n-1 do begin
    k1 = wx[i]+1 & k2 = wy[i]+1
    
    cp1=wx[i]+2 & cp2 = wy[i]+1
    
    x[wx[i],wy[i]] = median(x2[k1-1:k1+1,k2-1:k2+1]) ; median in 3x3  ; ORIGINAL
    
  ;  x[wx[i],wy[i]] = fix(robust_mean(x2[k1-2:k1+2,k2-2:k2+2], 2)) ; median in 5x5  ; NEW
    
    
    if wx[i] lt xdim-1 then x[wx[i]+1,wy[i]] = fix(robust_mean(x2[cp1-2:cp1+1,cp2-1:cp2+1], 2))
    
  end
  
  ;plot_image, bytscl(im,20,100)
  ;oplot,wx, wy, color=cgcolor('orange'), psym=2
  
;wait, 0.1
endofclean:

  return,x
end

