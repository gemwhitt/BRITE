function TwoDfit,x,y
; SMR Dec.2013
;
; two-dimensional LSQ fits
;     errors assuming same weights for all pixels
; input: 
;    X & Y: two images, same size, X fitted to Y
; output: 
;    5-element vector: 
;        [slope,intercept,err_slope,err_intcept,err_per_pix]
; use:
;    z = TwoDfit(x,y)
;
; execution time for 256x256 images: about 0.003 sec on my MacBookPro

   n = n_elements(x)

   tx = total(double(x))
   ty = total(double(y))
   tx2 = total(double(x^2))
   ty2 = total(double(y^2))
   txy = total(double(x)*double(y))

   d = n*tx2 - tx^2
   a = (tx2*ty - tx*txy)/d
   b = (n*txy - tx*ty)/d
;   z = [b,a]     ; b = slope, a = baseline/offeset

   sgm = total((y - a - b*x)^2)/(n-2)
   erra = sqrt(sgm*tx2/d)    ; error offset
   errb = sqrt(sgm*n/d)      ; error slope

   sgm = sqrt(sgm)      ; error for a single pixel

   z = [b,a,errb,erra,sgm]

return,z
end
