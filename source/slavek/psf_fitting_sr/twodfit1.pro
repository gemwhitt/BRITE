function TwoDfit1,x,y,wgt
; SMR Dec.2013
;
; two-dimensional LSQ fits: Y = a + b*X
;     errors assuming weights as in the array wgt 
;             typically: wgt ~ X
; for stellar flux determinations: flux ~ b
;     [flux proportional to the slope b]
;
; input: 
;    X & Y: two images, same size, X fitted to Y
;    wgt: array of weights, same size as the images
; output: 
;    5-element vector: 
;        [slope,intercept,err_slope,err_intcept,err_per_pix]
; use:
;    z = TwoDfit1(x,y,wgt)
;
; execution time for 256x256 images: 0.005 sec on my MacBookPr
; predicted fit:
;    y1 = z[1]+z[0]*x

;TIC                     ; used with TOC for timing
   n = n_elements(x)
   tw = total(wgt)
   tx = total(double(x)*wgt)
   ty = total(double(y)*wgt)
   tx2 = total(double(x^2*wgt))
   ty2 = total(double(y^2*wgt))
   txy = total(double(x)*double(y)*wgt)

   d = tw*tx2 - tx^2
   a = (tx2*ty - tx*txy)/d
   b = (tw*txy - tx*ty)/d

   sgm = total((y - a - b*x)^2)/(n-2) ; dispersion per pixel
   erra = sqrt(sgm*tx2/d)             ; error of the offset
   errb = sqrt(sgm*tw/d)              ; error of the slope
   sgm = sqrt(sgm)                    ; error per pixel, in units of Y

   z = [b,a,errb,erra,sgm]
;TOC

return,z
end
