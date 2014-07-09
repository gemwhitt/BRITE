pro read_bert

; read berts files - make a test file, do psf-fitting and compare results to gems data

Compile_opt idl2

indir='~/BRITE/data/UB/bert/fits/'

filesin=file_search(indir+'HD36486-5*', count=nf)

outdir='~/BRITE/data/UB/bert/sav/'

fileout=outdir+'Orion-CF1-2_HD36486_p3.sav'

data1=[]
jd=[]
ccd_temp=[]
p2_trend=[]

for i=0, nf-1 do begin
  
  ex0=mrdfits(filesin[i], 0, header, /silent)
  
  jd=[jd,sxpar(header, 'JD-OBS')]
  
  data=mrdfits(filesin[i], 1, header, /silent)
  
  data1=[[[data1]],[[data]]]
  
  telem=mrdfits(filesin[i], 3, header, /silent)
  
  temp0=telem.temp0
  temp1=telem.temp1
  temp2=telem.temp2
  temp3=telem.temp3
  
  ccd_temp=[ccd_temp,(total([temp0,temp1,temp2,temp3]))/4.]
  
  p2_trend=[p2_trend,0.]
  
endfor

; save info
save, filename=fileout, jd, data1, ccd_temp, p2_trend


print, 'end of program'

end


