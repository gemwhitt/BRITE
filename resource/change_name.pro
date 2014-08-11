pro change_name

indir='/Users/gemmawhittaker/BRITE/data/UB/p4/'

filesin=file_search(indir+'*CF1-*.sav', count=nsav)

fname=file_basename(filesin, '_p2.sav_p4_gem2.sav')

for ii=0, nsav-1 do begin
  
  newname=indir+fname[ii]+'_p4_gem2.sav'

  spawn, 'mv '+filesin[ii]+' '+newname
  
  
endfor
print, 'end of program'
end