FUNCTION mimic_ghost, x, p
  xshift = p[0]
  yshift = p[1]
;  o0 = p[2]
  o1 = p[2]
  o2 = p[3]
;  o3 = p[5]
;  o4 = p[6]

  new = fshift(x, xshift, yshift)

  return, o1*new + o2*new^2. ;+ o3*new^3. ;+ o4*new^4.

END


FUNCTION remove_ghost, file, ghost_yhi, xshift=xshift, yshift=yshift, scale1=scale1, scale2=scale2

im = mrdfits(file, 0, h)

dimen = size(im)

ghost_xlo = [183, 230, 270]
ghost_xhi = [203, 250, 290]

ghost_ylo = 0
;ghost_yhi = 220
;ghost_yhi = 120
;ghost_yhi = 105

source_xlo = 120
source_xhi = 160

subim = im[*,ghost_ylo:ghost_yhi]

backsub = fltarr(dimen[1], ghost_yhi-ghost_ylo+1)
aindex = [indgen(source_xlo), indgen(ghost_xlo[0]-source_xhi+1)+source_xhi, indgen(ghost_xlo[1]-ghost_xhi[0]+1)+ghost_xhi[0], indgen(ghost_xlo[2]-ghost_xhi[1]+1)+ghost_xhi[1], indgen(dimen[1]-ghost_xhi[2]+1)+ghost_xhi[2]] 

for i = 0, ghost_yhi-ghost_ylo do begin

  out = bspline_iterfit((findgen(dimen[1]))[index], subim[index, i], nord=3, bkspace=40, invvar=1/sqrt(subim[index, i]))
  b = bspline_valu(findgen(dimen[1]), out)
  subim[*,i] = subim[*,i] - b
endfor


test = subim

ghost0 = subim[ghost_xlo[0]:ghost_xhi[0], *]
ghost1 = subim[ghost_xlo[1]:ghost_xhi[1], *]
ghost2 = subim[ghost_xlo[2]:ghost_xhi[2], *]

start = [0.1, -10, 0.2, 0.01]
out1 = mpfitfun('mimic_ghost', ghost0, ghost1, sqrt(abs(ghost1)) > 1, start)
start = [0.1, -10, 0.2, 0.01]
;start = [xshift+out1[0], yshift+out1[1], scale1*out1[2], scale2*out1[3]]
;pi = replicate({fixed:1, limited:[0,0], limits:[0.D,0.D]},4)
out2 = mpfitfun('mimic_ghost', ghost0, ghost2, sqrt(abs(ghost2)) > 1, start, parinfo=pi)

yshift = 10
test[ghost_xlo[1]:ghost_xhi[1], 0:(ghost_yhi-ghost_ylo)-yshift] = subim[ghost_xlo[1]:ghost_xhi[1], 0:(ghost_yhi-ghost_ylo)-yshift] - subim[ghost_xlo[0]:ghost_xhi[0], yshift:*]*0.25


subghost1 = ghost1 - mimic_ghost(ghost0, out1)
subghost2 = ghost2 - mimic_ghost(ghost0, out2)
subghost3 = ghost2 - mimic_ghost(ghost0, [xshift+out1[0], yshift+out1[1], scale1*out1[2], scale2*out1[3]])


;writefits, 'temp.fits', test, h
writefits, 'temp.fits', ghost1, h
writefits, 'temp2.fits', subghost1, h
writefits, 'temp3.fits', ghost2, h
writefits, 'temp4.fits', subghost2, h
writefits, 'temp5.fits', subghost3, h

print, out2[0]-out1[0], out2[1]-out1[1], out2[2]/out1[2], out2[3]/out1[3]
return,[out1, out2]
END

;    P(0) =            0.0175038  
;    P(1) =             -9.86690  
;    P(2) =             0.243509  
;    P(3) =             0.114969  
;    P(4) =          0.000294144  

;    P(0) =            0.0308853  
;    P(1) =             -9.46994  
;    P(2) =           -0.0675198  
;    P(3) =             0.152524  
;    P(4) =          0.000227861  

