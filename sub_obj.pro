FUNCTION psf, x, p

return, p[0] + gauss1(x, [p[1], p[2], p[3]]) + gauss1(x, [p[4], p[5], p[6]]) + gauss1(x, [p[7], p[8], p[9]])

END

FUNCTION sub_obj, file

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

subim = im

index = [indgen(source_xlo), indgen(ghost_xlo[0]-source_xhi+1)+source_xhi, indgen(ghost_xlo[1]-ghost_xhi[0]+1)+ghost_xhi[0], indgen(ghost_xlo[2]-ghost_xhi[1]+1)+ghost_xhi[1], indgen(dimen[1]-ghost_xhi[2]+1)+ghost_xhi[2]] 

for i = 0, dimen[2]-1 do begin

  out = bspline_iterfit((findgen(dimen[1]))[index], subim[index, i], nord=3, bkspace=40, invvar=1/sqrt(subim[index, i]))
  b = bspline_valu(findgen(dimen[1]), out)
  subim[*,i] = subim[*,i] - b
endfor

writefits, 'backsub.fits', subim, h

obj = subim[source_xlo:source_xhi, *]
params = fltarr(dimen[2], 3)
arr = findgen(source_xhi-source_xlo+1)

for i = 0, dimen[2]-1 do begin
  start = [20., 4, 1000]
  out = mpfitfun('gauss1', arr, subim[source_xlo:source_xhi,i], sqrt(abs(subim[source_xlo:source_xhi,i])) > 1, start)
  params[i,*] = out
endfor

index = indgen(1700)+200

out1 = robust_poly_fit(indgen(1700), params[index,0], 3)
out2 = robust_poly_fit(indgen(1700), params[index,0], 4)

loadct,2
plot,indgen(2227),params[*,0],yr=[20,24] 
oplot,indgen(2227),poly(indgen(2227)-200,out1),col=23
oplot,indgen(2227),poly(indgen(2227)-200,out2),col=100

out3 = robust_poly_fit(indgen(1700), params[index,1], 3)
out4 = robust_poly_fit(indgen(1700), params[index,1], 4)

plot,indgen(2227),params[*,1],yr=[3,5] 
oplot,indgen(2227),poly(indgen(2227)-200,out3),col=23
oplot,indgen(2227),poly(indgen(2227)-200,out4),col=100


params2 = fltarr(dimen[2], 10)
subim2 = subim
for i = 0, dimen[2]-1 do begin
;for i = 1000, 1200 do begin
;  pi = replicate({fixed:0, limited:[1,1], limits:[0.D,0.D]},10)
  pi = replicate({fixed:0, limited:[1,1], limits:[0.D,0.D]},10)

  pi[0].limits = [-30,30]
  pi[1].limits = poly(i-200,out2)+[-2.5,2.5]
  pi[4].limits = poly(i-200,out2)+[-2.5,2.5]
  pi[7].limits = poly(i-200,out2)+[-2.5,2.5]
  pi[2].limits = params[i,1]*[0.4,0.9]
  pi[5].limits = params[i,1]*[0.9,1.5]
  pi[8].limits = params[i,1]*[1.5,3.5]

  pi[3].limited=[1,0]
  pi[6].limited=[1,0]
  pi[9].limited=[1,0]

  pi[3].limits[0] = [0]
  pi[6].limits[0] = [0]
  pi[9].limits[0] = [0]

  pi[3].limited=[0,0]
  pi[6].limited=[0,0]
  pi[9].limited=[0,0]

  start = [0, poly(i-200,out2), params[i,1]*0.5, params[i,2]*0.85, poly(i-200,out2), params[i,1]*2, params[i,2]*0.1, poly(i-200,out2), params[i,1]*5, params[i,2]*0.05]
  out = mpfitfun('psf', arr, subim[source_xlo:source_xhi,i], sqrt(abs(subim[source_xlo:source_xhi,i])) > 1, start, parinfo=pi)
  out = mpfitfun('psf', arr, subim[source_xlo:source_xhi,i], sqrt(abs(subim[source_xlo:source_xhi,i])) > 1, start)
  params2[i,*] = out
  subim2[*,i] = subim[*,i] - psf(findgen(dimen[1])-source_xlo, out)
endfor

writefits, 'sub.fits', subim2, h


out01 = robust_poly_fit(indgen(1700), params2[index,0], 3)
out02 = robust_poly_fit(indgen(1700), params2[index,0], 1)

loadct,2
plot,indgen(2227),params2[*,0],yr=[0,20] 
oplot,indgen(2227),poly(indgen(2227)-200,out01),col=23
oplot,indgen(2227),poly(indgen(2227)-200,out02),col=100


out11 = robust_poly_fit(indgen(1700), params2[index,1], 3)
out12 = robust_poly_fit(indgen(1700), params2[index,1], 4)

plot,indgen(2227),params2[*,1],yr=[15,25] 
oplot,indgen(2227),poly(indgen(2227)-200,out11),col=23
oplot,indgen(2227),poly(indgen(2227)-200,out12),col=100


out21 = robust_poly_fit(indgen(1700), params2[index,2], 3)
out22 = robust_poly_fit(indgen(1700), params2[index,2], 4)

plot,indgen(2227),params2[*,2],yr=[1,3] 
oplot,indgen(2227),poly(indgen(2227)-200,out21),col=23
oplot,indgen(2227),poly(indgen(2227)-200,out22),col=100


out41 = robust_poly_fit(indgen(1700), params2[index,4], 3)
out42 = robust_poly_fit(indgen(1700), params2[index,4], 4)

plot,indgen(2227),params2[*,4],yr=[15,25] 
oplot,indgen(2227),poly(indgen(2227)-200,out41),col=23
oplot,indgen(2227),poly(indgen(2227)-200,out42),col=100


out51 = robust_poly_fit(indgen(1700), params2[index,5], 2)
out52 = robust_poly_fit(indgen(1700), params2[index,5], 3)

plot,indgen(2227),params2[*,5],yr=[2,5] 
oplot,indgen(2227),poly(indgen(2227)-200,out51),col=23
oplot,indgen(2227),poly(indgen(2227)-200,out52),col=100


out71 = robust_poly_fit(indgen(1700), params2[index,7], 3)
out72 = robust_poly_fit(indgen(1700), params2[index,7], 4)

plot,indgen(2227),params2[*,7],yr=[15,25] 
oplot,indgen(2227),poly(indgen(2227)-200,out71),col=23
oplot,indgen(2227),poly(indgen(2227)-200,out72),col=100


out81 = robust_poly_fit(indgen(1700), params2[index,8], 3)
out82 = robust_poly_fit(indgen(1700), params2[index,8], 4)

plot,indgen(2227),params2[*,8],yr=[5,15] 
oplot,indgen(2227),poly(indgen(2227)-200,out81),col=23
oplot,indgen(2227),poly(indgen(2227)-200,out82),col=100




writefits, 'sub.fits', subim2, h


cut = 3
subim2 = subim

for i = 0, dimen[2]-1 do begin
  shift = round(poly(i-200,out11)+cut)

  pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},10)
  pi[[1,2,4,5,7,8]].fixed = 1

  start = [0, poly(i-200,out11), poly(i-200,out21), params2[i,3], poly(i-200,out41), poly(i-200,out51), params2[i,6], poly(i-200,out71), poly(i-200,out82), params2[i,9]]

  fit = mpfitfun('psf', arr[shift:*], subim[source_xlo+shift:source_xhi,i], sqrt(abs(subim[source_xlo+shift:source_xhi,i])) > 1, start, parinfo=pi)
;  params2[i,*] = out
  fit2 = fit
  fit2[1] = fit[1]+source_xlo
  fit2[4] = fit[1]+source_xlo
  fit2[7] = fit[1]+source_xlo
  subim2[*,i] = subim[*,i] - psf(findgen(dimen[1]), fit2)
endfor

writefits, 'sub2.fits', subim2, h


return, 0
END
