% straightforward calculation of the dipolar energy
% of a rotated matrix: rotate it, then calculate dipolar energy
function e = en_dr1(r, th)
  e = en_d0(rot_th(r,th));
end
