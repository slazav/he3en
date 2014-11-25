function r = rmatr_abt(a,b,t)
% convert alpha, beta, theta angles (rad) to the rotation matrix

  n=[sin(b)*cos(a) sin(b)*sin(a) cos(b)];
  r=rmatr_nt(n,t);
end

