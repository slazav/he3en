function r1 = rot_th0(r, th)
% rotate matrix r by angle th
  t=sqrt(sum(th.^2)); % angle
  n=th/t;             % axes
  r1 = rmatr_nt(n,t) * r;
end
