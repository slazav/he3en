function r = rmatr_nt(n,t)
% convert n and theta to the rotation matrix

  nn = [ n(1)^2    n(2)*n(1) n(3)*n(1)
         n(1)*n(2) n(2)^2    n(3)*n(2)
         n(1)*n(3) n(2)*n(3) n(3)^2 ];

  en = [    0   n(3) -n(2)
         -n(3)    0   n(1)
          n(2) -n(1)    0];

  dd = [1 0 0; 0 1 0; 0 0 1];

  r = dd*cos(t) + nn*(1-cos(t)) - en*sin(t);
end

