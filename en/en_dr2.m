function e = en_dr2(a,b,t, th)
  n(1) = sin(b)*cos(a);
  n(2) = sin(b)*sin(a);
  n(3) = cos(b);
  ct=cos(t);
  st=sin(t);
  nt=n*th'; tt=th*th';

  e = - 0.5 + 1/2*(4*ct+1)^2 ...
      - 4*(4*ct+1) * st * nt ...
      - (4*ct^2+5*ct+1)*tt ...
      + (9 + 3*ct - 12*ct^2) *nt^2;
end

