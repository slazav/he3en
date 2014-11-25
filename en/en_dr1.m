function e = en_dr1(r, th)
  e=0;
  et = [    0   th(3) -th(2)
         -th(3)    0   th(1)
          th(2) -th(1)    0];
  dd = [1 0 0
        0 1 0
        0 0 1];
  tt = th*th'; %'

  for j=1:3; for k=1:3;
    e = e + (r(j,j)*r(k,k) + r(j,k)*r(k,j))  *  (1-tt);
    for a=1:3;
      e = e + (r(a,j)*r(k,k) + r(a,k)*r(k,j))  *  (th(j)*th(a) - 2*et(j,a));
      for b=1:3;
         e = e + (r(a,j)*r(b,k) + r(a,k)*r(b,j))  *  et(j,a)*et(k,b);
      end;
    end;
  end; end
end
