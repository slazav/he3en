% dipolar energy
function e = en_d0(r)
  e=0;
  for j=1:3; for k=1:3;
    e = e + r(j,j)*r(k,k) + r(j,k)*r(k,j);
  end; end
end
