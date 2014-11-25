% strightforward way: calculate differences of the matrix
function [e1 e2 e3] = en_gr0(a,b,t, ga,gb,gt, dx)
  r = rmatr_abt(a, b, t);
  gr(:,:,1) = (rmatr_abt(a+ga(1)*dx, b+gb(1)*dx, t+gt(1)*dx) - r)/dx;
  gr(:,:,2) = (rmatr_abt(a+ga(2)*dx, b+gb(2)*dx, t+gt(2)*dx) - r)/dx;
  gr(:,:,3) = (rmatr_abt(a+ga(3)*dx, b+gb(3)*dx, t+gt(3)*dx) - r)/dx;
  e1=0; e2=0; e3=0;
  for k=1:3; for j=1:3; for a=1:3;
    e1 = e1 + gr(a,j,k) * gr(a,j,k);
    e2 = e2 + gr(a,j,k) * gr(a,k,j);
    e3 = e3 + gr(a,j,j) * gr(a,k,k);
  end; end; end
end
