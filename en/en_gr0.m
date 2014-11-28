function [e1 e2 e3] = en_gr0(a,b,t, ga,gb,gt, th, gth)

  % vector n and its gradient
  n(1) = sin(b)*cos(a);
  n(2) = sin(b)*sin(a);
  n(3) = cos(b);
  gn(1,:) = cos(b)*cos(a)*gb - sin(b)*sin(a)*ga;
  gn(2,:) = cos(b)*sin(a)*gb + sin(b)*cos(a)*ga;
  gn(3,:) = -sin(b)*gb;

  ee = zeros(3,3,3);
  ee(1,2,3) = 1;
  ee(2,3,1) = 1;
  ee(3,1,2) = 1;
  ee(3,2,1) = -1;
  ee(2,1,3) = -1;
  ee(1,3,2) = -1;

  dd=[1 0 0; 0 1 0; 0 0 1];

  ct=cos(t);
  st=sin(t);

  % undistorted matrix and its gradients
  r0 = rmatr_abt(a, b, t);
  gr0 = zeros(3,3,3);
  for l=1:3; for k=1:3; for j=1:3; for i=1:3;
    gr0(i,j,k) = gr0(i,j,k) + ...
      ((1-ct)*(dd(i,l)*n(j) + dd(j,l)*n(i)) - st*ee(i,j,l)) * gn(l,k) + ...
      (st*(n(i)*n(j) - dd(i,j))*dd(j,l) - ct*ee(i,j,l)*n(l))*gt(k);
  end; end; end; end;

  et = [    0   th(3) -th(2)
         -th(3)    0   th(1)
          th(2) -th(1)    0];
  tt = th*th'; %'

  % gradient of the distorted matrix
  % \nabla (R(th) R0) = R(th) (\nabla R0) + (\nabla R(th)) R0
  gr1 = zeros(3,3,3);
  for i=1:3; for j=1:3; for k=1:3; for l=1:3;
    gr1(i,k,l) = gr1(i,k,l) +...
      (dd(i,j) - et(i,j) + th(i)*th(j)/2.0 - dd(i,j)*tt/2.0) *gr0(j,k,l) +...
      1/2.0*(th(i)*gth(j,l) + th(j)*gth(i,l)) * r0(j,k);
    for m=1:3;
      gr1(i,k,l) = gr1(i,k,l) - ee(i,j,m)*gth(m,l)*r0(j,k);
      gr1(i,k,l) = gr1(i,k,l) - th(m)*gth(m,l)*dd(i,j)*r0(j,k);
    end
  end;end;end;end;

  e1=0; e2=0; e3=0;
  for k=1:3; for j=1:3; for i=1:3;
    e1 = e1 + gr1(i,j,k) * gr1(i,j,k);
    e2 = e2 + gr1(i,j,k) * gr1(i,k,j);
    e3 = e3 + gr1(i,j,j) * gr1(i,k,k);
  end; end; end
end
