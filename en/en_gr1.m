function [e1 e2 e3] = en_gr1(a,b,t, ga,gb,gt, th, gth)
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

  r0 = rmatr_abt(a, b, t);
  gr0 = zeros(3,3,3);
  for l=1:3; for k=1:3; for j=1:3; for i=1:3;
    gr0(i,j,k) = gr0(i,j,k) + ...
      ((1-ct)*(dd(i,l)*n(j) + dd(j,l)*n(i)) - st*ee(i,j,l)) * gn(l,k) + ...
      (st*(n(i)*n(j) - dd(i,j))*dd(j,l) - ct*ee(i,j,l)*n(l))*gt(k);
  end; end; end; end;

  e1=0; e2=0; e3=0;
%  e1a=0; e2a=0; e3a=0;
  for k=1:3; for j=1:3;
    for c=1:3;
%     e1 = e1 + gr0(c,j,k) * gr0(c,j,k);
%     e2 = e2 + gr0(c,j,k) * gr0(c,k,j);
%     e3 = e3 + gr0(c,j,j) * gr0(c,k,k);
    end
    for a=1:3; for b=1:3; for c=1:3;
%      xx = 2*ee(a,b,c)* gth(c,j);
      xx = 2*ee(a,b,c)* gth(c,j) +  dd(c,1)*th(b) * gth(a,j) - dd(c,1)*th(a)*gth(b,j);
      e1 = e1 + xx*r0(a,k)*gr0(b,k,j);
      e2 = e2 + xx*r0(a,k)*gr0(b,j,k);
      e3 = e3 + xx*r0(a,j)*gr0(b,k,k);
    end; end; end;
    for a=1:3; for b=1:3;
      e2 = e2 - r0(a,k)*r0(b,j)*gth(b,j)*gth(a,k);
      e3 = e3 - r0(a,j)*r0(b,k)*gth(b,j)*gth(a,k);
    end end
  end; end;
  for j=1:3; for a=1:3;
      e1 = e1 + 2*gth(a,j)*gth(a,j);
      e2 = e2 + gth(a,j)*gth(a,j);
      e3 = e3 + gth(a,j)*gth(a,j);
  end; end

end
