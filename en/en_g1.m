function [e1 e2 e3] = en_gr1(a,b,t, ga,gb,gt)
  % best for calculations

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

  gr = zeros(3,3,3);
  for a=1:3; for j=1:3; for k=1:3; for l=1:3;
    gr(a,j,k) = gr(a,j,k) + ...
      ((1-ct)*(dd(a,l)*n(j) + dd(j,l)*n(a)) - st*ee(a,j,l)) * gn(l,k) + ...
      (st*(n(a)*n(j) - dd(a,j))*dd(j,l) - ct*ee(a,j,l)*n(l)) * gt(k);
  end; end; end; end;

  e1=0; e2=0; e3=0;
  for a=1:3; for j=1:3; for k=1:3;
    e1 = e1 + gr(a,j,k) * gr(a,j,k);
    e2 = e2 + gr(a,j,k) * gr(a,k,j);
    e3 = e3 + gr(a,j,j) * gr(a,k,k);
  end; end; end

%% the same
%  for a=1:3; for j=1:3; for k=1:3; for l=1:3; for m=1:3;
%    e1 = e1 +...
%      (((1-ct)*(dd(a,l)*n(j) + dd(j,l)*n(a)) - st*ee(a,j,l)) * gn(l,k) + ...
%      (st*(n(a)*n(j) - dd(a,j))*dd(j,l) - ct*ee(a,j,l)*n(l)) * gt(k)) * ...
%      (((1-ct)*(dd(a,m)*n(j) + dd(j,m)*n(a)) - st*ee(a,j,m)) * gn(m,k) + ...
%      (st*(n(a)*n(j) - dd(a,j))*dd(j,m) - ct*ee(a,j,m)*n(m)) * gt(k));
%    e2 = e2 +...
%      (((1-ct)*(dd(a,l)*n(j) + dd(j,l)*n(a)) - st*ee(a,j,l)) * gn(l,k) + ...
%      (st*(n(a)*n(j) - dd(a,j))*dd(j,l) - ct*ee(a,j,l)*n(l)) * gt(k)) * ...
%      (((1-ct)*(dd(a,m)*n(k) + dd(k,m)*n(a)) - st*ee(a,k,m)) * gn(m,j) + ...
%      (st*(n(a)*n(k) - dd(a,k))*dd(k,m) - ct*ee(a,k,m)*n(m)) *gt(j));
%    e3 = e3 +...
%      (((1-ct)*(dd(a,l)*n(j) + dd(j,l)*n(a)) - st*ee(a,j,l)) * gn(l,j) + ...
%      (st*(n(a)*n(j) - dd(a,j))*dd(j,l) - ct*ee(a,j,l)*n(l)) * gt(j)) * ...
%      (((1-ct)*(dd(a,m)*n(k) + dd(k,m)*n(a)) - st*ee(a,k,m)) * gn(m,k) + ...
%      (st*(n(a)*n(k) - dd(a,k))*dd(k,m) - ct*ee(a,k,m)*n(m)) *gt(k));
%  end; end; end; end; end;

end
