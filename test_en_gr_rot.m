%% test a change of gradient energy after small rotation

function test_en_gr_rot()

  for i=1:1000;
    a=rand*pi;
    b=rand*pi;
    t=rand*pi;

    ga=(rand(1,3)-0.5)*2;
    gb=(rand(1,3)-0.5)*2;
    gt=(rand(1,3)-0.5)*2;

    th=(rand(1,3)-0.5)*1e-2;
%    th=[0 0 0.01];
    gth=(rand(3,3)-0.5)*1e-2;

    [e1a e2a e3a] = en_gr1(a,b,t,ga,gb,gt, th,gth);
    [e1b e2b e3b] = en_gr2(a,b,t,ga,gb,gt, th,gth);
%    if abs(e1a-e1b)>1e-12 error('DIFF1 %e %e  %e\n', e1a, e1b, abs(e1a-e1b)); end
%    if abs(e2a-e2b)>1e-12 error('DIFF1 %e %e  %e\n', e2a, e2b, abs(e2a-e2b)); end
%    if abs(e3a-e3b)>1e-12 error('DIFF1 %e %e  %e\n', e3a, e3b, abs(e3a-e3b)); end

    aa(i) = max(abs(th));
    gg(i) = max(max(abs(gth)));
    bb1(i) = abs(e1a-e1b);
    bb2(i) = abs(e2a-e2b);
    bb3(i) = abs(e3a-e3b);
  end
  figure; hold on;
  plot(aa,bb1./aa,'*r');
  plot(gg,bb1./gg,'*b');

  plot(aa,bb2./aa,'xr');
  plot(gg,bb2/gg,'xb');
  
  plot(aa,bb3./aa,'or');
  plot(gg,bb3/gg,'ob');

end

%%%
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

  r0 = abt2r(a, b, t);
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

  gr1 = gr0;
  for i=1:3; for j=1:3; for k=1:3; for l=1:3;
          gr1(i,k,l) = gr1(i,k,l) +...
          (-et(i,j) + th(i)*th(j)/2.0 - dd(i,j)*tt/2.0) *gr0(j,k,l);
          gr1(i,k,l) = gr1(i,k,l) +...
              1/2.0*(th(i)*gth(j,l) + th(j)*gth(i,l))*r0(j,k);
    for m=1:3;
      gr1(i,k,l) = gr1(i,k,l) - ee(i,j,m)*gth(m,l)*r0(j,k);
      gr1(i,k,l) = gr1(i,k,l) - th(m)*gth(m,l)*dd(i,j)*r0(j,k);
    end
  end;end;end;end;

  e1=0; e2=0; e3=0;
  for k=1:3; for j=1:3; for i=1:3;
    e1 = e1 + gr1(i,j,k) * gr1(i,j,k) - gr0(i,j,k) * gr0(i,j,k);
    e2 = e2 + gr1(i,j,k) * gr1(i,k,j) - gr0(i,j,k) * gr0(i,k,j);
    e3 = e3 + gr1(i,j,j) * gr1(i,k,k) - gr0(i,j,j) * gr0(i,k,k);
  end; end; end
end

function [e1 e2 e3] = en_gr2(a,b,t, ga,gb,gt, th, gth)
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

  r0 = abt2r(a, b, t);
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
