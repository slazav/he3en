%% test a change of dipolar energy after a small rotation th
function test_dr()
  addpath en
  for i=1:1000;
    a=rand*pi;
    b=rand*pi;
    t=rand*pi;
    th=(rand(1,3)-0.5)*1e-3;
    th(find(abs(th)<1e-6))=1e-6;

    r = rmatr_abt(a,b,t);
    E0 = en_d0(rot_th0(r,th));
    E1 = en_dr1(r,th);
    E2 = en_dr2(a,b,t,th);
    xx(i) = sqrt(sum(th.^2));
    yy1(i) = E1-E0;
    yy2(i) = E1-E2;
  end
  figure; hold on;
  plot(xx,yy1./xx.^3,'.r');
  plot(xx,yy2./xx.^3,'.b');
  xlabel('|ths|');
  ylabel('difference/|ths|^3');
end
