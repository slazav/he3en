%% Test dipolar energies en_d0, en_d1, en_d2

function test_d()
  addpath en
  for i=1:1000;
    a=rand*pi;
    b=rand*pi;
    t=rand*pi;
    r = abt2r(a,b,t);
    E0 = en_d0(r);
    E1 = en_d1(r);
    E2 = en_d2(t);
    xx(i) = t;
    yy1(i) = E1-E0;
    yy2(i) = E2-E0;
  end
  figure; hold on;
  plot(xx,yy1, '*r');
  plot(xx,yy2, '*b');
  xlabel('theta');
  ylabel('difference');
end
