function test_g()
  addpath en
  for i=1:1000;
    a=rand*pi;
    b=rand*pi;
    t=rand*pi;

    ga=(rand(1,3)-0.5)*2;
    gb=(rand(1,3)-0.5)*2;
    gt=(rand(1,3)-0.5)*2;

    dx(i) = 1e-8 + (rand)*1e-6;

    [e1a e2a e3a] = en_g0(a,b,t,ga,gb,gt, dx(i));
    [e1b e2b e3b] = en_g1(a,b,t,ga,gb,gt);
    d1(i) = e1b-e1a;
    d2(i) = e2b-e2a;
    d3(i) = e3b-e3a;
  end

  figure; hold on;
  plot(dx, d1, '.r');
  plot(dx, d2, '.b');
  plot(dx, d3, '.g');
  xlabel('dx')
  ylabel('difference')
end
