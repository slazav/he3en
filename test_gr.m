%% test a change of dipolar energy after a small rotation th
function test_gr()
  addpath en
  for i=1:1000;
    a=rand*pi;
    b=rand*pi;
    t=rand*pi;

    ga=(rand(1,3)-0.5)*2;
    gb=(rand(1,3)-0.5)*2;
    gt=(rand(1,3)-0.5)*2;

    th=(rand(1,3)-0.5)*1e-3;
    th(find(abs(th)<1e-6))=1e-6;
    gth=(rand(3,3)-0.5)*1e-3;
    gth(find(abs(th)<1e-6))=1e-6;

    [e1a e2a e3a] = en_gr0(a,b,t, ga,gb,gt, th,gth);
    [e1b e2b e3b] = en_gr1(a,b,t, ga,gb,gt, th,gth);

    x1(i) = sqrt(sum(th.^2) );
    x2(i) = sqrt(sum(sum(gth.^2)));

    y1(i) = e1b-e1a;
    y2(i) = e2b-e2a;
    y3(i) = e3b-e3a;

  end
  figure; hold on;
  plot(x1,y1./x1.^2,'.r');
  plot(x1,y2./x1.^2,'.g');
  plot(x1,y3./x1.^2,'.b');

  plot(x2,y1./x2.^3,'or');
  plot(x2,y2./x2.^3,'og');
  plot(x2,y3./x2.^3,'ob');

  xlabel('|ths|');
  ylabel('difference/|ths|^3');
end
