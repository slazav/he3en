%% Test small rotation rot_th0, rot_th1

function test_rot()
  addpath en
  for i=1:1000;
    a=rand*pi;
    b=rand*pi;
    t=rand*pi;
    th=(rand(1,3)-0.5)*1e-3;

    r = rmatr_abt(a,b,t);
    r0 = rot_th0(r, th);
    r1 = rot_th1(r, th);

    xx(i) = sqrt(sum(th.^2));
    yy(i) = sqrt(sum(sum((r1-r0).^2)));
  end
  figure; hold on;
  plot(xx,yy./xx.^2, '.r');
  xlabel('|ths|');
  ylabel('|difference|/ths^2');
end
