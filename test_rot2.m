%% Test double rotations
function test_rot2()
  addpath en
  for i=1:1000;
    a=rand*pi;
    b=rand*pi;
    t=rand*pi;
    th1=(rand(1,3)-0.5)*1e-3;
    th2=(rand(1,3)-0.5)*1e-3;

    tt = [ th1(2)*th2(3)-th1(3)*th2(2)
           th1(3)*th2(1)-th1(1)*th2(3)
           th1(1)*th2(2)-th1(2)*th2(1)]';

    r = rmatr_abt(a,b,t);
    r0 = rot_th0(rot_th0(r, th1),th2);
    r1 = rot_th0(r, th1+th2-tt/2);

    xx(i) = sqrt(sum(th1.^2) + sum(th2.^2));
    yy(i) = sqrt(sum(sum((r1-r0).^2)));
  end
  figure; hold on;
  plot(xx,yy./xx.^3, '.r');
  xlabel('|ths|');
  ylabel('|difference|/ths^3');
end
