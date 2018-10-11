% Take 1D spin current from test_torque.m
% Check its derivatives on n, theta, n', theta'
% (it is used in vmcw boundary conditions)

function test_scurr()

% Spin currents as a function of (n, th, Gn, Gth) in 1D.
% Taken from test_torque.m
function j = JA(vn, t, gn, gt)
  st = sin(t);
  ct = cos(t);
  j(1) = ...
     - 2*vn(1)*gt ...
     - 2*st*gn(1) ...
     - 2*(1-ct) *(vn(2)*gn(3) - vn(3)*gn(2)) ...
  ;
  j(2) = ...
     - 2*vn(2)*gt ...
     - 2*st*gn(2) ...
     - 2*(1-ct) *(vn(3)*gn(1) - vn(1)*gn(3)) ...
  ;
  j(3) = ...
     - 2*vn(3)*gt ...
     - 2*st*gn(3) ...
     - 2*(1-ct) *(vn(1)*gn(2)-vn(2)*gn(1))...
  ;
end

% dJa/dn
function dj = DJAn(vn, t, gn, gt)
  st = sin(t);
  ct = cos(t);
  dj(1,1) = - 2*gt;
  dj(1,2) = - 2*(1-ct)*gn(3);
  dj(1,3) = + 2*(1-ct)*gn(2);

  dj(2,1) = + 2*(1-ct)*gn(3);
  dj(2,2) = - 2*gt;
  dj(2,3) = - 2*(1-ct)*gn(1);

  dj(3,1) = - 2*(1-ct)*gn(2);
  dj(3,2) = + 2*(1-ct)*gn(1);
  dj(3,3) = - 2*gt;
end

% dJa/dgn
function dj = DJAgn(vn, t, gn, gt)
  st = sin(t);
  ct = cos(t);
  dj(1,1) = - 2*st;
  dj(1,2) = + 2*(1-ct)*vn(3);
  dj(1,3) = - 2*(1-ct)*vn(2);
  dj(2,1) = - 2*(1-ct)*vn(3);
  dj(2,2) = - 2*st;
  dj(2,3) = + 2*(1-ct)*vn(1);
  dj(3,1) = + 2*(1-ct)*vn(2);
  dj(3,2) = - 2*(1-ct)*vn(1);
  dj(3,3) = - 2*st;
end

function dj = DJAgt(vn, t, gn, gt)
  dj = - 2*vn;
end

function j = JB(vn, t, gn, gt)
  st = sin(t);
  ct = cos(t);
  j(1) = ...
   - (1 - (1-ct)*vn(3)^2) * vn(1)*gt ...
   + st * vn(2)*vn(3)*gt ...
   -     st*(ct + (1-ct)*vn(3)^2) *gn(1) ...
   + (1-ct)*(ct + (1-ct)*vn(3)^2) *vn(3)*gn(2) ...
   + (1-ct)*(ct - (1-ct)*vn(3)^2) *vn(2)*gn(3) ...
   + 2*st*(1-ct)*vn(1)*vn(3)*gn(3) ...
   ;

  j(2) = ...
   - (1 - (1-ct)*vn(3)^2) * vn(2)*gt ...
   - st * vn(1)*vn(3)*gt ...
   - (1-ct)*(ct+(1-ct)*vn(3)^2) *vn(3)*gn(1) ...
   -     st*(ct+(1-ct)*vn(3)^2) *gn(2) ...
   - (1-ct)*(ct-(1-ct)*vn(3)^2) *vn(1)*gn(3) ...
   + 2*st*(1-ct)*vn(2)*vn(3)*gn(3) ...
   ;

  j(3) = ...
   - (1-ct)*(1 - vn(3)^2) * vn(3)*gt ...
   + (st*st + (1-ct)^2 *vn(3)^2) *vn(2)*gn(1) ...
   - (st*st + (1-ct)^2 *vn(3)^2) *vn(1)*gn(2) ...
   - st*(1-ct) * (1-vn(3)^2) *gn(3) ...
   ;
end

function dj = DJBn(vn, t, gn, gt)
  st = sin(t);
  ct = cos(t);
  dj(1,1) = - (1 - (1-ct)*vn(3)^2)*gt + 2*st*(1-ct)*vn(3)*gn(3);
  dj(1,2) = + st*vn(3)*gt + (1-ct)*(ct - (1-ct)*vn(3)^2)*gn(3);
  dj(1,3) = ...
   + 2*(1-ct)*vn(3)*vn(1)*gt ...
   + st * vn(2)*gt ...
   +     2*st*(1-ct)*(vn(1)*gn(3)-vn(3)*gn(1)) ...
   + 2*(1-ct)*(1-ct)*vn(3)*(vn(3)*gn(2)-vn(2)*gn(3)) ...
   + (1-ct)*(ct + (1-ct)*vn(3)^2)*gn(2) ...
   ;
  dj(2,1) = - st*vn(3)*gt - (1-ct)*(ct-(1-ct)*vn(3)^2)*gn(3);
  dj(2,2) = - (1 - (1-ct)*vn(3)^2)*gt + 2*st*(1-ct)*vn(3)*gn(3);
  dj(2,3) = ...
   + 2*(1-ct)*vn(3)*vn(2)*gt ...
   - st*vn(1)*gt ...
   - (1-ct)*(ct+(1-ct)*vn(3)^2) *gn(1) ...
   +     2*st*(1-ct)*(vn(2)*gn(3)-vn(3)*gn(2)) ...
   + 2*(1-ct)*(1-ct)*vn(3) *(vn(1)*gn(3)-vn(3)*gn(1)) ...
   ;

  dj(3,1) = - (st*st + (1-ct)^2 *vn(3)^2)*gn(2);
  dj(3,2) = + (st*st + (1-ct)^2 *vn(3)^2)*gn(1);
  dj(3,3) = ...
   - (1-ct)*gt ...
   + 3*(1-ct)*vn(3)^2*gt ...
   + 2*(1-ct)^2 *vn(3) *(vn(2)*gn(1)-vn(1)*gn(2)) ...
   + 2*st*(1-ct)*vn(3)*gn(3) ...
  ;
end

function dj = DJBgn(vn, t, gn, gt)
  st = sin(t);
  ct = cos(t);
  dj(1,1) = - st*(ct + (1-ct)*vn(3)^2);
  dj(1,2) = + (1-ct)*(ct + (1-ct)*vn(3)^2)*vn(3);
  dj(1,3) = + (1-ct)*(ct - (1-ct)*vn(3)^2)*vn(2) ...
            + 2*st*(1-ct)*vn(1)*vn(3);
  dj(2,1) = - (1-ct)*(ct+(1-ct)*vn(3)^2) *vn(3);
  dj(2,2) = - st*(ct+(1-ct)*vn(3)^2);
  dj(2,3) = - (1-ct)*(ct-(1-ct)*vn(3)^2) *vn(1) ...
            + 2*st*(1-ct)*vn(2)*vn(3);
  dj(3,1) = + (st*st + (1-ct)^2 *vn(3)^2) *vn(2);
  dj(3,2) = - (st*st + (1-ct)^2 *vn(3)^2) *vn(1);
  dj(3,3) = - st*(1-ct) * (1-vn(3)^2);
end

function dj = DJBgt(vn, t, gn, gt)
  st = sin(t);
  ct = cos(t);
  dj(1) = - (1 - (1-ct)*vn(3)^2) * vn(1) ...
          + st * vn(2)*vn(3);
  dj(2) = - (1 - (1-ct)*vn(3)^2) * vn(2) ...
          - st * vn(1)*vn(3);
  dj(3) = - (1-ct)*(1 - vn(3)^2) * vn(3);
end


dx=1e-6;

# two values for n*th
nt1 = 2*pi*[rand, rand, rand] - pi;
gnt = (2*pi*[rand, rand, rand] - pi);
nt2 = nt1+gnt*dx;

t1 = sqrt(sum(nt1.^2));
n1 = nt1/t1;

t2 = sqrt(sum(nt2.^2));
n2 = nt2/t2;

gt = sum(nt1.*gnt)/t1;
gn = gnt/t1 - nt1*gt/t1^2;


r1 = (gt-(t2-t1)./dx)/gt;
r2 = (gn-(n2-n1)./dx);
r2 = sqrt(sum(r2.^2))/sqrt(sum(gn.^2));
fprintf ("Test grad: %e %e\n", r1, r2);

%%%%%%%%%%%

dd = ([rand, rand, rand] - 0.5)*dx;
Ja1 = JA(n1, t1, gn, gt);
Jb1 = JB(n1, t1, gn, gt);

Ja2 = JA(n1+dd, t1, gn, gt);
Jb2 = JB(n1+dd, t1, gn, gt);

d1 = DJAn(n1,t1,gn,gt);
r1 = Ja2-Ja1-[sum(d1(1,:).*dd) sum(d1(2,:).*dd) sum(d1(3,:).*dd)];
r1 = sqrt(sum(r1.^2))/sqrt(sum((Ja2-Ja1).^2));

d2 = DJBn(n1,t1,gn,gt);
r2 = Jb2-Jb1-[sum(d2(1,:).*dd) sum(d2(2,:).*dd) sum(d2(3,:).*dd)];
r2 = sqrt(sum(r2.^2))/sqrt(sum((Jb2-Jb1).^2));

fprintf ("Test DJAn:  %e %e\n", r1, r2);

%%%%%%%%%%%

dd = ([rand, rand, rand] - 0.5)*dx;
Ja1 = JA(n1, t1, gn, gt);
Jb1 = JB(n1, t1, gn, gt);

Ja2 = JA(n1, t1, gn+dd, gt);
Jb2 = JB(n1, t1, gn+dd, gt);

d1 = DJAgn(n1,t1,gn,gt);
r1 = Ja2-Ja1-[sum(d1(1,:).*dd) sum(d1(2,:).*dd) sum(d1(3,:).*dd)];
r1 = sqrt(sum(r1.^2))/sqrt(sum((Ja2-Ja1).^2));

d2 = DJBgn(n1,t1,gn,gt);
r2 = Jb2-Jb1-[sum(d2(1,:).*dd) sum(d2(2,:).*dd) sum(d2(3,:).*dd)];
r2 = sqrt(sum(r2.^2))/sqrt(sum((Jb2-Jb1).^2));


fprintf ("Test DJAgn: %e %e\n", r1, r2);

%%%%%%%%%%%

dd = (rand - 0.5)*dx;
Ja1 = JA(n1, t1, gn, gt);
Jb1 = JB(n1, t1, gn, gt);

Ja2 = JA(n1, t1, gn, gt+dd);
Jb2 = JB(n1, t1, gn, gt+dd);

r1 = Ja2-Ja1-DJAgt(n1,t1,gn,gt)*dd;
r1 = r1/(Ja2-Ja1);

r2 = Jb2-Jb1-DJBgt(n1,t1,gn,gt)*dd;
r2 = r2/(Jb2-Jb1);

fprintf ("Test DJAgt: %e %e\n", r1, r2);

%%%%%%%%%%%


end