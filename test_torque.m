% test gradient torque
% test for a single random matrix with gradient ~d
function test_torque()
  addpath en
  d = 0.00012

  % standard matrices
  ee = zeros(3,3,3);
  ee(1,2,3) = 1;
  ee(2,3,1) = 1;
  ee(3,1,2) = 1;
  ee(3,2,1) = -1;
  ee(2,1,3) = -1;
  ee(1,3,2) = -1;

  dd=[1 0 0; 0 1 0; 0 0 1];

  % 3x3 coordinate grid
  for i=1:3; for j=1:3, for k=1:3
    xx(i,j,k) = (i-2)*d;
    yy(i,j,k) = (j-2)*d;
    zz(i,j,k) = (k-2)*d;
  end; end; end

  ma = rand*pi + xx*(rand-0.5) + yy*(rand-0.5) + zz*(rand-0.5) ...
     + xx.*xx*2*(rand-0.5) + yy.*yy*2*(rand-0.5) + zz.*zz*2*(rand-0.5) ...
     + xx.*yy*(rand-0.5) + yy.*zz*(rand-0.5) + zz.*xx*(rand-0.5);

  mb = rand*pi + xx*(rand-0.5) + yy*(rand-0.5) + zz*(rand-0.5) ...
     + xx.*xx*2*(rand-0.5) + yy.*yy*2*(rand-0.5) + zz.*zz*2*(rand-0.5) ...
     + xx.*yy*(rand-0.5) + yy.*zz*(rand-0.5) + zz.*xx*(rand-0.5);

  mt = rand*pi + xx*(rand-0.5) + yy*(rand-0.5) + zz*(rand-0.5) ...
     + xx.*xx*2*(rand-0.5) + yy.*yy*2*(rand-0.5) + zz.*zz*2*(rand-0.5) ...
     + xx.*yy*(rand-0.5) + yy.*zz*(rand-0.5) + zz.*xx*(rand-0.5);

  % n vector and R matrix
  mn = zeros(3,3,3,3);
  mR = zeros(3,3,3,3,3);
  for i=1:3; for j=1:3, for k=1:3
    mn(i,j,k,:)=[sin(mb(i,j,k))*cos(ma(i,j,k)) sin(mb(i,j,k))*sin(ma(i,j,k)) cos(mb(i,j,k))];
    mR(i,j,k,:,:) = rmatr_nt(mn(i,j,k,:),mt(i,j,k));
  end; end; end

  % gradients: ggR(k,m,a,j) == nabla_k nabla_m R_aj
  for a=1:3; for j=1:3
    [vR(a,j), gR(:,a,j), ggR(:,:,a,j)] = vmatr2grad(squeeze(mR(:,:,:,a,j)), d);
  end; end
  for a=1:3;
    [vn(a), gn(:,a), ggn(:,:,a)] = vmatr2grad(squeeze(mn(:,:,:,a)), d);
  end
  [vt, gt, ggt] = vmatr2grad(mt, d);
  ct = cos(vt); st=sin(vt);


  %%%%%%%%%%%%%%%%%%%
  %test 1 -- R matrix
  testR = zeros(3,3);
  for a=1:3; for j=1:3;
    testR(a,j) = testR(a,j) + ct*dd(a,j) + (1-ct)*vn(a)*vn(j);
    for k=1:3
      testR(a,j) = testR(a,j) - st*ee(a,j,k)*vn(k);
  end; end; end
  r = sum(sum((testR-vR).^2));
  r0 = sum(sum(testR.^2));
  fprintf('test   R: %e / %e\n', r, r0);

  %%%%%%%%%%%%%%%%%%%
  %test 2 -- gR
  testgR = zeros(3,3,3);
  for k=1:3; for a=1:3; for j=1:3;
    testgR(k,a,j) = testgR(k,a,j) ...
       + st*(vn(a)*vn(j)-dd(a,j))*gt(k) ...
       + (1-ct)*(vn(j)*gn(k,a) + vn(a)*gn(k,j));

    for m=1:3
      testgR(k,a,j) = testgR(k,a,j) ...
       - ct*ee(a,j,m)*vn(m) * gt(k) ...
       - st*ee(a,j,m)*gn(k,m);
    end
  end; end; end

  r  = sum(sum(sum((testgR-gR).^2)));
  r0 = sum(sum(sum((testgR).^2)));
  fprintf('test  gR: %e / %e\n', r, r0);

  %%%%%%%%%%%%%%%%%%%
  %test 3 -- ggR
  testggR = zeros(3,3,3,3);
   for l=1:3; for k=1:3; for a=1:3; for j=1:3;
    testggR(l,k,a,j) = testggR(l,k,a,j) ...
       + ct*(vn(a)*vn(j)-dd(a,j))*gt(k)*gt(l) ...
       + st*(vn(a)*gn(l,j)+vn(j)*gn(l,a))*gt(k) ...
       + st*(vn(a)*vn(j)-dd(a,j))*ggt(k,l) ...
       + (1-ct)*(gn(l,j)*gn(k,a) + gn(l,a)*gn(k,j)) ...
       + (1-ct)*(vn(j)*ggn(k,l,a) + vn(a)*ggn(k,l,j)) ...
       + st*(vn(j)*gn(k,a) + vn(a)*gn(k,j))*gt(l);
    for m=1:3
      testggR(l,k,a,j) = testggR(l,k,a,j) ...
        + st * ee(a,j,m) * vn(m) * gt(l)*gt(k) ...
        - ct * ee(a,j,m) * gn(l,m) * gt(k) ...
        - ct * ee(a,j,m) * vn(m) * ggt(l,k) ...
        - ct * ee(a,j,m) * gn(k,m) * gt(l) ...
        - st * ee(a,j,m) * ggn(l,k,m);
    end
  end; end; end; end

  r = sum(sum(sum(sum((testggR-ggR).^2))));
  r0 = sum(sum(sum(sum(testggR.^2))));

  fprintf('test ggR: %e / %e\n', r, r0);


  %%%%%%%%%%%%%%%%%%%
  % torque-1
  T1a=zeros(3,1);
  T1b=zeros(3,1);
  for a=1:3; for b=1:3; for c=1:3; for j=1:3; for k=1:3;
    T1a(a) = T1a(a) + ee(a,b,c)*vR(c,j)*ggR(k,k,b,j);
    T1b(a) = T1b(a) + ee(a,b,c)*vR(c,j)*ggR(k,j,b,k);
  end; end; end; end; end;

  % torque-2
  T2a=zeros(3,1);
  T2b=zeros(3,1);
  for a=1:3; for b=1:3; for c=1:3; for j=1:3; for k=1:3; for m=1:3;
    T2a(b) = T2a(b) + ee(b,a,c)*vR(c,j)* (
       + dd(m,1) * ct*(vn(a)*vn(j)-dd(a,j))*gt(k)*gt(k) ...
       + dd(m,1) * st*(vn(a)*gn(k,j)+vn(j)*gn(k,a))*gt(k) ...
       + dd(m,1) * st*(vn(a)*vn(j)-dd(a,j))*ggt(k,k) ...
       + dd(m,1) * (1-ct)*(gn(k,j)*gn(k,a) + gn(k,a)*gn(k,j)) ...
       + dd(m,1) * (1-ct)*(vn(j)*ggn(k,k,a) + vn(a)*ggn(k,k,j)) ...
       + dd(m,1) * st*(vn(j)*gn(k,a) + vn(a)*gn(k,j))*gt(k) ...
       + st * ee(a,j,m) * vn(m) * gt(k)*gt(k) ...
       - ct * ee(a,j,m) * gn(k,m) * gt(k) ...
       - ct * ee(a,j,m) * vn(m) * ggt(k,k) ...
       - ct * ee(a,j,m) * gn(k,m) * gt(k) ...
       - st * ee(a,j,m) * ggn(k,k,m) ...
    );
    T2b(b) = T2b(b) + ee(b,a,c)*vR(c,j)* (
       + dd(m,1) * ct*(vn(a)*vn(k)-dd(a,k))*gt(j)*gt(k) ...
       + dd(m,1) * st*(vn(a)*gn(k,k)+vn(k)*gn(k,a))*gt(j) ...
       + dd(m,1) * st*(vn(a)*vn(k)-dd(a,k))*ggt(j,k) ...
       + dd(m,1) * (1-ct)*(gn(k,k)*gn(j,a) + gn(k,a)*gn(j,k)) ...
       + dd(m,1) * (1-ct)*(vn(k)*ggn(j,k,a) + vn(a)*ggn(j,k,k)) ...
       + dd(m,1) * st*(vn(k)*gn(j,a) + vn(a)*gn(j,k))*gt(k) ...
       + st * ee(a,k,m) * vn(m) * gt(k)*gt(j) ...
       - ct * ee(a,k,m) * gn(k,m) * gt(j) ...
       - ct * ee(a,k,m) * vn(m) * ggt(k,j) ...
       - ct * ee(a,k,m) * gn(j,m) * gt(k) ...
       - st * ee(a,k,m) * ggn(k,j,m) ...
    );
  end; end; end; end; end; end;

  r1 = sum((T1a-T2a).^2);
  r2 = sum((T1b-T2b).^2);
  fprintf('test torque 1. A: %e B: %e\n', r1,r2);

  % torque-3
  T2a=zeros(3,1);
  T2b=zeros(3,1);
  for a=1:3; for b=1:3; for j=1:3; for k=1:3; for m=1:3;
    T2a(b) = T2a(b) ...
... % gn gt
       - 2*dd(a,1)*dd(m,1)*dd(j,1) *(1+ct)  * gn(k,b)*gt(k) ...
       + 2*dd(a,1)                 *st   *ee(b,j,m) *vn(m)*gn(k,j)*gt(k) ...
... % ggn
       - dd(a,1)*dd(m,1)*dd(j,1) *2*st        *ggn(k,k,b) ...
       - dd(a,1)                 *2*(1-ct)    *ee(m,j,b) * vn(m)*ggn(k,k,j) ...
       + dd(a,1)*dd(m,1)         *2*st*(1-ct) *vn(b)*vn(j)*ggn(k,k,j) ...
... % ggt
       - dd(a,1)*dd(m,1)*dd(j,1) *2*vn(b)     *ggt(k,k) ...
... % gn gn
       + dd(a,1)*dd(m,1)         *2*(1-ct)*st *vn(b)*gn(k,j)^2;



    T2b(b) = T2b(b) ...
...% gt gt
       - dd(a,1)*dd(m,1)*dd(k,1) *st  *vn(j) * gt(b)*gt(j) ...
       + dd(a,1)*dd(m,1)         *st  *vn(b)*vn(j)*vn(k) *gt(k)*gt(j) ...
       + dd(m,1)                 *ct  *ee(b,j,k)    *vn(j)*vn(a) *gt(k)*gt(a) ...
...% gn gn
       + dd(j,1)*dd(m,1) *vn(b) *  st*(1-ct)*(gn(k,k)*gn(a,a) + gn(k,a)*gn(a,k)) ...
       + dd(m,1) *ct*(1-ct) *ee(b,a,j)*gn(k,k)*gn(j,a) ...
       + dd(m,1) *ct*(1-ct) *ee(b,a,j)*gn(k,a)*gn(j,k) ...
       + (1-ct)^2  *ee(b,a,m)*vn(m)*vn(j) *gn(k,k)*gn(j,a) ...
       + (1-ct)^2  *ee(b,a,m)*vn(m)*vn(j) *gn(k,a)*gn(j,k) ...
...% gn gt
       - dd(a,1)*dd(m,1)*dd(k,1) *2*ct^2      * gn(j,b)*gt(j) ...
       + dd(a,1)*dd(m,1)*dd(k,1) *(ct^2-st^2) * gn(b,j)*gt(j) ...
       + dd(a,1)*dd(m,1)*dd(k,1) *(ct^2-st^2) * gn(j,j)*gt(b) ...
       - dd(a,1)*dd(m,1) *(1-ct)*ct* vn(k)*vn(j)  *2*gn(k,b) * gt(j) ...
       + dd(a,1)*dd(k,1) *2*st*st* vn(b) * (vn(m)*gn(j,j) + vn(j)*gn(j,m))*gt(m) ...
       + dd(m,1)  *ct*st  *ee(b,a,j) *vn(a) * gn(k,k) * gt(j) ...
       + dd(m,1)  *ct*st  *ee(b,a,j) *vn(k) * gn(k,a) * gt(j) ...
       + dd(m,1)  *ct*st  *ee(b,a,j) *vn(k) * gn(j,a) * gt(k) ...
       + dd(m,1)  *ct*st  *ee(b,a,j) *vn(a) * gn(j,k) * gt(k) ...
       + dd(m,1)  *ct*st  *ee(a,k,j) *vn(a) * gn(b,j) * gt(k) ...
       + dd(j,1)  *ct*st  *ee(a,k,m) *vn(a) * gn(k,m) * gt(b) ...
       - dd(j,1)  *st*ct  *ee(a,k,m) *vn(b) * gn(k,m) * gt(a) ...
       - dd(j,1)  *st*ct  *ee(a,k,m) *vn(b) * gn(a,m) * gt(k) ...
       + 2*(1-ct)*st *ee(b,a,m)*vn(m)*vn(j)*vn(k)*gn(k,a)*gt(j) ...
...% ggt
       + dd(a,1)*dd(m,1)*dd(k,1) *ct *vn(j)*ggt(b,j) ...
       - dd(a,1)*dd(m,1)*dd(k,1)     *vn(b)*ggt(j,j) ...
       + dd(a,1)*dd(m,1) *(1-ct)    *vn(b)*vn(k)*vn(j) * ggt(k,j) ...
       + dd(k,1)  *st     *ee(b,a,j)*vn(a)*vn(m)*ggt(j,m) ...
       - dd(k,1)  *st*ct  *ee(a,j,m)*vn(b)*vn(m) * ggt(j,a) ...
...% ggn
       - dd(a,1)*dd(m,1)*dd(k,1) *ct * st     *ggn(j,j,b) ...
       + dd(a,1)*dd(m,1)*dd(k,1) *(2*ct-1)*st *ggn(b,j,j) ...
       - dd(a,1)*dd(m,1) *(1-ct)*st *vn(k)*vn(j) * ggn(k,j,b) ...
       + dd(j,1)*dd(m,1) *2*(1-ct)*st *vn(b)*vn(k) * ggn(a,k,a) ...
       + dd(m,1) *ct*(1-ct) *ee(b,a,j) *vn(k) * ggn(j,k,a) ...
       + dd(m,1) *ct*(1-ct) *ee(b,a,j) *vn(a) * ggn(j,k,k) ...
       + dd(m,1) *st*st     *ee(a,k,j) *vn(a) * ggn(k,b,j) ...
       - dd(m,1) *st*st     *ee(a,k,j) *vn(b) * ggn(k,a,j) ...
       + (1-ct)^2* ee(b,a,m)*vn(m)*vn(j)*vn(k)*ggn(j,k,a);

  end; end; end; end; end;

  r1 = sum((T1a-T2a).^2);
  r2 = sum((T1b-T2b).^2);
  fprintf('test torque 2. A: %e B: %e\n', r1,r2);

  %%%%%%%%%%%%%%%%%%%
  % spin current 1
  J1a=zeros(3,3);
  J1b=zeros(3,3);
  for a=1:3; for b=1:3; for c=1:3; for j=1:3; for k=1:3;
    J1a(a,k) = J1a(a,k) + ee(a,b,c)*vR(c,j)*gR(k,b,j);
    J1b(a,k) = J1b(a,k) + ee(a,b,c)*vR(c,j)*gR(j,b,k);
  end; end; end; end; end;

  %%%%%%%%%%%%%%%%%%%
  % spin current 2
  J2a=zeros(3,3);
  J2b=zeros(3,3);
  for a=1:3; for b=1:3; for c=1:3; for j=1:3; for k=1:3; for m=1:3;

    J2a(a,k) = J2a(a,k) + ee(a,b,c)*vR(c,j)*(
       + dd(m,1)*st*(vn(b)*vn(j)-dd(b,j))*gt(k) ...
       + dd(m,1)*(1-ct)*(vn(j)*gn(k,b) + vn(b)*gn(k,j)) ...
       - ct*ee(b,j,m)*vn(m) * gt(k) ...
       - st*ee(b,j,m)*gn(k,m) ...
    );
    J2b(a,k) = J2b(a,k) + ee(a,b,c)*vR(c,j)*(
       + dd(m,1)*st*(vn(b)*vn(k)-dd(b,k))*gt(j) ...
       + dd(m,1)*(1-ct)*(vn(k)*gn(j,b) + vn(b)*gn(j,k)) ...
       - ct*ee(b,k,m)*vn(m) * gt(j) ...
       - st*ee(b,k,m)*gn(j,m) ...
    );

  end; end; end; end; end; end;

  r1 = sum(sum((J1a-J2a).^2));
  r2 = sum(sum((J1b-J2b).^2));
  fprintf('test spin current 2. A: %e B: %e\n', r1,r2);

  %%%%%%%%%%%%%%%%%%%
  % spin current 3
  J3a=zeros(3,3);
  J3b=zeros(3,3);
  for a=1:3; for b=1:3; for c=1:3; for j=1:3; for k=1:3; for m=1:3;

    J3a(a,k) = J3a(a,k) + ee(a,b,c)*vR(c,j)*(
       + dd(m,1)*st*(vn(b)*vn(j)-dd(b,j))*gt(k) ...
       + dd(m,1)*(1-ct)*(vn(j)*gn(k,b) + vn(b)*gn(k,j)) ...
       - ct*ee(b,j,m)*vn(m) * gt(k) ...
       - st*ee(b,j,m)*gn(k,m) ...
    );
    J3b(a,k) = J3b(a,k) + ee(a,b,c)*vR(c,j)*(
       + dd(m,1)*st*(vn(b)*vn(k)-dd(b,k))*gt(j) ...
       + dd(m,1)*(1-ct)*(vn(k)*gn(j,b) + vn(b)*gn(j,k)) ...
       - ct*ee(b,k,m)*vn(m) * gt(j) ...
       - st*ee(b,k,m)*gn(j,m) ...
    );
  end; end; end; end; end; end;

  r1 = sum(sum((J1a-J3a).^2));
  r2 = sum(sum((J1b-J3b).^2));
  fprintf('test spin current 2. A: %e B: %e\n', r1,r2);


  %%%%%%%%%%%%%%%%%%%
  % spin current 1z
  J1za=zeros(3,1);
  J1zb=zeros(3,1);
  for a=1:3; for b=1:3; for c=1:3; for j=1:3;
    J1za(a) = J1za(a) + ee(a,b,c)*vR(c,j)*gR(3,b,j);
    J1zb(a) = J1zb(a) + dd(j,3)*ee(a,b,c)*vR(c,3)*gR(3,b,3);
  end; end; end; end;

  %%%%%%%%%%%%%%%%%%%
  % spin current 2z
  J2za=zeros(3,1);
  J2zb=zeros(3,1);

  for a=1:3; for b=1:3; for c=1:3; for j=1:3; for m=1:3;

    J2za(a) = J2za(a) + ee(a,b,c)*vR(c,j)*(
       + dd(m,1)*st*(vn(b)*vn(j)-dd(b,j))*gt(3) ...
       + dd(m,1)*(1-ct)*(vn(j)*gn(3,b) + vn(b)*gn(3,j)) ...
       - ct*ee(b,j,m)*vn(m) * gt(3) ...
       - st*ee(b,j,m)*gn(3,m) ...
    );
    J2zb(a) = J2zb(a) + ee(a,b,c)*vR(c,3)*dd(j,3)*(
       + dd(m,1)*st*(vn(b)*vn(3)-dd(b,3))*gt(3) ...
       + dd(m,1)*(1-ct)*(vn(3)*gn(3,b) + vn(b)*gn(3,3)) ...
       - ct*ee(b,3,m)*vn(m) * gt(3) ...
       - st*ee(b,3,m)*gn(3,m) ...
    );

  end; end; end; end; end;
  r1 = sum((J1za-J2za).^2);
  r2 = sum((J1zb-J2zb).^2);
  fprintf('test spin current 2z. A: %e B: %e\n', r1,r2);

  %%%%%%%%%%%%%%%%%%%
  % spin current 3z
  J3za=zeros(3,1);
  J3zb=zeros(3,1);

  for a=1:3; for b=1:3; for m=1:3;

    J3za(a) = J3za(a) + (...
     -    2      *dd(m,1)*dd(b,1)*(vn(a)*gt(3) + st*gn(3,a)) ...
     + 2*(1-ct)  *ee(a,m,b)*vn(b)*gn(3,m) ...
    );
    J3zb(a) = J3zb(a) + (
     + st       *dd(m,1) *ee(a,b,3) *vn(b)    *vn(3)*gt(3) ...
     + ct       *dd(m,1)*dd(b,1)*dd(a,3)      *vn(3)*gt(3) ...
     -           dd(b,1)*dd(m,1)              *vn(a)*gt(3) ...
     + (1-ct)   *dd(b,1)*dd(m,1) *vn(3)*vn(3) *vn(a)*gt(3) ...
%
     + ct*(1-ct) *dd(m,1) *ee(a,b,3) *vn(3)*gn(3,b) ...
     + ct*(1-ct) *dd(m,1) *ee(a,b,3) *vn(b)*gn(3,3) ...
     + st*st     *dd(a,3) *ee(m,3,b)*vn(m)*gn(3,b)
     + (1-ct)^2  *ee(a,b,m)*vn(m)*vn(3)*vn(3)*gn(3,b) ...
     + st*(1-ct)   *dd(b,1)*dd(m,1) * 2*vn(a)*vn(3)*gn(3,3) ...
     - st*(1-2*ct) *dd(m,1)*dd(b,1)*dd(a,3)      *gn(3,3) ...
     - ct*st       *dd(b,1)*dd(m,1)              *gn(3,a) ...
     - st*(1-ct)   *dd(b,1)*dd(m,1) *vn(3)*vn(3) *gn(3,a) ...
    );

  end; end; end;
  r1 = sum((J1za-J3za).^2);
  r2 = sum((J1zb-J3zb).^2);
  fprintf('test spin current 3z. A: %e B: %e\n', r1,r2);


  %%%%%%%%%%%%%%%%%%%
  % spin current 4z
  J4za=zeros(3,1);
  J4zb=zeros(3,1);

    J4za(1) = ...
       - 2*vn(1)*gt(3) ...
       - 2*st*gn(3,1) ...
       - 2*(1-ct) *(vn(2)*gn(3,3) - vn(3)*gn(3,2)) ...
    ;
    J4za(2) = ...
       - 2*vn(2)*gt(3) ...
       - 2*st*gn(3,2) ...
       - 2*(1-ct) *(vn(3)*gn(3,1) - vn(1)*gn(3,3)) ...
    ;
    J4za(3) = ...
       - 2*vn(3)*gt(3) ...
       - 2*st*gn(3,3) ...
       - 2*(1-ct) *(vn(1)*gn(3,2)-vn(2)*gn(3,1))...
    ;

    J4zb(1) = ...
     - (1 - (1-ct)*vn(3)^2) * vn(1)*gt(3) ...
     + st * vn(2)*vn(3)*gt(3) ...
%
     -     st*(ct + (1-ct)*vn(3)^2) *gn(3,1) ...
     + (1-ct)*(ct + (1-ct)*vn(3)^2) *vn(3)*gn(3,2) ...
     + (1-ct)*(ct - (1-ct)*vn(3)^2) *vn(2)*gn(3,3) ...
     + st*(1-ct)    * 2*vn(1)*vn(3)*gn(3,3) ...
     ;

    J4zb(2) = ...
     - (1 - (1-ct)*vn(3)^2) * vn(2)*gt(3) ...
     - st * vn(1)*vn(3)*gt(3) ...
%
     - (1-ct)*(ct+(1-ct)*vn(3)^2) *vn(3)*gn(3,1) ...
     -     st*(ct+(1-ct)*vn(3)^2) *gn(3,2) ...
     - (1-ct)*(ct-(1-ct)*vn(3)^2) *vn(1)*gn(3,3) ...
     + st*(1-ct)   * 2*vn(2)*vn(3)*gn(3,3) ...
     ;

    J4zb(3) = ...
     - (1-ct)*(1 - vn(3)^2) * vn(3)*gt(3) ...
%
     + (st*st + (1-ct)^2 *vn(3)^2) *vn(2)*gn(3,1) ...
     - (st*st + (1-ct)^2 *vn(3)^2) *vn(1)*gn(3,2) ...
     - st*(1-ct) * (1-vn(3)^2) *gn(3,3) ...
     ;


  r1 = sum((J1za-J4za).^2);
  r2 = sum((J1zb-J4zb).^2);
  fprintf('test spin current 4z. A: %e B: %e\n', r1,r2);


  %%%%%%%%%%%%%%%%%%%
  % spin current 5z (in Dmitriev's program)
  J5za=zeros(3,1);
  J5zb=zeros(3,1);

  DD45 = vn(1)*gn(3,2)-vn(2)*gn(3,1); % Nx Ny' - Ny Nx'
  FTN=(1-ct)*DD45 - st*gn(3,3) - gt(3)*vn(3);
  UJX = 2.0*(gt(3)*vn(1)+st*gn(3,1)+(1-ct)*(vn(2)*gn(3,3)-gn(3,2)*vn(3))) ...
     + ((1-ct)*vn(1)*vn(3)+vn(2)*st)*FTN;
  UJY=2.0D0*(gt(3)*vn(2)+st*gn(3,2)-(1-ct)*(vn(1)*gn(3,3)-gn(3,1)*vn(3))) ...
     +   ((1-ct)*vn(2)*vn(3)-vn(1)*st)*FTN;
  UJZ=2.0D0*(gt(3)*vn(3)+st*gn(3,3)+(1-ct)*(vn(1)*gn(3,2)-gn(3,1)*vn(2))) ...
     +   ((1-ct)*vn(3)**2+ct)*FTN;

  J5z = [UJX; UJY; UJZ];
  r1 = sum((J1za + 2*J1za - J5z).^2);
  fprintf('test spin current 5z. %e\n', r1);

end
