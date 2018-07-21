% test gradient torque
% create a random rotation matrix with a random gradient of

% test for a single random matrix with gradient ~d
function test_torque()
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

  % create value matrices
  ma = zeros(3,3,3);
  mb = zeros(3,3,3);
  mt = zeros(3,3,3);

  % central value
  ma(2,2,2)=rand*pi;
  mb(2,2,2)=rand*pi;
  mt(2,2,2)=rand*pi;

  % other values
  for i=1:3; for j=1:3, for k=1:3
    if i==2 && j==2 && k==2; continue; end
    ma(i,j,k)=ma(2,2,2)+(rand-0.5)*d*d;
    mb(i,j,k)=mb(2,2,2)+(rand-0.5)*d*d;
    mt(i,j,k)=mt(2,2,2)+(rand-0.5)*d*d;
  end; end; end

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
  fprintf('test R: %e\n', r);

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

  r = sum(sum(sum((testgR-gR).^2)));
  fprintf('test gR: %e\n', r);

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

  fprintf('test ggR: %e\n', r);


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
  ct = cos(vt); st=sin(vt);
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

  r = sum((T1a-T2a).^2);
  fprintf('test torque 2: %e\n', r);



end
