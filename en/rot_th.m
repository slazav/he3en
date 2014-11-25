function r1 = rot_th(r, th)
% rotate matrix r by small Euler angles th
  et = [    0   th(3) -th(2)
         -th(3)    0   th(1)
          th(2) -th(1)    0];

  ttn = [th(1)*th(1) th(2)*th(1) th(3)*th(1)
         th(1)*th(2) th(2)*th(2) th(3)*th(2)
         th(1)*th(3) th(2)*th(3) th(3)*th(3)];
  dd=[1 0 0; 0 1 0; 0 0 1];
  tt = th*th'; %'
  r1 = (dd - et - dd*tt*0.5 + ttn*0.5)*r;
end

