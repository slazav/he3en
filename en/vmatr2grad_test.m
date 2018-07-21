# Text vmatr2grad function using polynomial function

function vmatr2grad_test()

  d=0.012;
  F0 = rand;
  G  = [rand, rand, rand];
  GG = [rand, rand, rand; rand, rand, rand; rand, rand, rand];
  # we will use only 6 components:
  GG(2,1) = GG(1,2);
  GG(3,2) = GG(2,3);
  GG(1,3) = GG(3,1);

  for i=1:3; for j=1:3; for k=1:3;
    x = (i-2)*d;
    y = (j-2)*d;
    z = (k-2)*d;
    F(i,j,k) = F0 ...
       + G(1)*x + G(2)*y + G(3)*z ...
       + (GG(1,1)*x^2 + GG(2,2)*y^2 + GG(3,3)*z^2)/2.0 ...
       + GG(1,2)*x*y + GG(2,3)*y*z + GG(3,1)*z*x;
  end; end; end

  [f,g,gg] = vmatr2grad(F, d);
  res = (f-F0)^2 ...
      + sum((G-g).^2) ...
      + sum(sum((GG-gg).^2));

  fprintf('%e\n', res);

end
