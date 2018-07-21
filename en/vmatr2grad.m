# From a matrix of function values at
# 3x3x3 grid with spacing d calculate
# first and second derivatives
#
# arguments: F is a 3x3x3 matrix, d is a grid size
# return valuse:
#   f - function in the central node,
#   g - gradient vector
#   gg - second derivatives, 3x3 matrix


function [v g gg] = vmatr2grad(F, d)
  v = F(2,2,2);

  g(1) = F(3,2,2) - F(1,2,2);
  g(2) = F(2,3,2) - F(2,1,2);
  g(3) = F(2,2,3) - F(2,2,1);
  g=g/(2*d);

  gg(1,1) = F(3,2,2) - 2*F(2,2,2) + F(1,2,2);
  gg(2,2) = F(2,3,2) - 2*F(2,2,2) + F(2,1,2);
  gg(3,3) = F(2,2,3) - 2*F(2,2,2) + F(2,2,1);

  gg(1,2) = ( F(3,3,2)-F(3,1,2)-F(1,3,2)+F(1,1,2) ) / 4;
  gg(2,3) = ( F(2,3,3)-F(2,1,3)-F(2,3,1)+F(2,1,1) ) / 4;
  gg(3,1) = ( F(3,2,3)-F(3,2,1)-F(1,2,3)+F(1,2,1) ) / 4;

  gg(2,1) = gg(1,2);
  gg(3,2) = gg(2,3);
  gg(1,3) = gg(3,1);
  gg=gg/d^2;
end
