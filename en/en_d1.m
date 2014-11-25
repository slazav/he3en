% dipolar energy
function e = en_d1(r)
  e = (trace(r^2) + trace(r)^2);
end
