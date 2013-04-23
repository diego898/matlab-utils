function epred = implicit_obj(ab,x,ye)
a = ab(1);
b = ab(2);
n = length(ye);
epred = zeros(1,n);
estart = .1;
for i = 1:n
  fun = @(ei) (ye(i)-ei - a*((ye(i)-ei) - b*x(i)).^2);
  epred(i) = fzero(fun,estart);
end