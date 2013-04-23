function x = testnestfun(y)
  function res = nestfun(xi)
    res = erf(xi) - yi;
  end
x = zeros(size(y));
start = 0;
for i=1:prod(size(y))
  yi = y(i);
  x(i) = fzero(@nestfun,start,optimset('disp','off'));
end
end % testnestfun terminator
