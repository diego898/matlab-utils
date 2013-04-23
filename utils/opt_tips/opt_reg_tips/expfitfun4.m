function [obj,lincoef] = expfitfun4(coef,xdata,ydata)
% Partitioned least squares solution for the batched
% single exponential , used from lsqnonlin

[m,n] = size(xdata);
M = exp(xdata*spdiags(-coef,0,n,n));
% solve for the linear parameters in one operation
i = (1:(m*n))';
j = repmat(1:n,m,1);
j=j(:);
M = sparse(i,j,M(:),n*m,n);
lincoef = M\ydata(:);
obj = ydata(:) - M*lincoef;


