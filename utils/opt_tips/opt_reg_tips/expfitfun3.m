function [obj,lincoef] = expfitfun3(coef,xdata,ydata)
% Partitioned least squares solution for the single
% exponential , uses lsqnonlin

M = exp(-coef*xdata);
lincoef = M\ydata;

obj = ydata - lincoef*M;


