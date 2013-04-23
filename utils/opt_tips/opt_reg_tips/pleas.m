function [INLP,ILP] = pleas(funlist,NLPstart,xdata,ydata,options)
% pleas: partitioned nonlinear least squares estimation
% usage: [INLP,ILP] = pleas(funlist,NLPstart,xdata,ydata)
% usage: [INLP,ILP] = pleas(funlist,NLPstart,xdata,ydata,options)
%
% arguments: (input)
%  funlist - cell array of functions comprising the nonlinear parts
%            of each term in the model. Each independent function in
%            this list must transform xdata using a vector of intrinsicly
%            nonlinear parameters into an array of the same size and
%            shape as ydata. The arguments to each function will be in
%            the order (coef,xdata).
%
%            These functions may be
%             - scalar (double) constants (E.g., 1)
%             - anonymous functions
%             - inline functions
%             - character function names
%
%  NLPstart - vector of starting values for the intrinsicly nonlinear
%            parameters only.
%
%  xdata   - array of independent variables
%            
%  ydata   - array of dependent variable data
%
%  options - options structure appropriate for lsqnonlin
%
%
% arguments (output)
%  INLP - optimized list of intrinsicly nonlinear parameters
%
%  ILP  - optimized list of intrinsicly linear parameters 
%
%
% Example usage:
%  Fit a simple exponential model plus a constant term to data
%
%   x = rand(100,1);
%   y = 4 - 3*exp(2*x) + randn(size(x));
%
%   funlist = {1, @(xdata,coef) exp(xdata*coef)};
%   NLPstart = 1;
%   options = optimset('disp','iter');
%   [INLP,ILP] = pleas(funlist,NLPstart,x,y,options)
%
% Output:
%                                          Norm of      First-order 
%  Iteration  Func-count     f(x)          step          optimality   CG-iterations
%     0          2         116.796                          40.5
%     1          4         74.6378        1.00406            2.6            1
%     2          6         74.4382      0.0758513         0.0443            1
%     3          8         74.4381     0.00131311       0.000597            1
% Optimization terminated: relative function value
%  changing by less than OPTIONS.TolFun.
%
% INLP =
%    2.0812
%
% ILP =
%    3.6687
%   -2.7327

% Check the functions
if ~iscell(funlist)
  error 'funlist must be a cell array of functions, even if only one fun'
end
nfun=length(funlist);
for i = 1:nfun
  fi = funlist{i};
  
  % There are two cases where we need to turn the supplied
  % function into an executable function
  if isa(fi,'double')
    % a constant
    funlist{i} = @(xdata,coef) repmat(fi,size(ydata));
  elseif ischar(fi)
    % a character function name
    funlist{i} = str2func(fi);
  end
end

% were any options supplied?
if (nargin<5) || isempty(options)
  options = optimset('lsqnonlin');
end

% make sure that ydata is a column vector
ydata = ydata(:);
ny = length(ydata);

% ================================================
% =========== begin nested function ==============
% ================================================
function [res,ILP] = pleas_obj(INLP)
  % nested objective function for lsqnonlin, so all
  % the data and funs from pleas are visible to pleas_obj
  
  % loop over funlist
  A = zeros(ny,nfun);
  for i=1:nfun
    fi = funlist{i};
    term = fi(xdata,INLP);
    A(:,i) = term(:);
  end
  
  % do the linear regression using \
  ILP = A\ydata;
  
  % residuals for lsqnonlin
  res = A*ILP - ydata;
  
end % nested function termination
% ================================================
% ============= end nested function ==============
% ================================================

% call lsqnonlin, using a nested function handle
LB=[];
UB=[];
INLP = lsqnonlin(@pleas_obj,NLPstart,LB,UB,options);

% call one final time to get the final linear parameters
[junk,ILP] = pleas_obj(INLP);

end % main function terminator




