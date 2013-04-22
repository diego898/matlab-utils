function [approx, error] = taylorApprox(func, point, order, range_offset)
%taylorApprox   Construct, return and plot the taylor series approx.
%
%   taylorApprox(func, point, order, range_offset) will construct the
%   taylor series approx of the string in func to the specified order
%   evaluated at some range offset centered around point.
%
%   approx is the symbolic approximation
%
%   error is the error values in the range



%% Construct Taylor Series

% make symbolic function
funct = sym(func);

% in terms of x
syms x;

% calculate DC Term
dc_term = subs(funct,x,point);

% begin building the approx
approx = dc_term;

for n = 1:order
    curr_term = (subs(diff(funct,n), x, point)/factorial(n))*(x-point)^(n);
    approx = approx + curr_term;
end



%% Calculate Error
eval_range = (point - range_offset):0.01:(point + range_offset);
approx_vals = subs(approx, x, eval_range);
funct_vals = subs(funct,x,eval_range);
error = abs(approx_vals - funct_vals);



%% Plot Results
figure;

% plot func and approx
subplot(2,1,1)
plot(eval_range, funct_vals, 'b', eval_range, approx_vals, 'r--');
xlabel('X Values');
ylabel('Function Values');
title(sprintf('%d^{th} Order Taylor Series Expansion of %s about %f', order, func, point));
legend('func', 'approx')

% plot error
subplot(2,1,2)
plot(eval_range, error, 'k*');
xlabel('X Values');
ylabel('Error');
title(sprintf('Error of %d^{th} Order Taylor Series Expansion of %s about %f', order, func, point));



end     % end function