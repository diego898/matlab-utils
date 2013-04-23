function stop = optimplot(x, optimValues, state)
% plots the current point of a 2-d otimization
global optimemory
stop = false;
optimemory(end+1,:) = x;

