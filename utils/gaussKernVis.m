% Demo to visualize the mapping with a Gaussian Radial Basis Function,
% especially in the context of Support Vector Machines.
%
% When this script is executed, first a collection of red points can be
% clicked on the graph.
% After that, the blue points can be generated.
% Then the user must provide a gamma value.
% The final graph can be rotated to inspect the 3D-space (use the "turn"
% icon next the hand in the toolbar)
%
% Created by Roemer Vlasveld (roemer.vlasveld@gmail.com)
% 
% The blog post where this is used:
% http://rvlasveld.github.io/blog/2013/07/12/introduction-to-one-class-support-vector-machines/
%
% Please feel free to use this script to your need.
 
figure;
axis([-10 10 -10 10])
hold on
grid on;
% Initially, the list of points is empty.
red = [];
blue = [];
 
% Loop, picking up the points for the red class.
disp('---')
disp('Click in the graph for the red points, e.g. in a wide circular form')
disp('Left mouse button picks points.')
disp('Right mouse button picks last point.')
but = 1;
n = 0;
while but == 1
    [xi,yi,but] = ginput(1);
    plot(xi,yi,'ro')
    n = n+1;
    red(:,n) = [xi;yi];
end
 
disp('Finished collection red points')
disp('---')
 
% Loop again, picking up the points for the blue class
disp('Now click in the graph for the blue points, e.g. in a smaller circular form')
disp('Left mouse button picks points.')
disp('Right mouse button picks last point.')
but = 1;
n = 0;
while but == 1
    [xi,yi,but] = ginput(1);
    plot(xi,yi,'bo')
    n = n+1;
    blue(:,n) = [xi;yi];
end
 
disp('Finished collection blue points')
disp('---')
 
sigma = input('sigma = ? (default value: 1): ');
if isempty(sigma)
    sigma = 1;
end
 
project = @(data, sigma) sum(exp(-(squareform( pdist(data, 'euclidean') .^ 2) ./ ( 2*sigma^2))));
 
blue_z = project(blue', sigma);
red_z = project(red', sigma);
 
clf;
hold on;
grid on;
scatter3(red(1,:), red(2,:), red_z, 'r');
scatter3(blue(1,:), blue(2,:), blue_z, 'b')