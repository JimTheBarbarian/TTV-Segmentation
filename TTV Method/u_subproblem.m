function [u_new] = u_subproblem(c,f, beta, v,p)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% f here is our image
% v here is our auxillary
% beta is the parameter beta1 in our overall algorithm.

[row,col,regions] = size(v);
F = zeros(row,col,regions);
for i = 1:regions
    F(:,:,i) = f-c(i); % only way I know to add a vector component wise to different channels.
end
F = F.^2;

u_new = v - 1/beta * (F + p);
u_new = reshape(u_new, row*col, regions);
u_new = projsplx(u_new); % this paralellizes the simplex projection 
u_new = reshape(u_new, row, col, regions);


end