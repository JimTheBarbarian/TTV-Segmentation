function [c] = c_subproblem(f,u,c)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
% Here f is our image and u is the fuzzy map
n = size(c,2);
for i = 1: n
    c(i) = sum(sum(f.*u(:,:,i))) / (sum(sum(u(:,:,i)))+1e-10);
end
