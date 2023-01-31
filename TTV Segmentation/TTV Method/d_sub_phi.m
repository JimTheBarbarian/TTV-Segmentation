function [output] = d_sub_phi(t, lambda, a)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
output = acos(1 - (27 * lambda * a * (a+1))./(2*(a + abs(t)).^3));

end