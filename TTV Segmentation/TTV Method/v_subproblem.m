function [v] = v_subproblem(v,p, qx, qy,u, beta1, beta2, dx, dy,L)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


 %build kernel: use the fft algorithm (5-pt stencil)

denominator = beta1 + beta2 * L; % Precompute the denominator found in the equation for v_i. 

[~,~,n] = size(v);
for i = 1:n 
    Ax = qx(:,:,i) - beta2 * (dx(:,:,i));
    Ay = qy(:,:,i) - beta2 * (dy(:,:,i));
    numerator = fft2(p(:,:,i) + beta1 * u(:,:,i) - Dxt(Ax) - Dyt(Ay));
    v(:,:,i) = ifft2(numerator./denominator);
end