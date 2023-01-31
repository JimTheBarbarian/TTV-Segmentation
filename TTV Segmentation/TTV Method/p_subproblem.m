function [p] = p_subproblem(tau,u,theta,p)
%UNTITLED2 Performs fixed point iteration to find the optimal p 
%   Detailed explanation goes here

[~,~,n] = size(u);



for i = 1:100
    for j = 1:n
        A = Div(p(:,:,j)) - 1/theta * u(:,:,j);
        numerator = p(:,:,j) + tau * Dx(A) + Dy(A);
        denominator = 1 + tau * abs(Dx(A) + Dy(A));

        p(:,:,j) = numerator./denominator;
    end


end


end