function [u,c,result] = L1FC(f,n,lambda,method)
%UNTITLED Summary of this function goes here
%   method is which method will be used to minimize over u. It should be
%   equal to either 1, 2 or 3.
MAX_ITER = 1000;
rel_tol = 1e-4;
theta = 0.1;
nu = 1000;


[row,col,~] = size(f);
if n == 2
    u = zeros(row,col,n);
    u(round(row/4): round(3*row/4), round(col/4):round(3*col/4),1) = 1;
    u(:,:,2) = 1- u(:,:,1);
else
    u = rand(row,col,n); u = reshape(u,row*col,n);
    u = projsplx(u); u = reshape(u,row,col,n);
end

v= u;

p = zeros(row,col,n);
c = zeros(1,n);


for i = 1:MAX_ITER

% c update step
c_new = c_subproblem(f,u,c);

% v update step

p = p_subproblem(.1,u,theta,p);
for j = 1:n
    v(:,:,j) = u(:,:,j)- theta * Div(p(:,:,j));
end

% u update step
if method == 1
    F = zeros(row,col,n);
    for j = 1:n
        F(:,:,j) = (f-c_new(j)).^2;
    end
    num = nu * theta  * (sum(v - lambda * theta * F ,3)-1);
    denom = 1 + n * nu * theta;
    for j = 1:n
        u(:,:,j) = v(:,:,j) - lambda * theta * F(:,:,j) - num/denom;
        u(:,:,j) = min(max(u(:,:,j),0),1);
    end
    
elseif method == 2

elseif method == 3


end


% Check termination condition 

error = norm(c_new(:) - c(:));
if error < rel_tol
    break
else
    fprintf("relative error %d\n", error);
    c = c_new;
    for j = 1:n
        subplot(1,n,j), imshow(u(:,:,j) > .5);
    end
       pause(.01)

    

end

end


result = zeros(size(f));
for i =1:n
    result = result + (u(:,:,i) == max(u,[],3)).*c(i);
end


end