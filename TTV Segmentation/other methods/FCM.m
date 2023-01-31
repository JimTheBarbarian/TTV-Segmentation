function [u,C] = FCM(f,n)
%UNTITLED2 Summary of this function goes here
%   Original FCM code written in Python by Yoo Jeonghwa.
%   That code is found at https://github.com/jeongHwarr/various_FCM_segmentation/blob/master/FCM.py
[row,col,~] = size(f);
X = reshape(f,row*col,1);
u = zeros(row*col,1,n);

% initializing u
idx = 1:row*col;
for i = 1:n 
    idxii = (mod(idx,n) == i-1);
    u(idxii,i) = 1;
end


C = zeros(1,n);
MAX_ITER = 100;
rel_tol = 1e-4;

for i = 1:MAX_ITER

% Update C step 
    F = u.^2;
    for k = 1:n
        F_k = F(:,:,k);
        numerator = dot(X,F_k);
        denominator = sum(F_k);
        C(k) = numerator/denominator;
    end

    old_u  = u;
    
    %update u
    [c_mesh,x_mesh] = meshgrid(C,X);
    p2 = sum(1./(c_mesh - x_mesh).^2,2);
    for k = 1:n
        p1 = (X - C(k)).^2;
        u(:,:,k) = (1./(p1.*p2))';
    end
    
    
    d  = sum(abs(u - old_u),3);
    if d < rel_tol
        break
    end
end

% Defuzzify
%result = zeros(size(u));
%for i = 1:n
%    result(:,:,i) = (u(:,:,i) == max(u,[],3));
%end

end