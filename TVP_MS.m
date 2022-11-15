function [u,c,result] = TVP_MS(f,n,lambda, mu, eta, r1, r2,p)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
MAX_ITER = 1000;
rel_tol = 1e-4;
[row,col,~] = size(f);

u=zeros(row,col,n);
v = zeros(size(u));
qx = Dx(v);
qy = Dy(v);



c = zeros(1,n);


I = imsegkmeans(single(f),n);

for i = 1:n
    u(:,:,n) = (I == i);
end
%v = u;
%qx = Dx(v);
%qy = Dy(v);
Lambda1 = zeros(size(u));
Lambda2_x = zeros(size(qx));
Lambda2_y = zeros(size(qy));

 uker = zeros(row,col);
 uker(1,1) = 4;uker(1,2)=-1;uker(2,1)=-1;uker(row,1)=-1;uker(1,col)=-1;


L = fft2(uker); %this is -F(\Delta)
denominator_v = r1 + r2 * L; % Precompute the denominator found in the equation for v_i and b_i

for i = 1:MAX_ITER
    % c update step, because we are performing image segmentation
    for j = 1: n
        c(j) = sum(sum(f.*u(:,:,j))) / (sum(sum(u(:,:,j)))+1e-10);
    end

    % q update step
    Ax = Dx(v)-r2*Lambda2_x; Ay = Dy(v)-r2*Lambda2_y;
    for k = 1:n
        qx(:,:,k) = GST(qx(:,:,k) - Ax(:,:,k),r2,p,30);
        qy(:,:,k) = GST(qy(:,:,k) - Ay(:,:,k),r2,p,30);
    end


    % u update step
    z = v - r1*Lambda1 - r1*f;
    z = reshape(z, row*col,n);
    u_new = projsplx(z);
    u_new = reshape(u_new,row,col,n);


    %v update step
    for k = 1:n
        g1 = qx(:,:,k) - r2 * Lambda2_x(:,:,k); g2 = qy(:,:,k) - r2*Lambda2_y(:,:,k);
        numerator = r2* fft2(u(:,:,k) + r1* Lambda1(:,:,k)) - fft2(Dxt(g1) - Dyt(g2));
        v(:,:,k) = ifft2(numerator./denominator_v);
    end

    Lambda1 = Lambda1 - 1/r1 * (v - u);
    Lambda2_x = Lambda2_x - 1/r2 * (qx - Dx(v));
    Lambda2_y = Lambda2_y - 1/r2 * (qy-Dy(v));

    error = norm(u_new(:) - u(:)) / max([norm(u(:)), norm(u_new(:))]);
    if error < rel_tol
        break
    else
       fprintf("relative error %d\n", error);
       u = u_new;
       for j = 1:n
        subplot(1,n,j), imshow(u(:,:,j) > .5);
       end
       pause(.01)
    end
end
result = zeros(size(f));
for  i= 1:n
    result = result + (u(:,:,i) == max(u,[],3)).*c(i);
end
end