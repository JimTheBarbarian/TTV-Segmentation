function [u,c,result] = L1TV(f, n, lambda, beta1,beta2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% f is an image
% n is number of regions, 
% a is the transformation parameter
MAX_ITER = 1000;
rel_tol = 1e-4;



[row,col] = size(f);
%I = imsegkmeans(single(f),n);

%u = zeros(row, col, n); 
%for i = 1:n
%    u(:,:,n) = (I == i);
%end


[u,c] = FCM(f,n);
%u = reshape(u,row*col, n);
%u = projsplx(u);
u = reshape(u,row, col,n);
v = u;
p = zeros(size(u));
dx = Dx(v);
dy = Dy(v);
qx= zeros(size(dx));
qy = zeros(size(dy));
%c = max(f(:)) *rand(1,n);
% for i = 1:n
%     ui = u(:,:,i);
%     ci = u(:,:,i).*f;
%     c(i)= sum(ci(:))/sum(ui(:));
% end
for i = 1:MAX_ITER
    % u update step
    u_new = u_subproblem(c,f,beta1,v,p);

    % d update step
    for j = 1:n
        d_new = [Dx(v(:,:,j)),Dy(v(:,:,j))] + 1/beta2*[qx(:,:,j),qy(:,:,j)];
        d_new = reshape(d_new,row*col,2);
        [dx1, dx2] = shrink2(d_new(:,1), d_new(:,2), lambda/beta2);
        dx(:,:,j) = reshape(dx1, row,col);
        dy(:,:,j) =  reshape(dx2, row,col);
    end
    %dx = d_subproblem(v, lambda, beta2,a,qx,1);
    %dy = d_subproblem(v,lambda,beta2,a,qy,2);
    %v update step
    v = v_subproblem(v,p,qx, qy,u_new,beta1,beta2,dx,dy);

    % adjust c and Lagrange multipliers
    p = p + beta1 * (u_new-v);
    qx = qx + beta2 * (Dx(v)- dx);
    qy = qy + beta2 * (Dy(v)- dy);
    c = c_subproblem(f, u_new ,c);
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
for i =1:n
    result = result + (u(:,:,i) == max(u,[],3)).*c(i);
end


end