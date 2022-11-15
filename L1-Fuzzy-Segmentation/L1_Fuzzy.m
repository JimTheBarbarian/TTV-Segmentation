function [u,c,result] = L1_Fuzzy(f,n,lambda,beta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
MAX_ITER = 1000;
rel_tol = 1e-4;
[row,col,~] = size(f);
%I = imsegkmeans(single(f),n);

%u = zeros(row, col, n); 
%for i = 1:n
%    u(:,:,n) = (I == i);
%end
[u,c] = FCM(f,n);
u = reshape(u,row, col,n);


u_new = u;
F = zeros(row,col,n); % this is to be used in the w-subproblem
w = u;
dx = Dx(w);
dy = Dy(w);
%c = max(f(:))*rand(1,n);
Lambda_dx = zeros(size(dx));
Lambda_dy = zeros(size(dy));

Lambda_w = zeros(size(w));
for iter = 1:MAX_ITER
    % d update step
    dx = L1Shrink(u,Lambda_dx,beta,1);
    dy = L1Shrink(u,Lambda_dy,beta,2);


    % w update step
    for channel = 1 : n
        F(:,:,channel) = f - c(channel); F = abs(F);    
    end
    w = u + 1\beta * Lambda_w + - lambda/beta * F;
    w = reshape(w,row*col, n);
    w = projsplx(w);
    w = reshape(w, row, col, n);

    % c update step 
    I_vec = f(:); I_vec = sort(I_vec);
    for i = 1:n
        w_vec = w(:,:,i); w_vec = sort(w_vec(:));
        A = sum(w_vec); B = cumsum(w_vec);
        if isequal(ismember(B,A/2),zeros(size(w_vec)))
            j_1 = find(A-2*B < 0,1); 
            c(i) = I_vec(j_1);
        else
            j_1 = find(A-2*B < 0,1); j_2 = j_1 -1;
            c(i) = 1/2 * (I_vec(j_1) + I_vec(j_2));
        end
    end
    

    

    % u update step 
    %build kernel: use the fft algorithm (5-pt stencil)
    uker = zeros(row,col);
    uker(1,1) = 4;uker(1,2)=-1;uker(2,1)=-1;uker(row,1)=-1;uker(1,col)=-1;


    L = fft2(uker); %this is -F(\Delta)
    denominator = 1+ L; %Precompute the denominator found in the equation for v_i. 
    for j = 1:n
        Ax = dx(:,:,j) - 1/beta * Lambda_dx(:,:,j);
        Ay = dy(:,:,j) - 1/beta * Lambda_dy(:,:,j);
        numerator = fft2(Dyt(Ay)+Dxt(Ax)+ w(:,:,j) -1/beta * Lambda_w(:,:,j));
        u_new(:,:,j) = ifft2(numerator./denominator);
    end
    % update dual variables 
    Lambda_dx = Lambda_dx + beta * (Dx(u) - dx);
    Lambda_dy = Lambda_dy + beta * (Dy(u) - dy);
    Lambda_w = Lambda_w + beta * (u- w);

    error = norm(u_new(:) - u(:)) / max([norm(u(:)), norm(u_new(:))]);
    if error < rel_tol
        break
    else
       fprintf("relative error %d\n", error);
       u = u_new;
       subplot(1,n,1), imshow(u(:,:,1) > .5);
       subplot(1,n,2), imshow(u(:,:,2) > .5);
       %subplot(1,n,3), imshow(u(:,:,3) > .5);
       %subplot(1,n,4), imshow(u(:,:,4) > .5);
    end

end
result = zeros(size(f));
for i =1:n
    result = result + (u(:,:,i) == max(u,[],3)).*c(i);
end
end