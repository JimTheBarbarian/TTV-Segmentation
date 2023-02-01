I1 = image1;
noise_level = 0.01;
rng(51729);
x1 = (.01:.01:.05); %List of potential Lambdas

p1 = [.25,.25,100]; % non tunable parameters for TTV
p2 = [0.25,0.25, 1/3]; % non tunable parameters for TVP-MS
p3 = [0,.25,.25]; % non tunable parameters for L1TV


I1(I1 == 0) = 10; 

I1 = double(I1)/154; 
I1g = I1 + sqrt(noise_level)*randn(size(I1)); 


csf_dice1 = zeros(3, 1);
wm_dice1 = zeros(3,1);
gm_dice1 = zeros(3,1);
ssim_result1 = zeros(3,1);
time_result1 = zeros(3,1);


 %Brain 1
[u11,c11,result11,lambda11,ssim_result1(1),time_result1(1)] = ParameterTL1(@TTV,x1,p1,I1g,I1,4);
[u12,c12,result12,lambda12,ssim_result1(2),time_result1(2)] = ParameterTL1(@TVp,x1,p2,I1g,I1,4);
[u13,c13,result13,lambda13,ssim_result1(3),time_result1(3)] = ParameterTL1(@L1TV,x1,p3,I1g,I1,4);

[csf_dice1(1), wm_dice(1), gm_dice(1)] = compute_brain_dice(u11,image1);
[csf_dice1(2), wm_dice(2), gm_dice(2)] = compute_brain_dice(u12,image1);
[csf_dice1(3), wm_dice(3), gm_dice(3)] = compute_brain_dice(u13,image1);

%%L1-L2 method
l1ml2_lambda = [20, 25, 30, 35, 40, 45];

pm2.outer_iter = 10;
pm2.alpha = 0.5;

pm2.c = 1e-8;
pm2.inner_iter = 300;
pm2.tau = 1/8;
pm2.nu = 5;
pm2.beta = 1;

[u,c] = FCM(I1g,4);
[M,N] = size(image1);

u_initial{1} = reshape(u(:,:,1), M,N);
u_initial{2} = reshape(u(:,:,2), M,N);
u_initial{3} = reshape(u(:,:,3), M,N);
u_initial{4} = reshape(u(:,:,4), M,N);
best_l1ml2_ssim_val = 0;
for i = 1:6
    pm2.lambda = l1ml2_lambda(i);
    tic;
    [L1_L2_f, L1_L2_c] = fuzzy_L1L2(I1g, u_initial, pm2, 4);
    toc

    L1_L2_u = zeros(M,N,4);
    for j = 1:4
        L1_L2_u(:,:,j) = L1_L2_f{j};
    end

    result = zeros(size(I1g));
    for j =1:4
        result = result + (L1_L2_u(:,:,j) == max(L1_L2_u,[],3)).*L1_L2_c{j};
    end

    l1ml2_ssim = ssim(result, I1);
    if l1ml2_ssim > best_l1ml2_ssim_val
        best_l1ml2_ssim_val = l1ml2_ssim;
        best_l1ml2_val = pm2.lambda;
        best_l1ml2_u = L1_L2_u;
    end
end
ssim_result1(4) = best_l1ml2_ssim_val;
[csf_dice1(4), wm_dice(4), gm_dice(4)] = compute_brain_dice(best_l1ml2_u,image1);
