I = cdata; I = rescale(double(cdata));
rng(51729);
noise_level = 0.1;
x1 = (.05:0.05:.25); %List of potential Lambdas for TTV

p1 = [.25,.25,1]; % non tunable parameters for TTV
p2 = [.25,.25,.1]; % non tunable parameters for TVP-MS
p3 = [.25,.25,100]; % non tunable parameters for L1TV

I(I==0) = 10/255;

Ig = I + sqrt(noise_level)*randn(size(I));

ssim_result1 = zeros(3,1);
time_result1 = zeros(3,1);


[u11,c11,result11,lambda11,ssim_result1(1),time_result1(1)] = ParameterTL1(@TTV,x1,p1,Ig,I,2);
[u12,c12,result12,lambda12,ssim_result1(2),time_result1(2)] = ParameterTL1(@TTV,x1,p2,Ig,I,2);
[u13,c13,result13, lambda13, ssim_result1(3), time_result1(3)] = ParameterTL1(@TTV, x1,p3,Ig, I,2);


dice_val(1) = compute_retina_dice(u11(:,:,1), cdata);
dice_val(2) = compute_retina_dice(u12(:,:,1), cdata);
dice_val(3) = compute_retina_dice(u13(:,:,1), cdata);



%%l1-l2 fuzzy region method
[M,N] = size(Ig);
lambda_val = [0.5, 1, 3, 5, 7, 10];

pm2.outer_iter = 10;
pm2.alpha = 0.5;

pm2.c = 1e-8;
pm2.inner_iter = 300;
pm2.tau = 1/8;
pm2.nu = 5;
pm2.beta = 1;

[u,c] = FCM(Ig,2);

u_initial{1} = reshape(u(:,:,1), M,N);
u_initial{2} = ones(M,N)-u_initial{1};

best_L1_L2_dice = 0; 
for i = 1:6
    pm2.lambda = lambda_val(i);
    tic;
    L1_L2_f = fuzzy_L1L2(Ig, u_initial, pm2, 2);
    toc
    L1_L2_dice = max(dice(double(L1_L2_f{1}>0.5), double(cdata)/255), dice(double(L1_L2_f{1}<0.5), double(cdata)/255));
    if L1_L2_dice > best_L1_L2_dice
        best_L1_L2_dice = L1_L2_dice;
        best_L1_L2_f = L1_L2_f;
        best_L1_L2_lambda = pm2.lambda;
    end
end

figure; imagesc(labeloverlay(cdata, double(u11(:,:,1)>0.5), 'Colormap', 'autumn')); colormap gray;
figure; imagesc(labeloverlay(cdata, double(u12(:,:,1)>0.5), 'Colormap', 'autumn')); colormap gray;
figure; imagesc(labeloverlay(cdata, double(u13(:,:,1)>0.5), 'Colormap', 'autumn')); colormap gray;
figure; imagesc(labeloverlay(cdata, double(best_L1_L2_f{1}>0.5), 'Colormap', 'autumn')); colormap gray;
