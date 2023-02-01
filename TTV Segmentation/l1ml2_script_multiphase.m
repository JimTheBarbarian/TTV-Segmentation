[M,N] = size(I1g);

%% fuzzy competition models
pm2.outer_iter = 10;
pm2.alpha = 0.5;

pm2.lambda = 25; %20, 25, 30, 35, 40
pm2.c = 1e-8;
pm2.inner_iter = 300;
pm2.tau = 1/8;
pm2.nu = 5;
pm2.beta = 1;

[u,c] = FCM(I1g,4);

u_initial{1} = reshape(u(:,:,1), M,N);
u_initial{2} = reshape(u(:,:,2), M,N);
u_initial{3} = reshape(u(:,:,3), M,N);
u_initial{4} = reshape(u(:,:,4), M,N);

tic;
[L1_L2_f, L1_L2_c] = fuzzy_L1L2(I1g, u_initial, pm2, 4);
toc

L1_L2_u = zeros(M,N,4);
for i = 1:4
    L1_L2_u(:,:,i) = L1_L2_f{i};
end

result = zeros(size(I1g));
for i =1:4
    result = result + (L1_L2_u(:,:,i) == max(L1_L2_u,[],3)).*L1_L2_c{i};
end
ssim(result, I1)

[~,threshold_L1_L2_u] = max(L1_L2_u,[],3);

csf_dice1(4) = max([dice(double(threshold_L1_L2_u==1), csf_m1), dice(double(threshold_L1_L2_u==2), csf_m1), dice(double(threshold_L1_L2_u==3), csf_m1), dice(double(threshold_L1_L2_u==4), csf_m1)]);
wm_dice1(4) = max([dice(double(threshold_L1_L2_u==1), white_matter_m1), dice(double(threshold_L1_L2_u==2), white_matter_m1), dice(double(threshold_L1_L2_u==3), white_matter_m1), dice(double(threshold_L1_L2_u==4), white_matter_m1)]);
gm_dice1(4) = max([dice(double(threshold_L1_L2_u==1), grey_matter_m1), dice(double(threshold_L1_L2_u==2), grey_matter_m1), dice(double(threshold_L1_L2_u==3), grey_matter_m1), dice(double(threshold_L1_L2_u==4), grey_matter_m1)]);
