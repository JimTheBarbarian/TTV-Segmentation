[M,N] = size(Ig10);

%% fuzzy competition models
pm2.outer_iter = 10;
pm2.alpha = 0.5;

pm2.lambda = 5; %1, 5, 10
pm2.c = 1e-8;
pm2.inner_iter = 300;
pm2.tau = 1/8;
pm2.nu = 5;
pm2.beta = 1;

[u,c] = FCM(Ig10,2);

u_initial{1} = reshape(u(:,:,1), M,N);
u_initial{2} = ones(M,N)-u_initial{1};

tic;
L1_L2_f = fuzzy_L1L2(Ig10, u_initial, pm2, 2);
toc