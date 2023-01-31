%load('C:\Users\edaya\OneDrive\Documents\MATLAB\TTV Segmentation\Images\brain_image\result.mat'); 
I = image1;
I(I==0) = 10; 
csf_m = double(image1 == 48);
white_matter_m = double(image1==154);
grey_matter_m = double(image1==106); 
I = double(I)/154; 
rng(51729);
Ig = I + sqrt(.01)*randn(size(I)); 
csf_dice = zeros(5, 1);
wm_dice = zeros(5,1);
gm_dice = zeros(5,1);
ssim_result = zeros(5,1);
time_result = zeros(5,1);
% TTV segmentation, lambda = .05
tic;
[u1,c1,result1]= TTV(Ig,4,.01,.25,.25,1,2);
time_result(1,1) = toc
TTV_approxn1 = result1;
ssim_result(1,1) = ssim(TTV_approxn1,I);
csf(1,1) = max([dice(double(u1(:,:,1)==1), csf_m), dice(double(u1(:,:,2)==1), csf_m), dice(double(u1(:,:,3)==1), csf_m), dice(double(u1(:,:,4)==1), csf_m)])
wm(1,1) = max([dice(double(u1(:,:,1)==1), white_matter_m), dice(double(u1(:,:,2)==1), white_matter_m), dice(double(u1(:,:,3)==1), white_matter_m), dice(double(u1(:,:,4)==1), white_matter_m)])
gm(1,1) = max([dice(double(u1(:,:,1)==1), grey_matter_m), dice(double(u1(:,:,2)==1), grey_matter_m), dice(double(u1(:,:,3)==1), grey_matter_m), dice(double(u1(:,:,4)==1), grey_matter_m)])


% TTV segmentation, lambda = .075
tic;
[u1,c1,result1]= TTV(Ig,4,.05,.25,.25,1,2);
time2 = toc
TTV_approxn5 = result1;
time_result(2,1) = toc
TTV_approxn5 = result1;
ssim_result(2,1) = ssim(TTV_approxn1,I);
csf(2,1) = max([dice(double(u1(:,:,1)==1), csf_m), dice(double(u1(:,:,2)==1), csf_m), dice(double(u1(:,:,3)==1), csf_m), dice(double(u1(:,:,4)==1), csf_m)])
wm(2,1) = max([dice(double(u1(:,:,1)==1), white_matter_m), dice(double(u1(:,:,2)==1), white_matter_m), dice(double(u1(:,:,3)==1), white_matter_m), dice(double(u1(:,:,4)==1), white_matter_m)])
gm(2,1) = max([dice(double(u1(:,:,1)==1), grey_matter_m), dice(double(u1(:,:,2)==1), grey_matter_m), dice(double(u1(:,:,3)==1), grey_matter_m), dice(double(u1(:,:,4)==1), grey_matter_m)])


% TTV segmentation, lambda = .1
tic;
[u1,c1,result1] = TTV(Ig,4,.075,.25,.25,10,2);
TTV_approxn10 = result1;
time_result(3,1) = toc
TTV_approxn10 = result1;
ssim_result(3,1) = ssim(TTV_approxn1,I);
csf(3,1) = max([dice(double(u1(:,:,1)==1), csf_m), dice(double(u1(:,:,2)==1), csf_m), dice(double(u1(:,:,3)==1), csf_m), dice(double(u1(:,:,4)==1), csf_m)])
wm(3,1) = max([dice(double(u1(:,:,1)==1), white_matter_m), dice(double(u1(:,:,2)==1), white_matter_m), dice(double(u1(:,:,3)==1), white_matter_m), dice(double(u1(:,:,4)==1), white_matter_m)])
gm(3,1) = max([dice(double(u1(:,:,1)==1), grey_matter_m), dice(double(u1(:,:,2)==1), grey_matter_m), dice(double(u1(:,:,3)==1), grey_matter_m), dice(double(u1(:,:,4)==1), grey_matter_m)])


%L1 total variation fuzzy 
tic;
[u1,c1,result1] = L1TV(Ig,4,.05,.25,.25,2);
time_result(4,1) = toc
L1TV_approxn = result1;
ssim_result(4,1) = ssim(TTV_approxn1,I);
csf(4,1) = max([dice(double(u1(:,:,1)==1), csf_m), dice(double(u1(:,:,2)==1), csf_m), dice(double(u1(:,:,3)==1), csf_m), dice(double(u1(:,:,4)==1), csf_m)])
wm(4,1) = max([dice(double(u1(:,:,1)==1), white_matter_m), dice(double(u1(:,:,2)==1), white_matter_m), dice(double(u1(:,:,3)==1), white_matter_m), dice(double(u1(:,:,4)==1), white_matter_m)])
gm(4,1) = max([dice(double(u1(:,:,1)==1), grey_matter_m), dice(double(u1(:,:,2)==1), grey_matter_m), dice(double(u1(:,:,3)==1), grey_matter_m), dice(double(u1(:,:,4)==1), grey_matter_m)])


%TV-P MS
tic;
[u1,c1,result1] = TVP_MS(Ig,4,120,10,200,1/2,2);
time_result(5,1) = toc
TVP_approxn = result1;
ssim_result(5,1) = ssim(TTV_approxn1,I);
csf(5,1) = max([dice(double(u1(:,:,1)==1), csf_m), dice(double(u1(:,:,2)==1), csf_m), dice(double(u1(:,:,3)==1), csf_m), dice(double(u1(:,:,4)==1), csf_m)])
wm(5,1) = max([dice(double(u1(:,:,1)==1), white_matter_m), dice(double(u1(:,:,2)==1), white_matter_m), dice(double(u1(:,:,3)==1), white_matter_m), dice(double(u1(:,:,4)==1), white_matter_m)])
gm(5,1) = max([dice(double(u1(:,:,1)==1), grey_matter_m), dice(double(u1(:,:,2)==1), grey_matter_m), dice(double(u1(:,:,3)==1), grey_matter_m), dice(double(u1(:,:,4)==1), grey_matter_m)])




subplot(1,7,1), imshow(I);
subplot(1,7,2), imshow(Ig);
subplot(1,7,3), imshow(TTV_approxn1); 
subplot(1,7,4), imshow(TTV_approxn5);
subplot(1,7,5), imshow(TTV_approxn10);
subplot(1,7,6), imshow(L1TV_approxn);
subplot(1,7,7), imshow(TVP_approxn)