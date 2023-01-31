%load('C:\Users\edaya\OneDrive\Documents\MATLAB\TTV Segmentation\Images\brain_image\result.mat'); 
I = image1;
I(I==0) = 10; 
csf_m = double(image1 == 48);
white_matter_m = double(image1==154);
grey_matter_m = double(image1==106); 
I = double(I)/154; 
Ig = I + sqrt(.02)*randn(size(I)); 
% TTV segmentation, lambda = .05
tic;
[u1,c1,result1]= TTV(Ig,4,.01,.25,.25,1,2);
time1 = toc
TTV_approxn1 = result1;
score1 = ssim(TTV_approxn1,I);
dice_csf_ttv1 = max([dice(double(u1(:,:,1)==1), csf_m), dice(double(u1(:,:,2)==1), csf_m), dice(double(u1(:,:,3)==1), csf_m), dice(double(u1(:,:,4)==1), csf_m)])
dice_wm_ttv1 = max([dice(double(u1(:,:,1)==1), white_matter_m), dice(double(u1(:,:,2)==1), white_matter_m), dice(double(u1(:,:,3)==1), white_matter_m), dice(double(u1(:,:,4)==1), white_matter_m)])
dice_gm_ttv1 = max([dice(double(u1(:,:,1)==1), grey_matter_m), dice(double(u1(:,:,2)==1), grey_matter_m), dice(double(u1(:,:,3)==1), grey_matter_m), dice(double(u1(:,:,4)==1), grey_matter_m)])


% TTV segmentation, lambda = .075
tic;
[u1,c1,result1]= TTV(Ig,4,.05,.25,.25,1,2);
time2 = toc
TTV_approxn5 = result1;
score2 = ssim(TTV_approxn5,I);
dice_csf_ttv2 = max([dice(double(u1(:,:,1)==1), csf_m), dice(double(u1(:,:,2)==1), csf_m), dice(double(u1(:,:,3)==1), csf_m), dice(double(u1(:,:,4)==1), csf_m)])
dice_wm_ttv2 = max([dice(double(u1(:,:,1)==1), white_matter_m), dice(double(u1(:,:,2)==1), white_matter_m), dice(double(u1(:,:,3)==1), white_matter_m), dice(double(u1(:,:,4)==1), white_matter_m)])
dice_gm_ttv2 = max([dice(double(u1(:,:,1)==1), grey_matter_m), dice(double(u1(:,:,2)==1), grey_matter_m), dice(double(u1(:,:,3)==1), grey_matter_m), dice(double(u1(:,:,4)==1), grey_matter_m)])


% TTV segmentation, lambda = .1
tic;
[u1,c1,result1] = TTV(Ig,4,.075,.25,.25,10,2);
TTV_approxn10 = result1;
time3 = toc
score3 = ssim(TTV_approxn10,I);
dice_csf_ttv3 = max([dice(double(u1(:,:,1)==1), csf_m), dice(double(u1(:,:,2)==1), csf_m), dice(double(u1(:,:,3)==1), csf_m), dice(double(u1(:,:,4)==1), csf_m)])
dice_wm_ttv3 = max([dice(double(u1(:,:,1)==1), white_matter_m), dice(double(u1(:,:,2)==1), white_matter_m), dice(double(u1(:,:,3)==1), white_matter_m), dice(double(u1(:,:,4)==1), white_matter_m)])
dice_gm_ttv3 = max([dice(double(u1(:,:,1)==1), grey_matter_m), dice(double(u1(:,:,2)==1), grey_matter_m), dice(double(u1(:,:,3)==1), grey_matter_m), dice(double(u1(:,:,4)==1), grey_matter_m)])

%L1 total variation fuzzy 
tic;
[u1,c1,result1] = L1TV(Ig,4,.05,.25,.25,2);
time4 = toc
L1TV_approxn = result1;
score4 = ssim(L1TV_approxn,I);
dice_csf_ttv4 = max([dice(double(u1(:,:,1)==1), csf_m), dice(double(u1(:,:,2)==1), csf_m), dice(double(u1(:,:,3)==1), csf_m), dice(double(u1(:,:,4)==1), csf_m)])
dice_wm_ttv4 = max([dice(double(u1(:,:,1)==1), white_matter_m), dice(double(u1(:,:,2)==1), white_matter_m), dice(double(u1(:,:,3)==1), white_matter_m), dice(double(u1(:,:,4)==1), white_matter_m)])
dice_gm_ttv4 = max([dice(double(u1(:,:,1)==1), grey_matter_m), dice(double(u1(:,:,2)==1), grey_matter_m), dice(double(u1(:,:,3)==1), grey_matter_m), dice(double(u1(:,:,4)==1), grey_matter_m)])



subplot(1,6,1), imshow(I);
subplot(1,6,2), imshow(Ig);
subplot(1,6,3), imshow(TTV_approxn1); 
subplot(1,6,4), imshow(TTV_approxn5);
subplot(1,6,5), imshow(TTV_approxn10);
subplot(1,6,6), imshow(L1TV_approxn);