% I2 = image2; 
%  I3 = image3; 
%I4 = image4;
I1 = image1;
rng(51729);
x1 = (.001:.0025:.156); %List of potential Lambdas for TTV
x2 = [.01:.01:1 2:1:20 40:20:200]; % List of lambdas for TVP-MS
x3 = (.001:.0035:.141); % List of lambdas for L1TV

p1 = [.25,.25,5]; % non tunable parameters for TTV
p2 = [10,200,1/3]; % non tunable parameters for TVP-MS
p3 = [0,.25,.25]; % non tunable parameters for L1TV


I1(I1 == 0) = 10; %I2(I2 == 0) = 10; 
 %I3(I3 == 0) = 10; 
 %I4(I4 == 0) = 10;


csf_m1 = double(image1 == 48);
white_matter_m1 = double(image1==154);
grey_matter_m1 = double(image1==106); 
I1 = double(I1)/154; 
I1g = I1 + sqrt(.02)*randn(size(I1)); 

%csf_m2 = double(image2 == 48);
%white_matter_m2 = double(image2>0.554);
%grey_matter_m2 = double(image2>0.506); 
%I2 = double(I2)/154; 
%I2g = I2 + sqrt(.01)*randn(size(I2)); 

%csf_m3 = double(image3 == 48);
%white_matter_m3 = double(image3>0.554);
%grey_matter_m3 = double(image3>0.506); 
%I3 = double(I3)/154; 
%I3g = I3 + sqrt(.01)*randn(size(I3));

% csf_m4 = double(image4 == 48);
% white_matter_m4 = double(image4>0.554);
% grey_matter_m4 = double(image4>0.506); 
% I4 = double(I4)/154; 
% I4g = I4 + sqrt(.04)*randn(size(I4));

csf_dice1 = zeros(3, 1);
wm_dice1 = zeros(3,1);
gm_dice1 = zeros(3,1);
ssim_result1 = zeros(3,1);
time_result1 = zeros(3,1);

%csf_dice2 = zeros(3, 1);
%wm_dice2 = zeros(3,1);
%gm_dice2 = zeros(3,1);
%ssim_result2 = zeros(3,1);
%time_result2 = zeros(3,1);

%csf_dice3 = zeros(3, 1);
%wm_dice3 = zeros(3,1);
%gm_dice3 = zeros(3,1);
%ssim_result3 = zeros(3,1);
%time_result3 = zeros(3,1);

csf_dice4 = zeros(3, 1);
wm_dice4 = zeros(3,1);
gm_dice4 = zeros(3,1);
ssim_result4 = zeros(3,1);
time_result4 = zeros(3,1);

 %Brain 1
[u11,c11,result11,lambda11,ssim_result1(1),time_result1(1)] = ParameterTL1(@TTV,x1,p1,I1g,I1,4);
[u12,c12,result12,lambda12,ssim_result1(2),time_result1(2)] = ParameterTL1(@TVP_MS,x2,p2,I1g,I1,4);%
[u13,c13,result13,lambda13,ssim_result1(3),time_result1(3)] = ParameterTL1(@L1TV,x3,p3,I1g,I1,4);

[~,threshold_u11] = max(u11,[],3);

csf_dice1(1) = max([dice(double(threshold_u11==1), csf_m1), dice(double(threshold_u11==2), csf_m1), dice(double(threshold_u11==3), csf_m1), dice(double(threshold_u11==4), csf_m1)]);
wm_dice1(1) = max([dice(double(threshold_u11==1), white_matter_m1), dice(double(threshold_u11==2), white_matter_m1), dice(double(threshold_u11==3), white_matter_m1), dice(double(threshold_u11==4), white_matter_m1)]);
gm_dice1(1) = max([dice(double(threshold_u11==1), grey_matter_m1), dice(double(threshold_u11==2), grey_matter_m1), dice(double(threshold_u11==3), grey_matter_m1), dice(double(threshold_u11==4), grey_matter_m1)]);


csf_dice1(1) = max([dice(double(u11(:,:,1)>0.5), csf_m1), dice(double(u11(:,:,2)>0.5), csf_m1), dice(double(u11(:,:,3)>0.5), csf_m1), dice(double(u11(:,:,4)>0.5), csf_m1)]);
wm_dice1(1) = max([dice(double(u11(:,:,1)>0.5), white_matter_m1), dice(double(u11(:,:,2)>0.5), white_matter_m1), dice(double(u11(:,:,3)>0.5), white_matter_m1), dice(double(u11(:,:,4)>0.5), white_matter_m1)]);
gm_dice1(1) = max([dice(double(u11(:,:,1)>0.5), grey_matter_m1), dice(double(u11(:,:,2)>0.5), grey_matter_m1), dice(double(u11(:,:,3)>0.5), grey_matter_m1), dice(double(u11(:,:,4)>0.5), grey_matter_m1)]);

csf_dice1(2) = max([dice(double(u12(:,:,1)>0.5), csf_m1), dice(double(u12(:,:,2)>0.5), csf_m1), dice(double(u12(:,:,3)>0.5), csf_m1), dice(double(u12(:,:,4)>0.5), csf_m1)]);
wm_dice1(2) = max([dice(double(u12(:,:,1)>0.5), white_matter_m1), dice(double(u12(:,:,2)>0.5), white_matter_m1), dice(double(u12(:,:,3)>0.5), white_matter_m1), dice(double(u12(:,:,4)>0.5), white_matter_m1)]);
gm_dice1(2) = max([dice(double(u12(:,:,1)>0.5), grey_matter_m1), dice(double(u12(:,:,2)>0.5), grey_matter_m1), dice(double(u12(:,:,3)>0.5), grey_matter_m1), dice(double(u12(:,:,4)>0.5), grey_matter_m1)]);

csf_dice1(3) = max([dice(double(u13(:,:,1)>0.5), csf_m1), dice(double(u13(:,:,2)>0.5), csf_m1), dice(double(u13(:,:,3)>0.5), csf_m1), dice(double(u13(:,:,4)>0.5), csf_m1)]);
wm_dice1(3) = max([dice(double(u13(:,:,1)>0.5), white_matter_m1), dice(double(u13(:,:,2)>0.5), white_matter_m1), dice(double(u13(:,:,3)>0.5), white_matter_m1), dice(double(u13(:,:,4)>0.5), white_matter_m1)]);
gm_dice1(3) = max([dice(double(u13(:,:,1)>0.5), grey_matter_m1), dice(double(u13(:,:,2)>0.5), grey_matter_m1), dice(double(u13(:,:,3)>0.5), grey_matter_m1), dice(double(u13(:,:,4)>0.5), grey_matter_m1)]);



%Brain 2
%[u21,c21,result21,lambda21,ssim_result2(1),time_result2(1)] = ParameterTL1(@TTV,x1,p1,I2g,I2,4);
%[u22,c22,result22,lambda22,ssim_result2(2),time_result2(2)] = ParameterTL1(@TVP_MS,x2,p2,I2g,I2,4);
%[u23,c23,result23, lambda23, ssim_result2(3), time_result2(3)] = ParameterTL1(@L1TV, x3,p3,I2g, I2,4);







%Brain 3
%[u31,c31,result31,lambda31,ssim_result3(1),time_result3(1)] = ParameterTL1(@TTV,x1,p1,I3g,I3,4);
%[u32,c32,result32,lambda32,ssim_result3(2),time_result3(2)] = ParameterTL1(@TVP_MS,x2,p2,I3g,I3,4);
%[u33,c33,result33, lambda33, ssim_result3(3), time_result3(3)] = ParameterTL1(@L1TV, x3,p3,I3g, I3,4);


%Brain 4
[u41,c41,result41,lambda41,ssim_result4(1),time_result4(1)] = ParameterTL1(@TTV,x1,p1,I4g,I4,4);
[u42,c42,result42,lambda42,ssim_result4(2),time_result4(2)] = ParameterTL1(@TVP_MS,x2,p2,I4g,I4,4);
[u43,c43,result43, lambda43, ssim_result4(3), time_result4(3)] = ParameterTL1(@L1TV, x3,p3,I4g, I4,4);

csf_dice4(1) = max([dice(double(u41(:,:,1)>0.5), csf_m4), dice(double(u41(:,:,2)>0.5), csf_m4), dice(double(u41(:,:,3)>0.5), csf_m4), dice(double(u41(:,:,4)>0.5), csf_m4)]);
wm_dice4(1) = max([dice(double(u41(:,:,1)>0.5), white_matter_m4), dice(double(u41(:,:,2)>0.5), white_matter_m4), dice(double(u41(:,:,3)>0.5), white_matter_m4), dice(double(u41(:,:,4)>0.5), white_matter_m4)]);
gm_dice4(1) = max([dice(double(u41(:,:,1)>0.5), grey_matter_m4), dice(double(u41(:,:,2)>0.5), grey_matter_m4), dice(double(u41(:,:,3)>0.5), grey_matter_m4), dice(double(u41(:,:,4)>0.5), grey_matter_m4)]);

csf_dice4(2) = max([dice(double(u42(:,:,1)>0.5), csf_m4), dice(double(u42(:,:,2)>0.5), csf_m4), dice(double(u42(:,:,3)>0.5), csf_m4), dice(double(u42(:,:,4)>0.5), csf_m4)]);
wm_dice4(2) = max([dice(double(u42(:,:,1)>0.5), white_matter_m4), dice(double(u42(:,:,2)>0.5), white_matter_m4), dice(double(u42(:,:,3)>0.5), white_matter_m4), dice(double(u42(:,:,4)>0.5), white_matter_m4)]);
gm_dice4(2) = max([dice(double(u42(:,:,1)>0.5), grey_matter_m4), dice(double(u42(:,:,2)>0.5), grey_matter_m4), dice(double(u42(:,:,3)>0.5), grey_matter_m4), dice(double(u42(:,:,4)>0.5), grey_matter_m4)]);

csf_dice4(3) = max([dice(double(u43(:,:,1)>0.5), csf_m4), dice(double(u43(:,:,2)>0.5), csf_m4), dice(double(u43(:,:,3)>0.5), csf_m4), dice(double(u43(:,:,4)>0.5), csf_m4)]);
wm_dice4(3) = max([dice(double(u43(:,:,1)>0.5), white_matter_m4), dice(double(u43(:,:,2)>0.5), white_matter_m4), dice(double(u43(:,:,3)>0.5), white_matter_m4), dice(double(u43(:,:,4)>0.5), white_matter_m4)]);
gm_dice4(3) = max([dice(double(u43(:,:,1)>0.5), grey_matter_m4), dice(double(u43(:,:,2)>0.5), grey_matter_m4), dice(double(u43(:,:,3)>0.5), grey_matter_m4), dice(double(u43(:,:,4)>0.5), grey_matter_m4)]);




