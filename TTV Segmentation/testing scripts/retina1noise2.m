I = cdata; I = rescale(double(cdata));
rng(51729);
x1 = (.04:.0025:.11); %List of potential Lambdas for TTV
x2 = (.01:.0025:.075); % List of lambdas for TVP-MS
x3 = (.0925:.0035:.1795); % List of lambdas for L1TV

p1 = [.25,.25,1]; % non tunable parameters for TTV
p2 = [.25,.25,.1]; % non tunable parameters for TVP-MS
p3 = [.25,.25,100]; % non tunable parameters for L1TV

I(I==0) = 10/255;

Ig2 = I + sqrt(.02)*randn(size(I));

ssim_result1 = zeros(3,1);
time_result1 = zeros(3,1);


[u11,c11,result11,lambda11,ssim_result1(1),time_result1(1)] = ParameterTL1(@TTV,x1,p1,Ig2,I,2);
[u12,c12,result12,lambda12,ssim_result1(2),time_result1(2)] = ParameterTL1(@TTV,x2,p2,Ig2,I,2);
[u13,c13,result13, lambda13, ssim_result1(3), time_result1(3)] = ParameterTL1(@TTV, x3,p3,Ig2, I,2);


Ig5 = I + sqrt(.05)*randn(size(I));


ssim_result2 = zeros(3,1);
time_result2 = zeros(3,1);

[u21,c21,result21,lambda21,ssim_result2(1),time_result2(1)] = ParameterTL1(@TTV,x1,p1,Ig5,I,2);
[u22,c22,result22,lambda22,ssim_result2(2),time_result2(2)] = ParameterTL1(@TTV,x2,p2,Ig5,I,2);
[u23,c23,result23, lambda23, ssim_result2(3), time_result2(3)] = ParameterTL1(@TTV, x3,p3,Ig5, I,2);


Ig10 = I + sqrt(.1)*randn(size(I));
ssim_result3 = zeros(3,1);
time_result3 = zeros(3,1);

[u31,c31,result31,lambda31,ssim_result3(1),time_result3(1)] = ParameterTL1(@TTV,x1,p1,Ig10,I,2);
[u32,c32,result32,lambda32,ssim_result3(2),time_result3(2)] = ParameterTL1(@TTV,x2,p2,Ig10,I,2);
[u33,c33,result33, lambda33, ssim_result3(3), time_result3(3)] = ParameterTL1(@TTV, x3,p3,Ig10, I,2);
