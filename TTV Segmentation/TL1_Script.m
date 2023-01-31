I = imread('C:\Users\edaya\OneDrive\Documents\MATLAB\TTV Segmentation\Images\retina\21_manual1.gif'); 
I = rescale(I); I(I == 0) = 1e-10;
% TTV segmentation, a = 1
tic;
[~,~,result1]= TTV(I,2,.005,.25,.25,1,1);
time = toc
TTV_approx1 = result1;

% TTV segmentation, a = 5
tic;
[~,~,result1]= TTV(I,2,.005,.25,.25,5,1);
time = toc
TTV_approx5 = result1;


% TTV segmentation, a = 10
tic;
[~,~,result1] = TTV(I,2,.005,.25,.25,10,1);
TTV_approx10 = result1;

%L1 total variation fuzzy 
tic;
[~,~,result1] = L1TV(I,2,.005,.25,.25,1);
time = toc