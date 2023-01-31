function [u,c,result,lambda,SegmentationScore,time] = ParameterTL1(f,x,p,I,ground,n)
%UNTITLED2 Summary of this function goes here
%   I is the image we are using
%   f is the function we are performing parameter tuning on 
%   x is a vector whose entries are the different parameters we will
%   attempt to use 
%   c is a vector containing all the non-tunable parameters
lambda = 0;
SegmentationScore = 0;
time = 0;

for i = 1:length(x)
    tic;
    [u1,c1,result1] = f(I,n,x(i),p(1),p(2),p(3),2);
    t = toc;
    score = ssim(result1,ground)
    if score > SegmentationScore
        SegmentationScore = score;
        lambda = x(i);
        time = t;
        u = u1; c = c1; result = result1;
    end
end



end