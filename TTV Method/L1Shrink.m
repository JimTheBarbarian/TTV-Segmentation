function [d_new] = L1Shrink(u, q, beta,direction)
%UNTITLED2 Summary of this function goes here
if direction == 1
    d_new = Dx(u); 
else
    d_new = Dy(u); 
end
[~,~,regions] = size(u);
d_new = d_new + 1/beta * q;

for i = 1: regions
thresh = max(abs(d_new)- 1/beta, 0);
d_new(:,:,i) = (d_new(:,:,i) / norm(d_new(:,:,i))).*thresh(:,:,i);
end


end