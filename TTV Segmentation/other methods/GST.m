function [channel] = GST(w,r,p,J)
%GST Performs generalized soft thresholding on a vector w based upon r,p,
%and number of iteration J. 
%  w: The vector that will be soft thresholded. Can be a matrix. 
%  r: A parameter that controls the threshold. The threshold increases as
%  least linearly in r. 
%  p: as in ell^p
%  J: The number of iterations that our iterative solver will run for.
[row,col,~] = size(w);

tau = (2 * r * (1-p))^(1/(2-p)) + r * p * ( 2 * r * (1-p))^((p-1)/(2-p));
w = w(:);
X = abs(w(:));
channel = abs(X); % we will vectorize this subproblem

for i = 1:J
    channel = X - (r * p * channel.^(p-1));
end
channel(X < tau) = 0;
channel = sign(w(:)).*channel;
channel= reshape(channel,row,col,1);

end