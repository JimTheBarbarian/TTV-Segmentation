%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% project an n-dim vector y to the simplex Dn
% Dn = { x : x n-dim, 1 >= x >= 0, sum(x) = 1}
%
%The original algorithm is provided by
% (c) Xiaojing Ye
% xyex19@gmail.com
%
% Algorithm is explained as in the linked document
% http://arxiv.org/abs/1101.6081
% or
% http://ufdc.ufl.edu/IR00000353/
%
% Jan. 14, 2011.
%
%It has been modified so that it runs more quickly.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = projsplx(y)


m = size(y,2); 
n = size(y,1);
s = sort(y, 2, 'descend'); tmpsum = zeros(n,1);

ind = zeros(n,1);
new_ind = zeros(n,1);
tmax = zeros(n,1);
for ii = 1:m-1
    tmpsum = tmpsum + s(:,ii).*(1-ind);
    tmax = tmax.*ind+(1-ind).*(tmpsum - 1)/ii;

    new_ind = new_ind + (tmax>=s(:,ii+1));
    ind = new_ind>0;
    
end

tmax2 = (tmpsum + s(:,m) -1)/m;

tmax3 = ind.*tmax + tmax2.*(1-ind);

x = max(bsxfun(@minus, y, tmax3), 0);

return;
end