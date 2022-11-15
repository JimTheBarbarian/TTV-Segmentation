%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function performs the discrete gradient operator in y-direction.
%
%Input:
%   u: function whose Dy will be computed
%Output:
%   d: Dy of u
%Note: Code is from Yifei Lou.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = Dy(u)
%obtain matrix size of u
[rows,cols, channels] = size(u);

%preinitialize matrix
d = zeros(rows,cols, channels);

%compute Dy(u)
d(2:rows,:,:) = u(2:rows,:,:)-u(1:rows-1,:,:);
d(1,:,:) = u(1,:,:)-u(rows,:,:);
return