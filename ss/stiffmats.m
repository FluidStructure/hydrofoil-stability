function [KDeta,KDtheta,KReta,KRtheta] = stiffmats(EImat,Lmat)

%-----------------------
% Add some zeros to the end of the matrices (for later for loops)
EImatT = [0;EImat;0];
LmatT = [1;Lmat;1];

%-----------------------
% Construct the stiffness matrices
nels = length(EImat);

% Transverse displacment forcing in vertical direction 
KDeta = zeros(nels+2);
for rw = (1:nels)+1
    KDeta(rw,rw-1) = (12*EImatT(rw)/(LmatT(rw).^3));
    KDeta(rw,rw) = -1*((12*EImatT(rw)/(LmatT(rw).^3)) + (12*EImatT(rw+1)/(LmatT(rw+1).^3)));
    KDeta(rw,rw+1) = (12*EImatT(rw+1)/(LmatT(rw+1).^3));
end
KDeta(1,:) = [];
KDeta(:,1) = [];
KDeta(nels+1,:) = [];
KDeta(:,nels+1) = [];

% Axial rotation forcing in vertical direction
KDtheta = zeros(nels+2);
for rw = (1:nels)+1
    KDtheta(rw,rw-1) = (6*EImatT(rw)/(LmatT(rw).^2));
    KDtheta(rw,rw) = ((6*EImatT(rw)/(LmatT(rw).^2)) - (6*EImatT(rw+1)/(LmatT(rw+1).^2)));
    KDtheta(rw,rw+1) = (-6*EImatT(rw+1)/(LmatT(rw+1).^2));
end
KDtheta(1,:) = [];
KDtheta(:,1) = [];
KDtheta(nels+1,:) = [];
KDtheta(:,nels+1) = [];

% Vertical displacement forcing in axial rotation
KReta = zeros(nels+2);
for rw = (1:nels)+1
    KReta(rw,rw-1) = (-6*EImatT(rw)/(LmatT(rw).^2));
    KReta(rw,rw) = ((6*EImatT(rw)/(LmatT(rw).^2)) - (6*EImatT(rw+1)/(LmatT(rw+1).^2)));
    KReta(rw,rw+1) = (6*EImatT(rw+1)/(LmatT(rw+1).^2));
end
KReta(1,:) = [];
KReta(:,1) = [];
KReta(nels+1,:) = [];
KReta(:,nels+1) = [];

% Axis rotation forcing in axial rotation
KRtheta = zeros(nels+2);
for rw = (1:nels)+1
    KRtheta(rw,rw-1) = (-2*EImatT(rw)/(LmatT(rw)));
    KRtheta(rw,rw) = ((-4*EImatT(rw)/(LmatT(rw))) - (4*EImatT(rw+1)/(LmatT(rw+1))));
    KRtheta(rw,rw+1) = (-2*EImatT(rw+1)/(LmatT(rw+1)));
end
KRtheta(1,:) = [];
KRtheta(:,1) = [];
KRtheta(nels+1,:) = [];
KRtheta(:,nels+1) = [];