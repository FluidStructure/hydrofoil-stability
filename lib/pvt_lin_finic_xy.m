function [ua,ub,va,vb] = pvt_lin_finic_xy(elemat)

%
% This function works out the velocity at some
% evaluation points by the POINT VORTEX sheet 
% elements with properties given in elemat.
% 
% FORMAT:
%  [ua,ub,va,vb] = pvt_lin_finic_xy(elemat)
%
% WHERE:
%
% 
%

% Get the evaluation points
xi = elemat(:,1);
yi = elemat(:,2);

% Get information from "elemat"
xj = elemat(:,1);
yj = elemat(:,2);
eleangj = elemat(:,3);
elewidj = elemat(:,4);

% Put the vectors into a form where bulk addition and subtraction can occur
[xj,xi] = meshgrid(xj,xi);
[yj,yi] = meshgrid(yj,yi);
eleangj = meshgrid(eleangj,ones(size(xj,1),1));
elewidj = meshgrid(elewidj,ones(size(xj,1),1));

% Calculate the relative x and y distances from element "j"
rir = (((yi - yj)).*(cos(eleangj))-((xi - xj)).*(sin(eleangj)));
rur = (((xi - xj)).*(cos(eleangj))+((yi - yj)).*(sin(eleangj)));

tmp1 = ((abs(rir)<1e-10)&(abs(rur)<1e-10));
tmp2 = atan2(rir,(rur - (elewidj/2)));
tmp2 = tmp2 - atan2(rir,(rur + (elewidj/2)));
tmp2 = tmp2.*(1-tmp1);
tmp2 = tmp2 + tmp1.*pi;

% Calculate the ua velocity parts
tmp1 = log((rur - (elewidj/2)).^2 + rir.^2);
tmp1 = tmp1 - log((rur + (elewidj/2)).^2 + rir.^2);
tmp1 = tmp1./(4*pi*elewidj);

uar = -1*tmp1.*(rir);
ubr = tmp1.*(rir);
 
tmp1 = tmp2./(2*pi*elewidj);

uar = uar - tmp1.*(rur - (elewidj/2));
ubr = ubr + tmp1.*(rur + (elewidj/2));


% Calculate the va velocity parts
tmp1 = log((rur + (elewidj/2)).^2 + rir.^2);
tmp1 = tmp1 - log((rur - (elewidj/2)).^2 + rir.^2);
tmp1 = tmp1./(4*pi*elewidj);                    % NOTE: Not sure about the 2 or 4? - Correct as is

var = tmp1.*(rur - (elewidj/2));
vbr = -1*tmp1.*(rur + (elewidj/2));

%tmp1 = atan2(rir,(rur - (elewidj/2)));
%tmp1 = tmp1 - atan2(rir,(rur + (elewidj/2)));
tmp1 = tmp2.*(rir./elewidj);
tmp1 = (1/(2*pi)).*(1 - tmp1);                  % NOTE: Not sure about the 2 or 4? - Correct as is

var = var - tmp1;
vbr = vbr + tmp1;

ua = zeros(size(uar));
ub = ua;
va = ua;
vb = va;
for i = 1:size(elemat,1)
    % Calculate the actual global x-component of velocity at "i"
    ua(:,i) = uar(:,i).*cos(elemat(i,3)) - var(:,i).*sin(elemat(i,3));
    ub(:,i) = ubr(:,i).*cos(elemat(i,3)) - vbr(:,i).*sin(elemat(i,3));
    % Calculate the actual global y-component of velocity at "i"
    va(:,i) = var(:,i).*cos(elemat(i,3)) + uar(:,i).*sin(elemat(i,3));
    vb(:,i) = vbr(:,i).*cos(elemat(i,3)) + ubr(:,i).*sin(elemat(i,3));
end

