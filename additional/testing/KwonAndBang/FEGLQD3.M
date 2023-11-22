function [point3,weight3]=feglqd3(nglx,ngly,nglz)

%-------------------------------------------------------------------
%  Purpose:
%     determine the integration points and weighting coefficients
%     of Gauss-Legendre quadrature for three-dimensional integration
%
%  Synopsis:
%     [point3,weight3]=feglqd3(nglx,ngly,nglz)
%
%  Variable Description:
%     nglx - number of integration points in the x-axis
%     ngly - number of integration points in the y-axis
%     nglz - number of integration points in the z-axis
%     point3 - vector containing integration points   
%     weight3 - vector containing weighting coefficients 
%-------------------------------------------------------------------

%  determine the largest one between nglx and ngly

   if nglx > ngly
     if nglx > nglz
       ngl=nglx;
     else
       ngl=nglz;
     end
   else
     if ngly > nglz  
       ngl=ngly;
     else
       ngl=nglz;
     end
   end

%  initialization

   point3=zeros(ngl,3);
   weight3=zeros(ngl,3);

%  find corresponding integration points and weights

 [pointx,weightx]=feglqd1(nglx);     % quadrature rule for x-axis
 [pointy,weighty]=feglqd1(ngly);     % quadrature rule for y-axis
 [pointz,weightz]=feglqd1(nglz);     % quadrature rule for z-axis

%  quadrature for two-dimension

 for intx=1:nglx                     % quadrature in x-axis
   point3(intx,1)=pointx(intx);
   weight3(intx,1)=weightx(intx);
 end

 for inty=1:ngly                     % quadrature in y-axis
   point3(inty,2)=pointy(inty);
   weight3(inty,2)=weighty(inty);
 end
  
 for intz=1:nglz                     % quadrature in z-axis
   point3(intz,3)=pointz(intz);
   weight3(intz,3)=weightz(intz);
 end

