function [PANEND,PLENGT,CP]=NACAsec(NumU,NumL,Naca,Series)

% Old implementation:
% function [x,y]=NACAsec(NumU,NumL,Naca,Series)
%
% This function returns the x and y coordinates for a NACA section
%
% [PANEND,PLENGT,CP]=NACAsec(NumU,NumL,Naca,Series)
%
% NumU = number of panels on the upper sufrace
% NumL = number of panels on the lower surface
% Naca = the 4 or 5 digit NACA number
% Series = 4 or 5 depending on the NACA series
%
% PANEND = [x y] coordinates of node points
% PLENGT = length of each panel
% CP = [x y] coordinates of control points
%
% NOTE: this function does seem to have some peculiarities
% 1 - The panel coordinates definitions start from the trailing edge and
% then loop around, going first down along the bottom edge and then up over
% the top edge and back to the trailing edge again.
% 2 - If the NumU is set to be larger than about 80 then a little tiny
% kink forms at the final panel of the top surface

%Initialize variables
	%Series=4;
	XUpper=zeros(NumU,1);
   CamberUpper=XUpper;
   dCamdxUpper=XUpper;
   BetaUpper=XUpper;
	XLower=zeros(NumL,1);
	CamberLower=XLower;
   dCamdxLower=XLower;
   BetaLower=XLower;
   
	%Convert NACA 4 digit to fraction of chord values.
	MaxCamber=floor(Naca/1000);
	PtMaxCamber=floor(Naca/100)-MaxCamber*10;
	MaxThickness=Naca-MaxCamber*1000-PtMaxCamber*100;
	MaxCamber=MaxCamber*.01;
	PtMaxCamber=PtMaxCamber*.1;
	MaxThickness=MaxThickness*.01;
	if MaxCamber>=10
		PtMaxCamber=0.2025;
		MaxCamber=2.6595*PtMaxCamber^3;
		Series=5;
	end

	%Find the x and y coordinates of the airfoil.
	for i=1:NumU
		XUpper(i)=.5*(1-cos(pi*(i-1)/NumU));
	end
	for i=1:NumL
		XLower(i)=.5*(1+cos(pi*(i-1)/NumL));
    end
	

   %Evaluate thickness and camber
	%for naca 4- or 5-digit airfoil
	ThickUpper=5*MaxThickness*(.2969*sqrt(XUpper)-XUpper.*(.126+XUpper.*(.35372-XUpper.*(.2843-XUpper*.1015))));
	ThickLower=5*MaxThickness*(.2969*sqrt(XLower)-XLower.*(.126+XLower.*(.35372-XLower.*(.2843-XLower*.1015))));
	if Series==4
		for i=1:NumU
			if XUpper(i)<PtMaxCamber
				CamberUpper(i)=(MaxCamber/PtMaxCamber^2)*(2*PtMaxCamber-XUpper(i))*XUpper(i);
				dCamdxUpper(i)=(2*MaxCamber/PtMaxCamber^2)*(PtMaxCamber-XUpper(i));
			else
				CamberUpper(i)=(MaxCamber/(1-PtMaxCamber)^2)*(1+XUpper(i)-2*PtMaxCamber)*(1-XUpper(i));
				dCamdxUpper(i)=(2*MaxCamber/(1-PtMaxCamber)^2)*(PtMaxCamber-XUpper(i));
			end
			BetaUpper=atan(dCamdxUpper);
		end
		for i=1:NumL
			if XLower(i)<PtMaxCamber
				CamberLower(i)=(MaxCamber/PtMaxCamber^2)*(2*PtMaxCamber-XLower(i))*XLower(i);
				dCamdxLower(i)=(2*MaxCamber/PtMaxCamber^2)*(PtMaxCamber-XLower(i));
			else
				CamberLower(i)=(MaxCamber/(1-PtMaxCamber)^2)*(1+XLower(i)-2*PtMaxCamber)*(1-XLower(i));
				dCamdxLower(i)=(2*MaxCamber/(1-PtMaxCamber)^2)*(PtMaxCamber-XLower(i));
			end
			BetaLower=atan(dCamdxLower);
		end
	else
		for i=1:NumU
         if XUpper(i)<PtMaxCamber
				W=XUpper(i)/PtMaxCamber;
				CamberUpper(i)=MaxCamber*W*((W-3)*W+3-PtMaxCamber);
				dCamdxUpper(i)=MaxCamber*3*W*(1-W)/PtMaxCamber;
			else
				CamberUpper(i)=MaxCamber*(1-XUpper(i));
				dCamdxUpper(i)=-MaxCamber;
			end
			BetaUpper=atan(dCamdxUpper);
		end
		for i=1:NumL
			if XLower(i)<PtMaxCamber
				W=XLower(i)/PtMaxCamber;
				CamberLower(i)=MaxCamber*W*((W-3)*W+3-PtMaxCamber);
				dCamdxLower(i)=MaxCamber*3*W*(1-W)/PtMaxCamber;
			else
				CamberLower(i)=MaxCamber*(1-XLower(i));
				dCamdxLower(i)=-MaxCamber;
			end
			BetaLower=atan(dCamdxLower);
		end
   end
   x=XLower+ThickLower.*sin(BetaLower);
   y=CamberLower-ThickLower.*cos(BetaLower);
   x=[x;XUpper-ThickUpper.*sin(BetaUpper);x(1)];
   y=[y;CamberUpper+ThickUpper.*cos(BetaUpper);y(1)];
   
% Shift the variables so that the first panel is at the leading edge, not the trailing
% firsthalf = NumL+1:NumU+NumL;
% xt = x(firsthalf);yt = y(firsthalf);
% secondhalf = 1:NumL;
% xt = [xt;x(secondhalf)];yt = [yt;y(secondhalf)];
% xt = [xt;x(NumL+1)];yt = [yt;y(NumL+1)];
% x = xt;y = yt;

   PANEND = [x y];
   N = size(x,1);
   PLENGT = sqrt((x(2:N)-x(1:N-1)).^2 + (y(2:N)-y(1:N-1)).^2);
   CP = [0.5*(x(2:N)+x(1:N-1)) 0.5*(y(2:N)+y(1:N-1))];
   
%********************************************************** 
function [xmid,Cp,Cl,Cm,Cd]=FindPandC(x,y,alpha)
	%This function finds a coefficient matrix using the information
   %given by the x and y coordinates of the body shape., and finds
   %the solution matrix.  It then uses that to find the pressures
   %along the surface and the lift,moment,drag coefficients.
   %Input  x     - x values of contour
   %       y     - y values of contour
   %       alpha - angle of attack
   %Output sol   - solution to A*sol=b
   
   %Initialize
   nodetotal=length(x);
   A=zeros(nodetotal);
   b=zeros(nodetotal,1);
   c=A;
   
   %Set variables for the use of making the coefficient matrix.
   dx=x(2:nodetotal)-x(1:nodetotal-1);
   dy=y(2:nodetotal)-y(1:nodetotal-1);
   sintheta=dy./sqrt(dx.*dx+dy.*dy);
   costheta=dx./sqrt(dx.*dx+dy.*dy);
   xmid=(x(1:nodetotal-1)+x(2:nodetotal))/2;
   ymid=(y(1:nodetotal-1)+y(2:nodetotal))/2;
   
   dxj=xmid-x(1:nodetotal-1);
   dxjp=xmid-x(2:nodetotal);
   dyj=ymid-y(1:nodetotal-1);
   dyjp=ymid-y(2:nodetotal);
   flog=.5*log((dxjp.*dxjp+dyjp.*dyjp)./(dxj.*dxj+dyj.*dyj));
   ftan=atan2(dyjp.*dxj-dxjp.*dyj,dxjp.*dxj+dyjp.*dyj);
   
   %Create A matrix
   for i=1:nodetotal-1
      for j=1:nodetotal-1
      	flog=0;
      	ftan=pi;
      	if j~=i
		      dxj=xmid(i)-x(j);
		      dxjp=xmid(i)-x(j+1);
		      dyj=ymid(i)-y(j);
		      dyjp=ymid(i)-y(j+1);
		      flog=.5*log((dxjp*dxjp+dyjp*dyjp)/(dxj*dxj+dyj*dyj));
      		ftan=atan2(dyjp*dxj-dxjp*dyj,dxjp*dxj+dyjp*dyj);
      	end
      	ctimtj=costheta(i)*costheta(j)+sintheta(i)*sintheta(j);
      	stimtj=sintheta(i)*costheta(j)-sintheta(j)*costheta(i);
      	A(i,j)=(.5/pi)*(ftan*ctimtj+flog*stimtj);
      	c(i,j)=(.5/pi)*(flog*ctimtj-ftan*stimtj);
         A(i,nodetotal)=A(i,nodetotal)+c(i,j);
      	if i==1|i==nodetotal-1
			%If ith panel touches trailing edge, add contribution
			%to kutta condition
	      	A(nodetotal,j)=A(nodetotal,j)-c(i,j);
   	   	A(nodetotal,nodetotal)=A(nodetotal,nodetotal)+A(i,j);
			end
   	end
		%Fill in known sides
		b(i)=sintheta(i)*cos(alpha)-costheta(i)*sin(alpha);
	end
	b(nodetotal)=-(costheta(1)+costheta(nodetotal-1))*cos(alpha)-(sintheta(1)+sintheta(nodetotal-1))*sin(alpha);
   sol=inv(A)*b;
   
	%Use solution to find Pressures along airfoil
	Cp=zeros(nodetotal-1,1);
	for i=1:nodetotal-1
		vtang=cos(alpha)*costheta(i)+sin(alpha)*sintheta(i);
		for j=1:nodetotal-1
         vtang=vtang-c(i,j)*sol(j)+sol(nodetotal)*A(i,j);
		end
		Cp(i)=1-vtang^2;
	end
   
   %
   Cfx=Cp'*dy;
   Cfy=-Cp'*dx;
   Cm=Cp'*(dx.*xmid+dy.*ymid);
   Cd=Cfx*cos(alpha)+Cfy*sin(alpha);
   Cl=Cfy*cos(alpha)-Cfx*sin(alpha);
%    --Message-Boundary-8915
% Content-type: text/plain; charset=US-ASCII
% Content-disposition: inline
% Content-description: Attachment information.
% 
% The following section of this message contains a file attachment
% prepared for transmission using the Internet MIME message format.
% If you are using Pegasus Mail, or any another MIME-compliant system,
% you should be able to save it or view it from within your mailer.
% If you cannot, please ask your system administrator for assistance.
% 
%    ---- File information -----------
%      File:  thinairfoil.m
%      Date:  11 Oct 1999, 16:54
%      Size:  6179 bytes.
%      Type:  Text
