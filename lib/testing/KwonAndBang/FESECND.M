function [t_p, t_r, t_s, M_p]=fesecnd(zeta, w_n)
%------------------------------------------------------------------------------------------
%  Purpose:
%     The function subroutine fesecnd.m calculates dynamic characteristics of a
%     typical standard second order system. Transfer function :
%    			       w_n^2
%		H(s)= -------------------------
%			s^2+2*zeta*w_n*s+w_n^2   
%  Synopsis:
%     [t_p, t_r, M_p, t_s]=fsecond(zeta, w_n)
%
%  Variable Description:
%     Input parameters - zeta - system damping ratio,   
%                        w_n  - undamped natural frequency
%
%     Output parameters - t_p - peak time,    t_r  - rise time 
%	   	          t_s - settling time,  M_p - maximum overshoot,     
%------------------------------------------------------------------------------------------

w_d=sqrt(1-zeta^2)*w_n;				     % Calculate undamped natural frequency

t_p=pi/w_d;						% Calculate peak time

t_r=atan2(sqrt(1-zeta^2), -zeta);			% Calculate rise time

t_s=4/zeta/w_n;						% Calculate settling time

M_p=exp(-zeta*pi/sqrt(1-zeta^2));			% Calculate maximum overshoot

%------------------------------------------------------------------------------------------

			 