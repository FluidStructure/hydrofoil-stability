function [yfft,freq]=fefft(y,t)

%--------------------------------------------------------------------------
%  Purpose:
%     This function subroutine calculates Fast Fourier Transform (FFT)
%     using the time domain signal. The time domain data are provided with
%     coorseponding time interval.
%  
%  Synopsis:
%     [yf, freq]=fefft(y,t)
%   
%  Variable Description:
%     Input parameters - y : Time domain data n by 1
%                        t : Time interval for y  of n by 1 size
%     Output - yf : Absolute value of FFT of the time domain data y	
%              freq : Frequency axis values
%
%  Notes: 
%    The number of data points for y should be power of 2, and 
%    truncation is needed to achieve the requirement
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Compute number of data points and sampling time interval
%------------------------------------------------------------------------

ntime=max(size(t));
dt=(t(1,ntime)-t(1,1))/ntime;

%-----------------------------------------------------------------------
%  Extract data points at the power of 2. Truncate extra data points
%  so that the final number of data points is in the power of two and
%  also as close as possible to the given number of data points
%-----------------------------------------------------------------------

N=fix(log10(ntime)/log10(2))


%  Calculate FFT of the time domain data and take absolute values of the result

yfft=fft(y(1:2^N,:));

yfft=abs(yfft(1:2^N/2,:))*dt;

%---------------------------------------------------------------------
%  Set up the frequency scale from the given sampling interval.
%  Apply the Nyquist criterion to establish the maximum frequency.
%---------------------------------------------------------------------

freq0=0; 

freqf= (1/dt)/2;        %  Maximum or final frequency value

df=freqf/(2^N/2);           %  Frequency interval

freq=0:df:freqf-df;     %  Frequency axis values

%----------------------------------------------------------------------