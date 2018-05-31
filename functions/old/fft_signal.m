function [ k,yts,phi ] = fft_signal( y,Fs, figtitle,fig_on )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 4
  fig_on = 1;
end
N = length(y);
x = (1:N)/Fs;

k = (-N/2:N/2-1)/N*Fs;
yt = 2/N*fftshift(fft(y-mean(y)));
phi = unwrap( angle( fft(y)));
yts = abs(yt);

X = yt;
X2=X;%store the FFT results in another array
%detect noise (very small numbers (eps)) and ignore them
 
threshold = max(abs(X))/100; %tolerance threshold
X2(abs(X)<threshold) = 0; %maskout values that are below the threshold
phi=wrapToPi(atan2(imag(X2),real(X2))); %phase information





[V,I] = max(yts);

% display(['Max = ',num2str(V),', \phi(400)= ',num2str( phi(I)), ' rad'])
if fig_on == 1
% figure() 
    figure('Name',figtitle,'NumberTitle','off')
        subplot(131); plot(x,y); 
            xlabel('Time [s]');ylabel('\epsilon [-]')
        subplot(132); plot(k,yts,'-ko' ); 
            axis([0,1e3,0,max(abs(yts))*1.2])
                xlabel('frequency [Hz]');ylabel('\epsilon [-]')
        subplot(133); plot(k,phi,'-ko' ); 
%         min(phi)
            axis([0,1e3,-180,180])
                xlabel('frequency [Hz]');ylabel('phase \phi [deg]')
end



end

