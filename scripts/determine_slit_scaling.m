


% Match I and A 

clc;clear all;close all


r1 = 50
r2 = 150
slitWidth = 10;
h = (r2-r1);
b = slitWidth;
d = (r2+r1)/2;

I = pi/2*(r2^4 - r1^4)-...
    ( b* (h)^3 / 12   + b *h * d/2);

A = pi*(r2^2- r1^2) - b*h;


r3 = sqrt(   pi/(2*A) * ( I*2/pi- A^2/pi^2)   )
r4 = sqrt( A/pi+r3^2 )


A = A
A2 = pi*(r4^2- r3^2)
I = I
I2 = pi/2*(r4^4 - r3^4)



%% 
% D = 0.66;
% f = 40;
% mu = 0;
% % lambda = 0.16;
% 
% lambda = 2*D/f 