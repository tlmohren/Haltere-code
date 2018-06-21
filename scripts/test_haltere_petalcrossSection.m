

clc;clear all;close all

A_sawtooth = 0:0.1:0.3;
% A_sawtooth = 0.3;
n_sawtooth = 1;

n = 50;
theta = linspace(0,pi*0.5,n);


sawtooth = zeros(4,n);
for j = 1:n_sawtooth
    for k = 1:2:( (j)*2 )
        An = 4* (1-(-1)^k) / (pi^2*k^2);
        sawtooth(j,:) = sawtooth(j,:)+ An*cos( 4*k*theta);
    end
end

figure();
plot(sawtooth')

figure( 'Position',[100,400,1000,400]);

for j = 1:n_sawtooth
    subplot(n_sawtooth,1,j)
    hold on
    for k = 1:length(A_sawtooth)
        r = 1-A_sawtooth(k) + A_sawtooth(k)* sawtooth(j,:) ;
        x = cos(theta).*r;
        y = sin(theta).*r;

        plot(x+k*3,y,'k' )
        plot(-x+k*3,y,'k' )
        plot(x+k*3,-y,'k' )
        plot(-x+k*3,-y,'k' )
    end
    axis equal
    axis off
end

