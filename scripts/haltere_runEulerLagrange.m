%----------------------
% fun haltere Euler Lagrange Function
% Annika Eberle 2016? 
% 
% Modified: Thomas Mohren 2017 - October -30
%-----------------------

% Run haltere Ode

clc;clear all;close all
factor = 9;

addpath('functions')
%% 
omega1 = 10;  % out of plane rotation 
% omega1 = 0;  % out of plane rotation 
omega2 = 0;
omega3 = 0;  % flapping angle 
% m1 = 0.5;  m2 = 0.25; m3 = 0.25 ;
m1 = 0.9; m2 = 0.05; m3 = 0.05;
% m1 = 1; m2 = 0; m3 = 0;   % central mass
R = 0.01;
tic 
% [T,gammain,thetaout, phiout] = haltereODEsolver(factor,omega1,omega2,omega3,m1,m2,m3,R);
% [T,gammain,thetaout, phiout] = haltereODEsolver_AE(factor,omega1,omega2,omega3,m1,m2,m3,R);
[T,gammain,thetaout, phiout] = haltereODEsolver_TMorig(factor,omega1,omega2,omega3,m1,m2,m3,R);
 
runtime = toc 
 % gammain = flapping
 % thetaout = bending 
 % phiout = twisting 
 %% 
 
save(['data' filesep 'haltere_eulerLagrange_sidebulbs_Om' num2str(omega1) ' .mat'],'T','gammain', 'thetaout', 'phiout')
 
 
 %% Post-process results 
% gammain = pi/2*sin(2*pi*fflap*T);
% thetaout = Y(:,1); 
% phiout = Y(:,2); 
% %%

subT = floor(length(T)/2):length(T);

figure()
plot(gammain,thetaout,'LineWidth',1);hold on
plot(gammain(subT),thetaout(subT))
xlabel('flapping angle');ylabel('bending angle')
% axis([-2 2 -0.015 0.015])

figure()
plot(gammain,phiout,'LineWidth',1);hold on
plot(gammain(subT),phiout(subT))
xlabel('flapping angle');ylabel('twisting angle')
% axis([-2 2 -0.015 0.015])

subT = floor(length(T/2)):length(T);

figure()
    subplot(311)
        plot(T,gammain)
        xlabel('Time [s]')
        ylabel('flapping angle')
    subplot(312)
        plot(T,thetaout)
        xlabel('Time [s]')
        ylabel('bending angle')
    subplot(313)
        plot(T,phiout)
        xlabel('Time [s]')
        ylabel('twist angle')

%%
