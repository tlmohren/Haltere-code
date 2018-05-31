%----------------------
% fun haltere Euler Lagrange Function
% Annika Eberle 2016? 
% 
% Modified: Thomas Mohren 2017 - October -30
%-----------------------

% Run haltere Ode

clc;clear all;close all

addpath('functions')
%% 
omega = [10,0,0]; 
m = [1, 0, 0];
% m = [0.9, 0.05,0.05];


[T,dataOut] = haltereODEsolver_TM( omega ,m );
 
 % dataOut(:,1) = flapping
 % dataOut(:,2) = bending 
 % dataOut(:,3) = twisting 
 
 
 
 %% Post-process results 
% dataOut(:,1) = pi/2*sin(2*pi*fflap*T);
% dataOut(:,2) = Y(:,1); 
% dataOut(:,3) = Y(:,2); 
% %%

subT = floor(length(T)/2):length(T);

figure()
plot(dataOut(:,1),dataOut(:,2),'LineWidth',1);hold on
% plot(dataOut(:,1),dataOut(:,2)(subT))
xlabel('flapping angle');ylabel('bending angle')
% axis([-2 2 -0.015 0.015])

figure()
plot(dataOut(:,1),dataOut(:,3),'LineWidth',1);hold on
% plot(dataOut(subT,1)(),dataOut(:,3)(subT))
xlabel('flapping angle');ylabel('twisting angle')
% axis([-2 2 -0.015 0.015])

subT = floor(length(T/2)):length(T);

figure()
    subplot(311)
        plot(T,dataOut(:,1))
        xlabel('Time [s]')
        ylabel('flapping angle')
    subplot(312)
        plot(T,dataOut(:,2))
        xlabel('Time [s]')
        ylabel('bending angle')
    subplot(313)
        plot(T,dataOut(:,3))
        xlabel('Time [s]')
        ylabel('twist angle')


