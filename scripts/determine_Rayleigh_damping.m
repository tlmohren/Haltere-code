% Rayleigh damping parameters to set 
clc;clear all;close all

%Damping, and frequencies to set damping at 
dampValue = 0.66;
w1 = 40;
w2 = 200;

%% 
wMat = [1/w1 , w1; ...
        1/w2 , w2];
eta = [1;1]*dampValue;

dampConstants = inv(wMat)*2*eta 
% dampConstants(1) = 0
check_damping = 0.5*wMat * dampConstants


w = 1:1.2e3;
figure(); 
subplot(311); hold on
    plot(w,1./(2*w)*dampConstants(1))
    plot(w,w./(2)*dampConstants(2))
    plot(w,1./(2*w)*dampConstants(1)+ w./(2)*dampConstants(2))
    plot([0,2e3],[1,1]*eta(1),'--k')
    legend('mass damping','stiffness damping','mass + stiffness','desired damping','Location','NorthEastOutside')
    axis([10,1100,0,1])

subplot(312); hold on
    plot(w,1./(2*w)*dampConstants(1))
    plot(w,w./(2)*dampConstants(2))
    plot(w,1./(2*w)*dampConstants(1)+ w./(2)*dampConstants(2))
    plot([0,2e3],[1,1]*eta(1),'--k')
    legend('mass damping','stiffness damping','mass + stiffness','desired damping','Location','NorthEastOutside')
    axis([10,200,0,1])
    
    
subplot(313); hold on
    plot(w,1./(2*w)*dampConstants(1))
    plot(w,w./(2)*dampConstants(2))
    plot(w,1./(2*w)*dampConstants(1)+ w./(2)*dampConstants(2))
    plot([0,2e3],[1,1]*eta(1),'--k')
    legend('mass damping','stiffness damping','mass + stiffness','desired damping','Location','NorthEastOutside')
    axis([1000,1200,0,1])
    