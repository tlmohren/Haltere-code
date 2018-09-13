clc;clear all;close all


% calculations for haltere cross setion
% r =1;

r2 = 150;
r1 = 50;
beta = 0.312; % roughly

Jratio = 50; % rougly used right now 
% Jratio = 77; % rougly used right now 

C1 = 3*pi*(r2^4-r1^4);
C2 = 1/Jratio*pi/4*(r2^4-r1^4)/beta;
% C2 = 1/Jratio*pi/2*(r2^4-r1^4)/beta;
A = sqrt(C1/C2)

b = (C1/A^3)^(1/4)
h = A*b

% currently in play: 
display('in simulations, an aspect ratio of 12 is used')
A_used = 486.87/40.573 

% 
% A_used = 10; 
b = (3*pi*(r2^4-r1^4) / A_used^3)^(1/4)
h = b*A_used
J_ratio = pi/4*(r2^4-r1^4)/ (beta*A_used*b^4)
% 

%% 
for A = 1:20
   
    b = (3*pi*(r2^4-r1^4) / A^3)^(1/4);
    h = b*A_used;
    J_ratioVec(A) = pi/4*(r2^4-r1^4)/ (beta*A*b^4);
end

figure()
plot(1:length(J_ratioVec),J_ratioVec)
hold on
scatter(12,J_ratioVec(12) )
xlabel('A')
ylabel('J ratio')
legend('','A = 12')