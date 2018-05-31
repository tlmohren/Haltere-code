function [T,gammain,thetaout, phiout] = haltereODEsolver_TMorig(factor,omega1,omega2,omega3,m1,m2,m3,R)
%%
% Generates and solves Lagrange's equation to compute the bending and 
% torsion angles for a simple haltere model with three
% masses at the tip (1 is at the center; 2 and 3 are on either side) 
%
% Sample command to execute: [T,gammain,thetaout,phiout] = haltereODEsolver(1,0,0,10,0.8,0.1,0.1,0.1);
%
% INPUTS
%
% factor = multiplier for the resonant frequency (as in Thompson's
% model, the resonant frequency is assumed to be either one or two times 
% the flapping freq
%
% omega1, omega2, and omega3 are the angular rotation rate (assumed
% constant) about the x, y, and z axes, respectively
% 
% m1, m2, and m3 are the mass of each point mass (again, 1 is at the 
% center, 2 and 3 are on either size) 
%
% R = distance from the main mass to the tip masses (2 and 3) 
% 
% OUTPUTS 
%
% T = time vector
%
% gammain = input flapping angle
%
% thetaout = output bending angle
%
% phiout = output torsion angle 


%% Pre-process variables

fflap = 40; %flapping frequency
w = 2*pi*fflap*factor; %resonant frequency
zeta = 2; %damping factor

M = 1; %total mass (normalized)
r = 1; %normalized distance
nu = 0.5; %poisson's ratio 

if m1+m2+m3 ~= M
    error('myApp:argChk','mass not conserved')
end
 
omega = [omega1;omega2;omega3]; %angular velocity vector for constant angular velocity


ttest = 0:0.001:0.35;
    sigprop = [1;10;3];
    sigd = sigprop(1);
    sigc = sigprop(2);
    sign = sigprop(3);
    sigmoid = (sigd.*(2*pi* fflap *ttest).^sign) ./ (sigc+sigd.*(2*pi*fflap*ttest).^sign);
    
%     figure();
%     plot(ttest,sigmoid)
% setup symbols for equations 
syms t 
theta = sym('theta(t)');
phi = sym('phi(t)');

%     sigmoid = (sigd.*(2*pi*fflap*t).^sign) ./ (sigc+sigd.*(2*pi*fflap*t).^sign);
    sigmoid = 1;
gamma = pi/2*sin(2*pi*fflap*t)*sigmoid;


% generate rotation matrices 
rot1 = [cos(phi),sin(phi),0;-sin(phi), cos(phi),0;0 0 1];
rot2 = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
rot3 = [cos(gamma) 0 sin(gamma); 0 1 0; -sin(gamma) 0 cos(gamma)]; 

%% Generate postions and velocities

% compute positions
r1 = rot3*rot2*[0;0;r];
r2 = rot3*rot2*rot1*[0;R;r];
r3 = rot3*rot2*rot1*[0;-R;r]; 

% compute velocities 
dr1 = diff(r1,t); 
dr2 = diff(r2,t);
dr3 = diff(r3,t); 
dr1_i = dr1 + cross(omega,r1); 
dr2_i = dr2 + cross(omega,r2);
dr3_i = dr3 + cross(omega,r3); 

% evaluate distances
rmag2 = r1(1)^2+r1(2)^2+r1(3)^2;
r2mag2 = r2(1)^2+r2(2)^2+r2(3)^2; 
r3mag2 = r3(1)^2+r3(2)^2+r3(3)^2; 

%% Compute energies

% compute kinetic energy 
T = 0.5*m1*(dr1_i(1)^2+dr1_i(2)^2+dr1_i(3)^2) ...
    + 0.5*m2*(dr2_i(1)^2+dr2_i(2)^2+dr2_i(3)^2) ...
    + 0.5*m3*(dr3_i(1)^2+dr3_i(2)^2+dr3_i(3)^2);

% calculate potential energy 
kbend = r^2*M*w^2; 
ktwist = kbend/(3/(1+nu)); 

V = 0.5*kbend*r1(3)^2 +0.5*kbend*r1(2)^2 + 0.5*ktwist*phi^2; 

%% Evaluate Lagrange's equation 

% generate Lagrangian
L = T-V; 

% evaluate Lagrange's equation 
eqnt = diff(diffDepVar(L,diff(theta,t)),t)- diffDepVar(L,theta);

eqn2t = (-(eqnt-diffDepVar(eqnt,diff(theta,t,2))*diff(theta,t,2))-2*w*zeta*diff(theta,t))/diffDepVar(eqnt,diff(theta,t,2));

eqn3t = subs(eqn2t,diff(theta,t),'y(3)');
eqn4at = subs(eqn3t,theta,'y(1)');
eqn4bt = subs(eqn4at,diff(phi,t),'y(4)');
eqn4t = subs(eqn4bt,phi,'y(2)');
 
eqng = diff(diffDepVar(L,diff(phi,t)),t)- diffDepVar(L,phi);

eqn2g = (-(eqng-diffDepVar(eqng,diff(phi,t,2))*diff(phi,t,2))-2*w*zeta*diff(phi,t))/diffDepVar(eqng,diff(phi,t,2));

% substitute variables for differential equation 
eqn3g = subs(eqn2g,diff(phi,t),'y(4)');
eqn4ag = subs(eqn3g,phi,'y(2)');
eqn4bg = subs(eqn4ag,diff(theta,t),'y(3)');
eqn4g = subs(eqn4bg,theta,'y(1)');

% set equation 4 to zero if eccentric masses are zero 
if (m2 == 0) && (m3 == 0)
    eqn4g = 0; 
end 

%% Generate function with equations for ODEs
%Delete prior PlateODE.m file     
delete(fullfile(pwd,'myODE_twist.m'))
clear myODE_twist
fid = fopen('myODE_twist.m','w');
fprintf(fid,'function dy = myODE_twist(t,y)');
fprintf(fid,'\n dy(1,1) =y(3)');
fprintf(fid,';');
fprintf(fid,'\n dy(2,1) =y(4)');
fprintf(fid,';');
fprintf(fid,'\n dy(3,1) =');
fprintf(fid,strrep(char(eqn4t),'matrix',''));
fprintf(fid,';');
if (m2 ==0) && (m3 == 0) 
    fprintf(fid,'\n dy(4,1) = 0;');
else 
    fprintf(fid,'\n dy(4,1) =');
    fprintf(fid,strrep(char(eqn4g),'matrix',''));
    fprintf(fid,';');
end 
fprintf(fid,'\n t'); 
fclose(fid);

%make sure file saves before solving the ODE 
iter =1; 
while exist('myODE_twist.m', 'file') ~= 2 && iter<5
    pause(2)
    iter = iter+1; 
end 
%% Solve ODEs

% set parameters for ODE solver
% sampfreq = 8000; %sampling frequncy for ODE solver
% cycles = 40; %number of cycles to be solved

init = zeros(4,1); 
options = odeset('RelTol',1e-3);
% teval = 0:1/sampfreq:1/fflap*cycles;
teval = 0:0.001:0.02;
[T,Y] = ode45(@myODE_twist,teval,init,options); 

%% Post-process results 
gammain = pi/2*sin(2*pi*fflap*T);
% t = teval;
% gammain = eval(gamma);
thetaout = Y(:,1); 
phiout = Y(:,2); 
% %%
figure
plot(gammain,thetaout)
axis([-2 2 -0.015 0.015])

%%
figure
plot(gammain,phiout)
axis([-2 2 -0.015 0.015])