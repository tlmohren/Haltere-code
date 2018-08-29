function [T,dataOut] = haltereODEsolver_TM( omega, m )
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
% m(1), m(2), and m(3) are the mass of each point mass (again, 1 is at the 
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
run('config.m')

w = 2*pi*fflap*res_factor  %resonant frequency

if sum(m) ~= M
    error('myApp:argChk','mass not conserved')
end

% setup symbols for equations 
syms t  
dz = sym('dz(t)');
dx = sym('dx(t)');

% gamma = pi/2*sin(2*pi*fflap*t);
phix = pi/2*sin(2*pi*fflap*t);

% generate rotation matrices 
rot1 = [cos(dx),sin(dx),0;-sin(dx), cos(dx),0;0 0 1];
rot2 = [1 0 0; 0 cos(dz) -sin(dz); 0 sin(dz) cos(dz)];
rot3 = [cos(phix) 0 sin(phix); 0 1 0; -sin(phix) 0 cos(phix)]; 

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
% T = 0.5*m(1)*(dr1_i(1)^2+dr1_i(2)^2+dr1_i(3)^2) ...
%     + 0.5*m(2)*(dr2_i(1)^2+dr2_i(2)^2+dr2_i(3)^2) ...
%     + 0.5*m(3)*(dr3_i(1)^2+dr3_i(2)^2+dr3_i(3)^2);
m1 = m(1); m2 = m(2); m3 = m(3);
T = 0.5*m1*(dr1_i(1)^2+dr1_i(2)^2+dr1_i(3)^2) ...
    + 0.5*m2*(dr2_i(1)^2+dr2_i(2)^2+dr2_i(3)^2) ...
    + 0.5*m2*(dr3_i(1)^2+dr3_i(2)^2+dr3_i(3)^2);

% calculate potential energy 
kbend = r^2*M*w^2; 
ktwist = kbend/(3/(1+nu)); 

V = 0.5*kbend*r1(3)^2 +0.5*kbend*r1(2)^2 + 0.5*ktwist*dx^2; 

%% Evaluate Lagrange's equation 

% generate Lagrangian
L = T-V; 

% evaluate Lagrange's equation 
eqnt = diff(diffDepVar(L,diff(dz,t)),t)- diffDepVar(L,dz);

eqn2t = (-(eqnt-diffDepVar(eqnt,diff(dz,t,2))*diff(dz,t,2))-2*w*zeta*diff(dz,t))/diffDepVar(eqnt,diff(dz,t,2));

eqn3t = subs(eqn2t,diff(dz,t),'y(3)');
eqn4at = subs(eqn3t,dz,'y(1)');
eqn4bt = subs(eqn4at,diff(dx,t),'y(4)');
eqn4t = subs(eqn4bt,dx,'y(2)');
 
eqng = diff(diffDepVar(L,diff(dx,t)),t)- diffDepVar(L,dx);

eqn2g = (-(eqng-diffDepVar(eqng,diff(dx,t,2))*diff(dx,t,2))-2*w*zeta*diff(dx,t))/diffDepVar(eqng,diff(dx,t,2));

% substitute variables for differential equation 
eqn3g = subs(eqn2g,diff(dx,t),'y(4)');
eqn4ag = subs(eqn3g,dx,'y(2)');
eqn4bg = subs(eqn4ag,diff(dz,t),'y(3)');
eqn4g = subs(eqn4bg,dz,'y(1)');

% set equation 4 to zero if eccentric masses are zero 
if (m(2) == 0) && (m(3) == 0)
    eqn4g = 0; 
end 

%% Generate function with equations for ODEs

% Determine function location
if exist(['myODE_twist.m']) == 2
    display('PlateODE exists, deleting now ')
end
clear(['myODE_twist.m'])
pause(1)
%Delete prior PlateODE.m file     
% delete(fullfile(pwd,'myODE_twist.m'))
% clear myODE_twist

% Build ode file  
% fid = fopen(['functions' filesep 'myODE_twist.m'],'w')
fid = fopen(['myODE_twist.m'],'w')
fprintf(fid,'function dy = myODE_twist(t,y)');
fprintf(fid,'\n dy(1,1) =y(3)');
fprintf(fid,';');
fprintf(fid,'\n dy(2,1) =y(4)');
fprintf(fid,';');
fprintf(fid,'\n dy(3,1) =');
fprintf(fid,strrep(char(eqn4t),'matrix',''));
fprintf(fid,';');
if (m(2) ==0) && (m(3) == 0) 
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
while exist(['myODE_twist.m'], 'file') ~= 2 && iter<5
    pause(2)
    iter = iter+1; 
end 
%% Solve ODEs

% set parameters for ODE solver
sampfreq = 1000; %sampling frequncy for ODE solver
% cycles = 40; %number of cycles to be solved

init = zeros(4,1); 
options = odeset('RelTol',1e-4);
% teval = 0:1/sampfreq:1/fflap*cycles;
% teval = 0:1/sampfreq:1;
teval = 0:1/sampfreq:0.35;
[T,Y] = ode45(@myODE_twist,teval,init,options); 

%% Post-process results 
dataOut(:,1) = pi/2*sin(2*pi*fflap*T);
dataOut(:,2) = Y(:,1);
dataOut(:,3) = Y(:,2);
% gammain = pi/2*sin(2*pi*fflap*T);
% thetaout = Y(:,1); 
% phiout = Y(:,2); 