% euler angles test

clc;clear all;close all


L = 5;
r1 = 1;
m1 = [L,0,0];
m2 = [L,r1,0];
m3 = [0,r1,0];


a1 = pi/8;
a1 = 0;

a2 = pi/8;
% a2 = 0;

a3 = pi/8;
% a3 = 0;%

Xrot1 = [ 1 0 0 ; 0 cos(a1) -sin(a1) ; 0 sin(a1) cos(a1)];
Xrot2 = [ 1 0 0 ; 0 cos(a2) -sin(a2) ; 0 sin(a2) cos(a2)];
Xrot3 = [ 1 0 0 ; 0 cos(a3) -sin(a3) ; 0 sin(a3) cos(a3)];
Yrot1 = [ cos(a1) 0 -sin(a1); 0 1 0; sin(a1) 0 cos(a1) ] ;
Yrot2 = [ cos(a2) 0 -sin(a2); 0 1 0; sin(a2) 0 cos(a2) ] ;
Yrot3 = [ cos(a3) 0 -sin(a3); 0 1 0; sin(a3) 0 cos(a3) ] ;
Zrot1 = [cos(a1) -sin(a1) 0; sin(a1) cos(a1) 0; 0 0 1];
Zrot2 = [cos(a2) -sin(a2) 0; sin(a2) cos(a2) 0; 0 0 1];
Zrot3 = [cos(a3) -sin(a3) 0; sin(a3) cos(a3) 0; 0 0 1];

% R = [cos(a1)*cos(a3)+sin(a1)*sin(a2)*sin(a3) cos(a1)*sin(a3)-sin(a1)*sin(a2)*cos(a3)]


% good 
% R = Yrot1*Xrot2* Zrot3  
% R = Zrot1*Yrot2* Xrot3  
R = Xrot1*Yrot2* Zrot3  





d1 = R * m1';
d2 = R * m2';
d3 = R * m3';

d4 = d1+d3


figure('Position',[300,300,800,500]);
subplot(121) 
    hold on
    plot3( [0,m1(1),m2(1)], [0,m1(2),m2(2)], [0,m1(3),m2(3)],'b')
    scatter3( [m1(1),m2(1)],[m1(2),m2(2)],[m1(3),m2(3)],20,'b','filled')
    axis square;axis equal
    view(50,20)
    axis([-1,5.5,-5.5,5.5,-5.5,5.5])
    % 
    plot3( [0,d1(1),d2(1)], [0,d1(2),d2(2)], [0,d1(3),d2(3)],'r')
    plot3( [0,d3(1)], [0,d3(2)], [0,d3(3)],'r')
    scatter3( [d1(1),d2(1),d3(1),d4(1)],[d1(2),d2(2),d3(2),d4(2)],[d1(3),d2(3),d3(3),d4(3)],20,'r','filled')


subplot(122) 
    hold on
    plot3( [0,m1(1),m2(1)], [0,m1(2),m2(2)], [0,m1(3),m2(3)],'b')
    scatter3( [m1(1),m2(1)],[m1(2),m2(2)],[m1(3),m2(3)],20,'b','filled')
    axis square;axis equal
    view(90,0)
    axis([-1,5.5,-5.5,5.5,-5.5,5.5])
    % 
    plot3( [0,d1(1),d2(1)], [0,d1(2),d2(2)], [0,d1(3),d2(3)],'r')
    plot3( [0,d3(1)], [0,d3(2)], [0,d3(3)],'r')
    scatter3( [d1(1),d2(1),d3(1)],[d1(2),d2(2),d3(2)],[d1(3),d2(3),d3(3)],20,'r','filled')



radiusOr = norm(m2)
radiusTot = norm(d2)
radiusSum = norm(d4)
