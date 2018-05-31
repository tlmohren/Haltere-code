% plots effect of rotation on strain at haltere sides 
clc;clear all;close all
addpathFolderStructureHaltere()

% FEA(1).name = 'Haltere_CraneFly_EllipsoidVer_Om0';
FEA(1).name = 'Haltere_CraneFly_EllipsoidHor_Om0';
% FEA(1).name = 'Haltere_CraneFly_EllipsoidVerWithBulbBoundary_Om10';

for j =  1:length(FEA)
    tic
    [FEA(j).xyz, FEA(j).strain, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'eXX' });
    [FEA(j).xyz, FEA(j).deform, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'u2','v2','w2'});
    toc 
end

% Determine Circle locations
circleDistance = 300;               % distance from base to haltere 
circleRadius = 150;                 % radius of haltere   

% xMatch = find( FEA(1).xyz(:,1) == circleDistance );
% yMatch = find( abs( FEA(1).xyz(:,2) ) == circleRadius );
% zMatch = find( abs( FEA(1).xyz(:,3) ) == circleRadius );

% 
% rMatch = find(  sqrt( FEA(1).xyz(:,2).^2 + FEA(1).xyz(:,3).^2) == circleRadius );
% xMatch = find( FEA(1).xyz(:,1) == 0  );
%% 
% close all
% t_ind = 100;
% % xMatch = find( FEA(1).xyz(:,1) == 0  );
% rMatch = find(  round(sqrt( FEA(1).xyz(:,2).^2 + FEA(1).xyz(:,3).^2),1) == circleRadius );
% xMatch = find( round(FEA(1).xyz(:,1),3) == 150  );
% n_points = intersect( rMatch,xMatch);
% 
% yGoal = -cos( linspace(0,2*pi,17) )*150;
% zGoal = sin( linspace(0,2*pi,17) )*150; 
% xGoal = 0:150:3450;
% 
% [Xm,Ym] = meshgrid(xGoal,yGoal);
% [~,Zm] = meshgrid(xGoal,zGoal);
% Cm = zeros(size(Xm));
% for j = 1:length(xGoal)
%     for k = 1:length(yGoal)
%         x = round(xGoal(j),3);
%         y = round(yGoal(k),3);
%         z = round(zGoal(k),3);
%         xm = round(FEA(1).xyz(:,1),3) == x;
%         ym = round(FEA(1).xyz(:,2),3) == y;
%         zm = round(FEA(1).xyz(:,3),3) == z;
%         ind = find( xm & ym & zm);
% %         for l = 1:length(zGoal)
%         Cm(k,j) = FEA(1).strain(ind,t_ind);
% %         end
%     end
% end
% % Cm = (Cm - min(Cm(:))) / max(Cm(:))*254+1;
% % Cm=Zm
% figure()
% surf(Xm,Ym,Zm,Cm)
% shading flat
% % shading faceted
% axis equal
%%
% clear XmNew YmNew ZmNew CmNew CmNewOm10
% theta = linspace(0,2*pi,17);
% % zGoal = sin( linspace(0,2*pi,17) )*150; 
% % xGoal = 0:150:3450;
% for j = 1:length(theta)
%     yGoal = -cos( theta(j) )*150;
%     zGoal = sin( theta(j) )*150; 
% %         z = round(zGoal,1);
%         y = round(yGoal,1);
%         z = round(zGoal,1);
%         y = round(yGoal,1);
%         ym = round(FEA(1).xyz(:,2),1) == y;
%         zm = round(FEA(1).xyz(:,3),1) == z;
%         indm =  find( (ym & zm) & (FEA(1).xyz(:,1)<=4221 ) );
%         [V,indsort] = sort( FEA(1).xyz(indm,1),'descend' )
%         indm(indsort)
%         XmNew(j,:) =FEA(1).xyz(indm(indsort) ,1);
%         YmNew(j,:) =FEA(1).xyz(indm(indsort),2);
%         ZmNew(j,:) =FEA(1).xyz(indm(indsort),3);
%         CmNew(j,:) = FEA(1).strain( indm(indsort),t_ind);
% end
% %
redPurple = [158,1,66
213,62,79
244,109,67
253,174,97
254,224,139
255,255,191
230,245,152
171,221,164
102,194,165
50,136,189
94,79,16];
% 
strainScheme = colorSchemeInterp(redPurple/255, 500);
% 
% 
% figure()
%     subplot(211)
%     surf(XmNew,YmNew,ZmNew,CmNew)
%     axis equal
%     axis off
%     colormap(strainScheme)
%     colorbar
%     view(70,20)
% % subplot(212)
%     subplot(212)
%     surf(XmNew,YmNew,ZmNew,CmNewOm10-CmNew)
% %     shading flat; % shading faceted
% %     shading faceted
%     axis equal
%     axis off
%     colorbar
%     view(70,20)


%%
theta = linspace(0,2*pi,17);
theta(1) = 2*pi;

t_ind = 100;
FEA(1).xrtheta(:,1) = FEA(1).xyz(:,1);
FEA(1).xrtheta(:,2) = sqrt( FEA(1).xyz(:,2).^2  +  FEA(1).xyz(:,3).^2 );
FEA(1).xrtheta(:,3) = atan2( FEA(1).xyz(:,3), FEA(1).xyz(:,2) ) +pi ;
xDes = [0:150:4800];
for j=1:length(xDes)
    dx = abs(FEA(1).xrtheta(:,1)-xDes(j) );
    dr = abs(FEA(1).xrtheta(:,2)-150 );
%     dr = ( FEA(1).xrtheta(:,2) >= 140 );
	for k = 1:length(theta)
        theta(k)
        da = abs(FEA(1).xrtheta(:,3) - theta(k) );
        J = dx.^2 + dr.^2+ (da*150.^2);
        [V,I] = min(J);
        Xj(j,k) = FEA(1).xyz(I,1);
        Yj(j,k) = FEA(1).xyz(I,2);
        Zj(j,k) = FEA(1).xyz(I,3);
        Cj(j,k) = FEA(1).strain(I,t_ind);
    end
end
figure()
    subplot(211)
    surf(Xj,Yj,Zj,Cj)
    axis equal
    axis off
    colormap(strainScheme)
    colorbar
    view(70,20)
% 
%     
    %% 
theta = linspace(0,2*pi,17);
theta(1) = 2*pi;

t_ind = 100;
FEA(1).xrtheta(:,1) = FEA(1).xyz(:,1);
FEA(1).xrtheta(:,2) = sqrt( FEA(1).xyz(:,2).^2  +  FEA(1).xyz(:,3).^2 );
FEA(1).xrtheta(:,3) = atan2( FEA(1).xyz(:,3), FEA(1).xyz(:,2) ) +pi ;
    dx = FEA(1).xrtheta(:,1) > 4730;
    dr = abs(FEA(1).xrtheta(:,2)-150 );
    
figure(); hold on
%     subplot(211); hold on 
    surf(Xj,Yj,Zj,Cj)
%     shading flat; % shading faceted
    axis equal
    axis off
    colormap(strainScheme)
    colorbar
    view(70,20)
    xc = 5000;
    zc = 0;
    yc = 0; 
    xr = 300; 
    yr = 300;
    zr = 946; 
    n = 16;
[z,y,x] = ellipsoid(zc,yc,xc,zr,yr,xr,n);
    C = zeros(size(x));
surf(x,y,z, C)


%%

theta = linspace(0,2*pi,17);
theta(1) = 2*pi;

% deform_mult = 10;
t_ind = 100;
FEA(1).xrtheta(:,1) = FEA(1).xyz(:,1);
FEA(1).xrtheta(:,2) = sqrt( FEA(1).xyz(:,2).^2  +  FEA(1).xyz(:,3).^2 );
FEA(1).xrtheta(:,3) = atan2( FEA(1).xyz(:,3), FEA(1).xyz(:,2) ) +pi ;
xDes = [0:150:4800];
for j=1:length(xDes)
    dx = abs(FEA(1).xrtheta(:,1)-xDes(j) );
    dr = abs(FEA(1).xrtheta(:,2)-150 );
%     dr = ( FEA(1).xrtheta(:,2) >= 140 );
	for k = 1:length(theta)
        theta(k)
        da = abs(FEA(1).xrtheta(:,3) - theta(k) );
        J = dx.^2 + dr.^2+ (da*150.^2);
        [V,I] = min(J);
        Xj(j,k) = FEA(1).xyz(I,1); 
        Yj(j,k) = FEA(1).xyz(I,2) ;
        Zj(j,k) = FEA(1).xyz(I,3) ;
        Cj(j,k) = FEA(1).strain(I,t_ind);
    end
end
figure()
    subplot(211)
    surf(Xj,Yj,Zj,Cj)
    axis equal
    axis off
    colormap(strainScheme)
    colorbar
    view(70,20)
% 
theta = linspace(0,2*pi,17);
theta(1) = 2*pi;

t_ind = 100;
FEA(1).xrtheta(:,1) = FEA(1).xyz(:,1);
FEA(1).xrtheta(:,2) = sqrt( FEA(1).xyz(:,2).^2  +  FEA(1).xyz(:,3).^2 );
FEA(1).xrtheta(:,3) = atan2( FEA(1).xyz(:,3), FEA(1).xyz(:,2) ) +pi ;
    dx = FEA(1).xrtheta(:,1) > 4730;
    dr = abs(FEA(1).xrtheta(:,2)-150 );
    
figure(); hold on
%     subplot(211); hold on 
    surf(Xj,Yj,Zj,Cj)
%     shading flat; % shading faceted
    axis equal
    axis off
    colormap(strainScheme)
    colorbar
    view(70,20)
    xc = 5000;
    zc = 0;
    yc = 0; 
    xr = 300; 
    yr = 300;
    zr = 946; 
    n = 16;
[z,y,x] = ellipsoid(zc,yc,xc,zr,yr,xr,n);
    C = zeros(size(x));
surf(x,y,z, C)


    
%     shading faceted