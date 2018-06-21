clc;clear all;close all

addpathFolderStructureHaltere()
run('config_file.m')

%%
loadName = 'figure4_strainData';
saveName = 'figure4_strainData';

renew_data_load = false
% renew_data_load = true
if renew_data_load
    FEA(1).name = 'Haltere_CraneFly_Sphere_Om0';
    FEA(2).name = 'Haltere_CraneFly_Sphere_Om10';
    FEA(3).name = 'Haltere_CraneFly_ellipsoidHor_Om0';
    FEA(4).name = 'Haltere_CraneFly_ellipsoidHor_Om10';
    FEA(5).name = 'Haltere_CraneFly_ellipsoidVer_Om0';
    FEA(6).name = 'Haltere_CraneFly_ellipsoidVer_Om10';   
    for j =  1:length(FEA)
        tic
        [FEA(j).xyz, FEA(j).strain, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'eXX' });        toc 
    end
    % Determine Circle locations
    for j = 1:length(FEA)
        circleDistance = 300;               % distance from base to haltere 
        circleRadius = 150;                 % radius of haltere   
        mindist =  min( abs( FEA(j).xyz(:,1) - circleDistance) );
        xMatch = find(  abs(FEA(j).xyz(:,1) - circleDistance) <= (mindist+1) );
        
        
        rMatch = find( sqrt(FEA(j).xyz(:,2).^2 + FEA(j).xyz(:,3).^2)  >= circleRadius*0.99 );

        yMatch = find( round( abs( FEA(j).xyz(:,2) ), 7) == circleRadius );
        zMatch = find( round( abs( FEA(j).xyz(:,3) ), 7) == circleRadius );

        FEA(j).circleIndsUnsorted = intersect( xMatch,rMatch); 
        
        angle = atan2( FEA(j).xyz( FEA(j).circleIndsUnsorted,3), ...
            FEA(j).xyz( FEA(j).circleIndsUnsorted,2) );
        angleDeg = rad2deg(angle);
        angleDeg(angleDeg<0) = angleDeg(angleDeg<0)+360;
        [V,I_sort] = sort(angleDeg,'ascend');
% 
        FEA(j).circleInds= FEA(j).circleIndsUnsorted(I_sort);

        FEA(j).sideInds = intersect(xMatch,yMatch);
        FEA(j).topInds = intersect(xMatch,zMatch);
    end
    save(['data' filesep saveName],'FEA')
else
    load(['data' filesep loadName],'FEA')
end







%% 

FEA(1).xrtheta(:,1) = FEA(1).xyz(:,1);
FEA(1).xrtheta(:,2) = sqrt( FEA(1).xyz(:,2).^2  +  FEA(1).xyz(:,3).^2 );
FEA(1).xrtheta(:,3) = atan2( FEA(1).xyz(:,3), FEA(1).xyz(:,2) ) +pi ;

for k=1:length(xDes)
    dx = abs(FEA(1).xrtheta(:,1)-xDes(k) );
    dr = abs(FEA(1).xrtheta(:,2)-150 );
    for l = 1:length(theta)
        da = abs(FEA(1).xrtheta(:,3) - theta(l) );
        J = dx.^2 + dr.^2+ (da*150.^2);
        [V,I] = min(J);
        Xfront(k,l) = FEA(1).xyz(I,1) ; 
        Yfront(k,l) = FEA(1).xyz(I,2) ; 
        Zfront(k,l) = FEA(1).xyz(I,3) ; 

    end
end


xDes = [0:150:300];
for k=1:length(xDes)
    dx = abs(FEA(1).xrtheta(:,1)-xDes(k) );
    dr = abs(FEA(1).xrtheta(:,2)-150 );
    for l = 1:length(theta)
        da = abs(FEA(1).xrtheta(:,3) - theta(l) );
        J = dx.^2 + dr.^2+ (da*150.^2);
        [V,I] = min(J);
        Xback(k,l) = FEA(1).xyz(I,1) ; 
        Yback(k,l) = FEA(1).xyz(I,2) ; 
        Zback(k,l) = FEA(1).xyz(I,3) ; 

    end
end
   
surfParamBackground = {'FaceAlpha',0.5,'EdgeAlpha',0.2};
surfParamForeground = {'FaceAlpha',0.5,'EdgeAlpha',0.2};
surfParamPlane = {'FaceAlpha',0.2,'EdgeAlpha',0.4};
%     
% 
% Figure 2
fig1 = figure();
    width = 4;     % Width in inches,   find column width in paper 
    height = 4;    % Height in inches
    set(fig1, 'Position', [fig1.Position(1:2)-[width*100,0] width*100, height*100]); %<- Set size
    colormap(strainScheme)%     colorbar
% 
% % subplot 311, rotation around x 
xc = 0; yc = 0; zc = 0;
xr = 440; yr = 440;zr = 440; 
[z,y,x] = ellipsoid(zc,yc,xc,zr,yr,xr,n);
% 
Cb = zeros(size(Xback));
Cf = zeros(size(Xfront));
C = zeros(size(x));
% 
% subplot(311); 
    hold on 
	sback = surf(Xback,Yback,Zback,Cb);
        set(sback,surfParamBackground{:})
    splane = surf( [1,1;1,1]*300, [-1,1;-1,1]*300,[-1,-1;1,1]*300,[1,1;1,1]*0)
        set(splane,surfParamPlane{:})
	sfront = surf(Xfront,Yfront,Zfront,Cf);
        set(sfront,surfParamForeground{:})
%     s2 = surf(Xj,Yj,Zj,Cj);
%         set(s2,surfParamForeground{:})
    s1 = surf(x+5e3,y,z,C);
        set(s1,surfParamForeground{:})
% % 

    thetTemp = linspace(0,2*pi,50);
    r = 150; 
    plot3( 300*ones(length(thetTemp)), sin(thetTemp)*r, cos(thetTemp)*r,'k') 
    
    
selected_dots = 5:9;
for kl = 1:length(selected_dots)
    k = selected_dots(kl);
    scatter3( 300,...
    FEA(1).xyz( FEA(1).circleInds(k), 2),...
    FEA(1).xyz( FEA(1).circleInds(k), 3),...
    10,'filled','k')
end
%     scatter3( FEA(1).xyz( FEA(1).sideInds,1), FEA(1).xyz( FEA(1).sideInds,2), FEA(1).xyz( FEA(1).sideInds,3),30,'filled','r') 
%     scatter3( FEA(1).xyz( FEA(1).topInds,1), FEA(1).xyz( FEA(1).topInds,2), FEA(1).xyz( FEA(1).topInds,3),30,'filled','b')
        shading faceted
        axis tight;  axis off; axis equal
        xlabel('x');ylabel('y');zlabel('z')
        view(60,15)
        
       
%% Setting paper size for saving 

set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
tightfig;
set(fig1,'InvertHardcopy','on');
set(fig1,'PaperUnits', 'inches');
papersize = get(fig1, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(fig1, 'PaperPosition', myfiguresize);
print(fig1, ['figs' filesep 'Figure4_crossSectionmesh' ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig1, 'PaperPosition', myfiguresize);
print(fig1, ['figs' filesep 'Figure4_crossSectionmesh'], '-dsvg', '-r600');
        

