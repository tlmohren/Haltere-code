clc;clear all;close all

addpathFolderStructureHaltere()
run('config_file.m')

%%
% loadName = 'figure2_strainData';
% saveName = 'figure2_strainData';

loadName = 'figure3_crossDeform';
saveName = 'figure3_crossDeform';

% renew_data_load = false
renew_data_load = true
if renew_data_load
%     FEA(1).name = 'Haltere_CraneFlyLowDensity_sphereCrossStalk_Om0';
%     FEA(2).name = 'Haltere_CraneFlyLowDensity_sphereCrossStalk_Om10';
    FEA(1).name = 'Haltere_CraneFly_ellipsoidVerCrossStalk_Om0';
    FEA(2).name = 'Haltere_CraneFly_ellipsoidVerCrossStalk_Om10';
%     FEA(1).name = 'Haltere_CraneFly_Sphere_Om0';
%     FEA(2).name = 'Haltere_CraneFly_Sphere_Om10';
%     FEA(3).name = 'Haltere_CraneFly_ellipsoidHor_Om0';
%     FEA(4).name = 'Haltere_CraneFly_ellipsoidHor_Om10';
%     FEA(5).name = 'Haltere_CraneFly_ellipsoidVer_Om0';
%     FEA(6).name = 'Haltere_CraneFly_ellipsoidVer_Om10';   
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
FEA(1).xrtheta(:,3) = wrapTo2Pi(  atan2( FEA(1).xyz(:,3), FEA(1).xyz(:,2) ) +pi -0.02 )+0.02;


for k=1:length(xDes)
    dx = abs(FEA(1).xrtheta(:,1)-xDes(k) );
    dr = abs(FEA(1).xrtheta(:,2)-150 );
%         r_right = find( dr < 10);
    for l = 1:length(theta)
        
%         minx =  min( abs( FEA(1).xyz(:,1) - xDes(k) ) );
        da = abs(FEA(1).xrtheta(:,3) - theta(l) );
%         a_right = find( da < 0.1 );
        
%         I_circ = intersect(r_right,a_right)
        
        J = dx.^2 + dr.^2*1+ (da*150.^2);
        [V,I] = min(J);
        Xfront(k,l) = FEA(1).xyz(I,1) ; 
        Yfront(k,l) = FEA(1).xyz(I,2) ; 
        Zfront(k,l) = FEA(1).xyz(I,3) ; 

    end
end
%% 

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
    splane = surf( [1,1;1,1]*300, [-1,1;-1,1]*300,[-1,-1;1,1]*300,[1,1;1,1]*0);
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
    whichInds = [5:9];
    plot3( 300*ones(length(thetTemp)), sin(thetTemp)*r, cos(thetTemp)*r,'k') 
%     scatter3( FEA(1).xyz( FEA(1).sideInds,1), FEA(1).xyz( FEA(1).sideInds,2), FEA(1).xyz( FEA(1).sideInds,3),30,'filled','r') 
    scatter3( FEA(1).xyz( FEA(1).circleInds(whichInds) ,1), FEA(1).xyz( FEA(1).circleInds(whichInds),2), FEA(1).xyz( FEA(1).circleInds(whichInds),3),10,'filled','k')
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
print(fig1, ['figs' filesep 'Figure4_crossSectionMesh' ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig1, 'PaperPosition', myfiguresize);
print(fig1, ['figs' filesep 'Figure4_crossSectionMesh'], '-dsvg', '-r600');
        







%% now 2d cross sections 

% Figure 2
fig2 = figure();
    width = 4;     % Width in inches,   find column width in paper 
    height = 4;    % Height in inches
    set(fig2, 'Position', [fig2.Position(1:2)-[width*100,0] width*100, height*100]); %<- Set size
    colormap(strainScheme)%     colorbar
    
    thetTemp = linspace(0,2*pi,50);
    r = 150; 
%     subplot(211)
        hold on
        splane = fill(  [1,-1,-1,1]*300,[-1,-1,1,1]*300,[1,1,1,1]);
        plot(  sin(thetTemp)*r, cos(thetTemp)*r,'k') 
        scatter(  FEA(1).xyz( FEA(1).circleInds(whichInds),2), FEA(1).xyz( FEA(1).circleInds(whichInds),3),50,'filled','k') 

    %     scatter( FEA(1).xyz( FEA(1).topInds,2), FEA(1).xyz( FEA(1).topInds,3),30,'filled','b')
        shading faceted
        axis tight;  axis off; axis equal
        xlabel('x');ylabel('y');%zlabel('z')
%         view(60,15)

%     subplot(212)
%         hold on
%         splane = fill(  [1,-1,-1,1]*300,[-1,-1,1,1]*300,[1,1,1,1]);
%         plot(  sin(thetTemp)*r, cos(thetTemp)*r,'k') 
%     %     scatter(  FEA(1).xyz( FEA(1).sideInds,2), FEA(1).xyz( FEA(1).sideInds,3),30,'filled','r') 
%         scatter( FEA(1).xyz( FEA(1).topInds,2), FEA(1).xyz( FEA(1).topInds,3),30,'filled','b')
%         shading faceted
%         axis tight;  axis off; axis equal
%         xlabel('x');ylabel('y');%zlabel('z')
% %         view(60,15)








       
%% Setting paper size for saving 

set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
tightfig;
set(fig2,'InvertHardcopy','on');
set(fig2,'PaperUnits', 'inches');
papersize = get(fig2, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(fig2, 'PaperPosition', myfiguresize);
print(fig2, ['figs' filesep 'Figure4_crossSectionMeshDot' ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig2, 'PaperPosition', myfiguresize);
print(fig2, ['figs' filesep 'Figure4_crossSectionMeshDot'], '-dsvg', '-r600');
        
% 
