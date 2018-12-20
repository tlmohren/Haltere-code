clc;clear all;close all

addpathFolderStructureHaltere()
run('config_file.m')

loadName = 'FEA_processed_data';
load(['data' filesep loadName],'FEA')

%% 
sideL = 243.435; 
circleDistance = 300;               % distance from base to haltere 
er = 0.01;
for j = 1:length(FEA)
    xMatch = find(   abs(abs( FEA(j).xyz(:,1) ) - circleDistance)  <er );
    yMatch = find(   abs(abs( FEA(j).xyz(:,2) ) - sideL) <er & ...
                 abs( FEA(j).xyz(:,3) )   <er);
    zMatch = find(   abs(abs( FEA(j).xyz(:,3) ) - sideL)  <er & ...
                 abs( FEA(j).xyz(:,2) )   <er);
    FEA(j).sideInds = intersect(xMatch,yMatch);
    FEA(j).topInds = intersect(xMatch,zMatch);
end

sideL = 243.435; 
topL = 20.2865;
topV = linspace(-topL,topL,5);
sideV = linspace(topL,sideL,5);
        % y, z 
crossLoc = [  topV(1:5)'        , [1,1,1,1,1]'*sideL; 
         [1,1,1,1]'*topL    ,  flipud(sideV(1:4)') ; 
             sideV(2:5)'    , [1,1,1,1]'*topL  ;
         [1,1,1,1]'*sideL   ,  -topV(2:5)'    ;
         flipud(sideV(1:4)') , -[1,1,1,1]'*topL;
         [1,1,1,1]'*topL, -sideV(2:5)';
          -topV(2:5)', -[1,1,1,1]'*sideL;
          -[1,1,1,1]'*topL, -  flipud(sideV(1:4)')
           -sideV(2:5)',  -[1,1,1,1]'*topL
           -[1,1,1,1]'*sideL , topV(2:5)'
           flipud(-sideV(1:4)'), [1,1,1,1]'*topL 
            -[1,1,1,1]'*topL, sideV(2:5)'
          ];
          
xDes = 0:150:300;
for k= 1:length(xDes)
    for k2 = 1:length(crossLoc)
        
       Xback(k,k2) = xDes(k);
       Yback(k,k2) = crossLoc(k2,1); 
       Zback(k,k2) = crossLoc(k2,2); 
    end
end

xDes = 300:150:4600;
for k= 1:length(xDes)
    for k2 = 1:length(crossLoc)
        
       Xfront(k,k2) = xDes(k);
       Yfront(k,k2) = crossLoc(k2,1); 
       Zfront(k,k2) = crossLoc(k2,2); 
    end
end
 
surfParamBackground = {'FaceAlpha',0.5,'EdgeAlpha',0.2};
surfParamForeground = {'FaceAlpha',0.5,'EdgeAlpha',0.2};
surfParamPlane = {'FaceAlpha',0.2,'EdgeAlpha',0.4};
%      
%%  Figure 2
fig1 = figure();
    width = 4;     % Width in inches,   find column width in paper 
    height = 4;    % Height in inches
    set(fig1, 'Position', [fig1.Position(1:2)-[width*100,0] width*100, height*100]); %<- Set size
    colormap(strainScheme)%     colorbar
% 
% % subplot 311, rotation around x 
xc = 0; yc = 0; zc = 0;
% xr = 440; yr = 440;zr = 440; 
xr = 500; yr = 500;zr = 500; 
[z,y,x] = ellipsoid(zc,yc,xc,zr,yr,xr,n);
% 
Cb = zeros(size(Xback));
Cf = zeros(size(Xfront));
C = zeros(size(x)); 

hold on 
sback = surf(Xback,Yback,Zback,Cb);
    set(sback,surfParamBackground{:})
splane = surf( [1,1;1,1]*300, [-1,1;-1,1]*400,[-1,-1;1,1]*400,[1,1;1,1]*0)
    set(splane,surfParamPlane{:})
sfront = surf(Xfront,Yfront,Zfront,Cf);
    set(sfront,surfParamForeground{:})
s1 = surf(x+5e3,y,z,C);
    set(s1,surfParamForeground{:})

thetTemp = linspace(0,2*pi,50);
r = sideL; 
plot3( 300*ones(length(crossLoc(:,1)),1), crossLoc(:,1)', crossLoc(:,2)','k') 
scatter3( FEA(1).xyz( FEA(1).sideInds,1), FEA(1).xyz( FEA(1).sideInds,2), FEA(1).xyz( FEA(1).sideInds,3),30,'filled','r') 
scatter3( FEA(1).xyz( FEA(1).topInds,1), FEA(1).xyz( FEA(1).topInds,2), FEA(1).xyz( FEA(1).topInds,3),30,'filled','b')
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
print(fig1, ['figs' filesep 'Figure5_crossSectionmesh' ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig1, 'PaperPosition', myfiguresize);
print(fig1, ['figs' filesep 'Figure5_crossSectionmesh'], '-dsvg', '-r600');
        


%% 

fig2 = figure();
    width = 4;     % Width in inches,   find column width in paper 
    height = 4;    % Height in inches
    set(fig2, 'Position', [fig2.Position(1:2)-[width*100,0] width*100, height*100]); %<- Set size
    colormap(strainScheme)%     colorbar
% 
hold on
    splane = surf( [1,1;1,1]*300, [-1,1;-1,1]*300,[-1,-1;1,1]*300,[1,1;1,1]*0)
        set(splane,surfParamPlane{:})
% 	sfront = surf(Xfront,Yfront,Zfront,Cf);
%         set(sfront,surfParamForeground{:})
%     s1 = surf(x+5e3,y,z,C);
%         set(s1,surfParamForeground{:})

    thetTemp = linspace(0,2*pi,50);
    r = sideL; 
    plot3( 300*ones(length(crossLoc(:,1)),1), crossLoc(:,1)', crossLoc(:,2)','k') 
    scatter3( FEA(1).xyz( FEA(1).sideInds,1), FEA(1).xyz( FEA(1).sideInds,2), FEA(1).xyz( FEA(1).sideInds,3),30,'filled','r') 
    scatter3( FEA(1).xyz( FEA(1).topInds,1), FEA(1).xyz( FEA(1).topInds,2), FEA(1).xyz( FEA(1).topInds,3),30,'filled','b')
        shading faceted
        axis tight;  axis off; axis equal
        xlabel('x');ylabel('y');zlabel('z')
        view(90,0)
%% 
set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
tightfig;
set(fig2,'InvertHardcopy','on');
set(fig2,'PaperUnits', 'inches');
papersize = get(fig2, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(fig2, 'PaperPosition', myfiguresize);
print(fig2, ['figs' filesep 'Figure5_crossDots' ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig2, 'PaperPosition', myfiguresize);
print(fig2, ['figs' filesep 'Figure5_crossDots'], '-dsvg', '-r600');
        