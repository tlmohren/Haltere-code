clc;clear all;close all
addpathFolderStructureHaltere()




%% 

% load(['data' filesep 'haltere_eulerLagrange_sphere.mat'],'T','gammain', 'thetaout', 'phiout')
load(['data' filesep 'haltere_eulerLagrange_sidebulbs_Om10 .mat'],'T','gammain', 'thetaout', 'phiout')
% haltere_eulerLagrange_sidebulbs_Om

width = 2;     % Width in inches,   find column width in paper 
height = 3;    % Height in inches
fig1 = figure();
set(fig1, 'Position', [fig1.Position(1:2) width*100, height*100]); %<- Set size

t = 0:0.001:0.35;
labels = {'$\Delta \phi$','$\Delta \theta$','$\Delta \gamma$'};
axOpts1 = {'XGrid','On','XLim',[0,0.2],'XTick',[0:0.05:0.2]}; 
axOpts2 = {'XGrid','On','XLim',[0,0.2],'XTick',[0:0.05:0.2]}; 
axOpts3 = {'XGrid','On','XLim',[0,0.2],'XTick',[0:0.05:0.2]}; 

for j = 1%:length(FEA)/2
    subplot(3,1,1); hold on 
%         plot(t,FEA(j*2-1).yAngle )
%         plot(t,FEA(j*2).yAngle )
%         ylabel( labels{1} );
%         ax = gca();
%         set(ax,axOpts1{:})
    subplot(312); hold on 
        plot(T,thetaout )
%         plot(T, )
        ylabel( labels{2} );
        ax = gca();
        set(ax,axOpts2{:})
    subplot(313); hold on 
        plot(T,phiout)
%         plot(t,FEA(j*2-1).twistAngle)
%         plot(t,FEA(j*2).twistAngle )
        xlabel('Time (s)'); ylabel( labels{3} );
        ax = gca();
        set(ax,axOpts3{:})
end


%% Setting paper size for saving 

set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
set(fig1,'InvertHardcopy','on');
set(fig1,'PaperUnits', 'inches');
papersize = get(fig1, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(fig1, 'PaperPosition', myfiguresize);
print(fig1, ['figs' filesep 'Figure1_ELdeformPlot' ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig1, 'PaperPosition', myfiguresize);
print(fig1, ['figs' filesep 'Figure1_ELdeformPlot'], '-dsvg', '-r600');








%% Figure 2

surfParamBackground = {'FaceAlpha',1,'EdgeAlpha',0};
surfParamForeground = {'FaceAlpha',1,'EdgeAlpha',0};


colBulbFront = [253,128,0]/255;
colSides = [255, 191, 128]/255;
width = 2;     % Width in inches,   find column width in paper 
height = 3;    % Height in inches
fig2 = figure();
set(fig2, 'Position', [fig2.Position(1:2)-[width*100,0] width*100, height*100]); %<- Set size
% colormap(strainScheme)%     colorbar+


xc = 0;
zc = 0;
yc = 0; 
xr = 300; 
yr = 300;
zr = 300; 
n = 16;
[z,y,x] = ellipsoid(zc,yc,xc,zr,yr,xr,n);

%     Cb = zeros(size(Xb));
    CB = zeros(size(x));
    
CB(:,:,1) = ones(size(x))*colBulbFront(1); % red
CB(:,:,2) = ones(size(x))*colBulbFront(2); % green
CB(:,:,3) = ones(size(x))*colBulbFront(3); % blue
    
CS(:,:,1) = ones(size(x))*colSides(1); % red
CS(:,:,2) = ones(size(x))*colSides(2); % green
CS(:,:,3) = ones(size(x))*colSides(3); % blue
        
dy = 800; 

subplot(311) ; hold on
    xLine = [linspace(0,5000,30) , 5000 , 5000];
    yLine = [zeros(1,30), 0, 0 ];
    zLine = [zeros(1,30), -dy, dy];
    angles = [0,0.2,0];
    for j = 1:length(xLine) 
        angleTemp = angles*j/length(xLine);
        eul_1 = [ 1       0                     0;...
                    0,  cos( angleTemp(1) ),  - sin(  angleTemp(1) )  ; ...
                   0   sin(  angleTemp(1) ) cos(  angleTemp(1) )  ]^-1;
        eul_2 = [cos(  angleTemp(2))      0       sin( angleTemp(2)) ; ...
                    0              1       0 ;...
                    -sin( angleTemp(2))    0        cos(  angleTemp(2))]^-1;
        eul_3 = [cos(  angleTemp(3))   sin(  angleTemp(3))    0 ; ...
                     -sin(  angleTemp(3)) cos(  angleTemp(3))     0 ;...
                     0          0               1]^-1;

            xyzTemp = [xLine(j), yLine(j),  zLine(j)];
            xyz1(j,:) = eul_1*eul_2*eul_3*xyzTemp'; 
    end


    plot3([0,5000],[0,0],[0,0],'Color',[0.5,0.5,0.5,0.5])
    plot3([5000,5000],[0,0],[-dy,dy],'Color',[0.5,0.5,0.5,0.5])
    xl = 0;    zl = 0;    yl = 0;     xbr = 100;     ybr =100;    zbr = 100; 
    [zb,yb,xb] = ellipsoid(zl,yl,xl,zbr,ybr,xbr,n);
    s1 = surf( x + 5e3, y,   z,     CS);

        
    plot3(xyz1(1:10,1),xyz1(1:10,2),xyz1(1:10,3) ,'b')
    plot3(xyz1(21:30,1),xyz1(21:30,2),xyz1(21:30,3) ,'b')
    plot3(xyz1(31:32,1),xyz1(31:32,2),xyz1(31:32,3) ,'b')
        
    
    s2 = surf( xb + 5e3, yb ,   zb+ dy,     CS);
    s3 = surf( xb + 5e3, yb ,   zb- dy ,     CS);
    
    set(s1,surfParamBackground{:})
    set(s2,surfParamBackground{:})
    set(s3,surfParamBackground{:})
    
    s4 = surf( x+  xyz1(30,1) ,y+  xyz1(30,2) , z+  xyz1(30,3) ,     CB);
    
    s5 = surf( xb+  xyz1(31,1), yb+  xyz1(31,2) ,  zb +  xyz1(31,3),     CB);
    s6 = surf( xb+  xyz1(32,1), yb+  xyz1(32,2) ,   zb+  xyz1(32,3),     CB);
    
    
    
    set(s4,surfParamForeground{:})
    set(s5,surfParamForeground{:})
    set(s6,surfParamForeground{:})
%     
%     
        axis off
        axis equal
        axis tight;
        shading faceted
        view(40,40)
        
    
    
    
    
    
    
    
    
    
    
    
    
subplot(312); hold on
        
    xLine = [linspace(0,5000,30) , 5000 , 5000];
    yLine = [zeros(1,30), 0, 0 ];
    zLine = [zeros(1,30), -dy, dy];
    angles = [0,0,0.2];
    for j = 1:length(xLine) 
        angleTemp = angles*j/length(xLine);
        eul_1 = [ 1       0                     0;...
                    0,  cos( angleTemp(1) ),  - sin(  angleTemp(1) )  ; ...
                   0   sin(  angleTemp(1) ) cos(  angleTemp(1) )  ]^-1;
        eul_2 = [cos(  angleTemp(2))      0       sin( angleTemp(2)) ; ...
                    0              1       0 ;...
                    -sin( angleTemp(2))    0        cos(  angleTemp(2))]^-1;
        eul_3 = [cos(  angleTemp(3))   sin(  angleTemp(3))    0 ; ...
                     -sin(  angleTemp(3)) cos(  angleTemp(3))     0 ;...
                     0          0               1]^-1;

            xyzTemp = [xLine(j), yLine(j),  zLine(j)];
            xyz1(j,:) = eul_1*eul_2*eul_3*xyzTemp'; 
    end


    plot3([0,5000],[0,0],[0,0],'Color',[0.5,0.5,0.5,0.5])
    plot3([5000,5000],[0,0],[-dy,dy],'Color',[0.5,0.5,0.5,0.5])
    xl = 0;    zl = 0;    yl = 0;     xbr = 100;     ybr =100;    zbr = 100; 
    [zb,yb,xb] = ellipsoid(zl,yl,xl,zbr,ybr,xbr,n);
    s1 = surf( x + 5e3, y,   z,     CS);
    

        
    plot3(xyz1(1:10,1),xyz1(1:10,2),xyz1(1:10,3) ,'b')
    plot3(xyz1(21:30,1),xyz1(21:30,2),xyz1(21:30,3) ,'b')
    plot3(xyz1(31:32,1),xyz1(31:32,2),xyz1(31:32,3) ,'b')
    s4 = surf( x+  xyz1(30,1) ,y+  xyz1(30,2) , z+  xyz1(30,3) ,     CB);
        
    
    s2 = surf( xb + 5e3, yb ,   zb+ dy,     CS);
    s3 = surf( xb + 5e3, yb ,   zb- dy ,     CS);
    
    set(s1,surfParamBackground{:})
    set(s2,surfParamBackground{:})
    set(s3,surfParamBackground{:})
    
    
    s5 = surf( xb+  xyz1(31,1), yb+  xyz1(31,2) ,  zb+  xyz1(31,3),     CB);
    s6 = surf( xb+  xyz1(32,1), yb+  xyz1(32,2) ,   zb+  xyz1(32,3),     CB);
    
    
    
    set(s4,surfParamForeground{:})
    set(s5,surfParamForeground{:})
    set(s6,surfParamForeground{:})
    
        axis off
        axis equal
        axis tight;
        shading faceted
        view(40,40)
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
subplot(313) ; hold on
        
    xLine = [linspace(0,5000,30) , 5000 , 5000];
    yLine = [zeros(1,30), 0, 0 ];
    zLine = [zeros(1,30), -dy, dy];
    angles = [0.5,0,0];
    for j = 1:length(xLine) 
        angleTemp = angles*j/length(xLine);
        eul_1 = [ 1       0                     0;...
                    0,  cos( angleTemp(1) ),  - sin(  angleTemp(1) )  ; ...
                   0   sin(  angleTemp(1) ) cos(  angleTemp(1) )  ]^-1;
        eul_2 = [cos(  angleTemp(2))      0       sin( angleTemp(2)) ; ...
                    0              1       0 ;...
                    -sin( angleTemp(2))    0        cos(  angleTemp(2))]^-1;
        eul_3 = [cos(  angleTemp(3))   sin(  angleTemp(3))    0 ; ...
                     -sin(  angleTemp(3)) cos(  angleTemp(3))     0 ;...
                     0          0               1]^-1;

            xyzTemp = [xLine(j), yLine(j),  zLine(j)];
            xyz1(j,:) = eul_1*eul_2*eul_3*xyzTemp'; 
    end


    plot3([0,5000],[0,0],[0,0],'Color',[0.5,0.5,0.5,0.5])
    plot3([5000,5000],[0,0],[-dy,dy],'Color',[0.5,0.5,0.5,0.5])
    xl = 0;    zl = 0;    yl = 0;     xbr = 100;     ybr =100;    zbr = 100; 
    [zb,yb,xb] = ellipsoid(zl,yl,xl,zbr,ybr,xbr,n);
%     s1 = surf( x + 5e3, y,   z,     CS);
       
        
    plot3(xyz1(1:10,1),xyz1(1:10,2),xyz1(1:10,3) ,'b')
    plot3(xyz1(21:30,1),xyz1(21:30,2),xyz1(21:30,3) ,'b')
    plot3(xyz1(31:32,1),xyz1(31:32,2),xyz1(31:32,3) ,'b')
    s4 = surf( x+  xyz1(30,1) ,y+  xyz1(30,2) , z+  xyz1(30,3) ,     CB);
        
    
    s2 = surf( xb + 5e3, yb ,   zb+ dy,     CS);
    s3 = surf( xb + 5e3, yb ,   zb- dy ,     CS);
    
%     set(s1,surfParamBackground{:})
    set(s2,surfParamBackground{:})
    set(s3,surfParamBackground{:})
    
    
    s5 = surf( xb+  xyz1(31,1), yb+  xyz1(31,2) ,  zb+  xyz1(31,3),     CB);
    s6 = surf( xb+  xyz1(32,1), yb+  xyz1(32,2) ,   zb+  xyz1(32,3),     CB);
    
    
    
    set(s4,surfParamForeground{:})
    set(s5,surfParamForeground{:})
    set(s6,surfParamForeground{:})
    
        axis off
        axis equal
        axis tight;
        shading faceted
        view(40,40)
        
        
 

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
print(fig2, ['figs' filesep 'Figure1_ELdeform' ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig2, 'PaperPosition', myfiguresize);
print(fig2, ['figs' filesep 'Figure1_ELdeform'], '-dsvg', '-r600');
        
  
        
        
