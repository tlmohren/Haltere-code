clc;clear all;close all
addpathFolderStructureHaltere()

loadName = 'figure1_FEAaxis';
saveName = 'figure1_FEAaxis';
renew_data_load =  false
if renew_data_load
    FEA(1).name = 'Haltere_CraneFly_Sphere_Om0';
    FEA(2).name = 'Haltere_CraneFly_ellipsoidHor_Om0';
    FEA(3).name = 'Haltere_CraneFly_ellipsoidVer_Om0';
%     FEA(2).name = 'Haltere_CraneFly_ellipsoidHor_Om10';
%     parameters = { 'disp','u2','v2','w2','flapangle','theta_angle', 'X (µm)'}; % select which parameters to load 
    for j =  1:length(FEA)
        tic
%         [FEA(j).xyz, FEA(j).data, ~] = loadCSV( ['data' filesep  FEA(j).name], parameters);
        [FEA(j).xyz, FEA(j).strain, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'eXX' });
%         [FEA(j).xyz, FEA(j).deform, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'u2','v2','w2'});
%         [FEA(j).xyz, FEA(j).a, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'u2','v2','w2'});
%         [~, FEA(j).angles, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'flapangle','theta_angle'});
        toc 
    end
    % Determine Circle locations
    for j = 1:length(FEA)
        circleDistance = 3000;               % distance from base to haltere 
        circleRadius = 150;                 % radius of haltere   
        mindist =  min( abs( FEA(j).xyz(:,1) - circleDistance) );
        xMatch = find(  abs(FEA(j).xyz(:,1) - circleDistance) <= (mindist+1) );
        yMatch = find( round( abs( FEA(j).xyz(:,2) ), 7) == circleRadius );
        zMatch = find( round( abs( FEA(j).xyz(:,3) ), 7) == circleRadius );

        FEA(j).sideInds = intersect(xMatch,yMatch);
        FEA(j).topInds = intersect(xMatch,zMatch);
    end
    save(['data' filesep saveName],'FEA')
else
    load(['data' filesep loadName],'FEA')
end


%% 
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


%%
theta = linspace(0,2*pi,17);
theta(1) = 2*pi;

t_ind = 100;
deform_mult = 30;
OOP_mult = -300;
% xDes = [0:150:4800];

xDes = [3000:150:4800];


surfParamBackground = {'FaceAlpha',1,'EdgeAlpha',0.4};
surfParamBackgroundTwistStalk = {'FaceAlpha',1,'EdgeAlpha',0.4};
surfParamForeground = {'FaceAlpha',1,'EdgeAlpha',0.4};
    
%% 
% FEA(1).xrtheta(:,1) = FEA(1).xyz(:,1);
% FEA(1).xrtheta(:,2) = sqrt( FEA(1).xyz(:,2).^2  +  FEA(1).xyz(:,3).^2 );
% FEA(1).xrtheta(:,3) = atan2( FEA(1).xyz(:,3), FEA(1).xyz(:,2) ) +pi ;
% % FEA(2).xrtheta = FEA(1).xrtheta;
% 
for j = 2%:length(FEA)
    FEA(j).xrtheta(:,1) = FEA(j).xyz(:,1);
    FEA(j).xrtheta(:,2) = sqrt( FEA(j).xyz(:,2).^2  +  FEA(j).xyz(:,3).^2 );
    FEA(j).xrtheta(:,3) = atan2( FEA(j).xyz(:,3), FEA(j).xyz(:,2) ) +pi ;
% FEA(2).xrtheta = FEA(1).xrtheta;
%     FEA(j).diffPerPoint = squeeze(FEA(j).xyzHaltereFrame(:,t_ind,:)) - FEA(j).xyz;
    for k=1:length(xDes)
        dx = abs(FEA(j).xrtheta(:,1)-xDes(k) );
        dr = abs(FEA(j).xrtheta(:,2)-150 );
        for l = 1:length(theta)
    %         theta(k)
            da = abs(FEA(j).xrtheta(:,3) - theta(l) );
            J = dx.^2 + dr.^2+ (da*150.^2);
            [V,I] = min(J);
%             Xj(k,l) = FEA(j).xyz(I,1)  + FEA(j).diffPerPoint(I,1)*deform_mult; 
%             Yj(k,l) = FEA(j).xyz(I,2) + FEA(j).diffPerPoint(I,2)*deform_mult; 
%             Zj(k,l) = FEA(j).xyz(I,3) + FEA(j).diffPerPoint(I,3)*deform_mult; 
%             Cj(k,l) = FEA(j).strain(I,t_ind);
%             
            Xb(k,l) = FEA(j).xyz(I,1) ; 
            Yb(k,l) = FEA(j).xyz(I,2) ; 
            Zb(k,l) = FEA(j).xyz(I,3) ; 
        end
    end
end







%% Figure 2
width = 2;     % Width in inches,   find column width in paper 
height = 3;    % Height in inches
fig1 = figure();
set(fig1, 'Position', [fig1.Position(1:2)-[width*100,0] width*100, height*100]); %<- Set size
colormap(strainScheme)%     colorbar

%% subplot 311, rotation around x 
xc = 0;
zc = 0;
yc = 0; 
xr = 440; 
yr = 440;
zr = 440; 
n = 16;
[z,y,x] = ellipsoid(zc,yc,xc,zr,yr,xr,n);

Cb = zeros(size(Xb));

subplot(311); hold on 
	sb = surf(Xb,Yb,Zb,Cb);
        set(sb,surfParamBackground{:})
        c_ell = zeros(size(x));
    s1 = surf(x+5e3,y,z,c_ell);
        set(s1,surfParamBackground{:})
        
        axis off
        axis equal
        axis tight;
        shading faceted
        xlabel('x');ylabel('y');zlabel('z')
        view(40,40)
        %% 



xr = 300; 
yr = 948;
zr = 300; 
n = 16;
[z,y,x] = ellipsoid(zc,yc,xc,zr,yr,xr,n);



% [z,y,x] = ellipsoid(0,0,0,948,300,300,16);
subplot(312); hold on 
	sb = surf(Xb,Yb,Zb,Cb);
        set(sb,surfParamBackground{:})
    s1 = surf(x+5e3,y,z,c_ell);
        set(s1,surfParamBackground{:})

        axis tight;
        shading faceted
        axis off
        axis equal
        xlabel('x');ylabel('y');zlabel('z')
        view(40,40)
        
        
        
%% twist 


xr = 300; 
yr = 300;
zr = 948; 
n = 16;
[z,y,x] = ellipsoid(zc,yc,xc,zr,yr,xr,n);


Ctwist = ones(size(Cb));
Ctwist(1,1) = -1.5;
Ctwist(1,2) = 1.5;
% [z,y,x] = ellipsoid(0,0,0,948,300,300,16);
subplot(313); hold on   
%     plot stalks 
	sb = surf(Xb,Yb,Zb,Cb);
        set(sb,surfParamBackgroundTwistStalk{:})
    s1 = surf(x+5e3,y,z,c_ell);
        set(s1,surfParamBackground{:})

        axis tight;
        shading faceted
        axis off
        axis equal
        xlabel('x');ylabel('y');zlabel('z')
        view(40,40)
        
%% Setting paper size for saving 

set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
set(fig1,'InvertHardcopy','on');
set(fig1,'PaperUnits', 'inches');
papersize = get(fig1, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(fig1, 'PaperPosition', myfiguresize);
print(fig1, ['figs' filesep 'Figure1_FEAlegend' ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig1, 'PaperPosition', myfiguresize);
print(fig1, ['figs' filesep 'Figure1_FEAlegend'], '-dsvg', '-r600');


