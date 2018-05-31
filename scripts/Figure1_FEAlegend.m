clc;clear all;close all
addpathFolderStructureHaltere()
run('config_file.m')

%% load mesh position data
loadName = 'figure1_FEAaxis';
saveName = 'figure1_FEAaxis';
renew_data_load =  true
if renew_data_load
    FEA(1).name = 'Haltere_CraneFly_Sphere_Om0';
    FEA(2).name = 'Haltere_CraneFly_ellipsoidHor_Om0';
    FEA(3).name = 'Haltere_CraneFly_ellipsoidVer_Om0';  
    for j =  1:length(FEA)
        tic
        [FEA(j).xyz, ~, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'eXX' });
        toc 
    end
    save(['data' filesep saveName],'FEA')
else
    load(['data' filesep loadName],'FEA')
end
    
% Sort mesh locations for surface plots 
xDes = [3000:150:4800];
for j = 2
    FEA(j).xrtheta(:,1) = FEA(j).xyz(:,1);
    FEA(j).xrtheta(:,2) = sqrt( FEA(j).xyz(:,2).^2  +  FEA(j).xyz(:,3).^2 );
    FEA(j).xrtheta(:,3) = atan2( FEA(j).xyz(:,3), FEA(j).xyz(:,2) ) +pi ;
    for k=1:length(xDes)
        dx = abs(FEA(j).xrtheta(:,1)-xDes(k) );
        dr = abs(FEA(j).xrtheta(:,2)-150 );
        for l = 1:length(theta)
    %         theta(k)
            da = abs(FEA(j).xrtheta(:,3) - theta(l) );
            J = dx.^2 + dr.^2+ (da*150.^2);
            [V,I] = min(J);
            Xb(k,l) = FEA(j).xyz(I,1) ; 
            Yb(k,l) = FEA(j).xyz(I,2) ; 
            Zb(k,l) = FEA(j).xyz(I,3) ; 
        end
    end
end

%% Figure 2
fig1 = figure();
    width = 2;     % Width in inches,   find column width in paper 
    height = 3;    % Height in inches
    set(fig1, 'Position', [fig1.Position(1:2)-[width*100,0] width*100, height*100]); %<- Set size
    colormap(strainScheme)%     colorbar

% subplot 311, sphere 
xc = 0; yc = 0; zc = 0;
xr = 440; yr = 440;zr = 440; n = 16;
[z,y,x] = ellipsoid(zc,yc,xc,zr,yr,xr,n);

c_ell = zeros(size(x));
Cb = zeros(size(Xb));
subplot(311); hold on 
    hold on   
	sb = surf(Xb,Yb,Zb,Cb);
    s1 = surf(x+5e3,y,z,c_ell);
        set(sb,surfParamBackgroundTwistStalk{:})
        set(s1,surfParamBackground{:})
        shading faceted
        axis tight;  axis off; axis equal
        xlabel('x');ylabel('y');zlabel('z')
        view(40,40)
        
%subplot 312 horozontal ellipsoid 
xr = 300; yr = 948;zr = 300; n = 16;
[z,y,x] = ellipsoid(zc,yc,xc,zr,yr,xr,n);

subplot(312); hold on 
	sb = surf(Xb,Yb,Zb,Cb);
    s1 = surf(x+5e3,y,z,c_ell);
        set(sb,surfParamBackgroundTwistStalk{:})
        set(s1,surfParamBackground{:})
        shading faceted
        axis tight;  axis off; axis equal
        xlabel('x');ylabel('y');zlabel('z')
        view(40,40)
        
%subplot 313, vertical ellipsoid 
xr = 300; yr = 300;zr = 948; n = 16;
[z,y,x] = ellipsoid(zc,yc,xc,zr,yr,xr,n);
subplot(313);     hold on   
	sb = surf(Xb,Yb,Zb,Cb);
    s1 = surf(x+5e3,y,z,c_ell);
        set(sb,surfParamBackgroundTwistStalk{:})
        set(s1,surfParamBackground{:})
        shading faceted
        axis tight;  axis off; axis equal
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
