clc;clear all;close all

run('config_file.m')
addpathFolderStructureHaltere()

loadName = 'FEA_processed_data';
load(['data' filesep loadName],'FEA')

%% Figure 2

for j = 1%:length(FEA)
    FEA(j).xrtheta(:,1) = FEA(j).xyz(:,1);
    FEA(j).xrtheta(:,2) = sqrt( FEA(j).xyz(:,2).^2  +  FEA(j).xyz(:,3).^2 );
    FEA(j).xrtheta(:,3) = atan2( FEA(j).xyz(:,3), FEA(j).xyz(:,2) ) +pi ;
    FEA(1).xrtheta(:,3) = wrapTo2Pi(  atan2( FEA(1).xyz(:,3), FEA(1).xyz(:,2) ) +pi -0.02 )+0.02;
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
            if I == 49
               klR = [k,l]; 
            elseif I == 109
               klL = [k,l]; 
            end
        end
    end
end

angles = [pi/12,0,pi*11/24]; 
xc = 0; yc = 0; zc = 0;
xr = 440; yr = 440; zr = 440; 
n = 16;
[z,y,x] = ellipsoid(zc,yc,xc,zr,yr,xr,n);

% deform ellpisoid 
eul_1 = euler_angle('X',angles(1))^-1;
eul_2 = euler_angle('Y',angles(2))^-1;
eul_3 = euler_angle('Z',angles(3))^-1;
         
for j = 1:size(x,1)
    for k = 1:size(x,2)
        xyzTemp = [x(j,k), y(j,k), z(j,k) ];
        xyzT = eul_1*eul_2*eul_3*xyzTemp'; 
        xOMdiff(j,k) = xyzT(1);
        yOMdiff(j,k) = xyzT(2);
        zOMdiff(j,k) = xyzT(3);
    end
end 

% deform haltere stalk 
for j = 1:size(Xb,1)
    angleTemp = angles;
    eul_1 = euler_angle('X',angleTemp(1))^-1;
    eul_2 = euler_angle('Y',angleTemp(2))^-1;
    eul_3 = euler_angle('Z',angleTemp(3))^-1;         
    for k = 1:size(Xb,2)
        xyzTemp = [Xb(j,k), Yb(j,k), Zb(j,k) ];
        xyzT = eul_1*eul_2*eul_3*xyzTemp'; 
        Xtwist(j,k) = xyzT(1);
        Ytwist(j,k) = xyzT(2);
        Ztwist(j,k) = xyzT(3);
    end
end

fig2 = figure();
    width = 3;     % Width in inches,   find column width in paper 
    height = 3;    % Height in inches
    hold on
    set(fig2, 'Position', [fig2.Position(1:2)-[width*100,0] width*100, height*100]); %<- Set size
    colormap(strainScheme)%     colorbar

    plot3( [0,1]*3e3,[0,1]*0,[0,1]*0 ,'k') 
    plot3( [0,1]*0,[0,1]*6e3,[0,1]*0 ,'k') 
    plot3( [0,1]*0,[0,1]*0,[0,1]*3e3 ,'k') 
     text(3.3e3,0,0,'x')
     text(0,6.3e3,0,'y')
     text(0,0,3.3e3,'z')
    Cb = zeros(size(Xb));
	sb = surf(Xtwist,Ytwist,Ztwist,Cb);
        set(sb,surfParamBackground{:})
        c_ell = zeros(size(x));
    xEl = mean([ Xtwist(klR(1),klR(2)) ,Xtwist(klL(1),klL(2)) ]);
    yEl = mean([ Ytwist(klR(1),klR(2)) ,Ytwist(klL(1),klL(2)) ]);
    zEl = mean([ Ztwist(klR(1),klR(2)) ,Ztwist(klL(1),klL(2)) ]);
    s3 = surf( xOMdiff +xEl,...
            yOMdiff +yEl,...
            zOMdiff +zEl,...
            c_ell);
        axis off; axis equal; axis tight;
        shading faceted
        xlabel('x');ylabel('y');zlabel('z')
        view(140,20)
        
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
print(fig2, ['figs' filesep 'Figure1_FEAaxes' ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig2, 'PaperPosition', myfiguresize);
print(fig2, ['figs' filesep 'Figure1_FEAaxes'], '-dsvg', '-r600');
        