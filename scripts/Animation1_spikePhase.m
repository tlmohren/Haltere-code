clc;clear all;close all

addpathFolderStructureHaltere()
run('config_file.m')

animation_name = 'haltereSpiking_highQualityTest';
%%
loadName = 'figure4_strainData';
saveName = 'figure4_strainData';

renew_data_load = false
% renew_data_load = true
if renew_data_load
    FEA(1).name = 'Haltere_CraneFlyLowDensitywBulb_Sphere_Om0'; 
    FEA(2).name = 'Haltere_CraneFlyLowDensityt8u7wBulb_Sphere_Om10'
    for j =  1:length(FEA)
        tic
%         [FEA(j).xyz, FEA(j).strain, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'eXX' });        toc 
        
        [~, FEA(j).strain, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'eXX' });
        [FEA(j).xyz, FEA(j).deform, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'u2','v2','w2'});
        [~, FEA(j).angles, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'flapangle','theta_angle'});
        toc
    end
    % Determine Circle locations
    for j = 1:2 
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
        [~, FEA(j).phi, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'flapangle'});
        
        FEA(j).phi = FEA(j).angles(1,:,1) ;
        FEA(j).theta = -FEA(j).angles(1,:,2) ; 
        [ n_points,n_times, n_deform] = size( FEA(j).deform );
        % compute actual point locations 
        FEA(j).xyzPoints = FEA(j).deform  + ...       
            permute( repmat( FEA(j).xyz,1,1, n_times), [1,3,2] ) ;
        
        
    end
    save(['data' filesep saveName],'FEA')
else
    load(['data' filesep loadName],'FEA')
end

%% apply neural encoding

STAfreq = 0.5;
STAwidth = 5;
STAdelay = 5;
NLDgrad = 20;
NLDshift = 0.8;
[STAfun,NLDfun]=createNeuralFilters( STAfreq,STAwidth,STAdelay,NLDgrad,NLDshift );

fSamp = 1000;
subSamp =10;
STAt = linspace(-39,0,40*subSamp);
parse_extra = 5;
% calib_param_max = [0.00298239604088481,0.0212443192274280,0.0372780009360082,0.0479350861876651,0.0514646550697229,0.0475135494972593,0.0363905581996612,0.0198191699014800,0.00299958384555981,0.0201797042231435,0.0370669398494605,0.0483996355340012,0.0524251480251299,0.0487363345563017,0.0377914840661745,0.0213792050886826];
calib_param_max = [0.000705312124793717,0.00521574728317378,0.00916987624184632,0.0117979355154750,0.0126693024102941,0.0116192206073620,0.00884906250717318,0.00484166706412837,0.000716359983875505,0.00495047057398658,0.00903488798016781,0.0118425516863485,0.0128937828713950,0.0119841717461618,0.00928771453961158,0.00524523737029323];

    selected_dots = 5:9;
    len = 101;
    start = 35;
    It = start:(start+len-1);
    t_plot = (0:len-1)*0.001;

    axCircle= { 'XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)]}; 
    axZero= { 'XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)],'XTickLabel',{'','',''}}; 
    axCircleSpike= { 'YLim',[0,1.6],'XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)]}; 
    axZeroSpike= { 'YLim',[0,1.6],'XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)],'XTickLabel',{'','',''}}; 

calib_param_cross = [0.0047, 0.0803, 0.0047, 0.0803]; 
for jj = 1:length(FEA)/2 
    for j = [1,2]+(2*(jj-1))
        for k = 1:length(FEA(j).circleInds)
         tL = size(FEA(j).strain,2);
            strainTemp = FEA(j).strain( FEA(j).circleInds(k),:);
            t_plot = (1:tL)/fSamp;
            tNew = linspace(t_plot(1),t_plot(end), tL*subSamp ) ; 
            strainInterp =  interp1(t_plot,strainTemp ,tNew,'spline');

            STA = STAfun(STAt);
            strainConv = conv( [zeros(1,length(STA)-1),strainInterp], fliplr( STA), 'valid');
%             strainConv = conv( [zeros(1,length(STA)-1),strainInterp], ( STA), 'valid');
            calib_param(j,k) = max(  strainConv );
    %         FEA(j).pFire(k,:) = NLDfun( strainConv/calib_param(j,k) );
            if jj ~=7
                FEA(j).pFire(k,:) = NLDfun( strainConv/calib_param_max(k) );
                FEA(j).spikeInds{k} = findSpikes( FEA(j).pFire(k,:) ); 
            elseif jj == 7
                
                FEA(j).pFire(k,:) = NLDfun( strainConv/calib_param_cross(k) );
                FEA(j).spikeInds{k} = findSpikes( FEA(j).pFire(k,:) ); 
            end
        end
    end
    for k = 1:length(FEA(j).circleInds)
       FEA(jj*2).dI = FEA(jj*2).spikeInds{k}(end) - FEA(jj*2-1).spikeInds{k}(end);
       FEA(jj*2).dT(k) = abs(FEA(jj*2).dI/1e4 );
    end
    
end


%% 

% only used for phi 
len1 = 26;
start1 = 166;
It1 = start1:(start1+len1-1);
t_plot1 = linspace(0,25,len1);

% used for spike timing 
len10 = 251;
start10 = 1661;
It10 = start10:(start10+len10-1);
t_plot10 = linspace(0,25,len10);

spike_order = [5:13, 13:16, 1:5];

%%


FEA(1).xrtheta(:,1) = FEA(1).xyz(:,1);
FEA(1).xrtheta(:,2) = sqrt( FEA(1).xyz(:,2).^2  +  FEA(1).xyz(:,3).^2 ); 
FEA(1).xrtheta(:,3) = wrapTo2Pi(  atan2( FEA(1).xyz(:,3), FEA(1).xyz(:,2) ) +pi -0.02 )+0.02;


theta = linspace(0,2*pi,17);
theta(1) = 2*pi;
xDes = 0:150:4500; 
for j = 1:200
    for k = 1:length(xDes)
        dx = abs(FEA(1).xrtheta(:,1)-xDes(k) );
        dr = abs(FEA(1).xrtheta(:,2)-150 );
        for l = 1:length(theta)
    %         theta(k)
            da = abs(FEA(1).xrtheta(:,3) - theta(l) );
            J = dx.^2 + dr.^2+ (da*150.^2);
            [V,I] = min(J);
           Xb(k,l,j) = FEA(1).xyzPoints(I,j,1);
           Yb(k,l,j) = FEA(1).xyzPoints(I,j,2);
           Zb(k,l,j) = FEA(1).xyzPoints(I,j,3);
           Cb(k,l,j) = FEA(1).strain(I,j); 
        end
    end
end
% 1 = 1; 

% [z,y,xBulb] = ellipsoid(0,0,0,440,440,440,16);

colormap(strainScheme)%     colorbar

%% make intermediate haltere positions 
tOrig = 1:26;
phiOrig = FEA(1).phi(start1: (start1+len1-1));
tNew = 1:0.1:26; 
phiInterp = interp1(tOrig,phiOrig, tNew );

%% 


[zBulb,yBulb,xBulb] = ellipsoid(0,0,0,440,440,440,16);
xBulb = xBulb + 4800; 

FEA(1).xrtheta(:,1) = FEA(1).xyz(:,1);
FEA(1).xrtheta(:,2) = sqrt( FEA(1).xyz(:,2).^2  +  FEA(1).xyz(:,3).^2 ); 
FEA(1).xrtheta(:,3) = wrapTo2Pi(  atan2( FEA(1).xyz(:,3), FEA(1).xyz(:,2) ) +pi -0.02 )+0.02;


XbStalk = zeros(length(xDes), length(theta),200); 
YbStalk = XbStalk;
ZbStalk = XbStalk; 

for j = 1:length(phiInterp)
    j 
    for k=1:length(xDes)
        dx = abs(FEA(1).xrtheta(:,1)-xDes(k) );
        dr = abs(FEA(1).xrtheta(:,2)-150 );
        for l = 1:length(theta)
    %         theta(k)
            da = abs(FEA(1).xrtheta(:,3) - theta(l) );
            J = dx.^2 + dr.^2+ (da*150.^2);
            [V,I] = min(J); 
           
           
            xyzTemp = [FEA(1).xyz(I,1), FEA(1).xyz(I,2), FEA(1).xyz(I,3) ];

            eul_2 = euler_angle('Y',phiInterp(j) )^-1;
            
            xyzTranspose =eul_2*xyzTemp'; 
            
            XbStalk(k,l,j) = xyzTranspose(1) ; 
            YbStalk(k,l,j) = xyzTranspose(2) ; 
            ZbStalk(k,l,j) = xyzTranspose(3) ; 
            
           
        end
    end
end









































%% 
C = ones(size(xBulb))*0;
ImaxBottom = 155;
ImaxTop = 39;


spike_order = [5:13, fliplr([13:16, 1:5]) ]; 
surfParamForeground = {'EdgeAlpha',0.4 };

% surfParamForeground = {'EdgeAlpha',0.4};
rectif = @(x) (x<=1)*x ;
 

%% 
 
 
fsz = 14;
set(0,'DefaultAxesFontSize',fsz)% .
set(0,'DefaultLegendFontSize',fsz)% .

fig7= figure();
    width = 12;     % Width in inches,   find column width in paper 
    height = 10;    % Height in inches
%     set(fig7, 'Position', [fig7.Position(1:2)-[width,height]*100 width*100, height*100]); %<- Set size
    set(fig7, 'Position', [100,100 width*100, height*100]); %<- Set size

vid = VideoWriter(['figs' filesep animation_name ],'MPEG-4');
vid.Quality = 100;
vid.FrameRate = 60;
% v.FileFormat = 'mp4';

open(vid);

colormap(strainScheme)%     colorbar
% for j = 1:10:len10
for j = 1:1:len10
% for j = 1:1:30
    % -------------------------------------------------------------    
    subplot(4,2,[1:2:5])
    
    frame = start1 + round(j/10);
%     frame = start1+j
    
    plot3( [0,1]*4e3,[0,1]*0,[0,1]*0 ,'k') 
    hold on
    plot3( [0,1]*0,[0,-1]*4e3,[0,1]*0 ,'k') 
    plot3( [0,1]*0,[0,1]*0,[0,1]*4e3 ,'k') 
     text(4.3e3,0,0,'X')
     text(0,-5.3e3,0,'Y')
     text(-100,0,4.5e3,'Z')
    angles = [0,  phiInterp(j),0];
    
    for jj = 1:size(xBulb,1)
        for k = 1:size(xBulb,2)
            xyzTemp = [xBulb(jj,k), yBulb(jj,k), zBulb(jj,k) ];

                    eul_1 = [ 1       0                     0;...
                                0,  cos( angles(1) ),  - sin(  angles(1) )  ; ...
                               0   sin(  angles(1) ) cos(  angles(1) )  ]^-1;
                    eul_2 = [cos(  angles(2))      0       sin( angles(2)) ; ...
                                0              1       0 ;...
                                -sin( angles(2))    0        cos(  angles(2))]^-1;

                    eul_3 = [cos(  angles(3))   sin(  angles(3))    0 ; ...
                                 -sin(  angles(3)) cos(  angles(3))     0 ;...
                                 0          0               1]^-1;
            xyzT = eul_1*eul_2*eul_3*xyzTemp'; 
            xBulbTranspose(jj,k) = xyzT(1);
            yBulbTranspose(jj,k) = xyzT(2);
            zBulbTranspose(jj,k) = xyzT(3);
        end
    end

%     s2 = surf(Xb(:,:,frame),Yb(:,:,frame),Zb(:,:,frame),Cb(:,:,frame));
%     Cb = ones(size(Cb))*0.2;
%     s2 = surf( squeeze( Xb(:,:,frame)), ...
%         squeeze( Yb(:,:,frame)),...
%         squeeze( Zb(:,:,frame)),...
%         squeeze( Cb(:,:,frame)));
%     
    sStalk = surf( squeeze( XbStalk(:,:,j)), ...
        squeeze( YbStalk(:,:,j)),...
        squeeze( ZbStalk(:,:,j)),...
        squeeze( Cb(:,:,1)));
    
    
    
    
%         set(s2,surfParamForeground{:})
%     centerPoint = [( FEA(1 ).xyzPoints(ImaxBottom,frame,1) + FEA(1 ).xyzPoints(ImaxTop,frame,1))/2;
%                 ( FEA(1 ).xyzPoints(ImaxBottom,frame,2) + FEA(1 ).xyzPoints(ImaxTop,frame,2))/2;
%                 ( FEA(1 ).xyzPoints(ImaxBottom,frame,3) + FEA(1 ).xyzPoints(ImaxTop,frame,3))/2];
    s3 = surf( xBulbTranspose ,...
                yBulbTranspose  ,...
                zBulbTranspose  ,...
                C); 
        set(sStalk,surfParamForeground{:})
        set(s3,surfParamForeground{:})
        axis equal
        caxis([-1,1]*0.0015)
        axis([-300 5500 -5500 500 -5500 5500])
        shading faceted
        view(30,15)
        axis off
        hold off 



    
    % -------------------------------------------------------------    
    theta = 0:0.001:pi;
    subplot(4,2,[ 7])
        plot( -sin(theta)*100  +120 ,cos(theta)*100 ,'k')
    hold on 
        plot( sin(theta)*100 +120 ,cos(theta)*100 ,'k')
        plot( -sin(theta)*100-120   ,cos(theta)*100 ,'k')
        plot( sin(theta)*100-120,cos(theta)*100 ,'k')

      for k = 1:16
        spike_ind_temp = FEA(1).spikeInds{k};
        which_in_range = (spike_ind_temp<=It10(j)) & (spike_ind_temp>=It10(1));
        for kk = fliplr( 1:sum(nonzeros(which_in_range)) )
            which_ = which_in_range;
            nonz = find(which_in_range);
            if sum(nonzeros(which_in_range)) >1
            which_(nonz(kk) ) = 0;
            end
            try 
                if any(spike_ind_temp) && ( t_plot10(j) > t_plot10(  spike_ind_temp(which_) - start10 ) )   && ( t_plot10(j) < 5+ t_plot10(  spike_ind_temp(which_) - start10 ) ) 
                    dt = t_plot10(j) - t_plot10(  spike_ind_temp(which_) - start10 );
                     dots = scatter( 100/150*FEA(1).xyz( FEA(1).circleInds( k ),2 )-120  ,  100/150*FEA(1).xyz( FEA(1).circleInds( k ),3 ),20 ,'k','filled');
                    dots.MarkerFaceAlpha = rectif( 1-dt/5);
                end
            end
        end

        spike_ind_temp = FEA(2).spikeInds{k};
        which_in_range = (spike_ind_temp<=It10(j)) & (spike_ind_temp>=It10(1));
        for kk = fliplr( 1:sum(nonzeros(which_in_range)) )
            which_ = which_in_range;
            nonz = find(which_in_range);
            if sum(nonzeros(which_in_range)) >1
            which_(nonz(kk) ) = 0;
            end
            try 
                if any(spike_ind_temp) && ( t_plot10(j) > t_plot10(  spike_ind_temp(which_) - start10 ) )   && ( t_plot10(j) < 5+ t_plot10(  spike_ind_temp(which_) - start10 ) ) 
                    dt = t_plot10(j) - t_plot10(  spike_ind_temp(which_) - start10 );
                     dots = scatter( 100/150*FEA(1).xyz( FEA(1).circleInds( k ),2 )+120  ,  100/150*FEA(1).xyz( FEA(1).circleInds( k ),3 )  ,20,'r','filled');
                    dots.MarkerFaceAlpha = rectif( 1-dt/5);
                end
            end
        end


    end
        axis equal
        axis([-220,220,-100,100])
        axis off
        hold off
        
% -------------------------------------------------------------   
    subplot(4,2,2:2:8) ;
        plot(  t_plot1, FEA(1).phi(1,It1) +10 ,'k')
        hold on 
%         rectangle('Position',[1,2.9,38,4.8],'Curvature',0,'FaceColor',[1,1,1]*0.95)
%         rectangle('Position',[1,-2.6,38,4.8],'Curvature',0,'FaceColor',[1,1,1]*0.95)

        plot( [1,1]*t_plot10(j), [-3,12],'k') 

    for kl =1:length(spike_order)
        k = spike_order( kl ) ;
            spike_ind_temp = FEA(1).spikeInds{k};
                which_ = (spike_ind_temp<=It10(j)) & (spike_ind_temp>=It10(1));
        if any(spike_ind_temp) && kl <=9
           plot( [1,1]'*t_plot10(  spike_ind_temp(which_) - start10 ),...
               ([9 ; 8.6] - kl*0.5 -1)  * ones(1, length(  spike_ind_temp(which_) )),'k')
        elseif any(spike_ind_temp) 
           plot( [1,1]'*t_plot10(  spike_ind_temp(which_) - start10 ),...
               ([9 ; 8.6] - kl*0.5 - 2)  * ones(1, length(  spike_ind_temp(which_) )),'k')
        end

        scatter( 100/150*FEA(1).xyz( FEA(1).circleInds( [5:13, 13:16, 1:5] ),2 ) +120 ,  100/150*FEA(1).xyz( FEA(1).circleInds( [5:13, 13:16, 1:5]),3 ) +900 ,'k','filled')

        spike_ind_temp = FEA(2).spikeInds{k};
            which_ = (spike_ind_temp<=It10(j)) & (spike_ind_temp>=It10(1));

        if any(spike_ind_temp) && (kl <=9)
           plot( [1,1]'*t_plot10(  spike_ind_temp(which_) - start10 ),...
               ([9 ; 8.6] - kl*0.5 -1)  * ones(1, length(  spike_ind_temp(which_) )),'r')
        elseif any(spike_ind_temp) 
           plot( [1,1]'*t_plot10(  spike_ind_temp(which_) - start10 ),...
               ([9 ; 8.6] - kl*0.5 - 2)  * ones(1, length(  spike_ind_temp(which_) )),'r')
        end
    end

        hold off 
        axis([0,25,-3,12])
        xlabel('Time (ms)')
        ylab = ylabel('Flapping angle $\phi(t)$', 'Position', [ -2 , 10, 0]); 
        set(gca(),'YTick', [10-pi/2,10,10+pi/2] ,'YTickLabel',{'$\frac{\pi}{2}$',0,'$\frac{\pi}{2}$' } )
        title('Spikes along circumference')
        drawnow
%         set(ylab, 'Units', 'Normalized');
        set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
%         set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
   frame = getframe(gcf);
   writeVideo(vid,frame);

end
close(vid);
