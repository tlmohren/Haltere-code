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
    FEA(7).name = 'Haltere_CraneFly_SphereVerOffset_Om0';
    FEA(8).name = 'Haltere_CraneFly_SphereVerOffset_Om10';
    FEA(9).name = 'Haltere_CraneFly_ellipsoidVerOffset_Om0';
    FEA(10).name = 'Haltere_CraneFly_ellipsoidVerOffset_Om10';
    FEA(11).name = 'Haltere_CraneFly_ellipsoidHorOffset_Om0';
    FEA(12).name = 'Haltere_CraneFly_ellipsoidHorOffset_Om10';
    FEA(13).name = 'Haltere_CraneFly_ellipsoidVerCrossStalk_Om0';
    FEA(14).name = 'Haltere_CraneFly_ellipsoidVerCrossStalk_Om10';
    for j =  1:length(FEA)
        tic
        [FEA(j).xyz, FEA(j).strain, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'eXX' });        toc 
    end
    % Determine Circle locations
    for j = 1:length(FEA)-2
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
        
    end
    er = 0.01;
    for j = [13,14]
       outer_dimension = [ 243.435,243.435];
       
        sideL = outer_dimension(j-12); 
        xMatch = find(   abs(abs( FEA(j).xyz(:,1) ) - circleDistance)  <er );
        yMatch = find(   abs(abs( FEA(j).xyz(:,2) ) - sideL) <er & ...
                     abs( FEA(j).xyz(:,3) )   <er);
        zMatch = find(   abs(abs( FEA(j).xyz(:,3) ) - sideL)  <er & ...
                     abs( FEA(j).xyz(:,2) )   <er);
        FEA(j).sideInds = intersect(xMatch,yMatch);
        FEA(j).topInds = intersect(xMatch,zMatch);
        
        FEA(j).circleIndsUnsorted = [FEA(j).sideInds; 
                                    FEA(j).topInds];
        
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

calib_param_max = [0.00298239604088481,0.0212443192274280,0.0372780009360082,0.0479350861876651,0.0514646550697229,0.0475135494972593,0.0363905581996612,0.0198191699014800,0.00299958384555981,0.0201797042231435,0.0370669398494605,0.0483996355340012,0.0524251480251299,0.0487363345563017,0.0377914840661745,0.0213792050886826];

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
%     for j = 1:2%length(FEA)
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

    
    len1 = 26;
    start1 = 166;
    It1 = start1:(start1+len1-1);
%     t_plot1 = (0:len1-1)*0.001;
%     t_plot1 = linspace(0,1,len1);
    t_plot1 = linspace(0,40,len1);

    len10 = 261;
    start10 = 1651;
    It10 = start10:(start10+len10-1);
%     t_plot10 = (0:len10-1)*0.0001;
%     t_plot10 = linspace(0,1,len10);
    t_plot10 = linspace(0,40,len10);
        
    spike_order = [5:13, 13:16, 1:5];

    
    
%% 
   
%     len1 = 26;
%     start1 = 151;
%     It1 = start1:(start1+len1-1);
%     t_plot1 = (0:len1-1)*0.001;
% 
%     len10 = 261;
%     start10 = 1501;
%     It10 = start10:(start10+len10-1);
%     t_plot10 = (0:len10-1)*0.0001;
%         
        
    spike_order = [5:13, fliplr([13:16, 1:5]) ];

    fig7= figure();
        width = 5;     % Width in inches,   find column width in paper 
        height = 5;    % Height in inches
        set(fig7, 'Position', [fig7.Position(1:2)-[width,height]*100 width*100, height*100]); %<- Set size
        
    subplot(1,2,1); hold on 
        scatter( FEA(1).xyz( FEA(1).circleInds(5:13),2 ),  FEA(1).xyz( FEA(1).circleInds(5:13),3 ) +150+350,'k','filled')

        scatter( FEA(1).xyz( FEA(1).circleInds( [13:16, 1:5] ),2 ),  FEA(1).xyz( FEA(1).circleInds( [13:16, 1:5] ),3 ) +150 ,'k','filled')
        theta = 0:0.001:pi;
        plot( -sin(theta)*150,cos(theta)*150 +150+350,'k')
        plot( sin(theta)*150,cos(theta)*150 +150,'k')

        axis equal
        axis([-180,180,-10,900])
    
    subplot(122) ; hold on 
    
    
    
rectangle('Position',[1,2.9,38,4.8],'Curvature',0,'FaceColor',[1,1,1]*0.95)
rectangle('Position',[1,-2.6,38,4.8],'Curvature',0,'FaceColor',[1,1,1]*0.95)

% for j = 1:20
% rectangle('Position',[2*j,-3,1,10],'Curvature',0,'FaceColor',[1,1,1]*0.8,'EdgeColor','none')
% end
        plot(  t_plot1, FEA(1).phi(1,It1) +10 ,'k')
%     subplot(133) ; hold on 
%         plot(  t_plot1, FEA(1).phi(1,It1) +10 ,'k')
        
    for kl =1:length(spike_order)
%         kl
        k = spike_order( kl ) ;
%         hold on
        for j = [1]%[1,2]+(2*(jj-1))
%             subplot(1,32)
            spike_ind_temp = FEA(j).spikeInds{k};
                which_ = (spike_ind_temp<=It10(end)) & (spike_ind_temp>=It10(1));
                
            if any(spike_ind_temp) & kl <=9
               plot( [1,1]'*t_plot10(  spike_ind_temp(which_) - start10 ),...
                   ([9 ; 8.6] - kl*0.5 -1)  * ones(1, length(  spike_ind_temp(which_) )),'k')
            elseif any(spike_ind_temp) 
               plot( [1,1]'*t_plot10(  spike_ind_temp(which_) - start10 ),...
                   ([9 ; 8.6] - kl*0.5 - 2)  * ones(1, length(  spike_ind_temp(which_) )),'k')
            end
        end
        for j = [2]%[1,2]+(2*(jj-1))
%             subplot(1,32)
            spike_ind_temp = FEA(j).spikeInds{k};
                which_ = (spike_ind_temp<=It10(end)) & (spike_ind_temp>=It10(1));
                
            if any(spike_ind_temp) & kl <=9
               plot( [1,1]'*t_plot10(  spike_ind_temp(which_) - start10 ),...
                   ([9 ; 8.6] - kl*0.5 -1)  * ones(1, length(  spike_ind_temp(which_) )),'r')
            elseif any(spike_ind_temp) 
               plot( [1,1]'*t_plot10(  spike_ind_temp(which_) - start10 ),...
                   ([9 ; 8.6] - kl*0.5 - 2)  * ones(1, length(  spike_ind_temp(which_) )),'r')
            end
        end
        if kl == 1
               plot( [1,1]'*t_plot10(  spike_ind_temp(which_) - start10 ),...
                   ([-3,10])  * ones(1, length(  spike_ind_temp(which_) )),':k')
        elseif kl == 9
               plot( [1,1]'*t_plot10(  spike_ind_temp(which_) - start10 ),...
                   ([-3,10])  * ones(1, length(  spike_ind_temp(which_) )),':k')
        end
    end
    
%     for j = 
%       plot( [1,1]'*t_plot10(  spike_ind_temp(which_) - start10 ),...
%     axis([0,0.025,-5,12])
%     subplot(132)
%     axis([0,1,-3,12])
    axis([0,40,-3,12])
    xlabel('Time (ms)')
    ylabel('Flapping angle $\phi(t)$')
%     ax = gca();
    set(gca(),'YTick', [10-pi/2,10,10+pi/2] ,'YTickLabel',{'$\frac{\pi}{2}$',0,'$\frac{\pi}{2}$' } )
    title('Spikes along circumference')
    
    
    
    set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
    set(fig7,'InvertHardcopy','on');
    set(fig7,'PaperUnits', 'inches');
    papersize = get(fig7, 'PaperSize');
    left = (papersize(1)- width)/2;
    bottom = (papersize(2)- height)/2;
    myfiguresize = [left, bottom, width, height];
    set(fig7, 'PaperPosition', myfiguresize);
    print(fig7, ['figs' filesep 'Figure7_spikeSphere' ], '-dpng', '-r600');
    stupid_ratio = 15/16;
    myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
    set(fig7, 'PaperPosition', myfiguresize);
    print(fig7, ['figs' filesep 'Figure7_spikeSphere'  ], '-dsvg', '-r600');
