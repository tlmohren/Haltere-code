clc;clear all;close all

addpathFolderStructureHaltere()
run('config_file.m')

%%
loadName = 'figure6_strainData';
saveName = 'figure6_strainData';

renew_data_load = false
% renew_data_load = true
outer_dimension = [150,150,...
                    150,150,...
                    150,150,...
                    150,150,...
                    150,150,...
                    150,150,...
                    243.435,243.435];
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
        [FEA(j).xyz, FEA(j).strain, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'eXX','eYY','eZZ','eXY','eXZ','eYZ' });        toc 
    end
    circleDistance = 300;               % distance from base to haltere 
    er = 0.01;
    for j = 1:length(FEA)
        sideL = outer_dimension(j); 
        xMatch = find(   abs(abs( FEA(j).xyz(:,1) ) - circleDistance)  <er );
        yMatch = find(   abs(abs( FEA(j).xyz(:,2) ) - sideL) <er & ...
                     abs( FEA(j).xyz(:,3) )   <er);
        zMatch = find(   abs(abs( FEA(j).xyz(:,3) ) - sideL)  <er & ...
                     abs( FEA(j).xyz(:,2) )   <er);
        FEA(j).sideInds = intersect(xMatch,yMatch);
        FEA(j).topInds = intersect(xMatch,zMatch);
    end
    save(['data' filesep saveName],'FEA')
else
    load(['data' filesep loadName],'FEA')
end

%% find principle strains

for j = 1:length(FEA)
%     for k = 1:2
       exTop = FEA(j).strain( FEA(j).topInds(1), : , [1] );
       eyTop = FEA(j).strain( FEA(j).topInds(1), : , [2] );
       exyTop = FEA(j).strain( FEA(j).topInds(1), : , [4] );
        
%        FEA(j).epTop(1,:) = 0.5*(exTop+eyTop) + sqrt( (0.5*(exTop-eyTop)).^2 + (exyTop).^2 );
%        FEA(j).epTop(2,:) = 0.5*(exTop+eyTop) - sqrt( (0.5*(exTop-eyTop)).^2 + (exyTop).^2 );
%        FEA(j).thetapTop(1,:) = 0.5*atan( exyTop./(exTop-eyTop) );
%        
    for k = 1:length( exTop )
       FEA(j).epTop(1,k) = (exTop(k)+eyTop(k))/2 + sqrt( ((exTop(k)-eyTop(k))/2).^2 + (exyTop(k)/2).^2 );
       FEA(j).epTop(2,k) = (exTop(k)+eyTop(k))/2 - sqrt( ((exTop(k)-eyTop(k))/2).^2 + (exyTop(k)/2).^2 );
%        FEA(j).thetapTop(1,k) = atan( exyTop(k)/ (exTop(k)-eyTop(k)) )/2;
       FEA(j).thetapTop(1,k) = atan2( exyTop(k), (exTop(k)-eyTop(k)) )/2;
       
       if abs( FEA(j).epTop(2,k)) > abs(FEA(j).epTop(1,k)) 
           temp = FEA(j).epTop(1,k);
          FEA(j).epTop(1,k) = FEA(j).epTop(2,k);
          FEA(j).epTop(2,k) = temp;
          FEA(j).thetapTop(1,k) = FEA(j).thetapTop(1,k);
          FEA(j).thetapTop(1,k) = FEA(j).thetapTop(1,k)+ pi/2*sign(temp);
       end
       if FEA(j).thetapTop(1,k) >  pi*0.999 % 
           FEA(j).thetapTop(1,k) = FEA(j).thetapTop(1,k)-pi;
           
%             epTemp = FEA(j).epTop(1,k) ;
%              FEA(j).epTop(1,k) = FEA(j).epTop(2,k) ;
%              FEA(j).epTop(2,k)  = epTemp;
       end
    end
       
       
       
       
       
       exSide = FEA(j).strain( FEA(j).sideInds(1), : , [1] );
       ezSide = FEA(j).strain( FEA(j).sideInds(1), : , [3] );
       exzSide = FEA(j).strain( FEA(j).sideInds(1), : , [5] );
        
    for k = 1:length( exSide )
       FEA(j).epSide(1,k) = 0.5*(exSide(k)+ezSide(k)) + sqrt( (0.5*(exSide(k)-ezSide(k))).^2 + (exzSide(k)).^2 );
       FEA(j).epSide(2,k) = 0.5*(exSide(k)+ezSide(k)) - sqrt( (0.5*(exSide(k)-ezSide(k))).^2 + (exzSide(k)).^2 );
%        FEA(j).thetapSide(1,k) = atan( exzSide(k)/ (exSide(k)-ezSide(k)) )/2;
       FEA(j).thetapSide(1,k) = atan2( exzSide(k), (exSide(k)-ezSide(k)) )/2;
       
       
       if abs( FEA(j).epSide(2,k)) > abs(FEA(j).epSide(1,k)) 
           temp = FEA(j).epSide(1,k);
          FEA(j).epSide(1,k) = FEA(j).epSide(2,k);
          FEA(j).epSide(2,k) = temp;
%           FEA(j).thetapSide(1,k) = FEA(j).thetapSide(1,k);
          FEA(j).thetapSide(1,k) = FEA(j).thetapSide(1,k)+ pi/2*sign(temp);
       end
%        if FEA(j).thetapSide(1,k) <  0% 
%            FEA(j).thetapSide(1,k) = FEA(j).thetapSide(1,k)+pi/2;
%            
%             epTemp = FEA(j).epSide(1,k) ;
%              FEA(j).epSide(1,k) = FEA(j).epSide(2,k) ;
%              FEA(j).epSide(2,k)  = epTemp;
%        end
    end
end
 

 j = 2;

figure()
subplot(211); hold on
plot(FEA(j).epSide(1,It))
plot(FEA(j).epSide(2,It))
subplot(212)
plot(FEA(j).thetapSide(1,It) )


%% order of magnitude comparison 

len = 101;
start = 35;
It = start:(start+len-1);
t_plot = (0:len-1)*0.001;

axOptsStrainTop = {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)],'XTickLabel',{'','',''} ,'YLim',[-1,1]*2.5e-3,'YTick',[-1,0,1]*2e-3}; 
axOptsStrainSide = {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)] ,'YLim',[-1,1]*1.5e-4}; 

strainLabel = {'$\epsilon$'};
        cols = linspecer(4); 
for k = 1:length(FEA)/2

        fig1 = figure();
            width = 2;     % Width in inches,   find column width in paper 
            height = 2;    % Height in inches
            set(fig1, 'Position', [fig1.Position(1:2)-[0,300] width*100, height*100]); %<- Set size
    for kk = 1:2
        j = 2*(k-1) + kk;
%         kk
%         title( FEA(j).name)
        subplot(2,2,kk); hold on 
%         subplot(4,2,kk); hold on 
%             p1 = plot(t_plot, FEA(j).strain( FEA(j).topInds(1), It) );
%             p2 = plot(t_plot, FEA(j).strain( FEA(j).topInds(2), It) );
            p1 = plot(t_plot, FEA(j).epTop( 1, It) );
            p2 = plot(t_plot, FEA(j).epTop( 2, It) );
%             p2 = plot(t_plot, FEA(j).strain( FEA(j).topInds(2), It) );
%                 ylabel( '$\epsilon_\text{Principle~TOP}$','Rotation',0 );
                ylabel( '\epsilon_{Principle TOP}','Interpreter','tex');
                ax = gca();
                set(ax,axOptsStrainTop{:})
                set(p1,'Color',cols(1,:))
                set(p2,'Color',cols(2,:))
                
            set(gca,'XAxisLocation','top')
            xlabel( FEA(j).name(18:end),'Interpreter','none')
        subplot(2,2,kk+2); hold on 
%         subplot(4,2,kk+2); hold on 
            p3 = plot(t_plot, FEA(j).epSide( 1, It) );
            p4 = plot(t_plot, FEA(j).epSide( 2, It) );
%             p3 = plot(t_plot, FEA(j).strain( FEA(j).sideInds(1), It) );
%             p4 = plot(t_plot, FEA(j).strain( FEA(j).sideInds(2), It) );
                xlabel('Time (s)')
                ylabel( '\epsilon_{Principle SIDE}' ,'Interpreter','tex');
                ax = gca();
                set(ax,axOptsStrainSide{:})
                set(p3,'Color',cols(3,:))
                set(p4,'Color',cols(4,:))

                
%                 
%         subplot(4,2,kk+4); hold on 
% %             p1 = plot(t_plot, FEA(j).strain( FEA(j).topInds(1), It) );
% %             p2 = plot(t_plot, FEA(j).strain( FEA(j).topInds(2), It) );
%             p1 = plot(t_plot, FEA(j).thetapTop( 1, It) );
%             
%         subplot(4,2,kk+6); hold on 
%             
%             p1 = plot(t_plot, FEA(j).thetapSide( 1, It) );
% %         subplot(3,2,kk+4); hold on 
% %             p3 = plot(t_plot, FEA(j).thetapSide( 1, It) );
% % %             p4 = plot(t_plot, FEA(j).thetapSide( 2, It) );
% % %             p3 = plot(t_plot, FEA(j).strain( FEA(j).sideInds(1), It) );
% % %             p4 = plot(t_plot, FEA(j).strain( FEA(j).sideInds(2), It) );
% %                 xlabel('Time (s)')
% %                 ylabel( '\epsilon_{Principle SIDE}' ,'Interpreter','tex');
% %                 ax = gca();
% % %                 set(ax,axOptsStrainSide{:})
% %                 set(p3,'Color',cols(3,:))
% % %                 set(p4,'Color',cols(4,:))

    end
    % Setting paper size for saving 
    set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
    set(fig1,'InvertHardcopy','on');
    set(fig1,'PaperUnits', 'inches');
    papersize = get(fig1, 'PaperSize');
    left = (papersize(1)- width)/2;
    bottom = (papersize(2)- height)/2;
    myfiguresize = [left, bottom, width, height];
    set(fig1, 'PaperPosition', myfiguresize);
    print(fig1, ['figs' filesep 'Figure6_princpleStrain_' FEA(j).name(18:end-5) ,'OM0_OM10'], '-dpng', '-r600');
    stupid_ratio = 15/16;
    myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
    set(fig1, 'PaperPosition', myfiguresize);
    print(fig1, ['figs' filesep 'Figure6_princpleStrain_' FEA(j).name(18:end-5) ,'OM0_OM10'], '-dsvg', '-r600');
end
        
% % % % % % % 
% % % % % % 
% % % % % % 
% % % % % % 
% % % % % % 
% % % % % % strainLabel = {'$\epsilon$'};
% % % % % % axOptsStrainTop = {'XGrid','On','XLim',[0,0.15],'XTick',[0:0.05:t_plot(end)],'XTickLabel',{'','','',''},'YLim',[-1,1]*2e-3,'YTick',[-1,0,1]*2e-3}; 
% % % % % % axOptsStrainSide = {'XGrid','On','XLim',[0,0.15],'XTick',[0:0.05:0.15],'YLim',[-1,1]*1e-4}; 
% % % % % % 
% % % % % % cols = linspecer(4); 
% % % % % % for j = [1,3,5]
% % % % % %     subplot(211); hold on 
% % % % % %         p1 = plot(t_plot, FEA(j).strain( FEA(j).topInds(1), It) );
% % % % % %         p2 = plot(t_plot, FEA(j).strain( FEA(j).topInds(2), It) );
% % % % % %         ylabel( strainLabel,'Rotation',0 );
% % % % % %         ax = gca();
% % % % % %         set(ax,axOptsStrainTop{:})
% % % % % %     subplot(212); hold on 
% % % % % %         p3 = plot(t_plot, FEA(j).strain( FEA(j).sideInds(1), It) );
% % % % % %         p4 = plot(t_plot, FEA(j).strain( FEA(j).sideInds(2), It) );
% % % % % %         xlabel('Time (s)')
% % % % % %         ylabel( strainLabel,'Rotation',0);
% % % % % %         ax = gca();
% % % % % %         set(ax,axOptsStrainSide{:})
% % % % % % end
% % % % % % set(p1,'Color',cols(1,:))
% % % % % % set(p2,'Color',cols(2,:))
% % % % % % set(p3,'Color',cols(3,:))
% % % % % % set(p4,'Color',cols(4,:))

%% 


    fig1 = figure();
        width = 2;     % Width in inches,   find column width in paper 
        height = 3;    % Height in inches
        set(fig1, 'Position', [fig1.Position(1:2) width*100, height*100]); %<- Set size

for k = 4% 1:length(FEA)/2

%         fig1 = figure();
%             width = 3;     % Width in inches,   find column width in paper 
%             height = 3;    % Height in inches
%             set(fig1, 'Position', [fig1.Position(1:2) width*100, height*100]); %<- Set size
    for kk = 1:2
        j = 2*(k-1) + kk;
            subplot(3,2,kk); hold on 
        %             p1 = plot(t_plot, FEA(j).strain( FEA(j).topInds(1), It) );
        %             p2 = plot(t_plot, FEA(j).strain( FEA(j).topInds(2), It) );
                p1 = plot(t_plot, FEA(j).strain( FEA(j).sideInds(1), It, 1) );
        %             p2 = plot(t_plot, FEA(j).epTop( 2, It) );
        %             p2 = plot(t_plot, FEA(j).strain( FEA(j).topInds(2), It) );
                    ylabel( strainLabel,'Rotation',0 );
        %                 ax = gca();
        %                 set(ax,axOptsStrainTop{:})
        %                 set(p1,'Color',cols(1,:))
        %                 set(p2,'Color',cols(2,:))
            subplot(3,2,kk+2); hold on 
                p1 = plot(t_plot, FEA(j).strain( FEA(j).sideInds(1), It, 3) );
        %             p3 = plot(t_plot, FEA(j).epSide( 1, It) );
        %             p4 = plot(t_plot, FEA(j).epSide( 2, It) );
        %             p3 = plot(t_plot, FEA(j).strain( FEA(j).sideInds(1), It) );
        %             p4 = plot(t_plot, FEA(j).strain( FEA(j).sideInds(2), It) );
                    xlabel('Time (s)')
                    ylabel( strainLabel,'Rotation',0);
            subplot(3,2,kk+4); hold on 
                p1 = plot(t_plot, FEA(j).strain( FEA(j).sideInds(1), It, 5) );
%                 ax = gca();
%                 set(ax,axOptsStrainSide{:})
%                 set(p3,'Color',cols(3,:))
%                 set(p4,'Color',cols(4,:))
    end
end

%% 


% for j = 1:4
% %         exTop = FEA(j).strain( FEA(j).topInds(1), : , [1] );
% %        eyTop = FEA(j).strain( FEA(j).topInds(1), : , [2] );
% %        exyTop = FEA(j).strain( FEA(j).topInds(1), : , [4] );
% %         
% %        FEA(j).epTop(1,:) = 0.5*(exTop+eyTop) + sqrt( (0.5*(exTop-eyTop)).^2 + (exyTop).^2 );
% %        FEA(j).epTop(2,:) = 0.5*(exTop+eyTop) - sqrt( (0.5*(exTop-eyTop)).^2 + (exyTop).^2 );
% %        FEA(j).thetapTop(1,:) = 0.5*atan( exyTop./(exTop-eyTop) );
%        
%        exSide = FEA(j).strain( FEA(j).sideInds(1), : , [1] );
%        ezSide = FEA(j).strain( FEA(j).sideInds(1), : , [3] );
%        exzSide = FEA(j).strain( FEA(j).sideInds(1), : , [5] );
% % %        FEA(j).thetapSideNoShear(k,:) = 0.5*atan( FEA(j).strain( FEA(3).sideInds(k), : , [5] )./(exSide-eySide) );
% 
%     for k = 1:length( exSide )
%        FEA(j).epSide(1,k) = (exSide(k)+ezSide(k))/2 + sqrt( ((exSide(k)-ezSide(k))/2).^2 + (exzSide(k)/2).^2 );
%        FEA(j).epSide(2,k) = (exSide(k)+ezSide(k))/2 - sqrt( ((exSide(k)-ezSide(k))/2).^2 + (exzSide(k)/2).^2 );
%        
%        FEA(j).thetapSide(1,k) = atan2( exzSide(k), (exSide(k)-ezSide(k)) )/2;
% %        FEA(j).thetapSide(1,k) = atan( exzSide(k) /  (exSide(k)-ezSide(k)) )/2;
% %        if FEA(j).thetapSide(1,k) <  0% 0.001%pi*1/1000
% % %            display(k)
% %            FEA(j).thetapSide(1,k) = FEA(j).thetapSide(1,k)+pi/2;
% %            
% %             epTemp = FEA(j).epSide(1,k) ;
% %              FEA(j).epSide(1,k) = FEA(j).epSide(2,k) ;
% %              FEA(j).epSide(2,k)  = epTemp;
% %        end
%     end
% end
%
%% 


% %  maxE = max( [ FEA(1).strain( FEA(j).sideInds(1), : , [1] ) , FEA(2).strain( FEA(j).sideInds(1), : , [1] )]);
% % figure('Position',[100,100,500,900])
% % for kk =100:200%length(exSide)
% % 
% %     
% %     for j = 1:2
% %        exSide = FEA(j).strain( FEA(j).sideInds(1), : , [1] );
% %        ezSide = FEA(j).strain( FEA(j).sideInds(1), : , [3] );
% %        exzSide = FEA(j).strain( FEA(j).sideInds(1), : , [5] );
% %     
% %     subplot(2,1,j)
% %     title( ['Time is = ', num2str(kk/1e3) ])
% %     hold off
% %     plot( [0  exSide(kk)] ,[0,0], 'r')
% %     hold on 
% %     plot( [0,0],[0  ezSide(kk)] , 'b')
% %     plot( [0, 1]*exzSide(kk),[0 1]*exzSide(kk) ,'g' )
% %     
% %     plot(  [0, cos(FEA(j).thetapSide(1,kk) )]*FEA(j).epSide(1,kk)  ,  [0,sin(FEA(j).thetapSide(1,kk) )]*FEA(j).epSide(1,kk) ,'k') 
% %     plot(  [0, -cos( FEA(j).thetapSide(1,kk) +pi/2 )]*FEA(j).epSide(1,kk)  ,  [0, sin( FEA(j).thetapSide(1,kk) + pi/2  )]*FEA(j).epSide(2,kk) ,'--k') 
% %     axis( [-1,1,-1,1]*maxE)
% %     axis square
% %    
% %     
% % %     subplot(211)
% % %     title( ['Time is = ', num2str(kk/1e3) ])
% % %     hold off
% % %     plot( [0  exSide(kk)] ,[0,0], 'r')
% % %     hold on 
% % %     plot( [0,0],[0  ezSide(kk)] , 'b')
% % %     plot( [0, 1]*exzSide(kk),[0 1]*exzSide(kk) ,'g' )
% % %     
% % %     plot(  [0,-cos(FEA(j).thetapSide(1,kk) )]*FEA(j).epSide(1,kk)  ,  [0,sin(FEA(j).thetapSide(1,kk) )]*FEA(j).epSide(1,kk) ,'k') 
% % %     plot(  [0,  cos( FEA(j).thetapSide(1,kk) +pi/2 )]*FEA(j).epSide(1,kk)  ,  [0,sin( FEA(j).thetapSide(1,kk) +pi/2  )]*FEA(j).epSide(2,kk) ,'k') 
% % %     axis( [-1,1,-1,1]*max(exSide))
% % %    
% %     
% %     
% %     drawnow
% %     
% %     pause(0.05)
% %     end
% % end
% % %
% % subplot(211)
% %      j = 1
% %     scatter( [cos(FEA(j).thetapSide(1,:) )].*FEA(j).epSide(1,:)  , [sin(FEA(j).thetapSide(1,:) )].*FEA(j).epSide(1,:) ,'ok') 
% %     scatter( [-cos(FEA(j).thetapSide(1,:) +pi/2  )].*FEA(j).epSide(2,:)  , [sin(FEA(j).thetapSide(1,:) +pi/2 )].*FEA(j).epSide(2,:) ,'or') 
% % 
% % subplot(212)
% %      j = 2
% %     scatter( [cos(FEA(j).thetapSide(1,:) )].*FEA(j).epSide(1,:)  , [sin(FEA(j).thetapSide(1,:) )].*FEA(j).epSide(1,:) ,'ok') 
% %     scatter( [-cos(FEA(j).thetapSide(1,:) +pi/2  )].*FEA(j).epSide(2,:)  , [sin(FEA(j).thetapSide(1,:) +pi/2 )].*FEA(j).epSide(2,:) ,'or')
% % %     plot(  [0,  cos( FEA(j).thetapSide(1,kk) +pi/2 )]*FEA(j).epSide(1,kk)  ,  [0,sin( FEA(j).thetapSide(1,kk) +pi/2  )]*FEA(j).epSide(2,kk) ,'k') 
% % 
