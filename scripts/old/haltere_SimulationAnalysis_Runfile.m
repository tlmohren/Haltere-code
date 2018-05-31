%-------------------------------
% TMohren
% load comsol data and plot resulting strains
% 2017-08-09
%------------------------------
clc; clear all; close all
addpathFolderStructureHaltere()

sim.name = 'Haltere_CraneFly_Sphere';
load(['data' filesep sim.name filesep sim.name '_allData'])

% Simulation parameters
circleDistance = 300;               % distance from base to haltere 
circleRadius = 150;                 % radius of haltere
strainPoints = [circleDistance, 0,   circleRadius; ...
                circleDistance, 0,   - circleRadius; ...
               circleDistance,  circleRadius, 0;...
               circleDistance, -circleRadius, 0];

%% Analyze strain
circleIndices = [];

% pointIndices= findPointIndices( round(sim.strainXYZ,7) , strainPoints );
circleIndices= findCircleIndices( round(sim.strainXYZ,7) , circleDistance,circleRadius);
xyz = sim.strainXYZ;
angle = atan2( xyz(3,circleIndices), xyz(2,circleIndices));
angleDeg = rad2deg(angle)-180;
angleDeg(angleDeg<0) = angleDeg(angleDeg<0)+360;
[V,I_sort] = sort(angleDeg,'ascend');

Ind = circleIndices(I_sort);
sidePoints = Ind( find( mod(V,90) == 0));


%% 
dotStyle = {'Marker','.'};
markerStyle = {'Marker','o','MarkerFaceColor','red','CData',[1,0,0]};
% 
figure()
    plotMesh( sim.strainXYZ ,dotStyle) 
    plotMesh( sim.strainXYZ(:,circleIndices) ,markerStyle) 
    
%  % plot strain at pointIndices
figure('position',[100,500,1000,600]);

subplot(231); hold on
    plot( sim.Om0.strain(:, sidePoints([2,4]),1) )
    plot( sim.Om10.strain(:, sidePoints([2,4]),1),'--' )
    xlabel('Time (ms)'); ylabel('Strain [\Delta L/L]')
    legend('Top F','Bottom F', 'Top FR','Bottom FR' )
subplot(232); hold on
    difference = sim.Om0.strain(:, sidePoints([2,4]),1) -...
                sim.Om10.strain(:, sidePoints([2,4]),1);
    plot(difference)
    xlabel('Time (ms)'); ylabel('Strain [\Delta L/L]')
    legend('Top F-FR','Bottom F-FR' )
subplot(233); hold on
    difference = sim.Om0.strain(:, sidePoints([2]),1) -...
                sim.Om0.strain(:, sidePoints([4]),1);
    plot(difference)
    difference = sim.Om10.strain(:, sidePoints([2]),1) -...
                sim.Om10.strain(:, sidePoints([4]),1);
    plot(difference,'--')
    xlabel('Time (ms)'); ylabel('Strain [\Delta L/L]')
    legend('Top-Bottom F','Top-Bottom FR' )
subplot(234); hold on
    plot( sim.Om0.strain(:, sidePoints([1,3]),1) )
    plot( sim.Om10.strain(:, sidePoints([1,3]),1),'--' )
    xlabel('Time (ms)'); ylabel('Strain [\Delta L/L]')
    legend('Left F','Right F', 'Left FR','Right FR' )
subplot(235); hold on
    difference = sim.Om0.strain(:, sidePoints([1,3]),1) -...
                sim.Om10.strain(:, sidePoints([1,3]),1);
    plot(difference)
    xlabel('Time (ms)'); ylabel('Strain [\Delta L/L]')
    legend('Left F-FR','Right F-FR' )
subplot(236); hold on
    difference = sim.Om0.strain(:, sidePoints([1]),1) -...
                sim.Om0.strain(:, sidePoints([1]),1);
    plot(difference)
    difference = sim.Om10.strain(:, sidePoints([1]),1) -...
                sim.Om10.strain(:, sidePoints([3]),1);
    plot(difference,'--')
    xlabel('Time (ms)'); ylabel('Strain [\Delta L/L]')
    legend('Left-Right F','Left-Right FR' )
    
% 
% %% Create neural encoding functions      
% STAfreq = 0.5;
% STAwidth = 5;
% STAdelay = 5;
% NLDgrad = 20;
% NLDshift = 0.8;
% [STAfun,NLDfun]=createNeuralFilters( STAfreq,STAwidth,STAdelay,NLDgrad,NLDshift );
% 
% %% Apply neural encoding 
% fSamp = 100;
% subSamp =10;
% t = (0:size(strainData(1).strain,1)-1)/fSamp;
% tNew = linspace(0,t(end), size(strainData(1).strain,1) *subSamp ); 
% STAt = linspace(-39,0,40*subSamp);
% 
% 
% 
% 
% parse_extra = 5;
% for j = 1:length(simName)
%     for k = 1:length(pointIndices(j).Ind)
%         % add tail of data for convolution 
%         strainTemp = strainPrinciple(j).strain(k,:);
%         strainTemp_preConv = [strainTemp, strainTemp(end-24:end), strainTemp(end-24:end-25+parse_extra)]; 
%         
%         % Subsampling, create time vector 
%         t = (0:size(strainData(1).strain,1) + 24+parse_extra)/fSamp;
%         tNew = linspace(0,t(end), (size(strainData(1).strain,1)+25+parse_extra) *subSamp ); 
%         
%         % Subsampaling, interpolate strain data 
%         data(j).strain(k,:) = interp1(t,strainTemp_preConv,tNew,'spline');
%         
%         % 1) Neural encoding, apply STA
%         data(j).strainConv(k,:) = conv( data(j).strain(k,:)  , STAfun(STAt),'valid') ;
%         
%         % 2) Neural encoding, apply NLD
%         calib = max(  data(j).strainConv(k,:)  );
%         data(j).pFire(k,:) = NLDfun( data(j).strainConv(k,:) / calib ); 
%         
%        % determine location of max probability of firing 
%        [V,I] = max(  data(j).pFire(k,end-24*subSamp:end)  );
%        data(j).spike(k) = I;
%     end
% end
% 
% %%
% % % 
% % % t_conv = tNew(length(STAt):end);
% % % t_spike = t_conv( end-24*subSamp:end);
% % % 
% % % for j = 4%1:4
% % % k = circleIndices(j).subInd(end-2:2:end);
% % % 
% % % figure()
% % %     subplot(311)
% %     
% %         plot( tNew , data(j).strain(k,:)' )
% % %             axis([tNew(length(STAt)), tNew(end), [-1,1]*2e-3 ] )
% %         ylabel('Strain')
% %     subplot(312)
% %         plot( t_conv, data(j).strainConv(k,:)' )
% % %             axis([tNew(length(STAt)), tNew(end), [-1,1]*5e-2 ] )
% %         ylabel('Projection')
% %     subplot(313)
% %         plot(t_conv,  data(j).pFire(k,:)' )
% % %             axis([tNew(length(STAt)), tNew(end), [0,1] ] )
% %         ylabel('Probability of Firing')
% %         xlabel('time [?s]');
% %     hold on
% %         plot( t_spike( data(j).spike(k) ),1,'bd')
% % end
% % 
% % %% difference in firing rate time 
% % for j = 1:length(simMat)
% %    data(j).dt =  abs( data(j).spike(2) -   data(j).spike(3))/subSamp;
% %    data(j).dt = data(j).dt; 
% % end
% % % 
% % symbolList = {'*','*','d','d'};
% % figure();hold on
% % for j = 1:length(simMat)
% %     plot(simMat{j,2},data(j).dt,symbolList{j})
% %     legendList{j} = [simMat{j,1}, ' ',num2str(simMat{j,2})];
% % end 
% % 
% % legend(legendList,'Location','NorthEastOutside')
% % xlabel('Rotation rate [rad/s]')
% % ylabel('Delta t [ms]')
% % axis([-1,11,-1,8])
% % 
% % %% Firing rate along circumference of haltere 
% % 
% % sim_run = [1,2];
% % sub_ind = reshape(1:32,[4,8])';
% % circle_ind = reshape( sub_ind(:,[1,3]),[16,1]);
% % plot_ind = reshape( sub_ind(:,[2,4]),[16,1]);
% % 
% % t_non = (1:length(data(sim_run(1)).pFire(k(1),end-24*subSamp:end)))/length(data(sim_run(1)).pFire(k(1),end-24*subSamp:end));
% % symb = {'d','+'};
% % col = {'r','r';'b','b'};
% % figure('Position',[100,100,800,1000])
% % for jj = 1:length(sim_run)
% %     k = 1:16;
% %     xyz = strainData(jj).xyzStrain(:,pointIndices(jj).Ind(:));
% %     for j = 1:length(k)
% %         subplot( 8,4,circle_ind(j) ); hold on
% %             plot( xyz(2,[1:end,1]),xyz(3,[1:end,1]),'b')
% %             scatter( xyz(2,k(j)),xyz(3,k(j)),'o','filled')
% %             axis square
% %             axis off
% % 
% % t_non = (1:length(data(sim_run(jj)).pFire(k(j),end-24*subSamp:end)))/length(data(sim_run(jj)).pFire(k(j),end-24*subSamp:end));
% % 
% %         subplot( 8,4,plot_ind(j) ); hold on
% %             plot( t_non,data(sim_run(jj)).pFire(k(j),end-24*subSamp:end))
% %             plot( data(sim_run(jj)).spike(k(j)) /241,1.2,symb{jj},'LineWidth',1.5)
% % %         axis off
% %             xlabel('time');ylabel('pFire')
% %     end
% % end
% % legend({'F','Fspike','F&R','F&R spike'},'Location','NorthEastOutside')
% % 
% % 
% % %% Difference in spike timing compared
% % figure();
% % subplot(311)
% %     plot(data(1).spike / subSamp)
% %     hold on 
% %     plot(data(2).spike / subSamp)
% %     ylabel('spike timing'); xlabel('location')
% %     ax = gca;
% %     ax.XTick = [1:4:13];
% %     ax.XTickLabel= {'left','bottom','right','top'};
% %     legend('Flapping','Flapping & Rot','Location','Best')
% %     axis([0,17,0,20])
% % subplot(312); hold on
% %     plot( abs( data(1).spike- data(2).spike)  / subSamp)
% %     plot([0,17],[1,1]*0.2,'--k','LineWidth',0.5)
% %     ylabel('spike timing difference'); xlabel('location')
% %     ax = gca;
% %     ax.XTick = [1:4:13];
% %     ax.XTickLabel= {'left','bottom','right','top'};
% %     axis([0,17,0,10])
% % subplot(313); hold on 
% %     plot( abs( data(1).spike- data(2).spike)  / subSamp)
% %     hold on 
% %     plot([0,17],[1,1]*0.2,'--k','LineWidth',0.5)
% %     ylabel('spike timing difference'); xlabel('location')
% %     axis([0,17,0,1])
% %     ax = gca;
% %     ax.XTick = [1:4:13];
% %     ax.XTickLabel= {'left','bottom','right','top'};
% % 
% % detectable =  find( (abs( data(1).spike- data(2).spike)  / subSamp) >= 0.2);
% % % detectable = detectable(1:5)
% % %% Plot locations with spike timing difference larger than jitter of 0.2 ms
% % figure(); 
% %     hold on
% %     plot( xyz(2,[1:end,1]),xyz(3,[1:end,1]),'b')
% %     scatter( xyz(2,detectable),xyz(3,detectable),'ro','filled')
% %     xlabel('cross section [y]'); ylabel('cross section [z]')
% %     legend('circumference','locations with large $\Delta t$','Location','Best')
% %     
    
    
    

%% Determine principal components 
% run('haltereSimulationAnalysis_principleComponents')
% 


%% Analyze deformations 
% run('haltereSimulationAnalysis_Deformation')
