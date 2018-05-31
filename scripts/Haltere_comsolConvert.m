% -----------------------------
% TMOHREN
% import COMSOL strain data
% -----------------------------
clc; clear all; close all
addpathFolderStructureHaltere()
% define location of data 
data_base = 'data';

% names of FEAulations to load 
FEA(1).name = 'Haltere_CraneFly_Sphere_Om0';
% FEA(2).name = 'Haltere_CraneFly_Sphere_Om10';
% FEA(3).name = 'Haltere_CraneFly_ellipsoidHor_Om0';
% FEA(4).name = 'Haltere_CraneFly_ellipsoidHor_Om10';
% FEA(5).name = 'Haltere_CraneFly_ellipsoidVer_Om0';
% FEA(6).name = 'Haltere_CraneFly_ellipsoidVer_Om10';
% FEA(7).name = 'Haltere_CraneFly_Sphere_eigenfrequency';
% FEA(8).name = 'Haltere_CraneFly_ellipsoidHor_eigenfrequency';
% FEA(9).name = 'Haltere_CraneFly_ellipsoidVer_eigenfrequency';

% outputName = ['data' filesep 'Cranefly_Sphere_FEAdeformtesting'];
% load( ['data' filesep 'Cranefly_SphereVsEllipsoid_FEAresults'],'FEA')

%% 

for j =  1:length(FEA)
    tic
    
    data = loadCSV( ['data' filesep  FEA(j).name], { 'eYY' });
%     nameParts = strsplit(FEA(j).name,'_'); 
%     FEA(j).shape = nameParts{3};
%     if  nameParts{4}(1:2) == 'Om'
%         FEA(j).omega = str2num( nameParts{4}(3:end) );
%         FEA(j).analysis = 'time dependent';
%         [ FEA(j).XYZ, FEA(j).strain, FEA(j).strainOrder ] = loadStrain( nameParts );
%         [ ~,FEA(j).deform, FEA(j).deformOrder ] = loadDeform( nameParts);
%         [ FEA(j).phi,FEA(j).theta ]  = loadAngles( nameParts ); 
%  
%     elseif nameParts{4} == 'eigenfrequency'
%         FEA(j).analysis = 'eigenfrequency';
%         [ ~,FEA(j).strain, FEA(j).strainOrder ] = loadStrain( nameParts);
%         [ FEA(j).XYZ,FEA(j).deform, FEA(j).deformOrder ] = loadDeform( nameParts);
%     end
    toc 
end

%% 
% save(outputName,'FEA')

%% test 
% % j = 3;
% % figure();
% % scatter3( FEA(j).XYZ(:,1), FEA(j).XYZ(:,2), FEA(j).XYZ(:,3)  )
% % axis equal
% % 
% % figure();
% %     subplot(311); plot( FEA(1).deform(1,:,1))
% %     subplot(312); plot( FEA(1).deform(1,:,2))
% %     subplot(313); plot( FEA(1).deform(1,:,3))
% %     %
% % figure();
% %     subplot(211); plot( FEA(1).phi) 
% %     subplot(212); plot( FEA(1).theta)