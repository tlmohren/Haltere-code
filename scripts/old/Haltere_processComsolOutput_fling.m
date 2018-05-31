% -----------------------------
% TMOHREN
% import COMSOL strain data
% -----------------------------
clc; clear all; close all
addpathFolderStructureHaltere()

% define location of data 
data_base = ['data' filesep 'Haltere_CraneFly_Sphere_longRunK0_00001'];

% names of simulations to load 
% sim(1).name = 'Haltere_CraneFly_Sphere_longRun_Om0';
% sim(2).name = 'Haltere_CraneFly_Sphere_longRun_Om10';

% sim(1).name = 'Haltere_CraneFly_Sphere_longRunK0_00001_Om0';
% sim(2).name = 'Haltere_CraneFly_Sphere_longRunK0_00001_Om10';

sim(1).name = 'Haltere_CraneFly_Sphere_fling';

outputName = ['data' filesep 'Cranefly_Sphere_fling_strainDeform'];

%% 

for j = 1:length( sim  )
    tic
    % find subfolder and final data directory
    nameParts = strsplit(sim(j).name,'_Om'); 
    dataDirect =  [data_base  filesep sim(j).name ];
%     sim(j).thetaDot = nameParts{2};
    
    [sim(j).XYZ ,  sim(j).strain    , sim(j).strainOrder  ] = loadComsolFlingStrain( [dataDirect  '_strain'] , {'YY'} );
    [~          ,  sim(j).deform    , sim(j).deformOrder  ] = loadComsolFlingDeform( [dataDirect  '_deform'] , {'disp',' u',' v',' w'} );
%     [sim(j).phi ,  sim(j).theta  ]                          = loadExcitationAngles( [dataDirect '_excitationAngles'] ); 
    toc 
end

%% 
save(outputName,'sim')