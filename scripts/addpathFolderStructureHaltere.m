function addpathFolderStructureHaltere()

% Find location of this script 
    scriptLocation = fileparts(fileparts(mfilename('fullpath') ));
    [folderLocation,baseName] = fileparts(scriptLocation);
    dataFolderLocation = [folderLocation, filesep, baseName,'Data'];
    
    cd(scriptLocation );

    addpath([scriptLocation filesep 'data'])
    addpath([scriptLocation filesep 'functions'])
    addpath([scriptLocation filesep 'scripts'])

