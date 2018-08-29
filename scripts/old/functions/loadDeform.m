function [ xyz, D,  which_disp ] = loadDeform( nameParts )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


if nameParts{4}(1:2) == 'Om'
    which_disp = {' u',' v',' w'};
    file = ['data' filesep strjoin(nameParts(1:3),'_') filesep strjoin(nameParts,'_') '_deform.txt'];
    X = importdata( file , 't');
    testline = char(X(9));
    testline_split = strsplit(testline,'@');
    for j = 1:length( which_disp)
        matches = strfind(testline_split, which_disp{j} );
        line_inds(j,:) = find(~cellfun(@isempty,matches)) +3 ;
    end
elseif nameParts{4} == 'eigenfrequency'
    which_disp = {'disp', ' u',' v',' w'};
    file = ['data' filesep strjoin(nameParts(1:3),'_') filesep strjoin(nameParts,'_') '_deform.txt'];
    X = importdata( file , 't');
    testline = char(X(9));
    testline_split = strsplit(testline,'@');
    for j = 1:length( which_disp)
        matches = strfind(testline_split, which_disp{j} );
        line_inds(j,:) = find(~cellfun(@isempty,matches));
    end
end

header = X(1:9);
nodes = str2num(header{5}(end-7:end));
xyz = zeros(nodes,3);

for k = 1:nodes
    X_j  = strsplit( X{k+9} , ' ') ;
    xyz(k,:) = cellfun(@str2num, X_j(1:3) );
end
n_times = (length(X_j)-3)/4;% 
D = zeros(nodes,n_times, length(which_disp) );
for j = 1:length( which_disp)
    for k = 1:nodes
%         mask = ~cellfun(@isempty,matches);
%         mask = [0 , 0 , mask];
        X_j  = strsplit( X{k+9} , ' ') ;
        D(k,:,j) = cellfun(@str2num, X_j( line_inds(j,:) ) );
    end
end