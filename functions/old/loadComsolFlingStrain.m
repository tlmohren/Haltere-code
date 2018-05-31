function [ z,E,order ] = loadComsolFlingStrain( full_dir , which_strain )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    % import data into X (do I need if statement here? )
    if strcmp(full_dir(end-3:end),'.txt')
        X = importdata( full_dir, 't');
    else
        X = importdata( [full_dir,'.txt'] , 't');
    end

    % Check if file contains header, 
    first_line = char(X(1));
    header = all(first_line(3:7) == 'Model');
    
    % Divide string into columns 
    testline = char(X(9));
    testline_split = strsplit(testline,'  ');
    matches = strfind(testline_split, which_strain{1} );
    
    
    testline_splitNums = strsplit(testline,'@');
    matches2 = strfind( testline_splitNums, which_strain{1}  );
    inds =  find(~cellfun(@isempty,matches2)) + 2;
%         for j = 1:length(matches)
%            matches{j} 
%         end
%     for j = 1:length( matches{4} )
% %         matches = strfind(testline_split, which_strain{j})
%         line_inds(j,:) = find(~cellfun(@isempty,matches))
%     end
%     length( matches{4} )
    n_eigs = length( matches{4} );
%    
    if header
        % take out first 9 rows 
        n_points = length(X)-9;
        % initialize strain matrix 
        E =zeros(n_points,n_eigs);
        
        for j = 1:n_points
            
            % split character line 
            y11 = char( X(j+9));
            y12  = strsplit( y11,'  ');
% %             y12
%             y12{inds}
            ntot = length(y12);
            % xyz locations
            z(1,j) = str2double(y12{1});
            z(2,j) = str2double(y12{2}); 
            z(3,j) = str2double(y12{3});

            for m = 1:n_eigs
                E(j,m) = str2double( y12{   inds(m)   }); 
            end
        end
    else
        display('Table not right format?')
    end
%     z = round(z,14);
    order = which_strain;
end