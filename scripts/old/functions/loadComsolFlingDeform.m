function [ z,D, order ] = LoadComsolFlingDeform( full_dir, which_disp )
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
%     testline = char(X(9));
%     testline_split = strsplit(testline,'t');
%     for j = 1:length( which_disp)
%         matches = strfind(testline_split, which_disp{j});
%         line_inds(j,:) = find(~cellfun(@isempty,matches));
%     end


    testline = char(X(9));
    testline_split = strsplit(testline,'  ');
%     for j = 1:length(which_disp) 
    matches = strfind(testline_split, which_disp{1} );
%     end
    
    testline_splitNums = strsplit(testline,'@');
    for j = 1:length(which_disp)
        matches2 = strfind( testline_splitNums, which_disp{j}  )
        inds(j,:) =  find(~cellfun(@isempty,matches2)) + 3;
        
    end
    
    
    [n_disps, n_eigs] = size(inds);
%     n_eigs = length( matches{4} );
%     [n_disps, n_times ] = size(line_inds);
%     
    if header
        % take out first 9 rows 
        n_points = length(X)-9;
%         D =zeros( n_times, n_points, n_disps);
        D =zeros(n_points,n_eigs,n_disps);
%         
        for j = 1:n_points
            % split character line 
            y11 = char( X(j+9));
            y12  = strsplit( y11,'  ');
            
            % xyz locations
            z(1,j) = str2double(y12{1});
            z(2,j) = str2double(y12{2}); 
            z(3,j) = str2double(y12{3});
% 
            % strain allocations
            ntot = length(y12);
            % xyz locations
            z(1,j) = str2double(y12{1});
            z(2,j) = str2double(y12{2}); 
            z(3,j) = str2double(y12{3});

            for k = 1:n_disps 
                for m = 1:n_eigs
%                     inds(k,m)  
                    D(j,m,k) = str2double( y12{   inds(k,m)   }); 
                end
            end
        end
    else
        display('Table not right format?')
    end
%     z = round(z,14);
    order = which_disp;
end

