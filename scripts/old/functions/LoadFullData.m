function [ z,E,order ] = LoadFullData( full_dir )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%     n_strains = 6;

    if strcmp(full_dir(end-3:end),'.txt')
        X = importdata( full_dir, 't');
    else

        X = importdata( [full_dir,'.txt'] , 't');
    end

    first_line = char(X(1));
    header = all(first_line(3:7) == 'Model');
    
  
   
    testline = char(X(7));
    strains = {'XX','YY','ZZ','XY','YZ','XZ',...
        'First principal strain','direction 1, X','direction 1, Y','direction 1, Z',...
        'Second principal strain','direction 2, X','direction 2, Y','direction 2, Z',...
        'Third principal strain','direction 3, X','direction 3, Y','direction 3, Z'};
    
    for j = 1:length(strains)
        matches = strfind(testline, strains{j});
        if matches >=0
            binar(j) = matches;
        else 
            binar(j) = 0;
        end
    end
    binar;
    n_strains = sum(binar>0);
    c = find(binar);
    [~,b]=sort(binar(c));
    
    for kk = 1:length(b)
        order{kk} = strains{c(b(kk))};
    end

    if header
        n_points = length(X)-9;
        t_steps =  (  length(strsplit(char(X(10))))-3 ) / sum(binar>0);
        
%         t_steps
%         n_points
%         n_strains
        
        E =zeros(t_steps,n_points,n_strains);
        
        for j = 1:n_points
            y11 = char( X(j+9));
            y12 = strsplit(y11);
            z(1,j) = str2double(y12{1});
            z(2,j) = str2double(y12{2}); 
            z(3,j) = str2double(y12{3});

            for k = 1:t_steps
                for m = 1:n_strains
                    E(k,j,m) = str2double(y12{   3+  k*n_strains + m -n_strains      }); 
                end
            end
        end
    else
        display('Table not right format?')
    end

    z = round(z,6);
    
end

