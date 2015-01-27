function [ train_data_norm, mean_std_condition, tdata_tmp ] = normalizedata( train_data );
% normalize data with genome-wide mean and standard deviation 

C = length(train_data);

train_data_norm = cell( C, 1 );
mean_std_condition = cell(C, 2);
T_species = zeros(C, 1);
tdata_tmp = cell(C, 1);
for cc=1:C
        M = cell2mat( train_data{cc} );
        T_species(cc) = size(M, 1);
        mean_global = mean( M, 1 );
        std_global = std( M, 1 );
        mean_std_condition{cc,1} = mean_global;
        mean_std_condition{cc,2} = std_global;

        M_norm = ( M-repmat( mean_global, size(M, 1),1 ) )./ repmat( std_global, size(M,1),1 );
        tdata_tmp{cc} = M_norm;
        R = zeros( length(train_data{cc} ), 1);
        for jj=1:length(train_data{cc} )
                R(jj) = size(train_data{cc}{jj},1 );
        end
        M_norm_c = mat2cell( M_norm, R );
        for jj=1:length(train_data{cc} )
                train_data_norm{cc}{jj} = M_norm_c{jj}';
        end
end

