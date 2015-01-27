%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% train_hiHmm.m 
%% by Kyung-Ah Sohn
%% January 2014 
%% Input 
%% - list_fname = list of text filenames ; 
%% - list_spcs = list of species id (should start from 1) ; 
%% - nproc = 4; % number of parallel processors to run 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ mapS, stats ] = train_hiHmm( list_fnames, finfo, param, odir )

	stream = RandStream('mt19937ar','seed',201);
	RandStream.setGlobalStream(stream);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% input check 
	if ( isfield( param, 'b_condition_specific_emission' ) == 0 )
	    param.b_condition_specific_emission = 1;
	end
	
	if ( isfield( param, 'parallel' ) == 0 )
		nproc = 1;
	else
		nproc = param.nproc;
	end
	
	if exist('odir','var') == 0;
	    odir = './training/';
	end
	if ( exist(odir)== 0 )
	    mkdir( odir )
	end
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% 1. load train data 
	[ train_data, locs, markers ] = read_data( list_fnames, finfo );
	C = length(train_data);
	
	% 2. normalize data 
	[ train_data_norm, mean_std_condition, tdata_tmp ] = normalizedata( train_data );
	save( [ odir '/train_mean_std.mat'], 'mean_std_condition' );
	
	T_species = zeros(C,1);
	for cc=1:C
	   T_species(cc) = size( tdata_tmp{cc}, 1 );
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% 2. set hyper-parameters 
	T = sum( T_species );
	hypers.alpha0_a = 4*0.01*T/C;  % 
	hypers.alpha0_b = 2;
	hypers.gamma_a = 3; % 3
	hypers.gamma_b = 6;
	hypers.mu_0 = 0;
	hypers.sigma2 = param.sigma2;
	hypers.sigma2_0 = param.sigma2_0;
	
	hypers.p0 = 0.5*ones(1, C);
	hypers.nproc =  nproc;
	
	tic
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% 3. get initial state sequence by K-means 
	K0 = param.K0;  % initial number of states 
	opts = statset( 'Display', 'iter' , 'MaxIter', 500  );
	[ cidx, ctrs ] = kmeans( cell2mat(tdata_tmp), K0, 'Options', opts , 'emptyaction', 'singleton' );
	S0 = cell(C,1);
	S00 = mat2cell( cidx, T_species );
	for cc=1:C
	      Tc = zeros( length(train_data_norm{cc}),1 );
	      for j = 1:length( train_data_norm{cc} )
		Tc(j) = size( train_data_norm{cc}{j}, 2 );
	      end
	      S0{cc} = mat2cell( S00{cc}, Tc );
	end
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% 4. training by hiHMM (run parallel if computing nodes available) 
	if ( param.b_condition_specific_emission ) 
		b_markzerocov = 0;  % mark zero-coverage state in the emission matrix 
		[S, stats, mapS] = hiHmmp_model1( train_data_norm, hypers, param.numb, param.nums, param.numi, S0, odir, ...
				markers, b_markzerocov );
	else 
		[S, stats, mapS] = hiHmmp_model2( train_data_norm, hypers, param.numb, param.nums, param.numi, S0, odir, markers );
	
	end
	
	toc
	
	%exit

end 

