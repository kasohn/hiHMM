%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a sample code to run hiHmm on a given toy dataset  
% 

% for parallel run setting 
nproc = 1;    % number of parallel processing nodes 
%matlabpool(nproc);  %% reserve 'nproc' workers. Remove % to do parallel processing 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Part 1. train the model 
%%% INPUT:  one text file per condition per chromosome 
%  -------------------- EXAMPLE FILE ------------------ 
%	H3K4me3 H3K4me1 H3K27ac H3K27me3  
% 100     1.273   1.845   1.255   6.678   
% 300     1.263   2.094   1.314   6.704   
% 500     -0.1197 -0.3921 -0.3642 1.862  
% ....
%  --------------------------- ------------------------- 
%  list_fnames_files: input file names 
%     { 'fname_condition1_chr1', ..., 'fname_condition1_chrX1', ... 
%  	'fname_condition1_chr2', ..., 'fname_condition2_chrX2', ... 
%	... 
%	'fname_conditionM_chr1', ..., 'fname_conditionM_chrXM' 
%	}
% 

%%% set the directory where the input and output file folders are located for the current analysis
analysis_folder = './sample_analysis/';
indir = 'Mapped_input_files/'; % input directory

%%% data parameters  
condition = { 'worm', 'fly' };   % condition labels, e.g. for species name, timepoint, celllines, etc
list_tnames_files = { 'worm_chrI.txt', 'worm_chrII.txt', 'fly_chr2L.txt' }; % input files
list_tnames = strcat(analysis_folder, indir, list_tnames_files);
tinfo.condition = [ 1 1 2 ];	 % integer labels for condition (e.g. species/stages/time etc. per file) 
tinfo.chrs = { 'chrI' 'chrII' 'chr2L' };   % chromosome numbers (for bedfiles)
bin_size = 200;		% data bin size: used when writing a bedfile after Viterbi decoding  

odir_name = 'output/';	% output directory name
odir = [analysis_folder odir_name]; 		% output directory path
if ( exist(odir) == 0 )
	mkdir(odir);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% % program parameters 
% for model selection between 1 and 2  
param.b_condition_specific_emission = 0;	%  1 for model 1, 0 for model 2
%%%%%%%%%%%%%%%%%%%%%%%
param.K0 = 7; 	% initial number of states 
% for emission model 
param.sigma2 = 1;   % variance in the emission model 
param.sigma2_0 = 4;  % variance in Gaussian for prior distribution of each cell in the emission matrix. 
			  % This parameter should be large enough as the prior mean of each cell in the emission matrix is zero
% for MCMC iteration 
param.numb = 200; % the max number of burnin iterations
param.nums = 10; % the number of samples to output
param.numi = 1; % the number of sampling, iterations between two samples
% parallel run setting 
param.parallel = 0; 	% currently not used 
param.nproc = nproc;

%% train the model and save the results under odir 
train_hiHmm( list_tnames, tinfo, param, odir  );
% OUTPUT: odir/
% - train-hihmm-model#.emission.csv
% - train-hihmm-model#.transition.csv
% - train-hihmm-model#.mat   % contains MAP solution from the training. If the model is trained with the whole genome, use Maximum a Posterior solution.  
% - train_mean_std.mat % mean and std of training data, per condition 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Part 2. Viterbi decoding (when training data and testing data are different)

% load genome-wide mean and sd of the training data, per condition, for normalization of test data
load( [ odir '/train_mean_std.mat'], 'mean_std_condition' );

model = 2 -  param.b_condition_specific_emission;   % model 1 or 2 

% load hmm parameters 
mu0 = csvread( [odir 'train-hihmm-model' num2str(model) '.emission.csv'], 1, 0 );
pi0 = csvread( [odir 'train-hihmm-model' num2str(model) '.transition.csv'], 1, 0 );
M = size(mu0, 2);
K = size(pi0, 2);
C = size(pi0, 1) / K;
if ( model == 1 )
	muc = permute( reshape( mu0, [ K, C, M ] ), [1 3 2] );
else
	muc = repmat(  mu0, [1, 1, C] );
end
pic = permute( reshape( pi0, [ K, C, K ] ), [1 3 2] );

% load color map for bed files 
if ( exist('statecolmap.txt') )
	col = load( sprintf('statecolmap.txt'));
	assert( size(col,1) >= K );
else
	xx = 4*((K-1):-1:0)/(K-1);
	col = [ min(xx-1.5, -xx+4.5); min(xx-0.5, -xx+3.5); min(xx+0.5, -xx+2.5) ];
	col( find(col(:)>1) ) = 1;
	col( find(col(:)<0) ) = 0;
	col = floor( col'*255 );
end

% reord states if necessary 
ord = [1:K];

% list of filenames to decode 
list_fnames_files = { 'worm_chrI.txt', 'worm_chrII.txt', 'fly_chr2L.txt' };
list_fnames = strcat(analysis_folder, indir, list_fnames_files);
finfo.condition = [ 1 1 2 ];
finfo.chrs = { 'chrI' 'chrII' 'chr2L' };

% load data 
[ decode_data, locs, markers, chrs ] = read_data( list_fnames, finfo );

tic;
fprintf( 'Viterbi decoding..\n');
% run viterbi on each file  of data  
for c=1:length(decode_data)  % for different conditions  
     S = cell( length(decode_data{c}), 1 );
     parfor j=1:length(decode_data{c})    % for each chromosome 
	T = size( decode_data{c}{j}, 1 );   
	% normalization by genome-wide mean and variance as in training 
	data_norm = ( decode_data{c}{j} - repmat( mean_std_condition{c,1}, T, 1 )  ) ./ ... 
			 ( repmat( mean_std_condition{c,2}, T, 1) )   ;

	S{j} = Viterbi( data_norm', pic(ord,ord,c), muc(ord,:,c)' );
     end	
     write_bedfile( S, locs{c}, chrs{c}, col, model, K, odir, condition{c}, bin_size );
end
toc;

% make sure to close matlabpools if any error occurs  (release workers)
%matlabpool close;   


