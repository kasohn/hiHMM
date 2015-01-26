function [S, stats, mapS] = hiHmmp_model1( Y, hypers, numb, nums, numi, S0, odir, markers, b_postproc  )

% HIHMMNORMALSAMPLEBEAM Samples states from the hierarchically linked iHMM with normal output 
% using the Beam sampler.
%
% Note: modified from the code for Infinite HMM by Jurgen Van Gael (http://mloss.org/software/view/205/)
%
% [S, stats] = hiHmmNormalSampleBeam(Y, hypers, numb, nums, numi, S0) uses 
% the beam sampling training algorithm for the infinite HMM.
%
%   Input Parameters:
%   - cY: multi-species training sequence of arbitrary length (cell type),
%   - hypers: a structure that describes the hyperparameters for the beam
%             sampler. If this structure contains alpha0 and gamma, it will
%             not resample these during the algorithm. If these are not
%             specified, one needs to specify hyperparameters for alpha0
%             and gamma (alpha0_a, alpha0_b, gamma_a, gamma_b). hypers
%             should also contain the normal parameters for the mean prior
%             mu_0, sigma2_0 and the output variance sigma2.
%   - numb: the number of burnin iterations,
%   - nums: the number of samples to output,
%   - numi: the number of sampling, iterations between two samples,
%   - cS0: is the initial assignment to the sequence.
%
%   Output Parameters:
%   - S: is a cell array of sample structures where each sample contains the
%        hidden state sequence S, the number of states K, the Beta, Pi,
%        Phi's used for that sample.
%   - stats: is a structure that contains a variety of statistics for every
%            iteration of the sampler: K, alpha0, gamma, the size of the
%            trellis and the marginal likelihood.

	head = sprintf( 'train-hihmm-model1' );
	if ( exist( 'b_postproc', 'var' ) == 0  ) 
		b_postproc = 1;
	end
	
	% Initialize the sampler.
	C = length( Y  ); 	% # of species 
	K = 1;
	for cc=1:C
	    Ts(cc) = 0;
	    for j=1:length(Y{cc})
	        Ts(cc) = Ts(cc) + size(Y{cc}{j},2 );
		K = max( K, max( S0{cc}{j} ) );
	    end
	end
	T = sum(Ts);                      % # of time-steps T
	M = size(Y{1}{1},1);
	K0 = K;
	
	
	maxjll = -inf;
	
	sample.S = S0;
	sample.K = K;
	
	% Setup structures to store the output.
	S = {};
	stats.K = zeros(1,(numb + (nums-1)*numi));
	stats.alpha0 = zeros(1,(numb + (nums-1)*numi));
	stats.gamma = zeros(1,(numb + (nums-1)*numi));
	stats.jml = zeros(1,(numb + (nums-1)*numi));
	stats.trellis = zeros(1,(numb + (nums-1)*numi));
	
	Ln = -1;
	np = 5; % convergence test using the last np samples
	% Initialize hypers; resample a few times as our inital guess might be off.
	if isfield(hypers, 'alpha0')
	    sample.alpha0 = hypers.alpha0;
	else
	    sample.alpha0 = gamrnd(hypers.alpha0_a, 1.0 / hypers.alpha0_b);
	end
	if isfield(hypers, 'gamma')
	    sample.gamma = hypers.gamma;
	else
	    sample.gamma = gamrnd(hypers.gamma_a, 1.0 / hypers.gamma_b);
	end
	for i=1:5
	    sample.Beta = ones(1, sample.K+1) / (sample.K+1);
	    [sample.Beta, sample.alpha0, sample.gamma] = hiHmmHyperSample1(sample.S, sample.Beta, sample.alpha0, sample.gamma, hypers, 20 );
	end
	
	% Sample the emission and transition probabilities.
	if isfield(hypers, 'Mus')
		sample.Mus = hypers.Mus;
	else
		sample.Mus = hiHmmSampleMultNormalMeans1( sample.S, Y, sample.K, hypers.sigma2, hypers.mu_0, hypers.sigma2_0 );
	end
	if isfield(hypers, 'Pi')
		sample.Pi = hypers.Pi;
	else
		sample.Pi = hiHmmSampleTransitionMatrix( sample.S, sample.alpha0*sample.Beta, 1 );
		sample.Pi(sample.K+1,:,:) = [];
	end
	
	iter = 1;
	
	alpha_g = 0.2 * (hypers.p0 .* Ts); 
	beta_g = 0.2 * Ts - alpha_g;
	
	if isfield(hypers, 'G0')
		sample.G = hypers.G0;
		sample.p0 = hypers.p0;
	else
		[ sample.G, sample.p0 ] = SampleAuxVariable( sample.S, hypers.p0, sample.Pi );
	end
	
	minseg = hypers.nproc*ones(1,C);
	while iter <= (numb + (nums-1)*numi)
	    
	    % Reset the trellis size count in case the previous iteration didn't
	    % return a samplable path.
	  stats.trellis(iter) = 0;
	    
	    % Sample the auxilary variables and extend Pi and Phi if necessary.
	  rid = randperm(C);
	  for c1 = rid
	    jid = randperm(length(Y{c1}));
	    nseg = min( length(jid), minseg(c1) );
	    [sample, uc] = BreakSticks( sample.S{c1}(jid(1:nseg) ), sample, hypers, c1 );     
	
	    Sc = sample.S{c1}(jid(1:nseg));
	    Gc = sample.G{c1}(jid(1:nseg));
	   
	   parfor j = 1:nseg   % jid
		[ Snew, Gnew ] = Sample_hiddenState( Y{c1}{jid(j)}, sample.S{c1}{jid(j)}, sample.G{c1}{jid(j)}, uc{j},  sample, hypers, c1, iter );
		Sc{j} = Snew;
		Gc{j} = Gnew;
	   end  %%% end of j 
	   sample.S{c1}(jid(1:nseg)) = Sc;
	   sample.G{c1}(jid(1:nseg)) = Gc;
	 
	        % Resample Beta given the transition probabilities.
	        [sample.Beta, sample.alpha0, sample.gamma] = hiHmmHyperSample1(sample.S, sample.Beta, sample.alpha0, sample.gamma, hypers, 20, sample.G);
	        
	        % Resample the Phi's given the new state sequences.
	        sample.Mus = hiHmmSampleMultNormalMeans1( sample.S, Y, sample.K, hypers.sigma2, hypers.mu_0, hypers.sigma2_0 );
	       
	        % Resample the transition probabilities.
	        sample.Pi = hiHmmSampleTransitionMatrix(sample.S, sample.alpha0*sample.Beta, 0,  sample.G );
	       	sample.Pi(sample.K+1,:,:) = [];
	
	        % Resample the self-transition probabilities and auxiliary variables.
	       for c=1:C
	        ng = sum( sample.G{c}{1} );
	        for j=2:length(sample.G{c}) 
	            ng = ng + sum( sample.G{c}{j} );
	       end
	        sample.p0(c) = betarnd( ng + alpha_g(c),  Ts(c) - ng + beta_g(c)  ) ;
	       end 
	  
	      % Cleanup our state space by removing redundant states.
		sset = [];
	        for c=1:C
		  for k=1:length(sample.S{c})
	                sset = union( sset, unique( sample.S{c}{k} ));
	          end
		end
	        zind = sort(setdiff(1:sample.K, sset ));
	        %zind = find(sum(N,1) == 0);
	        %zind = setdiff(zind, 1);
		if ( length( zind ) > 0 ) 
	          for i = length(zind):-1:1
	            sample.Beta(end) = sample.Beta(end) + sample.Beta(zind(i));
	            sample.Beta(zind(i)) = [];
	            sample.Pi(:,zind(i),:) = [];
	            sample.Pi(zind(i),:,:) = [];
	            sample.Mus(:,zind(i),:) = [];
		    for cc=1:C 
		      for k=1:length(sample.S{cc})
			idx = find( sample.S{cc}{k}  > zind(i) ) ;
	                sample.S{cc}{k}(idx) = sample.S{cc}{k}(idx) - 1;
		      end
		    end
	          end
		end
	        sample.K = size(sample.Pi,1);
	        % Safety checks
	        assert(size(sample.Pi,1) == sample.K);
	        assert(size(sample.Pi,2) == sample.K+1);
	        assert(sample.K == length(sample.Beta) - 1);
	        %assert(min(sample.Beta) > 0);
	        assert(min(sample.Pi(:)) >= 0);
	%        assert(sample.K == max( cell2mat(sample.S)));
	   end   %%% end of c1 
	      
	        % Prepare next iteration.
	        stats.alpha0(iter) = sample.alpha0;
	        stats.gamma(iter) = sample.gamma;
	        stats.K(iter) = sample.K;
	        stats.jll(iter) = hiHmmMultNormalJointLogLikelihood1(sample.S, Y, sample.Beta, ...
	           sample.alpha0, hypers.mu_0, hypers.sigma2_0, hypers.sigma2);
	
	        if ( iter > np+1  ) 
	            Ln2 = mean(abs(stats.jll(iter-np:iter)));   
	            Ln1= mean(abs(stats.jll(iter-np-1:iter-1)));   
	            Ln = abs( (Ln2-Ln1)/Ln1 );
	            if ( Ln < 1.e-5 & iter < numb  )
	                numb = iter;
	            end
	        end
		
	        if ( iter >= numb & stats.jll(iter) > maxjll )
	                mapS = sample;
	        end
	
	        fprintf('Iteration: %d: JL=%.4f, L=%.6f\n', iter, stats.jll(iter), Ln );
	        if iter >= numb && mod(iter-numb, numi) == 0
	            S{end+1} = sample;
	        end
	        
	        iter = iter + 1;
	end
	
	
	stats.K = stats.K(1:numb+(nums-1)*numi);
	stats.alpha0 = stats.alpha0(1:numb+(nums-1)*numi);
	stats.gamma = stats.gamma(1:numb+(nums-1)*numi);
	stats.jml = stats.jml(1:numb+(nums-1)*numi);
	stats.trellis = stats.trellis(1:numb+(nums-1)*numi);
	stats.jll = stats.jll(1:numb+(nums-1)*numi);
	
	%%%% check coverage of each state and mark zero coverage state  
	if ( b_postproc) 
	   cnt = zeros( 3, mapS.K );
	   for cc=1:length(mapS.S)
		[ nn ] = hist( cell2mat( mapS.S{cc} ), [1:mapS.K] );
		cnt(cc,:) = cnt(cc,:) + nn;
	   end
	
	   %% mark low coverage state per species 
	   for kk=1:size(cnt,2)
		for cc=1:C 
			if ( cnt(cc,kk) == 0 )
				mapS.Mus(:,kk,cc) = 9999;	
			end	
		end
	   end 
	end
	
	save( [ odir head '.mat'] , 'S', 'stats', 'mapS' );
	writehmm( mapS.Pi, mapS.Mus, mapS.p0, head, odir, markers );

end

%%%%%%%%%%%%%%%%%%%%%%%
function [ sample, uc ] = BreakSticks( Sc, sample, hypers, c1 )
    C = size( sample.Mus, 3) ;
    M = size( sample.Mus, 1 );
    n = length(Sc);
    uc = cell(n,1);
    minu = 1;
    for j=1:n 
	      T1 = length( Sc{j} );
	      uc{j} = zeros(1,T1);
	      rn = rand(T1,1);
	      uc{j}(1) = rn(1) * sample.Pi(1, Sc{j}(1), c1);
	      for t=2:T1
	        pikk = sample.Pi(Sc{j}(t-1), Sc{j}(t), c1);
	        uc{j}(t) = rn(t) * ( sample.p0(c1)*pikk + ( 1-sample.p0(c1) )* ( Sc{j}(t-1)==Sc{j}(t) ));
	      end
	      minu = min( minu, min(uc{j}) );
    end
    nit = 0;
    while ( max( sample.Pi(:, end, c1) ) > minu & nit < 5 )     % Break the Pi{k} stick some more.
%    while max( max(sample.Pi(:, end, :)) ) > min( cell2mat(u))     % Break the Pi{k} stick some more.
	nit = nit + 1;

        pl = size(sample.Pi, 2);
        bl = length(sample.Beta);

        % Safety check.
        assert(bl == pl);

        % Add row to transition matrix.
        muv = randn(M,1) * sqrt(hypers.sigma2_0) + hypers.mu_0;
        for cc=1:C 
                sample.Pi(bl,:,cc) = sample.p0(cc)*dirichlet_sample(sample.alpha0 * sample.Beta);
                sample.Pi(bl,bl,cc) = 1-sample.p0(cc);
                sample.Mus(:,bl,cc) = muv; %randn(M,1) * sqrt(hypers.sigma2_0) + hypers.mu_0;
        end
        % Break beta stick.
        be = sample.Beta(end);
        bg = betarnd(1, sample.gamma);
        sample.Beta(bl) = bg * be;
        sample.Beta(bl+1) = (1-bg) * be;

        for cc=1:C 
                pe = sample.Pi(:, end, cc);
                a = repmat(sample.alpha0 * sample.Beta(end-1), bl, 1);
                b = sample.alpha0 * (1 - sum(sample.Beta(1:end-1)));
                pg = betarnd( a, b );
                if isnan(sum(pg))                   % This is an approximation when a or b are really small.
                        pg = binornd(1, a./(a+b));
                end
                sample.Pi(:, pl, cc) = pg .* pe;
                sample.Pi(:, pl+1, cc) = (1-pg) .* pe;
        end
      end	% end of while 
      sample.K = size(sample.Pi, 1);
% Safety check.
    assert(sample.K == length(sample.Beta) - 1);
    assert(sample.K == size(sample.Mus, 2));

end
 
function [ Snew, Gnew ] = Sample_hiddenState( Yj, Sj, Gj, uc, sample, hypers, c1, iter ) 
	%      Sj = sample.S{c1}{j};
	%      Yj = Y{c1}{j};
	      C = size( sample.Pi, 3 );
	      M = size( Yj, 1 );
	      T1 = length(Sj);
	    % Resample the hidden state sequence.
	%  for cc=1:C
	    dyn_prog = zeros(sample.K, T1);
	   
	    dyn_prog(:,1) = sample.Pi(1,1:sample.K, c1) > uc(1);
	%    stats.trellis(iter) = stats.trellis(iter) + sum(sum(dyn_prog(:,1)));
	    for k=1:sample.K
	        dyn_prog(k,1) = exp(-0.5*(Yj(:,1) - sample.Mus(:,k,c1))'*(Yj(:,1) - sample.Mus(:,k,c1))/hypers.sigma2) * dyn_prog(k,1);
	    end
	    dyn_prog(:,1) = dyn_prog(:,1) / sum(dyn_prog(:,1));
	    
	    A0 = sample.p0(c1) * sample.Pi(1:sample.K, 1:sample.K,c1) + ( 1-sample.p0(c1) ) * eye(sample.K);
	    for t=2:T1
	        A =  A0 > uc(t);
		%%%% check 
	        dyn_prog(:,t) = A' * dyn_prog(:,t-1);
	%        stats.trellis(iter) = stats.trellis(iter) + sum(sum(A));
	        for k=1:sample.K
	            dyn_prog(k,t) = exp(-0.5*(Yj(:,t) - sample.Mus(:,k,c1))'*(Yj(:,t) - sample.Mus(:,k,c1))/hypers.sigma2) * dyn_prog(k,t);
	        end
	        dyn_prog(:,t) = dyn_prog(:,t) / sum(dyn_prog(:,t));
	    end
	  
	    Snew = Sj; 
	    Gnew = Gj;
	    % Backtrack to sample a path through the HMM.
	    if sum(dyn_prog(:,T1)) ~= 0.0 && isfinite(sum(dyn_prog(:,T1)))
	        rn = rand(T1, 1);
	        rn1 = rand(T1,1);
	        s1 = 1 + sum(rn(T1) > cumsum(dyn_prog(:,T1)));
	        Snew(T1) = s1;
	        for t=T1-1:-1:1
	            pkk = sample.p0(c1) * sample.Pi(:, Snew( t+1 ), c1) ;
	            pkk(s1) = pkk(s1) + ( 1-sample.p0(c1) );
	            r = dyn_prog(:,t) .* ( pkk  > uc( t+1 ) );
	            r = r ./ sum(r);
	            s2 = s1;
	            s1 = 1 + sum(rn(t) > cumsum(r));
	            Snew( t ) = s1;
	            if ( s1 ~= s2 )
	                Gnew(t+1) = 1;
	            else
	                prob1 = sample.p0(c1) * sample.Pi( s1, s2 );
	                prob2 = ( 1 - sample.p0(c1) );
	                if ( rn1(t) < prob1/(prob1+prob2) )
	                    Gnew(t+1) = 1;
	                else
	                    Gnew(t+1) = 0;
	                end
	            end
	        end
	        
	        % Safety check.
	        assert(~isnan(sum(Snew(t))));
	    else
	        fprintf('Wasted computation as there were no paths through the iHMM.\n');
	    end
	
	%    Snew = sample.S;
	%    Gnew = sample.G;

end

