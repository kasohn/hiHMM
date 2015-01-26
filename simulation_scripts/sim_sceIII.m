%% scenario III - zero to two species-specific site 

stream = RandStream('mt19937ar','seed',20); %for rr=11:100
RandStream.setGlobalStream(stream);

sce0 = 'sceIII';

C = 3;
K = 10;
M = 8;	% num marker 
vT = [2000 2000 2000]; % sequence length per species 
sigma0 = 1.0;
nset = 50;

sf = 0.1 ;  %% noise variance across species 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% sceIII-1, III-2, III-3 
for n_specific = 0:2    % number of species-specific states 

  %%%%%% repeat 
  for rr=1:nset 

	% self-transition 
	p0 = [ 0.9  0.9  0.9 ];

	T0 = 0.01*eye(K) + abs( sprand( K, K, 0.3 ) );
	% transition 
	for cc=1:C 
		Tc = full( T0 + abs( sprand(K, K, sf ) ) );
		Tc = Tc - diag(diag(Tc)) + 0.2*diag(rand(K,1)); 
		Tc(1,:) = rand(1, size(Tc,2) );
		ssum = find( sum( Tc - diag( diag( Tc ) ), 2 ) == 0 )';
		if ( length(ssum) > 0 )
			for kk=1:length(ssum)
				Tc(ssum(kk),max(1, ssum(kk)-1)) = rand();
			end
		end
		Tc = Tc ./ ( sum( Tc, 2 )*ones(1,K) );
		pi{cc} = p0(cc)*eye(K) + (1-p0(cc)) * Tc  ;
	end

    	% mean emission
    	mu0 = zeros(K,M);
    	is = 1;
    	rid = [ randperm(M) randperm(M) randperm(M) ];
    	for kk=1:K
        	nf = ceil( rand()*2.2 );   %% number of enriched mark 1 ~ 3  
        	ie = is+nf-1;
        	if ( kk < K )
                	mu0(kk,rid(is:ie)) = abs(1 + 3.0*rand(1,nf));
        	else
                	mu0(kk,rid(is:ie)) = abs(0.05 + 0.95*rand(1,nf));
        	end
        	is = ie+1;
    	end

    	% emission 
    	zr = find( mu0(:) > 0 );
    	mu{1} = mu0;
    	for cc=2:C
       	 mu{cc} = mu0;
       	 sgv = 1+(-1)^(cc-1)*sf  ;
       	 mu{cc} = mu{cc} .* sgv ;
    	end
    	if  ( n_specific > 0 ) 
    		[sr, sc] = ind2sub( size(mu0), zr(1) );
    		[sr2, sc2] = ind2sub( size(mu0), zr(2) );
    		if ( length( find( mu0(sr,:) > 0 ) ) > 1 )
       		 	mu{2}(zr(1)) = 0;
    			mu{3}(zr(1)) = 0;
 	   	elseif ( length( find( mu0(sr2,:) > 0 ) ) > 1 )
        		mu{2}(zr(2)) = 0;
    			mu{3}(zr(2)) = 0;
   	 	end
		if ( n_specific == 2 ) 
    			[sr, sc] = ind2sub( size(mu0), zr(3) );
    			zeromark = find( mu{2}(sr,:) == 0 );
 	   		mu{2}(sr,zeromark(1)) = abs( 1+3*rand(1,1));
		end	
    	end
   	 % observation and hidden sequence 
    	for cc=1:C
		S{cc}{1}(1) = 1 + sum(rand() > cumsum( pi{cc}( 1, : )));
		j = S{cc}{1}(1);
		Y{cc}{1}(:,1) = mu{cc}(j,:)' + sigma0*randn( M, 1 );  
		for tt=2:vT(cc)
			k = 1 + sum(rand() > cumsum( pi{cc}( S{cc}{1}(tt-1),: )));
			S{cc}{1}(tt) = k;
			Y{cc}{1}(:,tt) = mu{cc}(k,:)' + sigma0*randn( M, 1 );  
		end
    	end
    	param.p0 = p0;
   	param.Mus = mu;
  	param.Pi = pi;
	Strue = S;
	fprintf(1, 'simulation_III_mat/%s%d_sf_%.1f_sigma%.1f_%d.mat\n', sce0, n_specific+1, sf, sigma0, rr );
	save( sprintf('simulation_III_mat/%s%d_sf_%.1f_sigma%.1f_%d.mat', sce0, n_specific+1, sf, sigma0, rr),'Y', 'Strue', 'param' );

  end  %% end of rr 

end   %%% end of n-specific 
