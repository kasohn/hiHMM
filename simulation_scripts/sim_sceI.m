%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scenario I: genome size effect 

stream = RandStream('mt19937ar','seed',20); 
RandStream.setGlobalStream(stream);

sce0 = 'sceI';

% sequence length per species 
vT_all = [ 2000 2000 2000;
	   2000	5000 10000;
	   2000	5000 10000];

% self-transition  probability
p0_all = [  ...
	0.9 0.9 0.9; 
	0.9 0.9 0.9; 
	0.7 0.8 0.9 ];

C = 3;
K = 10;
M = 8;	% number of markers
sigma0 =  1; %% noise variance  
nset = 50;

%%%%%%%%%%% scenario I-1, I-2, I-3 %%% 
for ss= 1:3 

% sequence length per species 
vT = vT_all(ss,:); 

% self-transition  probability
p0 = p0_all(ss,:);

%% generate 50 datasets  
for rr=1:nset

  for sf = [ 0.1 ]	%% noise variance across species 

    T0 = 0.01*eye(K) + abs( sprand( K, K, 0.3 ) );
    % generate transition matrix
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

    % emission
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
    [sr, sc] = ind2sub( size(mu0), zr(1) );
    [sr2, sc2] = ind2sub( size(mu0), zr(2) );
    % one species-specific state 
    if ( length( find( mu0(sr,:) > 0 ) ) > 1 )
        mu{2}(zr(1)) = 0;
    	mu{3}(zr(1)) = 0;
    elseif ( length( find( mu0(sr2,:) > 0 ) ) > 1 )
        mu{2}(zr(2)) = 0;
    	mu{3}(zr(2)) = 0;
    end
    % observation and hiddne sequences 
    for cc=1:C
	S{cc}{1}(1) = 1 + sum(rand() > cumsum( pi{cc}( 1, : )));
	j = S{cc}{1}(1);
	Y{cc}{1}(:,1) = abs(mu{cc}(j,:)' + sigma0*randn( M, 1 ));  
	for tt=2:vT(cc)
		k = 1 + sum(rand() > cumsum( pi{cc}( S{cc}{1}(tt-1),: )));
		S{cc}{1}(tt) = k;
		Y{cc}{1}(:,tt) = abs(mu{cc}(k,:)' + sigma0*randn( M, 1 ));  
	end
    end
    param.p0 = p0;
    param.Mus = mu;
    param.Pi = pi;
    Strue = S;
    fprintf(1, 'simulation_I_mat/%s%d_sf_%.1f_sigma%.1f_%d.mat\n', sce0, ss, sf, sigma0, rr );
    save( sprintf('simulation_I_mat/%s%d_sf_%.1f_sigma%.1f_%d.mat', sce0, ss, sf, sigma0, rr),'Y', 'Strue', 'param' );

  end

end

end 	% end of ss (scenario)
