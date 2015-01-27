function [S] = Viterbi( Y, Pi, Mus, sigma2 ) % Viterbi in log-scale 
%

K = size( Pi, 1 );
logPi = log( Pi + 1.e-20 );
logminf = -50;
if ( nargin < 4 )
	sigma2 = 0.8;
end
%sigma2 = 1.5;

T = size(Y, 2);
S = ones( 1, T );

len = 100000;
vte = [ len:len:T-ceil(len/2) T ];
vts = [ 1 vte(1:end-1) ];
for ss=1:length(vts)
    ts = vts(ss);
    te = vte(ss);
    fprintf( '[ %d, %d ] / %d \n', ts,te, T );
    T1 = te-ts+1;
    dyn_prog = zeros(K, T1);
    psi = zeros(K, T1);
    if ( ss== 1 )
    	dyn_prog(:,1) = logPi(1,1:K) ;
    else
	dyn_prog(:,1) = logminf;
        dyn_prog( S(ts), 1 ) = 0;
    end
    for k=1:K
        dyn_prog(k,1) = (-0.5*(Y(:,ts) - Mus(:,k))'*(Y(:,ts) - Mus(:,k))/sigma2) + dyn_prog(k,1);
    end
    dyn_prog(:,1) = dyn_prog(:,1) - min(dyn_prog(:,1));
    for t=2:T1
        for k=1:K
        	[ dyn_prog(k,t), psi(k,t) ] = max( dyn_prog(:,t-1) + logPi(:,k) );
            	dyn_prog(k,t) = (-0.5*(Y(:,t+ts-1) - Mus(:,k))'*(Y(:,t+ts-1) - Mus(:,k))/sigma2) + dyn_prog(k,t);
        end
        dyn_prog(:,t) = dyn_prog(:,t) - min(dyn_prog(:,t));
    end
  
    [ p, S(te) ] = max( dyn_prog(:,T1) ); 
    
%    disp('Backtrack to sample a path through the HMM.');
    for t=T1-1:-1:1
	S(ts+t-1) = psi( S(ts+t), t+1 );
    end
end


