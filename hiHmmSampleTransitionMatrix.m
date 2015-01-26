function [ Pi ] = SampleTransitionMatrixhihmmc( S, H, binit, G  )
%SAMPLETRANSITIONMATRIXHIHMM Samples a transition matrix from a state sequence S, species label L
% and Dirichlet prior H (row vector).

if ( nargin < 3 )
	binit= 0;
end

C = length(S);
K = size(H,2);
Pi = zeros(K,K,C);

N = ones(K,K,C);
for cc=1:C 
 for j=1:length(S{cc})
  N(1,S{cc}{j}(1),cc) = 1;
  for t=2:length(S{cc}{j})
    	N(S{cc}{j}(t-1), S{cc}{j}(t),cc) = N( S{cc}{j}(t-1), S{cc}{j}(t),cc ) + 1;
  end
 end
end

for cc=1:C
  for k=1:K
    Pi(k, :, cc) = dirichlet_sample( N(k,:,cc) + H  );
  end
  if ( binit ) 
    Pi(:,:,cc ) = eye(K) * 0.8 + 0.2 * Pi(:,:,cc);
  end
end



