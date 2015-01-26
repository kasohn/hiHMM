function [ Mus ] = hiHmmSampleMultNormalMeans2( S, Y, K, sigma2, mu_0, sigma2_0 )
%% SAMPLENORMALMEANS
% Samples a mean for a normal distribution with variance
% sigma2 for every state K given a sequence S, a corresponding observation
% vector Y and a prior normal with mean mu_0 and variance sigma2_0.

C = length(S);
M = size(Y{1}{1},1);
Mus = zeros(M,K);

for k=1:K
    Ys = [];
    N = 0;
    for cc=1:C 
     for j=1:length(S{cc})
	idxK = find( S{cc}{j} == k );
	Ys = [Ys Y{cc}{j}(:,idxK)];
	N = N + length(idxK);
     end
    end
    if N == 0
        Mus(:,k) = randn(M,1) * sqrt(sigma2_0) + mu_0;
    else
        Mu_ml = sum(Ys,2) / N;                                % The ML mean.
        Mu_n = (sigma2 * mu_0 + N * sigma2_0 * Mu_ml) ...   % The posterior mean.
               / (N * sigma2_0 + sigma2);
        S2_n = 1 / (1 / sigma2_0 + N / sigma2);             % The posterior variance.
        Mus(:,k) = randn(M,1) * sqrt(S2_n) + Mu_n;
    end
end


