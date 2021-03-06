function [ logp ] = hiHmmMultNormalJointLogLikelihood2( S, Y, Beta, alpha0, mu_0, sigma2_0, sigma2)
%HIHMMNORMALJOINTLOGLIKELIHOOD Computes the joint log likelihood of 
% generating a particular hidden state sequence and emission string.
        
C = length(Y);
M = size(Y{1}{1},1);
K = length(Beta)-1;
N = zeros(K,K,C);
E = zeros(K,3);

% Compute all transition counts and emission sufficient statistics.
for cc=1:C 
 for j=1:length(S{cc})
  Sc = S{cc}{j};
  N(1,Sc(1),cc) = 1;
  E(Sc(1), 1) = E(Sc(1), 1) + sum(Y{cc}{j}(:,1));
  E(Sc(1), 2) = E(Sc(1), 2) + 1*M;
  E(Sc(1), 3) = E(Sc(1), 3) + sum(Y{cc}{j}(:,1).^2);
  for t=2:length(Sc)
    N(Sc(t-1), Sc(t), cc) = N(Sc(t-1), Sc(t), cc) + 1;
    E(Sc(t), 1) = E(Sc(t), 1) + sum(Y{cc}{j}(:,t));
    E(Sc(t), 2) = E(Sc(t), 2) + 1*M;
    E(Sc(t), 3) = E(Sc(t), 3) + sum(Y{cc}{j}(:,t).^2);
  end
 end
end

% Compute the log likelihood.
logp = 0;
for k=1:K
  for cc=1:C 
    R = [N(k,:,cc) 0] + alpha0 * Beta;
    ab = alpha0 * Beta;
    nzind = find(R ~= 0);
    % Add transition likelihood.
    logp = logp + gammaln(alpha0) ...
            - gammaln(sum([N(k,:,cc) 0]) + alpha0) ...
            + sum(gammaln( R(nzind)  )) ...
            - sum(gammaln( ab(nzind) ));
  end
  % Add emission likelihood.
  sigma2_n = 1 / (1 / sigma2_0 + E(k, 2) / sigma2);
  mu_n = (mu_0 / sigma2_0 + E(k,1) / sigma2) * sigma2_n;
  logp = logp + 0.5 * log(sigma2_n) ...
                - 0.5 * (log(sigma2_0) + E(k,2)/2 * log(2*pi*sigma2)) ...
                - 0.5 * (E(k,3)/sigma2 + mu_0^2 / sigma2_0 - mu_n^2 / sigma2_n);
end


