function [ ll ] = loglik( Data, nbStates, Priors, Mu, Sigma)
% Log-likelihood of data modeled with GMM

%Compute the new probability p(data|i)
for i=1:nbStates
    Pxi(:,i) = gaussPDF(Data, Mu(:,i), Sigma(:,:,i));
end

% Compute the log-likelihood
F = Pxi*Priors';
F(find(F<realmin)) = realmin;
ll = sum(log(F))/size(Data,2);

end

