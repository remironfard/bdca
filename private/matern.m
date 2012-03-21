function covf = matern(nu,sigma,ddd)
% Matern
tmp  = sqrt(2*nu)*ddd;
covf = exp(log(besselk(nu, tmp, 1)) - tmp + nu*log(tmp) + log(sigma)-(nu-1)*log(2)-gammaln(nu));
covf(1:(size(covf,1)+1):end) = sigma;
