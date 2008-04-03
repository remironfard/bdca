function [Ey,pot] = bilinlogistregmultigp_run(X3,w0,a,b)
% [Ey,pot] = bilinlogistregmultigp_run(x,w0,a,b)
%
% Author: Mads Dyrholm
[I,J,N] = size(X3);
R = size(a,2);

betas_l = a;
betas_r = b;
alpha  = w0;

[I,J,N] = size(X3);
X = reshape(X3,[I,J*N]);
bxb = 0;
for k = 1:R
  bxbplus = betas_l(:,k)'*X;
  bxbplus = (reshape(bxbplus,[J,N])'*betas_r(:,k));
  bxb = bxb + bxbplus;
end
bxb = bxb';
pot = +alpha+bxb;
Ey = 1./(1+exp(-pot));
