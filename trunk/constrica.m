function [G,alpha,S,Ci,Co,info] = constrica(X,U,V,alpha)
% CONSTRICA
%
%   [G,alpha,S,Ci,Co,info] = constrica(X,U,V)
%   [G,alpha,S,Ci,Co,info] = constrica(X,U,V,alphafix)
%   [G,alpha,S,Ci,Co,info] = constrica(X,U,V,1)
%
% Author: Mads Dyrholm
iter = 10000;
[D,T,N] = size(X);
R = size(U,2);

init = eye(R);
%init = randn(R); 
init=init(:);

if nargin<4
  alphainit = 100;
  init = [init;alphainit];
  [G,info] = dampnewton(@cost_constraintica_alphafix_outerhess, init,[nan nan nan iter], X,U,V);
  %  [G,info] = ucminf('cost_constraintica',init(:),[nan nan nan 1000],[],X,U,V); 
  alpha = G(end);
  G = reshape(G(1:R*R),[R R]);
else
  [G,info] = dampnewton(@cost_constraintica_alphafix_outerhess, init,[nan nan nan iter], X,U,V,alpha);
  %  [G,info] = ucminf('cost_constraintica_alpha',[init(:);alphainit],[],[],X,U,V); 
  G = reshape(G(1:R*R),[R R]);
end

if nargout>1
  if nargout>3, 
    S = pinv(kr(V,U))*reshape(X,[D*T,N]);
    Ci = corrcoef(S');
  end
  S = pinv(kr(V*inv(G),U*G'))*reshape(X,[D*T,N]);
  Co = corrcoef(S');
end
%for k =1:R
% S(k,:) = S(k,:)-mean(S(k,:));
%end
%for k =1:R
%  S(k,:) = S(k,:)/std(S(k,:))/alpha^2;
%end

function [f,J,H] = cost_constraintica_alphafix_outerhess(G,X,U,V,alpha)
% 
[D,T,N] = size(X);
R = size(U,2);
if nargin<5
  alpha = max(G(end),0);
  alphafix = 0;
else
  alphafix = 1;
end
%if alpha<0, f=inf; J = G(:);return, end
G = reshape(G(1:R*R),[R,R]);

UGt = U*G';
Ginv = inv(G);
VGinv = V*Ginv;
W = kr(VGinv,UGt);
Wpinv = pinv(W);
WWinv = inv(W'*W);
logp = inline('-log(cosh(x))-log(pi)');
psi = inline('-tanh(x)');
reshapeX = reshape(X,[D*T,N]);

S = Wpinv*reshapeX/(alpha^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Null space hack
for k =1:R
  S(k,:) = S(k,:)-mean(S(k,:));
end
for k =1:R
  S(k,:) = S(k,:)/std(S(k,:))/alpha^2;
end

loglik = ... % -2*N*R*log(alpha)
    -0.5*N*log(det(W'*W)) + sum(sum(logp(S)));
fprintf('Log likelihood: %f\n',loglik);

J = zeros(size(G));

if ~alphafix
  sqrtH_alpha = sum(psi(S) .* S)*(-2*alpha^(-3)) -2*R/alpha;
  dl_dalpha = sum(sqrtH_alpha);
end

sqrtH_G = zeros(R,R,N);

for i = 1:R
  for j = 1:R
    %
    ELij = zeros(R); ELij(i,j)=1;
    dW = kr(-V*(Ginv*ELij*Ginv),UGt) + kr(VGinv,U*ELij');
    d1 = -N*sum(sum(Wpinv .* dW'));
    %
    dWWinv = -WWinv * 2*sum(sum(dW.*W)) *WWinv;
    ds = (dWWinv * W' + WWinv * dW') * reshapeX;
    
    sqrtH_G(i,j,:) = sum(psi(S) .* ds) - sum(sum(Wpinv .* dW'));
    
    d2 = sum(sqrtH_G(i,j,:));
    %
    J(i,j) = d2;
  end
end

sqrtH = -[reshape(sqrtH_G,[R*R,N])]; % G
if ~alphafix
  sqrtH = [sqrtH; -sqrtH_alpha]; % alpha
  J = [J(:); dl_dalpha];
end
H = sqrtH*sqrtH';
J = -J(:);
f = -loglik;
