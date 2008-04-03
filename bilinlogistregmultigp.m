function [w0,a,b,Kl,Kr] = bilinlogistregmultigp(x_train,y_train,R,sigw0,gpa,gpb,spaceunits,iter,initguess)
% [w0,a,b,Ka,Kb] = bilinlogistregmultigp(x_train,y_train,R,sigw0,gpa,gpb,spaceunits)
%
% synopsis:
%
%    [w0,a,b,Ka,Kb] = bilinlogistregmultigp(x_train,y_train,R,sigw0)
%    [w0,a,b,Ka,Kb] = bilinlogistregmultigp(x_train,y_train,R,sigw0,gpa,gpb,spaceunits)
%
% inputs:
%
%     x_train  : D-by-T-by-N data
%     y_train  : 1-by-N, binary labels in {0,1}
%           R  : number of components
%       sigw0  : prior stddev of w0
%     gpa,gpb  : prior covariance parameters for a and b resp.
%                : if left out, defaults to sigw0.
%                : otherwise, a scalar gives the prior stddev.
%                : otherwise, a vector [sigma l] defines a GP with SE cov. function.
%                : otherwise, a vector [sigma l nu] defines a GP with Matern cov. function
%  spaceunits  : (optional) spatial coords in euclidian basis [x y z; x y z; etc].
%
%  
%   outputs:
%
%        w0,a,b : you know.
%         Ka,Kb : covariancematrices for a and b respectively.  
%
% Dependencies: OGP, IMMOPTIBOX
%
% Author: Mads Dyrholm

[I,J,N] = size(x_train);

if nargin<8, iter=inf; end
opts = [1e-6 1e-4 1e-8 min(iter,10000)];

Ka = sigw0;
fprintf('Prior std.dev. for w0: %f\n',sigw0);

if nargin<5
  gpa = sigw0;
  gpb = sigw0;
end

if (nargin>=5) & ~isempty(spaceunits)
  A = 0;
  for c=1:3
    A = A + (repmat(spaceunits(:,c),[1,I]) - repmat(spaceunits(:,c),[1,I])').^2; 
  end
  A = sqrt(A);
else
  A = I;
end

if length(gpa)==1
  Kl = gpa*eye(I);
  fprintf('Prior std.dev. for a: %f\n',gpa);
elseif length(gpa)==2
  Kl = gpa(1)*sqexp(A,gpa(2));
  fprintf('Prior std.dev. for a: %f, SE cov. function l=%f\n',gpa(1),gpa(2));
elseif length(gpa)==3
  Kl = gpa(1)*Kmatern(A,gpa(2),gpa(3));
  fprintf('Prior std.dev. for a: %f, Matern cov. function l=%f, nu=%f\n',gpa(1),gpa(2),gpa(3));
end

if length(gpb)==1
  Kr = gpb*eye(J);
  fprintf('Prior std.dev. for b: %f\n',gpb);
elseif length(gpb)==2
  Kr = gpb(1)*sqexp(J,gpb(2));
  fprintf('Prior std.dev. for b: %f, SE cov. function l=%f\n',gpb(1),gpb(2));
elseif length(gpb)==3
  Kr = gpb(1)*Kmatern(J,gpb(2),gpb(3));
  fprintf('Prior std.dev. for b: %f, Matern cov. function l=%f, nu=%f\n',gpb(1),gpb(2),gpb(3));
end
fprintf('************\n');
%for retr=1:20
%  if nargin<9, initguess = randn((I+J)*R+1,1)/100; end
%  if 0
%    [X, info] = dampnewton('cost_logistreg_bilinmultigp', initguess(:),opts,x_train,y_train,Ka,Kl,Kr);
%  else

if nargin<9, initguess = randn((I+J)*R+1,1)/100; end
x0 = initguess;
if R==1, REGO = 1; else REGO = R+1; end
for rego = 1:REGO
  rego
  for comp = 1:R
    fprintf('Alternate %i\n',comp);
    subidx = [1 , 1+((1+(comp-1)*(I+J)):(comp*(I+J)))];
    for retr=1:20
      [X, info] = dampnewton(@cost_logistreg_bilinmultigp_alt, ...
			     x0(subidx),opts,x_train,y_train,...
			     Ka,...
			     pinv(Kl), sum(2*log(diag(chol(Kl)))),...
			     pinv(Kr), sum(2*log(diag(chol(Kr)))),...
			     comp,x0,subidx);
      if ~isnan(info(1)), break; end
      fprintf('Convergence problem detected!! ');
      if retr==20, 
	fprintf('Giving up now. :(\n'); 
      else 
	fprintf('Trying again...\n'); 
	if nargin<9, initguess = randn((I+J)*R+1,1)/100; end
	x0(subidx) = initguess(subidx);
      end
    end
    x0(subidx) = X; 
  end
  X = x0;
end
%opts = [1  1e-4  1e-8  10000];
%[X, info] = ucminf('cost_logistreg_bilinmulti', randn((I+J)*R+1,1)/100,opts,[],x_train,y_train);

info
w0 = X(1);
a = [];
b = [];
X = X(2:end);
for r=1:R
  a = [a X((r-1)*(I+J)+(1:I))];
  b = [b X((r-1)*(I+J)+(I+1:I+J))];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K = sqexp(dim,l)
if sum(size(dim))==2
  K = repmat(1:dim,[dim 1]);
  K = abs(K-K')/l;
else
  K = dim/l;
end
K = exp(-K.^2 ./ (2*(l+eps)^2));

function K = Kmatern(dim,l,nu)
if sum(size(dim))==2
  K = repmat(1:dim,[dim 1]);
  K = abs(K-K')/l;
else
  K = dim/l;
end
K = matern(nu,1,K);

function [f,G,H0,sal] = cost_logistreg_bilinmultigp_alt(x0sub,X3,y,Ka,invKl,logdetKl,invKr,logdetKr,comp,x0,subidx)
% COST_LOGISTREG_BILINMULTI Bilinear Logistic Regression costfunction for IMMOPTIBOX
%   w. alternating
%
%    [f,G,H] = cost_bilinlogistreg_alt(x0sub,X3,y,Ka,invKl,logdetKl,invKr,logdetKr,comp,x0,subidx)
%
%        x0 = [alpha;beta_l;beta_r];
%
% Author: Mads Dyrholm
% COPYRIGHT! Mads Dyrholm
x0(subidx) = x0sub;
[I,J,N] = size(X3);
R = (length(x0)-1)/(I+J);
H = zeros((I+J)*1+1); % ALT - zeros((I+J)*R+1); % hessian
betas_l = [];
betas_r = [];
for idx = 2:(I+J):length(x0)-1
  betas_l = [betas_l x0(idx+(0:I-1),1)];
  betas_r = [betas_r x0(idx+(I:I+J-1),1)];
end
alpha  = x0(1);
y = y(:);

%imagesc(Kl), drawnow
% kernel %%%%%%

%invKl = pinv(Kl);
%invKr = pinv(Kr);
%logdetKl = sum(2*log(diag(chol(Kl))));
%logdetKr = sum(2*log(diag(chol(Kl))));

Ckl = -0.5*I*log(2*pi) -0.5*logdetKl;
Ckr = -0.5*J*log(2*pi) -0.5*logdetKr;
Cka = -0.5*1*log(2*pi) -0.5*log(Ka);

%%%%%%%%%%%%%%%

[pies_it,bxb] = pies(X3,alpha,betas_l,betas_r);
loglik = sum( y.*(alpha+bxb) - log(1+exp(alpha+bxb)) );
logPl = 0;
logPr = 0;
logPa = 0;
for k = 1:R
  logPl  = logPl + Ckl - 0.5*betas_l(:,k)'*(invKl*betas_l(:,k));
  logPr  = logPr + Ckr - 0.5*betas_r(:,k)'*(invKr*betas_r(:,k));
end
logPa  = logPa + Cka - 0.5*(alpha^2)/Ka;
f = -loglik -logPl -logPr -logPa;
fprintf('*** log likelihood: %f,\tlog posterior: %f\n',loglik,-f);

% gradient alpha
dl_da = sum(y - pies_it);
dl_da = dl_da -alpha/Ka; %% prior
jaco = dl_da;
for k = comp % ALT - 1:R
  % gradient beta_l  
  xibr = reshape(reshape(permute(X3,[1 3 2]),[I*N,J])*betas_r(:,k),[I,N]);  % I,N
  dl_dbl = xibr * (y-pies_it);
  dl_dbl = dl_dbl -invKl*betas_l(:,k); %% prior
  
  % gradient beta_r
  blxi = reshape(betas_l(:,k)' * reshape(X3,[I,J*N]),[J,N]);   % J,N
  dl_dbr = blxi * (y-pies_it);
  dl_dbr = dl_dbr -invKr*betas_r(:,k); %% prior
  
  % output gradient
  jaco  = [jaco;dl_dbl;dl_dbr];
  
  %%%%%%%%%%%%%%%%%%   HESSIAN  %%%%%%%%%%%%%%%%%%%
  if nargout>2    
    IDX = 0; % ALT - (k-1)*(I+J);
    % alpha alpha
    H(1) = -pies_it'*(1-pies_it);
    % alpha beta_l
    H(IDX+(2:I+1),1) = - xibr * (pies_it.*(1-pies_it));
    H(1,IDX+(2:I+1)) = H(IDX+(2:I+1),1)';
    % alpha beta_r
    H(IDX+((I+2):(I+1+J)),1) = -  blxi * (pies_it.*(1-pies_it));
    H(1,IDX+(I+2:I+1+J)) = H(IDX+(I+2:I+1+J),1)';
    for kmark = comp % ALT - 1:R
      xibrmark = reshape(reshape(permute(X3,[1 3 2]),[I*N,J])*betas_r(:,kmark),[I,N]);  % I,N
      blximark = reshape(betas_l(:,kmark)' * reshape(X3,[I,J*N]),[J,N]);   % J,N
      IDXmark = 0; % ALT - (kmark-1)*(I+J);
      % beta_l beta_l
      H(IDX+(2:I+1),IDXmark+(2:I+1))    = - xibr * ( repmat(pies_it.*(1-pies_it),[1,I]) .* xibrmark' );
      if k==kmark
	H(IDXmark+(2:I+1),IDX+(2:I+1)) = H(IDX+(2:I+1),IDXmark+(2:I+1)) - invKl; %% prior 
      else
	H(IDXmark+(2:I+1),IDX+(2:I+1)) = H(IDX+(2:I+1),IDXmark+(2:I+1));
      end
      % beta_r beta_r
      H(IDX+(I+2:I+1+J),IDXmark+(I+2:I+1+J))= - blxi * ( repmat(pies_it.*(1-pies_it),[1,J]) .* blximark' );
      if k==kmark
	H(IDXmark+(I+2:I+1+J),IDX+(I+2:I+1+J)) = H(IDX+(I+2:I+1+J),IDXmark+(I+2:I+1+J)) - invKr; %% prior
      else
	H(IDXmark+(I+2:I+1+J),IDX+(I+2:I+1+J)) = H(IDX+(I+2:I+1+J),IDXmark+(I+2:I+1+J));
      end
      % beta_r beta_l
      ww = reshape(X3,[I*J,N]) * (y-pies_it); % I*J,1 --- i er snap
      
      %qq = reshape(permute(X3,[1 3 2]),[I*N,J]) * betas_r(:,k);  % I*N,1 --- i er snap
      %qq = reshape(qq,[I,N]);                              % I,N
      qq = repmat(xibr,[J,1]);                                % I*J,N --- i er snap
      %%%qq = repmat(qq,[J,1]);                               % I*J,N --- i er snap
      
      %ee = betas_l(:,kmark)' * reshape(X3,[I,J*N]);        % 1,J*N
      %ee = reshape(ee,[J,N]);                    % J,N
      ee = repmat(blximark,[1 1 I]);                   % J,N,I
      %%%ee = repmat(ee,[1 1 I]);                   % J,N,I
      ee = permute(ee,[3 1 2]);                  % I,J,N
      ee = reshape(ee,[I*J,N]);                  % I*J,N
      
      qq = qq .* ee;               % 
      
      qq = qq * (pies_it.*(1-pies_it)); %( exp(-alpha-bxb) ./ (1+exp(-alpha-bxb)).^2 ); % I*J,1 
      
      H(IDX+(2:I+1),IDXmark+(I+2:I+1+J)) = reshape(ww*(k==kmark) - qq,[I,J]);
      H(IDXmark+(I+2:I+1+J),IDX+(2:I+1)) = H(IDX+(2:I+1),IDXmark+(I+2:I+1+J))';
    end
    %  [R,p] = chol(-H);
    %  if p==0, disp('H IS NEGATIVE DEFINITE'); end   
    %  
    H0 = -H;
    
  end
  
end
G = -jaco;

% symmetric trick !!!!
H0 = (H0+H0')/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pie,bxb] = pies(X,alpha,betas_l,betas_r)
R = size(betas_l,2);
[I,J,N] = size(X);
X = reshape(X,[I,J*N]);
bxb = 0;
for k=1:R
  bxbplus = betas_l(:,k)'*X;
  bxbplus = (reshape(bxbplus,[J,N])'*betas_r(:,k));
  bxb = bxb + bxbplus;
end
pie = 1./(1+exp(-alpha-bxb));
