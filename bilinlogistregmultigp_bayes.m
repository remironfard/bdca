function [BIC,Laplace] = bilinlogistregmultigp_bayes(w0,a,b,X3,y,sigw0,Kl,Kr)
% [BIC,Laplace] = bilinlogistregmultigp_bayes(w0,a,b,X3,y,sigw0,Kl,Kr)
%
% Author: Mads Dyrhholm
[f,G,H0] = cost_logistreg_bilinmultigp(w0,a,b,X3,y,sigw0,Kl,Kr);
logdetH0 = 2*sum(log(diag(chol(H0))));
dim_theta = length(w0) + length(a(:)) + length(b(:));
Laplace = -f + 0.5 * dim_theta * log(2*pi) - 0.5 * logdetH0;
BIC = -f - 0.5 * dim_theta * log(size(X3,3));


function [f,G,H0] = cost_logistreg_bilinmultigp(w0,a,b,X3,y,Ka,Kl,Kr)
% COST_LOGISTREG_BILINMULTI Bilinear Logistic Regression costfunction for IMMOPTIBOX
%
%    [f,G,H] = cost_bilinlogistreg(w0,a,b,X3,y,Ka,Kl,Kr)
%
% Author: Mads Dyrholm
[I,J,N] = size(X3);
R = size(a,2);
H = zeros((I+J)*R+1); % hessian
beta_l = a;
beta_r = b;
alpha  = w0; %x0(1);
y = y(:);

%imagesc(Kl), drawnow
% kernel %%%%%%
invKl = pinv(Kl);
invKr = pinv(Kr);
logdetKl = sum(2*log(diag(chol(Kl))));
logdetKr = sum(2*log(diag(chol(Kl))));

Ckl = -0.5*I*log(2*pi) -0.5*logdetKl;
Ckr = -0.5*J*log(2*pi) -0.5*logdetKr;
Cka = -0.5*1*log(2*pi) -0.5*log(det(Ka));

%%%%%%%%%%%%%%%

[pies_it,bxb] = pies(X3,alpha,betas_l,betas_r);
loglik = sum( y.*(alpha+bxb) - log(1+exp(alpha+bxb)) );
logPl = 0;
logPr = 0;
logPa = 0;
for k=1:R
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
for k = 1:R
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
    IDX = (k-1)*(I+J);
    % alpha alpha
    H(1) = -pies_it'*(1-pies_it);
    % alpha beta_l
    H(IDX+(2:I+1),1) = - xibr * (pies_it.*(1-pies_it));
    H(1,IDX+(2:I+1)) = H(IDX+(2:I+1),1)';
    % alpha beta_r
    H(IDX+((I+2):(I+1+J)),1) = -  blxi * (pies_it.*(1-pies_it));
    H(1,IDX+(I+2:I+1+J)) = H(IDX+(I+2:I+1+J),1)';
    for kmark = 1:R
      xibrmark = reshape(reshape(permute(X3,[1 3 2]),[I*N,J])*betas_r(:,kmark),[I,N]);  % I,N
      blximark = reshape(betas_l(:,kmark)' * reshape(X3,[I,J*N]),[J,N]);   % J,N
      IDXmark = (kmark-1)*(I+J);
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
    
    % saliency
%    sal = 0.5*diag(H0).*(x0.^2);
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
