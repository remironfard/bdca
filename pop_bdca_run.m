function EEG = pop_bdca_run(EEG)
% Synopsis:
%
%    EEG = pop_bdca_run(EEG)
%
% Output:
%
%    EEG.bdca.Ey    : Ey(n)  = E[y(n)]
%    EEG.bdca.pot   : pot(n) = pi(Xn)
%
%    if G is present (i.e. if you ran POP_BDCA_CONSTRICA)
%
%       EEG.bdca.ica.S : Independent "source" across trials. 
%
% Author: Mads Dyrholm

supportframes  = EEG.bdca.cht.supportframes;

u = EEG.bdca.cht.u;
t = EEG.bdca.cht.t;
b = EEG.bdca.cht.b;

c1_idx = find(EEG.bdca.labels==1);
c0_idx = find(EEG.bdca.labels==0);

idx = sort([c1_idx,c0_idx]);

[Ey,pot] = bilinlogistregmultigp_run(EEG.data(:,supportframes,idx), b,u,t);
   
EEG.bdca.Ey  = Ey;
EEG.bdca.pot = pot;

try
  EEG.bdca.ica.S = bilin_sources(EEG.bdca.cht.u,EEG.bdca.cht.t,double(EEG.data(:,EEG.bdca.cht.supportframes,:)),EEG.bdca.ica.G);
end

% aux
function S = bilin_sources(U,V,X,G)
%  S = bilin_sources(U,V,X,G)

[D,T,N] = size(X);
if nargin<4, G = eye(size(U,2)); end

S = pinv(kr(V*inv(G),U*G'))*reshape(X,[D*T,N]);
