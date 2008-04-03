function EEG = pop_bdca_train(EEG,R,sigma,lt,nut,ls,nus,supportframes,folds)
% EEG = pop_bdca_train(EEG,R,sigma,lt,nut,ls,nus,supportframes,folds)

% process commandline
if nargin < 9
  prompt={'Number of components','sigma','lt (ms)','nut','ls (electrode space)','nus','support frames, []=all','folds'};
  name='Input for pop_bilin';
  numlines=1;
  defaultanswer={'1','0.5','50','2.5','0.1','100','','5'};
  try 
    defaultanswer{1} = sprintf('%i',EEG.bdca.cht.R);
    defaultanswer{2} = sprintf('%f',EEG.bdca.cht.sigma);
    defaultanswer{3} = sprintf('%f',EEG.bdca.cht.lt);
    defaultanswer{4} = sprintf('%f',EEG.bdca.cht.nut);
    defaultanswer{5} = sprintf('%f',EEG.bdca.cht.ls);
    defaultanswer{6} = sprintf('%f',EEG.bdca.cht.nus);
    defaultanswer{7} = sprintf('%i ',EEG.bdca.cht.supportframes);
    defaultanswer{8} = sprintf('%i',EEG.bdca.cht.folds);
  end
  try 
    defaultanswer{1} = sprintf('%i',R);
    defaultanswer{2} = sprintf('%f',sigma);
    defaultanswer{3} = sprintf('%f',lt);
    defaultanswer{4} = sprintf('%f',nut);
    defaultanswer{5} = sprintf('%f',ls);
    defaultanswer{6} = sprintf('%f',nus);
    defaultanswer{7} = sprintf('%i ',supportframes);
    defaultanswer{8} = sprintf('%i',folds);
  end
  answer=inputdlg(prompt,name,numlines,defaultanswer);
  if isempty(answer), return, end
  R              = eval(sprintf('[%s]',answer{1}));  EEG.bdca.cht.R             = R;
  sigma          = eval(sprintf('[%s]',answer{2}));  EEG.bdca.cht.sigma         = sigma;
  lt             = eval(sprintf('[%s]',answer{3}));  EEG.bdca.cht.lt            = lt;
  nut            = eval(sprintf('[%s]',answer{4}));  EEG.bdca.cht.nut           = nut;
  ls             = eval(sprintf('[%s]',answer{5}));  EEG.bdca.cht.ls            = ls;
  nus            = eval(sprintf('[%s]',answer{6}));  EEG.bdca.cht.nus           = nus;
  supportframes  = eval(sprintf('[%s ]',answer{7}));  EEG.bdca.cht.supportframes = supportframes;
  folds          = eval(sprintf('[%s]',answer{8}));  EEG.bdca.cht.folds         = folds;
end


% train
spaceunits = [];
if ~isempty(EEG.chanlocs)
  for d = 1:length(EEG.chanlocs)
    spaceunits(d,1:3) = [EEG.chanlocs(d).X EEG.chanlocs(d).Y EEG.chanlocs(d).Z];
  end
end
  
c1_idx = find(EEG.bdca.labels==1);
c0_idx = find(EEG.bdca.labels==0);

if isempty(supportframes), supportframes = 1:EEG.pnts; end

if folds<2 % train
%  [u,t,b] = bilin_train(cat(3,EEG.data(:,:,c0_idx),EEG.data(:,:,c1_idx)),cat(2,EEG.bdca.labels(c0_idx),EEG.bdca.labels(c1_idx)),sigma,ls,nus,EEG.srate * lt / 1000,nut,R,spaceunits);
 
  sigw0 = 5;
  gpa = [sigma ls nus];
  gpb = [sigma (EEG.srate * lt / 1000) nut];
  
  [b,u,t] = bilinlogistregmultigp(double(cat(3,EEG.data(:,supportframes,c0_idx),EEG.data(:,supportframes,c1_idx))),...
                                  cat(2,EEG.bdca.labels(c0_idx),EEG.bdca.labels(c1_idx)),...
                                  R,...
                                  sigw0,...
                                  gpa,...
                                  gpb,...
                                  spaceunits);
  
  EEG.bdca.cht.u = u;
  EEG.bdca.cht.t = t;
  EEG.bdca.cht.b = b;
  
else % cv
  
  huskeL=[];
  huskefi=[];
  huskepot = [];
  
  fold = folds;
  for lo = 1:fold  , lo
    
    [c0_idx_train,c0_idx_validate] = loo(c0_idx,fold,lo);
    [c1_idx_train,c1_idx_validate] = loo(c1_idx,fold,lo);
    
    % train
    
    sigw0 = 5;
    gpa = [sigma ls nus];
    gpb = [sigma (EEG.srate * lt / 1000) nut];
    
    [b,u,t] = bilinlogistregmultigp(double(cat(3,EEG.data(:,supportframes,c0_idx_train),EEG.data(:,supportframes,c1_idx_train))),...
			  cat(2,EEG.bdca.labels(c0_idx_train),EEG.bdca.labels(c1_idx_train)),...
			  R,...
			  sigw0,...
			  gpa,...
			  gpb,...
			  spaceunits);
    
    EEG.bdca.cvstats.cht{lo} = {b,u,t,c0_idx_train,c0_idx_validate,c1_idx_train,c1_idx_validate};
    
%    [u,t,b] = bilin_train(double(cat(3,EEG.data(:,supportframes,c0_idx_train),EEG.data(:,supportframes,c1_idx_train))),...
%			  cat(2,EEG.bdca.labels(c0_idx_train),EEG.bdca.labels(c1_idx_train)),...
%			  sigma,ls,nus,EEG.srate * lt / 1000,nut,R,spaceunits);
    
    % validate
    
    [Ey,pot] = bilinlogistregmultigp_run(cat(3,EEG.data(:,supportframes,c0_idx_validate),EEG.data(:,supportframes,c1_idx_validate)), b,u,t);
    
    huskeL = [huskeL, cat(2,EEG.bdca.labels(c0_idx_validate),EEG.bdca.labels(c1_idx_validate))];
    huskefi= [huskefi, Ey];
    huskepot = [huskepot, pot];
  end
  size(huskeL)
  size(huskefi)
  A = auc(huskeL,huskefi)
  loglik = sum( huskeL.*huskepot - log(1+exp(huskepot)) )
  EEG.bdca.cvstats.Az = A;
  EEG.bdca.cvstats.loglik = loglik;
  EEG.bdca.cvstats.huskeL = huskeL;
  EEG.bdca.cvstats.huskeEy = huskefi;
  EEG.bdca.cvstats.huskepot = huskepot;
    
%  figure,rocpoints(huskeL,huskefi);
end



% huske
EEG.bdca.cht.R             = R;  
EEG.bdca.cht.sigma         = sigma;
EEG.bdca.cht.lt            = lt;
EEG.bdca.cht.nut           = nut;
EEG.bdca.cht.ls            = ls;
EEG.bdca.cht.nus           = nus;
EEG.bdca.cht.supportframes = supportframes;