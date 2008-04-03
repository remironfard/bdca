function EEG = pop_bdca_constrica(EEG,alpha)
% synopsis:
%
%   EEG = pop_bdca_constrica(EEG,alpha)
%
% input:
%
%   alpha  :  Set it so that the algorithm converges smoothly.
%             try e.g. alpha = 1000 
%
% output:
%
%   EEG.bdca.ica.G      : ambiguity matrix
%   EEG.bdca.ica.alpha  : remember the alpha used
%
% Author: Mads Dyrholm

% process commandline
if nargin < 2
  prompt={'Alpha'};
  name='Input for sidekick';
  numlines=1;
  defaultanswer={'1000'};
  try 
    defaultanswer{1} = sprintf('%i',EEG.bdca.ica.alpha);
  end
  try 
    defaultanswer{1} = sprintf('%i',alpha);
  end
  answer=inputdlg(prompt,name,numlines,defaultanswer);
  if isempty(answer), return, end
  alpha              = eval(sprintf('[%s]',answer{1}));
end

if ~isempty(alpha)
  [G,alpha,S,Ci,Co,info] = constrica(double(EEG.data(:,EEG.bdca.cht.supportframes,:)),EEG.bdca.cht.u,EEG.bdca.cht.t,alpha);
else
  [G,alpha,S,Ci,Co,info] = constrica(double(EEG.data(:,EEG.bdca.cht.supportframes,:)),EEG.bdca.cht.u,EEG.bdca.cht.t);
end
EEG.bdca.ica.G = G;


% huske
EEG.bdca.ica.alpha             = alpha;  
