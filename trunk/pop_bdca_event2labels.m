function [EEG,LASTCOM] = pop_bdca_event2labels(EEG,contrast,versus)
% [EEG,LASTCOM] = pop_bdca_event2labels(EEG,contrast,versus)

% process commandline
if nargin < 3
  % 'Contrast event types','Versus event types',
  prompt={'Contrast events','Versus'};
  name='Input for pop_bdca_event2labels';
  numlines=1;
  defaultanswer={'',''};
  try 
    defaultanswer{1} = sprintf('%i',contrast);
    defaultanswer{2} = sprintf('%i',versus);
  end
  answer=inputdlg(prompt,name,numlines,defaultanswer);
  if isempty(answer), return, end
  contrast  = eval(sprintf('[%s]',answer{1}));
  versus     = eval(sprintf('[%s]',answer{2}));
end


% find labels etc.
labels = nan * zeros(1,EEG.trials);
for epo = 1:EEG.trials
  if iscell(EEG.epoch(epo).eventlatency)
    idx = find([EEG.epoch(epo).eventlatency{:}] == 0);
    tmp = [EEG.epoch(epo).eventtype{idx}];
  else
    idx = find([EEG.epoch(epo).eventlatency] == 0); 
    tmp = [EEG.epoch(epo).eventtype(idx)];
  end
  
  if ismember(tmp,contrast)
    labels(epo) = 1;
  end
  if ismember(tmp,versus)
    labels(epo) = 0;
  end
end


EEG.bdca.labels = labels;

try
  LASTCOM = sprintf('EEG = pop_bdca_event2labels(EEG,[%s],[%s]);',answer{1},answer{2});
catch
  LASTCOM = '';
end
