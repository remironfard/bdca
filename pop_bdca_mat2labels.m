function EEG = pop_bdca_mat2labels(EEG,labels)
% EEG = pop_bdca_mat2labels(EEG,labels)
%
% input
%
%    labels : can be the name of a variable, or a numeric vector
%             with the labels in {0,1}
%
% output
%
%    EEG.bdca.labels  :  the labels
% 

% process commandline
if nargin < 2
  % 
  prompt={'Event vector variable name'};
  name='Input for pop_bilin_mat2labels';
  numlines=1;
  defaultanswer={''};
  answer=inputdlg(prompt,name,numlines,defaultanswer);
  if isempty(answer), return, end
  labels = answer{1};
end

if isnumeric(labels)
 EEG.bdca.labels = labels;
else
  try
    EEG.bdca.labels = evalin('base', labels);
  catch
    error('Variable not found');
  end
end
EEG.bdca.labels = EEG.bdca.labels(:)'; % make row
