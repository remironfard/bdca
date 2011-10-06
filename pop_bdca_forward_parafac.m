function EEG = pop_bdca_forward_parafac(EEG)
% synopsis:
%
%   EEG = pop_bdca_forward_parafac(EEG)
%
% output:
%
%   EEG.bdca.fwd.parafac  :  the 'Factors' output from PARAFAC ( >> help parafac )
%
% dependencies:
%
%   NWAY toolbox
%
% Author: Mads Dyrholm

try 
  datoract = EEG.bdca.datoract;
catch
  datoract = 1;
end
if datoract==1
  DATORACT=EEG.data;
else
  DATORACT=EEG.icaact;
end

Options(5) = 1
OldLoad = {EEG.bdca.cht.u,EEG.bdca.cht.t,EEG.bdca.ica.S};
Factors = parafac(double(DATORACT(:,EEG.bdca.cht.supportframes,:)),...
		  EEG.bdca.cht.R,...
		  Options,... %Options
		  zeros(EEG.bdca.cht.R,1),... %const
		  OldLoad,...
		  [0 0 1]'); %,...%FixMode,
	          %Weights);
% icaact?		  
if datoract==0
  tmp = [];
  for ff=1:size(Factors{1},2)
    tmp = [tmp, EEG.icawinv * Factors{1}(:,ff)];
  end
  EEG.bdca.fwd.subspace = Factors{1};
  Factors{1} = tmp;
end

% output
EEG.bdca.fwd.parafac = Factors;			     
