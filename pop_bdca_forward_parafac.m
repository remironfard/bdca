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

OldLoad = {EEG.bdca.cht.u,EEG.bdca.cht.t,EEG.bdca.ica.S};
Factors = parafac(double(EEG.data(:,EEG.bdca.cht.supportframes,:)),...
		  EEG.bdca.cht.R,...
		  [],... %Options
		  zeros(EEG.bdca.cht.R,1),... %const
		  OldLoad,...
		  [0 0 1]'); %,...%FixMode,
	          %Weights);
		  
EEG.bdca.fwd.parafac = Factors;			     
