function [xr,xl] = loo(x,fold,leave)
% [xr,xl] = loo(x,fold,leave)

% Author: Mads Dyrholm
if (mod(size(x,2),fold)~=0), warning('Fold does not match size(x,2)');, end
idxr = [];
jump = floor(size(x,2)/fold);
for f=1:fold
  if (f==leave)
    idxl = (f-1)*jump + (1:jump);
  else
    idxr = [idxr ((f-1)*jump+(1:jump))];
  end
end
xr = x(:,idxr);
xl = x(:,idxl);
