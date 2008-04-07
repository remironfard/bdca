function M = kr(X,Y)
%

[I,K]=size(X);

M=zeros(I*size(Y,1),K);
for k=1:K
  M(:,k)=vec(Y(:,k)*X(:,k)');
end

function vv = vec(xx)
%
vv = xx(:);
