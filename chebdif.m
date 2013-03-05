function [x,D,W] = chebdif(N)

x = -cos(pi*(0:N)/N)';
c = [2; ones(N-1,1); 2].*(-1).^(0:N)';

X = repmat(x,1,N+1);
dX = X - X';

D = (c*(1./c)')./(dX + eye(N+1));
D = D - diag(sum(D,2));

W = (pi/N)*cos(repmat((1:N+1)',1,N+1).*acos(X'));

end