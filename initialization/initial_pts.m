function [xU,xE,Ain, bin, A, b] = initial_pts(bnd1,bnd2)
global n
% Input the equality constraints
Ain=[eye(n);-eye(n)];
bin=[bnd2 ;-bnd1];
% Calculate the initial points
xU=bounds(bnd1,bnd2, n);

xE=[0.25;0.5;0.75];
delta0=0.15;
for ii=1:n
    e=zeros(n,1); e(ii)=1;
    xE(:,ii+1)=xE(:,1)+delta0*e;
end

% extra constraints for nonconvex constraint 
% nonequality ||x_i - x_j || > 0.1   : i,j=1,2,3.

 A=[1 -1 0; 0 1 -1; 1 0 -1];
 b=0.1*ones(n,1);

end