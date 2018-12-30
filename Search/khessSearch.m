% kth iteration hessian for the search function
function [ H ] = khessSearch(x,inter_par,xc,R2)
global K
% p=interpolate_val(x,inter_par);
Hp=interpolate_hessian(x,inter_par);
% gp=interpolate_grad(x,inter_par);
% ge=-2*(x-xc);
% e=R2-norm(x-xc)^2;
% Constant K hessian for the search function 
H = Hp + K .* 2*eye(size(x,1));
%H=Hp/e-(gp*ge.'+ge*gp.')/e^2+(p-y0)*(2*ge*ge.'/e^3+2*eye(n)/e^2);
% adaptive K
% H=2*eye(n)/(p-y0)+(gp*ge.'+ge*gp.')/(p-y0)^2-e*(2*gp*gp.'/(p-y0)^3-2*Hp/(p-y0)^2); 
end