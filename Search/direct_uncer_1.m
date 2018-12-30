function [e,T,s] = direct_uncer(x,xi,inter_par_p,inter_par_g,tri)
% direct calculation of the uncertainty function (uneffient, just for
% debug)
global y0 n
% Uncertainty function
e=0;
% keyboard
for ind=1:size(tri)
 [xc,R2]=circhyp(xi(:,tri(ind,:)), n);
 e=max([e,R2-norm(x-xc)^2]);
end
% search function
T=max([interpolate_val(x,inter_par_p)-y0,interpolate_val(x,inter_par_g{1})]);
s=-e/T;
end
