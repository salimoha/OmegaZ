% kth iteration gradient for the search function
function [ y ] = kgradSearch(x,inter_par,xc,R2)
global K
%K=Search.constant;
% p=interpolate_val(x,inter_par);
gp=interpolate_grad(x,inter_par);
ge=-2*(x-xc);
% e=R2-(x-xc)'*(x-xc);
% constant K gradient search function 
y = gp - K .* ge;
%y=gp/e-(p-y0)*ge/e^2;
% adaptive K gradient search function 
% y=-ge/(p-y0)+e*gp/(p-y0)^2;
end