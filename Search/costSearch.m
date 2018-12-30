function [M]=costSearch(x,inter_par,xc,R2)
global K
% constant K search function 
M=  interpolate_val(x,inter_par) - K .* (R2-norm(x-xc)^2); 
% adaptive K search function 
% M=-(R2-norm(x-xc)^2)/(interpolate_val(x,inter_par)-y0); 
% if interpolate_val(x,inter_par)<y0
%     M=-inf;
% end
end