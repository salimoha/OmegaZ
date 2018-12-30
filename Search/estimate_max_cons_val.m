function [cons]= estimate_max_cons_val(x,inter_par_p,inter_par_g,y0,m)
% Calculate the estimated value of maximum constraints violations at point x   
cons=(interpolate_val(x,inter_par_p)-y0);
   for ii=1:m
     cons1=interpolate_val(x,inter_par_g{ii});
     cons=max(cons,cons1);
   end
end