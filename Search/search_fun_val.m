function [y,indm]=search_fun_val(x, inter_par,inter_Par,xi,tri)
% calculate the adaptive search function at x in the orginal format.
%  y: adaptive search function = max { f_i} i = 1, ...,m+1
% and
% indm: the index of deluany triangualtion simplex assocoieted with point x
% ""warning!!: the search function we are using is
% -e/search_fun_val(x,inter_par,inter_Par,xi,tri)!!""

y = MaxSearchFunc(x,inter_par,inter_Par);  %calculates max {p-y0/e , g_i/e} , i=1,...m
[e,indm]=Uncertainty_quad(x,xi,tri);
% calculates max {p-y0/e , g_i/e} , i=1,...m
y=y/e;

end

function [e,indm] = Uncertainty_quad(x,xi,tri)
% calculates the uncertainty function at x and the index of triangulation which x is located at, 
% w.r.t. the existing vertices (x_i)_ and
% the delauny triangluations (tri)

% output:
% indm: tells us the which simplex the point x located at.
% e: e(x) at that simplex
% e^k(x) = max { e^k_j(x) }, j: different simplices
e=0;
for ii=1:size(tri,1)
    [xc,R2]=circhyp(xi(:,tri(ii,:)), length(x));
    if R2-norm(x-xc)^2>e
        e=R2-norm(x-xc)^2;
        indm=ii;
    end
end
end

function CS = MaxSearchFunc(x,inter_par,inter_Par)
% Search function 1st format
% calculates max {p-y0/e , g_i/e} , i=1,...m
global ms y0 
 CS = interpolate_val(x,inter_par)-y0;
for jj=1:ms
CS=max([CS, interpolate_val(x,inter_Par{jj})]);
end
    
end