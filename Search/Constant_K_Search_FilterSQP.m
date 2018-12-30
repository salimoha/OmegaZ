function [x y CS]=Constant_K_Search_FilterSQP(x,inter_par_, inter_Par_,xc_,R2_,Search)
% This funciton finds the minimizer of the search function i a simplex 
%              minimize s^k_i(x) = p^k(x) - K e^k_i(x)
%             subject to CS^k_(l,i) = g^k_l(x) - K e^k_i(x) <= 0 
%                        A x <= b 
% input:
% x: interpolating initial point for search function 
% inter_par: interpolation parameters
% xc: circumcenter for the simplex that x is located at
% R2: square of circumradiaus for the simplex that x is located at
% output:
% x: the minimizer of search function 
% y: minimum of search fucntion at x
% cse: 
% created by: Shahrouz Alimo & Pooriya Beyhaghi
% last modification: Oct/7/2015
%
% inter_Par = inter_par.inter_par_C;
% inter_par = inter_par.inter_par_S;
global n Ain bin ms K Eps inter_par inter_Par xc R2
% inter_par_, inter_Par_,xc_,R2_
inter_par = inter_par_;
inter_Par =inter_Par_;
xc =xc_; R2 = R2_;
% parameters of backtracking
%CS: CASE:
% CS == 1: infeasible domain  
% CS == 0: feasible domain
%
% First find a feasible point
K=Search.constant;
Eps = 1e-2;
cc=0.01;
rho=0.9; 
% keyboard
% Initialize the point in the simplex
%% Find initial feasible point 
% Calculate an initail feasible point
[x,y]=Constrained_feasible_point_finder(x,inter_Par,xc,R2);
% y here is metric of comparison for different senarios
% here is sum norm2 of c_i's
% keyboard
if y>Eps
    CS=1; % infeasible point
else
CS=0; % feasible point
% [x,y]= spFilterSPQ(fun, gradfun, hessfun, consfun, consgrad, x, Ain,bin,ms)
 [x,y]= spFilterSPQ(@fun, @gradfun, @hessfun, @consfun, @consgrad, x, Ain,bin,ms)
end
end
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
%%
%% objective functions
function [f] = fun(x)
global inter_par xc R2 
    f = costSearch(x,inter_par,xc,R2);
end
% gard obj
function [df] = gradfun(x)
 global inter_par xc R2 
    df = kgradSearch(x,inter_par,xc,R2);  
end
% hessian objg
function d2f = hessfun(x)
   global inter_par xc R2 
    d2f = khessSearch(x,inter_par,xc,R2);  
end

%% constraint functions
function [C] = consfun(x)
   global inter_Par xc R2 ms
   for jj = 1:ms
   C{jj} = costSearch(x,inter_Par{jj},xc,R2); % Inequality constraints
   end
end

function [DC] = consgrad(x)
   global inter_Par xc R2 ms
   for jj = 1:ms
   DC{jj} = kgradSearch(x,inter_Par{jj},xc,R2); % Inequality constraints
   end
end


%%
% 
% function [y]= meritfunction(x,rho,inter_par,inter_Par,xc,R2)
% global ms
% %keyboard
% rho = rho / (1+rho);
% y=(1-rho)*costSearch(x,inter_par,xc,R2);
% for jj=1:ms
%     y=y+(rho)*max(0,costSearch(x,inter_Par{jj},xc,R2));
% end
% end
% 
% 
% function [M]=feasible_cost(x,inter_par,xc,R2)
% global K ms
% % keyboard
% % constant K search function
% M=0;
% for jj=1:ms
% M=  M+0.5*max([0,interpolate_val(x,inter_par{jj}) - K .* (R2-norm(x-xc)^2)])^2;
% end
% 
% end
% function [M]=feasible_grad(x,inter_par,xc,R2)
% global K ms
% % constant K search function
% M=zeros(length(x),1);
% for jj=1:ms
% cmi=max([0,interpolate_val(x,inter_par{jj}) - K .* (R2-norm(x-xc)^2)]);
% g=interpolate_grad(x,inter_par{jj})+2*K*(x-xc);
% M=  M+cmi*g;
% end
% 
% end
% function [M]=feasible_hessian(x,inter_par,xc,R2)
% global K ms
% % constant K search function
% M=zeros(length(x),length(x));
% % constant K search function 
% for jj=1:ms
% cmi=max([0,interpolate_val(x,inter_par{jj}) - K .* (R2-norm(x-xc)^2)]);
% g=interpolate_grad(x,inter_par{jj})+2*K*(x-xc);
% H=interpolate_hessian(x,inter_par{jj})+2*K*eye(length(x));
% if cmi>0
% M=  M+cmi*H+g*g';
% end
% end
% 
% end
% % % search function 
% function [M]=costSearch(x,inter_par,xc,R2)
% global K
% constant K search function 
% M=  interpolate_val(x,inter_par) - K .* (R2-norm(x-xc)^2); 
% adaptive K search function 
% M=-(R2-norm(x-xc)^2)/(interpolate_val(x,inter_par)-y0); 
% if interpolate_val(x,inter_par)<y0
%     M=-inf;
% end
% end
% % kth iteration gradient for the search function
% function [ y ] = kgradSearch(x,inter_par,xc,R2)
% global  K
% % p=interpolate_val(x,inter_par);
% gp=interpolate_grad(x,inter_par);
% ge=-2*(x-xc);
% % e=R2-(x-xc)'*(x-xc);
% % constant K gradient search function 
% y = gp - K .* ge;
% %y=gp/e-(p-y0)*ge/e^2;
% % adaptive K gradient search function 
% % y=-ge/(p-y0)+e*gp/(p-y0)^2;
% end
% % kth iteration hessian for the search function
% function [ H ] = khessSearch(x,inter_par,xc,R2)
% global K n
% % p=interpolate_val(x,inter_par);
% Hp=interpolate_hessian(x,inter_par);
% % gp=interpolate_grad(x,inter_par);
% % ge=-2*(x-xc);
% % e=R2-norm(x-xc)^2;
% % Constant K hessian for the search function 
% H = Hp + K .* 2*eye(size(x,1));
% %H=Hp/e-(gp*ge.'+ge*gp.')/e^2+(p-y0)*(2*ge*ge.'/e^3+2*eye(n)/e^2);
% % adaptive K
% % H=2*eye(n)/(p-y0)+(gp*ge.'+ge*gp.')/(p-y0)^2-e*(2*gp*gp.'/(p-y0)^3-2*Hp/(p-y0)^2); 
% end
% %% Constrained functions
% function [Mc]=lkConstraint(x,inter_Par,xc,R2)
% global K m
% % constant K constrained function 
% Mc = [];
% for l = 1:m
% Mc = vertcat(Mc, [xm ym indm]= tringulation_search_constantK_constraints(inter_par_p,inter_par_g,xi,tri); - K .* (R2-norm(x-xc)^2)); 
% end 
% end
% % kth iteration gradient for the search function
% function [ yc ] = lkGradConstraint(x,inter_Par,xc,R2)
% global  K m
% yc =[];
% for l=1:m
% gp=interpolate_grad(x,inter_Par{l});
% ge=-2*(x-xc);
% % e=R2-(x-xc)'*(x-xc);
% % constant K gradient constrained function 
% yc = vertcat(yc, gp - K .* ge);
% end
% end
% % kth iteration hessian for the search function
% function [ Hc ] = lkHessConstraint(x,inter_Par,xc,R2)
% global K m n
% Hc =[];
% for l=1:m
% Hp=interpolate_hessian(x,inter_Par{l});
% % Constant K hessian for the constrained function 
% Hc = vertcat(Hc, Hp + K .* 2*eye(n));
% end
% end