function [x M CS]=Constraint_Adaptive_K_Search(x,inter_par,inter_Par, xc,R2, Search)
% This funciton finds the minimizer of the search function in a simplex
% Using adaptive K algorithm
% input:
% x: interpolating initial point for search function 
% inter_par: interpolation parameters
% xc: circumcenter for the simplex that x is located at
% R2: square of circumradiaus for the simplex that x is located at
% output:
% x: the minimizer of search function 
% y: minimum of search fucntion at x
% CS: case senarios for feasibility of points
% if CS==1: % there is NO point that pk < y0 && gk < 0
% if CS==0:  there is a point that pk < y0 && gk < 0
% created by: Shahrouz Alimo & Pooriya Beyhaghi
% last modification: Jan/1/2016
%
% keyboard
global n Ain bin ms y0 xi tri
y0 = Search.constant;
% parameters of backtracking
gamma=0.9; 
% Initialize the point in the simplex
%[xc,R2]=circhyp(xiT(:,tri(index,:)), n); 

% Calculate the Newton direction   
% search function 
% y=cost(x,inter_par,xc, R2);
CS=1; % there is no point that pk < y0 && gk < 0
rho=1;
for iterrr=1:20 %for iterrr=1:25
iter=1;
% solve for given rho
x_pre=x;
% while iter<250
% keyboard
while iter<50
% keyboard
% Calculate the exponential fitting search function and its derivatives 
M = MaxSearchFunc(x,inter_par,inter_Par,xc,R2);
if M == -inf 
    CS=0; 
    break
end
% % % What is it for???
% rho=max(rho,-3*M/(R2-norm(x-xc)^2));
% rho=min(rho,(R2-norm(x-xc)^2)/M);
% rho = min(rho,1/M);
 rho = min(rho,100/(M*(R2-norm(x-xc)^2)));
%  keyboard
% disp(strcat(' iterrr is ...' , num2str(iterrr), '  and  iter is ... ' , num2str(iter) )  )
% if iterrr ==15
%      keyboard
% end
M = costSearch_A(x,inter_par,inter_Par,xc,R2,rho);
DM= kgradSearch_A(x,inter_par,inter_Par,xc,R2,rho);
D2M= khessianSearch_A(x,inter_par,inter_Par,xc,R2, rho);
D2M=modichol(D2M,0.1,20);
D2M=(D2M+D2M')/2;

% inital value for search resulted in p-y0/e
%% new added
% if abs(M)==inf
%     Search.method = 1;
%     Search.constant = 0;
%     [x M]=Constant_K_Search_FilterSQP(x,inter_par, inter_Par,xc,R2,Search);
%     iter=iter+1;
% %      continue
%   CS =0;
% return
% end
%%
% quadprog
options=optimoptions('quadprog','Display','none');
p=quadprog(D2M,DM,Ain ,bin - Ain*x,[],[],[],[],zeros(n,1),options);
% backtracking
a=1;
while 1
    x1=x+a*p;
        M1=costSearch_A(x1,inter_par,inter_Par,xc,R2, rho);
        if (M-M1)>0
        x=x1; M=M1;
             break
        else
        a=a*gamma;
           if norm(a*p)<1e-4
             break
           end
        end
end
% keyboard
%%%%%%
if norm(a*p)<1e-2
    break
end

iter=iter+1;
end
%%%%%
% impro=norm(x_pre-x);
% keyboard
%EXPFIT_RES = log(ms+1)*errEXP / search_fun_val_AdaptiveK2(x, inter_par,inter_Par,xi,tri);
% % if impro<1e-3 
% %     break
% % end
rho=rho*10;
% rho=rho*2;
end

    % inital value for search resulted in p-y0/e
if M==-inf
    Search.method = 1;
    Search.constant = 0;
    [x M]=Constant_K_Search_FilterSQP(x,inter_par, inter_Par,xc,R2,Search);
end

% keyboard
end


% search function 
function [M]=costSearch_A(x,inter_par,inter_Par,xc,R2,rho)
% global y0
% constant K search function 
% M=  interpolate_val(x,inter_par) - K .* (R2-norm(x-xc)^2); 
% adaptive K search function
e=R2-norm(x-xc)^2;
p=costSmoothing_A(x,inter_par,inter_Par,xc,R2,rho);
% M=-e/p; 
M = p/e;
if e < 0
    M = inf;
elseif p < 0 
    M = -inf;
end
% if interpolate_val(x,inter_par)<y0
%      M=-inf;
% end
end
% kth iteration gradient for the search function
function [ y ] = kgradSearch_A(x,inter_par,inter_Par,xc,R2,rho)
% global  y0
p=costSmoothing_A(x,inter_par,inter_Par,xc,R2,rho);
gp=gradSmoothing_A(x,inter_par,inter_Par,xc,R2,rho);
ge=-2*(x-xc);
e=R2-norm(x-xc)^2;
y=gp/e-(p)*ge/e^2;
end
% kth iteration hessian for the search function
function [ H ] = khessianSearch_A(x,inter_par,inter_Par,xc,R2,rho)
global n 
p=costSmoothing_A(x,inter_par,inter_Par,xc,R2,rho);
gp=gradSmoothing_A(x,inter_par,inter_Par,xc,R2,rho);
Hp=hessianSmoothing_A(x,inter_par,inter_Par,xc,R2,rho);
ge=-2*(x-xc);
e=R2-norm(x-xc)^2;
% Constant K hessian for the search function 
% H = Hp + K .* 2*eye(size(x,1));
H=Hp/e-(gp*ge.'+ge*gp.')/e^2+(p)*(2*ge*ge.'/e^3+2*eye(n)/e^2);
% adaptive K
% H=2*eye(n)/p+(gp*ge.'+ge*gp.')/p^2-e*(2*gp*gp.'/p^3-2*Hp/p^2); 
end

% search function smoothing
function [M]=costSmoothing_A(x,inter_par,inter_Par,xc,R2,rho)
global y0 ms FUN_MAX
e=R2-norm(x-xc)^2;
CS = MaxSearchFunc(x,inter_par,inter_Par,xc,R2)*e;
if CS ==-inf
    M = CS;
    return
end
% keyboard
% adaptive K search function
objFun_inter = interpolate_val(x,inter_par);
% if objFun_inter > FUN_MAX
%     objFun_inter = FUN_MAX;
% end
M =  exp(rho*(objFun_inter-y0));

for jj=1:ms
   M = M + exp(rho*interpolate_val(x,inter_Par{jj}));
end
M=log(M)/rho;
% M=log(M/rho;

end
% kth iteration gradient for the search function
function [ DM ] = gradSmoothing_A(x,inter_par,inter_Par,xc,R2,rho)
global  y0 ms FUN_MAX
% adaptive K gradient search function smoothing
% objFun_inter = interpolate_val(x,inter_par);
% if objFun_inter > FUN_MAX
%     objFun_inter = FUN_MAX;
% end
[M]=costSmoothing_A(x,inter_par,inter_Par,xc,R2,rho);

DM = exp(rho*(interpolate_val(x,inter_par)-y0)) / exp(rho*M) * interpolate_grad(x,inter_par);
% DM = exp(rho*( objFun_inter -y0)) / exp(rho*M) * interpolate_grad(x,inter_par);
for jj = 1:ms
DM = DM + exp(rho*interpolate_val(x,inter_Par{jj})) / exp(rho*M) * interpolate_grad(x,inter_Par{jj});
end


end
% kth iteration hessian for the search function
function [ D2M ] = hessianSmoothing_A(x,inter_par,inter_Par,xc,R2,rho)
global y0 ms  FUN_MAX
% objFun_inter = interpolate_val(x,inter_par);
% if objFun_inter > FUN_MAX
%     objFun_inter = FUN_MAX;
% end
M =costSmoothing_A(x,inter_par,inter_Par,xc,R2, rho);
DM = gradSmoothing_A(x,inter_par,inter_Par,xc,R2, rho);

D2M = -rho * DM*DM';
% D2M = D2M + exp(rho*(objFun_inter-y0)) / exp(rho*M) * ...
%     (interpolate_hessian(x,inter_par) + rho*interpolate_grad(x,inter_par) * interpolate_grad(x,inter_par)') ;

D2M = D2M + exp(rho*(interpolate_val(x,inter_par)-y0)) / exp(rho*M) * ...
    (interpolate_hessian(x,inter_par) + rho*interpolate_grad(x,inter_par) * interpolate_grad(x,inter_par)') ;

for jj = 1:ms
D2M = D2M + exp(rho*interpolate_val(x,inter_Par{jj})) / exp(rho*M) * ...
    (interpolate_hessian(x,inter_Par{jj}) + rho* interpolate_grad(x,inter_Par{jj}) * interpolate_grad(x,inter_Par{jj})' );
end


end


function CS = MaxSearchFunc(x,inter_par,inter_Par,xc,R2)
global ms y0
% inital value for search resulted in p-y0/e
e = R2 - (norm(x-xc))^2; 
% CS_1 = -e/costSearch_A(x,inter_par,xc,R2,y0);
 CS = interpolate_val(x,inter_par)-y0;

for jj=1:ms
CS=max([CS, interpolate_val(x,inter_Par{jj})]);
end
if CS < 0 
    CS = -inf;
else
       CS=CS/e;
%     CS=-e/CS;
end
    
end