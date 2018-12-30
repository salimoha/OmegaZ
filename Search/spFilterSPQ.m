function [x,y,Xtraj,Ytraj,Ctraj]= spFilterSPQ(fun, gradfun, hessfun, consfun, consgrad, x, Ain,bin,m)
% Basic Filter SQP method
% Inputs:
% fun: the objective function. 
% gradfun: gradient of objective function.
% hessfun: Hessian of objective function.
% consfun: constraint functions. This is a cell.
% consgrad: Constraint gradient function. Thi is a cell
% m: number of nonlinear constraints
% Initialize the filter set
% Decreading factor for trust region radius
rho=0.4;
% % IterMAX = 25;
IterMAX = 200;
% keyboard
x0=x; y=fun(x);% assume that first point is feasible 
n = length(x);
filterpoint=x;
filterobj=fun(x);
filterconst=0;
cc = consfun(x);
d=ones(n,1);
for jj=1:m
filterconst=filterconst+max(cc{jj},0);
end
Xtraj=[]; Ytraj =[]; Ctraj=[];
iter=0;
 while 1 & iter <= IterMAX
%  keyboard
iter = iter +1;
filter_size = length(filterobj);
% Objective grad and hessian
g=gradfun(x); H=hessfun(x);
H=modichol(H,0.1,20);
H=(H+H')/2;
% Linear constraint
CSl =[]; Jl=[]; %Hc=[];
CC = consfun(x);
gC = consgrad(x);
CSl = cell2mat(CC)';
for l=1:m
% CSl=vertcat(CSl,lkConstraint(x,inter_Par{l},xc, R2));
Jl=[Jl; gC{l}'];
end
% linearized constrained for sqp
A = vertcat(Jl,Ain,eye(n),-eye(n));
b = vertcat(-1.*CSl,bin - Ain*x,d,d);
% Solve the SQP:
options=optimoptions('quadprog','Display','none');
[p,py,exitflag]=quadprog(H,g,A,b,[],[],[],[],zeros(n,1),options);
%%
if exitflag~=-2
% % %     if norm(p) < 1e-2
    if norm(p) < 1e-5
    break
    end
    
x1=x+p; y1=fun(x1);
c1=0;
Cl= consfun(x1);
for jj=1:m
c1=c1+max(Cl{jj},0);
end
ind=1:filter_size;
indo=ind(filterobj<1.1*y1); indc=ind(filterconst<1.1*c1);
% indo=ind(filterobj<1.01*y1); indc=ind(filterconst<1.01*c1);
     if isempty(intersect(indo,indc))
         indo=ind(filterobj>y1); indc=ind(filterconst>c1);
         indu=intersect(indo,indc);
         filterobj(indu)=[]; filterconst(indu)=[]; filterpoint(:,indu)=[];
         filterobj=[filterobj y1]; filterconst=[ filterconst c1];
         x = x1; y=y1;
         filterpoint = [ filterpoint, x];
%          d=d/rho;
     else
         d=rho*d;
     end
%% x is not feasible 
else
  x = x0 + rho*(x-x0);   
end
% keyboard
Ytraj = vertcat(Ytraj,y);
Xtraj = horzcat(Xtraj,x);
Ctraj = vertcat(Ctraj,max(CSl));
 end