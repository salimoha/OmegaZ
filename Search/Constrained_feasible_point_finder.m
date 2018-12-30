function [x,y]=Constrained_feasible_point_finder(x,inter_par,xc,R2)
% find an initial feasible point for x s.t. g^k-K e^k(x)<0

global n Ain bin eps

rho=0.9;
eps=1e-3;
while 1
y=feasible_cost(x,inter_par,xc,R2);
if y==0
    break
end
% keyboard
g=feasible_grad(x,inter_par,xc,R2);
H=feasible_hessian(x,inter_par,xc,R2);
H=(H+H')/2;
H=modichol(H,0.1,20);

% quadprog
options=optimoptions('quadprog','Display','none');
p=quadprog(H,g,Ain,bin - Ain*x,[],[],[],[],zeros(n,1),options);

% backtracking
a=1;
while 1
    x1=x+a*p;
        y1=feasible_cost(x1,inter_par,xc,R2);
        if (y-y1)>0
        x=x1; y=y1;
        break
    else
        a=a*rho;
        if norm(a*p)<1e-4
            break
        end
    end
end
if norm(a*p)<1e-4
    break
end
end


end

function [M]=feasible_cost(x,inter_par,xc,R2)
global K ms eps
% constant K search function
M=0;

for jj=1:ms
M=  M+0.5*max([0,interpolate_val(x,inter_par{jj}) - K .* (R2-norm(x-xc)^2)+eps])^2;
end

end
function [M]=feasible_grad(x,inter_par,xc,R2)
global K ms eps
% constant K search function 
M=zeros(length(x),1);
for jj=1:ms
cmi=max([0,interpolate_val(x,inter_par{jj}) - K .* (R2-norm(x-xc)^2)+eps]);
g=interpolate_grad(x,inter_par{jj})+2*K*(x-xc);
M=  M+cmi*g;
end

end
function [M]=feasible_hessian(x,inter_par,xc,R2)
global K ms eps
M=zeros(length(x),length(x));
% constant K search function 
for jj=1:ms
cmi=max([0,interpolate_val(x,inter_par{jj}) - K .* (R2-norm(x-xc)^2)+eps]);
g=interpolate_grad(x,inter_par{jj})+2*K*(x-xc);
H=interpolate_hessian(x,inter_par{jj})+2*K*eye(length(x));
if cmi>0
M=  M+cmi*H+g*g';
end
end

end