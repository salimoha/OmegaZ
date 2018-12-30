function [x]=min_projection_control_feas1(x0,lb,ub)

% find the minimum distance from a domain:
% consfun: R^n \rightarrow R^m
fun=@(x) fdistance(x,x0);
nonlcon=@ fcigar;
keyboard
options = optimoptions('fmincon','algorithm','sqp','GradObj','on','GradConstr','on','display','iter-detailed');
x = fmincon(fun,x0,[],[],[],[],lb,ub,nonlcon,options);

end

function [y,Dy]=fdistance(x,x0)
  y=norm(x-x0)^2;
  Dy=2*(x-x0);
end

function [y,ye,Dy,Dye]=fcigar(x)
  A=[1 -1 0; 0 1 -1; 1 0 -1];
  b=0.1*ones(3,1);
  y=b-abs(A*x);
  Dy=-diag(sign(A*x+eps))*A;
  ye=[]; Dye=[];
end