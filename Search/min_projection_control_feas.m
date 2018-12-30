 function [x,A1,b1]=min_projection_control_feas(x0,lb,ub,A,b)

% find the minimum distance from a domain:
% consfun: R^n \rightarrow R^m
%
% 
%           minimize     (x-x0)^2  
%           subject to   | A*x |  > b
%
%           minimize        x'x -2x_0'x
%           subject to      -|A*x| <= -b
%  
% 
% make initial point feasible with lb,ub
x=max(x0,lb); x=min(x,ub);

A1=-diag(sign(A*x+eps))*A; b1=-b;
%keyboard
options = optimoptions('quadprog','algorithm','active-set','display','none');
x=quadprog(eye(length(x0)),-x0,A1,b1,[],[],lb,ub,x,options);   % H = I; f = -x0
    



end


% function xf=initial_feaspoint(x,A,b)
% [tt,ind]=sort(x);
% delta=0.1;
% tt(2)=max(delta,tt(2));
% tt(2)=min(1-delta,tt(2));
% tt(1)=tt(2)-delta; tt(3)=tt(2)+delta;
% xf(ind)=tt; xf=xf';
% end

%function A1=fcigar(x,A)
  %A=[1 -1 0; 0 1 -1; 1 0 -1];
  %delta=0.1;
  %b=-delta*ones(3,1);
  %A1=-diag(sign(A*x+eps))*A;
%end