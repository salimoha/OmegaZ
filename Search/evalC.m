function [yy,xx,C] = evalC(x1,x2,alpha,CON_NUM)
%  x2 =[.5, .52, .45];
% x1 = [.5, .49999, .49];
%  alpha = 0:0.1:1;
%  CON_NUM =5;
FUN=@(x)IMEXRK_Solver3(x);
for ii =1:numel(alpha)
   xx(ii,:) = x1+alpha(ii)*(x2-x1); 
   [con,ff] = IMEXRK_Solver3( xx(ii,:));
   for j=1:numel(con)
   C{j}(ii) =con{j};
   end
   C{numel(con)+1}(ii) = ff.val;
   yy(ii) = con{CON_NUM};
   disp([num2str(ii/numel(alpha)*100), ' % Completed'])
end


end