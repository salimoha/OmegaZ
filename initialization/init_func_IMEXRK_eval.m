% evaluting the initial func and cosnt evalution on the initial points xi
% at the vertices
function [yi, C] = init_func_IMEXRK_eval(xi,Ain, bin)
global acon m ms DX
disp('evaluating the function evalution on the vertices...')
for ii=1:size(xi,2)
     % [xE(:,ii),xT(:,ii)] = mapping(xi(:,ii)');
      [con , FUN] = IMEXRK_Solver3(xi(:,ii)',DX);
       yi(ii)=real(FUN.val);
      if imag(FUN.val) ~= 0 
         warning(strcat('The obj function:   ', ' has imaginary part!!!!!'))
     end
    for jj=1:ms
%         C{jj}(ii)=con{jj}.val(xi(:,ii));
     C{jj}(ii) = real(con{jj});
     if imag(con{jj}) ~= 0 
         warning(strcat('The constraint:   ' , num2str(jj), ' has imaginary part!!!!!'))
     end
    end    
    
    
end

disp('finished the initial function evalutions...')
disp('=============================================')
save init_point xi yi acon C DX  
end