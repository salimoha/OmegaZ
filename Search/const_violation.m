function [  ViolCon ] = const_violation( C, ConBound )
% calculates const. violation for IMEXRK problem
global ms
 ViolCon=C;
for kk = 1:numel(C{1})
ViolCon{1}(kk) = C{1}(kk) - ConBound(1);
ViolCon{2}(kk) = C{2}(kk) + ConBound(2);
ViolCon{3}(kk) = C{3}(kk) - ConBound(3);
ViolCon{4}(kk) = C{4}(kk) + ConBound(4);
ViolCon{5}(kk) = C{5}(kk) - ConBound(5); % Eq4  
ViolCon{6}(kk) = C{6}(kk) - ConBound(6); % Eq9
if ms ==7
ViolCon{7}(kk) = C{7}(kk) + ConBound(7); 
end
   
    
end


end

