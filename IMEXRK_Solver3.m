 function [con , FUN, CS , CON, SolC, Be, Bi ] = IMEXRK_Solver3(SolC,Eps)
% function [con , FUN, CS , CON, SolC ] = IMEXRK_Solver3(SolC,Eps)
% IMAXRK 3rd order accuracy coefficient 
% The challenging IMEXRK problems we are looking at compel
% us to think closely about how we do nonconvex optimization with
% constraints in an efficient manner.
% Author:
% Shahrouz Alimo & Daniele Cavaglieri
% 03/06/2015
% output
% con: constraints
% FUN: obj. function val
%CS: constraitn violation. if 1: infeasible else 0 is feasible point for i
%    constraint related to con{i}
global FUN_MAX CON_MAX 
%% smoothing function
conf = readConf('./conf/IMEXRK.conf',0);
smooth_fun = conf.saturation_function;
% if smooth_fun == 'atanh' 
% fun2=@(y,ymax) (tanh(y./ymax).*ymax);
% S12_34= 10; % scaling C12 wrt C34 is 10
S12_34= 1;
% fun2=@(y,ymax) (tanh(y).*ymax);
% else
fun2=@(x,xmax) (  (1+sign(x))*log(1+x)  + (sign(x)-1)*log(1-x)  )./2;
% end
%   fun2=@(y,ymax)(y/(1+y/ymax))*sign(y); 
%        con{ii} = min(fun1(yy(ii), CON_MAX), abs(con{ii})*sign(real(con{ii})));
%        con{ii} = max(-fun1(con{ii}, CON_MAX), con{ii});   
%%
if nargin <2
Eps = 5e-3;
end
if nargin <1
%     SolC = [1/2, 9/10,7/10];
    SolC = [ 14/25; 4/5; 7/10]'; % x_star
%     %% change of corrdinate for c2,c3,c4
% SolC = mapping(SolC)
elseif length(SolC) > 1
%% change of corrdinate for c2,c3,c4
% SolC = mapping(SolC);    
% if size(SolC,1) > 1
%     % mathematic vector to programing vector
%     SolC = SolC';
% end
% %     if SolC(1) == SolC(2) || SolC(2) == SolC(3) || SolC(1) == 0 || SolC(2) == 0 
%     if SolC(1) == SolC(2) || SolC(1) == 0 
%         if SolC(1) == 1 
%         SolC(1) = SolC(1) - Eps;
%         else SolC(1) = SolC(1) + Eps; end
%         if SolC(2) == 1 
%         SolC(2) = SolC(2) - 2*Eps; end
%       
%     end
%     
%      if SolC(2) == 1 
%         SolC(2) = SolC(2) - Eps;
%      end
%     if SolC(3) == SolC(2) || SolC(2) == 0 
%          if SolC(2) == 1 
%         SolC(2) = SolC(2) - Eps;
%          else SolC(2) = SolC(2) + Eps*1.5; end
%         if SolC(3) == 1 
%         SolC(3) = SolC(3) - Eps; end
%        
%    end   
%   
%      if SolC(3) == SolC(2) & SolC(1) == SolC(2) 
%       SolC(3) = SolC(2) + Eps*1.5; 
%       SolC(2) = SolC(1) + Eps*3;
%       SolC(1) = SolC(1) + Eps;
%      end
   
SolC =  Eps + (1-2*Eps)*SolC;  
end
%% Define variables
syms bI1 bI2 bI3 bI4 bI5
syms bE1 bE2 bE3 bE4
syms     c2  c3  c4

bi = [bI1, bI2, bI3, bI4, bI5];
be = [bE1, bE2, bE3, bE4, 0];

Ai = [0, 0, 0, 0, 0; bI1, c2-bI1, 0, 0, 0;
      bI1, bI2, c3-bI1-bI2, 0, 0; bI1, bI2, bI3, c4-bI1-bI2-bI3, 0;
      bI1, bI2, bI3, bI4, bI5];

Ae = [0, 0, 0, 0, 0; c2, 0, 0, 0, 0; bE1, c3-bE1, 0, 0, 0; 
      bE1, bE2, c4-bE1-bE2, 0, 0; bE1, bE2, bE3, bE4, 0];

c = diag([0 c2 c3 c4 1]);

e = ones(5,1);


%% Define constraints

t11i = bi*e - 1;
t11e = be*e - 1; % Eq1

t21i = bi*c*e - 1/2;
t21e = be*c*e - 1/2; % Eq2

% t31i = bi*c*c*e/2 - 1/6;
t31e = be*c*c*e/2 - 1/6; % Eq3

t32ee = be*Ae*c*e - 1/6; % Eq4 -> 2nd order
t32ei = be*Ai*c*e - 1/6; % Eq
t32ie = bi*Ae*c*e - 1/6;
t32ii = bi*Ai*c*e - 1/6; % Eq9 --> bI5


%% Define objective functions

t44eee = be*Ae*Ae*c*e - 1/24;
Lstab = - (bI1 * (bI1 + bI2 - c2) * (bI1 + bI2 + bI3 - c3) * (bI1 + bI2 + bI3 + bI4 - c4)) / ...
        (bI5 * (bI1 - c2) * (bI1 + bI2 - c3) * (bI1 + bI2 + bI3 - c4));

t42ei = be*c*Ai*c*e - 3/24;
t42ee = be*c*Ae*c*e - 3/24;
t43ie = bi*Ae*c*c*e/2 - 1/24;
t43ee = be*Ae*c*c*e/2 - 1/24;
t44iii = bi*Ai*Ai*c*e - 1/24;
t44iie = bi*Ai*Ae*c*e - 1/24;
t44iei = bi*Ae*Ai*c*e - 1/24;
t44iee = bi*Ae*Ae*c*e - 1/24;
t44eii = be*Ai*Ai*c*e - 1/24;
t44eie = be*Ai*Ae*c*e - 1/24;
t44eei = be*Ae*Ai*c*e - 1/24;

tau4 = sqrt(t42ei^2 + t42ee^2 + t43ie^2 + t43ee^2 + t44iii^2 + t44iie^2 + ...
       t44iei^2 + t44iee^2 + t44eii^2 + t44eie^2 + t44eei^2 + t44eee^2);

%% Solve nonlinear system



t11eC = subs(t11e, [c2, c3, c4], SolC);
t21eC = subs(t21e, [c2, c3, c4], SolC);
t31eC = subs(t31e, [c2, c3, c4], SolC);
SolE1 = solve(t11eC, t21eC, t31eC, bE1, bE2, bE3);

t32eeEC = subs(t32ee, [c2, c3, c4], SolC);
t32eeEC = subs(t32eeEC, [bE1, bE2, bE3], [SolE1.bE1, SolE1.bE2, SolE1.bE3]);
SolE2 = solve(t32eeEC, bE4);

t11iC = subs(t11i, [c2, c3, c4], SolC);
t21iC = subs(t21i, [c2, c3, c4], SolC);
% t31iC = subs(t31i, [c2, c3, c4], SolC);
t32ieEC = subs(t32ie, [c2, c3, c4], SolC);
t32ieEC = subs(t32ieEC, [bE1, bE2, bE3], [SolE1.bE1, SolE1.bE2, SolE1.bE3]);
t32eiEC = subs(t32ei, [c2, c3, c4], SolC);
t32eiEC = subs(t32eiEC, [bE1, bE2, bE3], [SolE1.bE1, SolE1.bE2, SolE1.bE3]);
SolI1 = solve(t11iC, t21iC, t32ieEC, t32eiEC, bI1, bI2, bI3, bI4);

t32iiIEC = subs(t32ii, [c2, c3, c4], SolC);
t32iiIEC = subs(t32iiIEC, [bI1, bI2, bI3, bI4], [SolI1.bI1, SolI1.bI2, SolI1.bI3, SolI1.bI4]);
SolI2 = solve(t32iiIEC, bI5);


%% Calculate radicands

SolE1rad = solve(t11e, t21e, t31e, bE1, bE2, bE3);
t32eeE = subs(t32ee, [bE1, bE2, bE3], [SolE1rad.bE1, SolE1rad.bE2, SolE1rad.bE3]);
[cE, tE] = coeffs(t32eeE, bE4);
radE = simplify(cE(2) ^ 2 - 4 * cE(1) * cE(3)); % Radicand formula b^2-4*a*c

t32ieE = subs(t32ie, [bE1, bE2, bE3], [SolE1rad.bE1, SolE1rad.bE2, SolE1rad.bE3]);
t32eiE = subs(t32ei, [bE1, bE2, bE3], [SolE1rad.bE1, SolE1rad.bE2, SolE1rad.bE3]);
SolI1rad = solve(t11i, t21i, t32ieE, t32eiE, bI1, bI2, bI3, bI4);
t32iiIErad = subs(t32ii, [bI1, bI2, bI3, bI4], [SolI1rad.bI1, SolI1rad.bI2, SolI1rad.bI3, SolI1rad.bI4]);
[cI, tI] = coeffs(t32iiIErad, bI5);
radI = cI(2) ^ 2 - 4 * cI(1) * cI(3); % Radicand formula b^2-4*a*c


%% Evaluate objective functions
nI2 = 1;
nE2 = 2;

t44eeeSOL = subs(t44eee, [c2, c3, c4], SolC);
t44eeeSOL = subs(t44eeeSOL, [bE1, bE2, bE3], [SolE1.bE1, SolE1.bE2, SolE1.bE3]);
if length(SolE2) ==1
t44eeeSOL = subs(t44eeeSOL, bE4, SolE2);    
else
t44eeeSOL = subs(t44eeeSOL, bE4, SolE2(nE2));
end
vpa(t44eeeSOL); % Must satisfy: -6/1000 < t44eeeSOL < 0

LstabSOL = subs(Lstab, [c2, c3, c4], SolC);
LstabSOL = subs(LstabSOL, [bI1, bI2, bI3, bI4], [SolI1.bI1, SolI1.bI2, SolI1.bI3, SolI1.bI4]);
LstabSOL = subs(LstabSOL, bI5, SolI2(nI2)); % TODO change the indexes
if length(SolE2) ==1
LstabSOL = subs(LstabSOL, bE4, SolE2);    
else
LstabSOL = subs(LstabSOL, bE4, SolE2(nE2));
end
% vpa(LstabSOL); % Must satisfy: -1/10 < LstabSOL < 1/10

tau4SOL = subs(tau4, [c2, c3, c4], SolC);
tau4SOL = subs(tau4SOL, [bE1, bE2, bE3], [SolE1.bE1, SolE1.bE2, SolE1.bE3]);
tau4SOL = subs(tau4SOL, [bI1, bI2, bI3, bI4], [SolI1.bI1, SolI1.bI2, SolI1.bI3, SolI1.bI4]);
tau4SOL = subs(tau4SOL, bI5, SolI2(nI2));
if length(SolE2) ==1
tau4SOL = subs(tau4SOL, bE4, SolE2);    
else
tau4SOL = subs(tau4SOL, bE4, SolE2(nE2));
end
vpa(tau4SOL); % Must satisfy: tau4SOL < 1/10


%% eq 4 and 9 quadratic ones
% radI 
radISOL = subs(radI, [c2, c3, c4], SolC);
% radISOL = subs(radISOL , [bE1, bE2, bE3], [SolE1.bE1, SolE1.bE2, SolE1.bE3]);
if length(SolE2) ==1
radISOL = subs(radISOL, bE4, SolE2);    
else
radISOL  = subs(radISOL , bE4, SolE2(nE2)); % Must satisfy: radESOL <= 0
end
% radE
radESOL = subs(radE, [c2, c3, c4], SolC); % Must satisfy: radESOL <= 0


%% generate constraints functions
% keyboard
CON{1}.field = strcat('LstabSOL');
CON{1}.val = double(vpa(LstabSOL));
CON{1}.lowb = 0;
CON{1}.upb = 0;
CON{2}.field = strcat('t44eeeSOL');
CON{2}.val = double(vpa(t44eeeSOL));
CON{2}.lowb = 0;
CON{2}.upb = 0; % % CON{2}.upb = -0.00001;

con{1} = (CON{1}.val - CON{1}.upb)*S12_34;
con{2} = (-CON{1}.val + CON{1}.lowb)*S12_34;
con{3} = CON{2}.val - CON{2}.upb; % con{3} = (CON{2}.val - CON{2}.upb)./100*1e-40;
con{4} = -CON{2}.val + CON{2}.lowb;
con{5} = -double(vpa(radESOL )) ; % Eq4  
con{6} = -double(vpa(radISOL )) ; % Eq9


con{7} = -abs(SolC(2) - SolC(3)); % > 1/10 % con{7}  = double(vpa(tau4SOL)) - 1/10 ;

% handing imaginary and large constraint violation
[con, CS] = constraint_smoothing(con,fun2, S12_34);
% generate obj. functions
[FUN] = cost_smoothing(tau4SOL,fun2);




%% beta coeffiecients
% % disp(strcat('real Function value is:  ' , num2str(real(FUN.val)), ' at xi = [ ' ,num2str(real(SolC)) , ' ] ' ))
Be=zeros(4,1);
Bi = zeros(5,1);
if nargout > 4
    if length(SolE2) ==2
 Be(1) = eval(subs(SolE1.bE1 ,bE4 , SolE2(nE2)));
 Be(2) = eval(subs(SolE1.bE2 ,bE4 , SolE2(nE2)));
 Be(3) = eval(subs(SolE1.bE3 ,bE4 , SolE2(nE2)));
 Be(4) = eval( SolE2(nE2));
 Bi(1) =  double(vpa(subs(subs(SolI1.bI1 ,bI5 , SolI2(nI2)) , bE4, SolE2(nE2))));
 Bi(2) =  double(vpa(subs(subs(SolI1.bI2 ,bI5 , SolI2(nI2)) , bE4, SolE2(nE2))));
 Bi(3) =  double(vpa(subs(subs(SolI1.bI3 ,bI5 , SolI2(nI2)) , bE4, SolE2(nE2))));
 Bi(4) =  double(vpa(subs(subs(SolI1.bI4 ,bI5 , SolI2(nI2)) , bE4, SolE2(nE2))));
 Bi(5) =  eval(subs( SolI2(nI2) , bE4, SolE2(nE2)));
    else
        nE2=1;
        Be(1) = eval(subs(SolE1.bE1 ,bE4 , SolE2(nE2)));
        Be(2) = eval(subs(SolE1.bE2 ,bE4 , SolE2(nE2)));
        Be(3) = eval(subs(SolE1.bE3 ,bE4 , SolE2(nE2)));
        Be(4) = eval( SolE2(nE2)); 
         Bi(1) =  double(vpa(subs(subs(SolI1.bI1 ,bI5 , SolI2(nI2)) , bE4, SolE2(nE2))));
         Bi(2) =  double(vpa(subs(subs(SolI1.bI2 ,bI5 , SolI2(nI2)) , bE4, SolE2(nE2))));
         Bi(3) =  double(vpa(subs(subs(SolI1.bI3 ,bI5 , SolI2(nI2)) , bE4, SolE2(nE2))));
         Bi(4) =  double(vpa(subs(subs(SolI1.bI4 ,bI5 , SolI2(nI2)) , bE4, SolE2(nE2))));
         Bi(5) =  eval(subs( SolI2(nI2) , bE4, SolE2(nE2))); 
    end 
end
end
 
%% saturation functions for cost and constrians
function [FUN] = cost_smoothing(tau4SOL,fun2)
global FUN_MAX 
FUN.field = strcat('func' , '_tau4SOL');
%  FUN.val  = double(vpa(tau4SOL))' .* double(vpa(tau4SOL));
 FUN.val  = double(vpa(tau4SOL));
 if imag(FUN.val) ~= 0
%  FUN.val  = abs(double(vpa(tau4SOL))'.*double(vpa(tau4SOL))) ;
      FUN.val  = abs(double(vpa(tau4SOL))') ;
 end
ff =abs(FUN.val)*sign(real(FUN.val));
FUN.val = fun2(ff,FUN_MAX);
FUN.upb = 1/10;           
end
%  
function [con_s, CS] = constraint_smoothing(con,fun2, S12_34)   
global CON_MAX
con_s = con;
CS = zeros(1,length(con));
for ii =1:numel(con)
yy(ii) = abs(con{ii})*sign(real(con{ii}));
  if con{ii} > 0 || imag(con{ii}) ~= 0
       CS(ii) = 1;
  end
% 
  if ii <=2
         con_s{ii} = fun2(yy(ii), CON_MAX*S12_34);        
%          con_s{ii} = fun2(yy(ii), CON_MAX*10);
% 
%   elseif ii ==4 || ii==3
%      con_s{ii} = fun2(yy(ii), CON_MAX);
% 
  else      
%       con_s{ii} =  yy(ii);
     con_s{ii} = fun2(yy(ii), CON_MAX);
  end
end

end
