% function [k] = consraints_delta_dogs(Var, Method, fun, con,lob,upb)
% The Delta_Dogs algorithm with nonlinear costraints.
% Simple bound costraints for the control parameter is imposed.
% bnd1<x<bnd2
% fun: The objective function
% con{i}= The costraint function
% n:dimension of the problem, ms: number of state constraints, m: number of
% control constraints
% Search.Method=1 constant_K
% Search.Method=2 Adaptive_K
clear all; close all ; clc
% run('../addTopath.m')
% run('C:\Users\shahrouz\Documents\Bewley_lab\Delta-dogs-nonconvex\addTopath.m')
% parameters
x_star = [ 14/25; 4/5; 7/10]; % x_star = [0.5600; 0.8000; 0.7000];
% DX = 5e-3; % DX = 1e-1;
global n m ms bnd1 bnd2  Search r ConBound RES FUN_MAX CON_MAX DX MESH_SIZE iter_max
opt_init_IMEXRK()
% MESH_SIZE=50; % grid size% % Mss=20;
%quantization=@(x,bnd1,bnd2,Nm) bnd1+(bnd2 - bnd1).*round((x-bnd1)./(bnd2 - bnd1).*Nm)/Nm;
% interpolaion strategy
global  Ain bin
inter_method=1;
% Input the equality constraints
% Calculate the function evaluation at initial points
[xU,xE, Ain, bin, A_kir, b_kir] = initial_pts(bnd1,bnd2);
% calculates acon , yi, C
% DX = 0.1;
% ms=6;
[yE, C] = init_func_IMEXRK_eval(xE,Ain, bin);
% xi(:,end+1) = round(x_mid*MESH_SIZE)./MESH_SIZE;
% [ym, Cm] = init_func_IMEXRK_eval(xi(:,end) ,Ain, bin);
% yi(end+1) = ym;
% for i=1:ms
% C{i}(end+1) = Cm{i};
% end
%%
 iter_max = 30;
 y0=Search.constant;
for k=1:iter_max
    %             save surr_pts_2 yi xi C tri acon
    %     disp([num2str(k/iter_max*100), ' % Completed'])
    vm1 = min(max([C{1};C{2};C{3};C{4};C{5};C{6};yE-Search.constant]));
    disp( [' started iteration ' ,num2str(k), ' ... '] )
    if mod(k,25) ==0
        time = datestr(now,'yyyymmddHHMMSS');
        save(strcat('./run_new/tanh7_iter_', num2str(k,'%02d'),'.mat'));
    end
    
    % Modify the interpolations
    inter_par_p= interpolateparametarization(xE,yE,inter_method);
    for jj=1:ms
        inter_par_g{jj}= interpolateparametarization(xE,C{jj},inter_method);
    end
    
    % Calculate Discrete search function
    yup=zeros(1,size(xU,2));
    for ii=1:size(xU,2)
         %yup(ii)=estimate_max_cons_val(xU(:,ii),inter_par_p,inter_par_g,y0,ms)/mindis(xU(:,ii),xE);
         yup(ii)=inf;
    end
    % Perform the search
    while 1
    [xm ym(k) Cs(k)]= tringulation_search_constraints(inter_par_p,inter_par_g,[xE xU],A_kir,b_kir);
    %keyboard
    xm=round((xm-bnd1)./(bnd2-bnd1).*MESH_SIZE)./MESH_SIZE.*(bnd2-bnd1)+bnd1;
    [xm,xE,xU,newadd,success]=points_neighbers_find(xm,xE,xU);
    if success==1
         break
    else
    yup=[yup estimate_max_cons_val(xm,inter_par_p,inter_par_g,y0,ms)/mindis(xm,xE)];
    end
    end
    
    if mindis(xm,[xE xU])<1e-6
        break
    end
    
    if (estimate_max_cons_val(xm,inter_par_p,inter_par_g,y0,ms)/mindis(xm,xE)>min(yup) || mindis(xm,xU)<1e-6)
       [t,ind]=min(yup);
       xm=xU(:,ind); xU(:,ind)=[];
    end
    % feasible constraint projection
    if mindis(xm,[xE xU])<1e-6
        break
    end 
    [con , FUN, CS , CON, Xm] = IMEXRK_Solver3(xm',DX);
    %% Stoping Criteria
    %impro_p=mindis(xmp,xE);
    %impro_q=mindis(xmq,xE);
    %
    [ViolCon ] = const_violation( con, ConBound );
    %
    conV = cell2mat(ViolCon)>0;
    CX = [con sqrt(FUN.val)];
    %     keyboard
    ConStop = max(0, max(cell2mat(ViolCon)));
    %  if impro<RES|| k == iter_max +1 || ConStop ==0
    xE=[xE xm];
    %if isx ~=1
     %   xU=[xU xm];
    %end
    yE = [yE FUN.val];
    for jj=1:ms
        C{jj}=[C{jj} real(con{jj})];
    end
    % figure(1)
    % plot(yi)
    % drawnow
    %% ploting constraints
    % plot_preview(xi,yi,C )
    %     prettyplot_IMEXRK( xi, yi, C )
    [  CV_S ] = const_violation( C, ConBound );
    prettyplot_IMEXRK( xE, yE, CV_S)
    drawnow
    disp([num2str(k/iter_max*100), ' % Completed'])
end

%%

for ii=1:length(yE)
    maxviol(ii)=yE(ii);
    for jj=1:ms
        maxviol(ii)=max( maxviol(ii), C{jj}(ii));
    end
end

%[mi, idm] = min(max([C{1};C{2};C{3};C{4};C{5};C{6};yE-Search.constant]))
%SolC = xE(:,idm); SolC = round(SolC*MESH_SIZE)./MESH_SIZE
% [con , FUN, CS , CON, SolC, Be, Bi ]=IMEXRK_Solver3(SolC)
