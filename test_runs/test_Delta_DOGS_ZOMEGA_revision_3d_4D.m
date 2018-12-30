% Search.Method=1 constant_K
% Search.Method=2 Adaptive_K

% For 3 and 4 dimension problem 2018 December

clear all; close all; clc






% global n m ms bnd1 bnd2 Ain bin acon Search r xi tri
global n m ms bnd1 bnd2  Search  Ain bin tri MESH_SIZE iter_max
% Search.constant=0;
%
% n=2;
% for InitNum = ['one', 'tre']

for n = 4
    % n= 3;
    ms=1;
    RES = 5e-2;
    x_star=ones(n,1) *0.153;
    % x_star=ones(n,1) *0.46;
    
    %
    Method = 2;
    %
    x0=zeros(n,1); KCf = 1; KC2 = 1;
    
    fun=@(x)  sum((x-x0).^2)*4; Var = 4*0.024*n;
    
    con{1}=@(x) (4*rastriginn2(x,n)-n )*KC2;
    
    
    
    
    %
    MESH_SIZE=8; % grid size% % Mss=20;
    % interpolaion strategy
    inter_method=1;
    
    
    
    for iiii = 1:5
        
        clearvars xE yE inter_par xU xi inter_par_g inter_par_p
        
        % Calculate the initial points
        
        % InitNum = 'one'; % bad initializaiton
        % InitNum = 'two'; %
        % InitNum = 'tre'; % good initialization
        % % if InitNum == 'one'
        %      InitNum = 'tre';
        
        
        InitNum = 'rdm';
        
        if InitNum == 'one'
            xE=ones(n ,1)*0.45;
            delta0=0.15;
        elseif InitNum == 'two'
            xE=ones(n,1)*0.6;
            delta0=0.15;
        elseif InitNum == 'tre'
            xE=ones(n,1)*0.3;
            delta0=0.15;
        else
            xE = rand(n,1);
            delta0=0.15;
            
        end
        
        for ii=1:n
            e=zeros(n,1); e(ii)=1;
            xE(:,ii+1)=xE(:,1)+delta0*e;
        end
        % calculates acon , yi, C
        lob=zeros(n,1); upb=ones(n,1);
        Search.method = Method;
        Search.constant = Var;
        bnd1 = lob;
        bnd2 = upb;
        xU=bounds(bnd1,bnd2, n);
        % Input the equality constraints
        Ain=[eye(n);-eye(n)];
        bin=[bnd2 ;-bnd1];
        % Calculate the function evaluation at initial points
        for ii=1:ms
            acon{ii}=[];
        end
        for ii=1:size(xE,2)
            ind=1:2*n;
            ind=ind(Ain*xE(:,ii)-bin>-0.01);
            for jj=1:length(ind)
                acon{ind(jj)}=[acon{ind(jj)} ii];
            end
            yE(ii)=fun(xE(:,ii));
            for jj=1:ms
                C{jj}(ii)=con{jj}(xE(:,ii));
            end
        end
        x_prev = xE(:,1);
        y_prev = yE(:,1);
        delta_tol=0.2;
        
        %
        iter_max = 300;
        y0=Search.constant;
        
        for kkk=1:9
            for k=1:iter_max
                %             save surr_pts_2 yi xi C tri acon
                disp( [' started iteration ' ,num2str(k), ' ... ', 'with number func. eval = ',num2str(size(xE,2))] )
                tri=delaunayn([xE xU].');
                xi=[xE xU];
                
                %%%%%%%%%%% Modify the interpolations %%%%%%%%%%%%%%%%%%%
                inter_par_p= interpolateparametarization(xE,yE,inter_method);
                for jj=1:ms
                    inter_par_g{jj}= interpolateparametarization(xE,C{jj},inter_method);
                end
                
                %%%%%%%%% Calculate Discrete search function %%%%%%%%%%%
                yup=zeros(1,size(xU,2));
                for ii=1:size(xU,2)
                    yup(ii)=estimate_max_cons_val(xU(:,ii),inter_par_p,inter_par_g,y0,ms)/mindis(xU(:,ii),xE);
                    %             yup(ii)=inf;
                end
                
                % Perform the search
                while 1
                    %              keyboard
                    [xm ym(k) Cs(k) indm2]= ...
                        tringulation_search_constraints(inter_par_p,inter_par_g,[xE xU], tri)
                    %keyboard
                    xm=round((xm-bnd1)./(bnd2-bnd1).*MESH_SIZE)./MESH_SIZE.*(bnd2-bnd1)+bnd1;
                    
                    [xm,xE,xU,newadd,success]=...
                        points_neighbers_find(xm,xE,xU);
                    
                    if success==1
                        break
                    else
                        yup=[yup estimate_max_cons_val(xm,inter_par_p,inter_par_g,y0,ms)/mindis(xm,xE)];
                    end
                end
                
                %                        % stopping criteria
                %         if mindis(xm,[xE xU])<1e-6
                %             break
                %         end
                
                
                %         evaluating step
                if (estimate_max_cons_val(xm,inter_par_p,inter_par_g,y0,ms)/mindis(xm,xE)>min(yup) )
                    [t,ind]=min(yup);
                    xm=xU(:,ind); xU(:,ind)=[];
                    %         end
                    %         evaluating step
                elseif  (mindis(xm,xU)<1e-6 && mindis(xm,xU)~=Inf)
                    [t,ind]=min(yup);
                    xm=xU(:,ind); xU(:,ind)=[];
                end
                
                
                
                
                %         % feasible constraint projection
                %         if mindis(xm,[xE xU])<1e-6
                %             break
                %         end
                
                
                %        identifying step
                if mindis(xm,xE)<1e-6
                    break
                end
                
                % Perform the function evalutions
                %     Evaluated set
                xE=[xE xm];
                yE = [yE fun(xm)];
                for jj=1:ms
                    C{jj}=[C{jj} con{jj}(xm)];
                end
                
                %
                if Search.method ==2
                    figure(11);
                    subplot(2,1,1)
                    plot(yE-Var)
                    subplot(2,1,2)
                    plot(max(yE-Var, C{1}))
                    
                end
                
                if min(max(yE-Var, C{1}))<= RES
                    
                    %                             if norm(xm-x_star)<RES
                    resFile = strcat('results_OmegaZ_paper_init_', InitNum, '_n_', num2str(n),'_new.mat');
                    save(resFile)
                           
    
                    break
                end
                
                %                      % stopping criteria
                %         if mindis(xm,[xE xU])<1e-6
                %             break
                %         end
                
                iplot_check=0;
                %% stoping
                disp([num2str(k/iter_max*100), ' % Completed'])
   
            end
            
            %     keyboard
            MESH_SIZE=MESH_SIZE*2;
        end
        
        stat.test = InitNum;
        stat.bndpt = sum(min(xE)<0.01 | max(xE)>0.99);
        stat.support = size(xU,2);
        stat.n = n;
        stat.meshrefinement = kkk;
        stat.NumfEval = size(xE,2);
        
        resFile = strcat('results_OmegaZ_stat_', InitNum, '_n_', num2str(n),'_new2.mat');
        save(resFile)
        % end
        % end
        
        %%
        %%
        % save results_OmegaZ_paper4D.mat
        disp('initial point &  \# of Func. Eval.  &  \# of support points &  \# of boundary points \n')
        fprintf(' %s & %d &  & %d  & %d & %d   \n ',InitNum, n, stat.NumfEval,  stat.support,  stat.bndpt);
        %
        
        
    end
    
    if n==3
        
        test3{iiii} = struct('xE',xE,'xU',xU, 'yE',yE, 'inter_par', inter_par, 'stat', stat, 'Nm', Nm );
    elseif n==4
        test4{iiii} = struct('xE',xE,'xU',xU, 'yE',yE, 'inter_par', inter_par, 'stat', stat, 'Nm', Nm );
        keyboard 
        
    end
end