function [xm ym CSm indm] = tringulation_search_constraints(inter_par_p,inter_par_g,xi,tri)
% This function solves the optimization problem for search functions within
% the feasible domain for the search fucntion subject to the constraints.
global n bnd1 bnd2 Search Eps ms MAXRAD
% keyboard
ym=inf; CSm=1;MAXRAD=2e2;
%   keyboard

% [ind0] = find_min_simplex(xi,tri, inter_par_p, inter_par_g);
% solve the problem in each simplex
for ind=1:size(tri,1)
%      keyboard
    [xc,R2]=circhyp(xi(:,tri(ind,:)), n);
    if sqrt(R2)<MAXRAD % bounded circumradius criteria for convergance
        x=xi(:,tri(ind,:))*ones(n+1,1)/(n+1); % midpoint 

        
        if Search.method ==1
            %                The constant K mehtod
            scheme = 'snopt';
            switch scheme
                case 'SQP_filter'
                    [x y CS]=Constant_K_Search_FilterSQP(x,inter_par_p, inter_par_g,xc,R2,Search);
                case 'SQP'
                    [x y CS]=Constant_K_Search_sqp(x,inter_par_p,inter_par_g,xc,R2, Search);
                case 'primal_dual'
                    [x y CS]=Constant_K_Search_2(x,inter_par_p,inter_par_g,xc,R2,Search);
                case 'snopt'
                    [x y CS]=Constant_K_Search_snopt_2(x,inter_par_p,inter_par_g,xc,R2,Search);
                case 'fmincon'
                    [x y CS]=Constant_K_Search_fmincon(x,inter_par_p, inter_par_g,xc,R2,Search);
%                 case 'interior'
%                     [x y CS]=Constant_K_Search_interior(x,inter_par_p, inter_par_g,xc,R2,Search);
                    %     case 'newton'
                    %     [x y CS]=Constant_K_Search_backup(x,inter_par_p, inter_par_g,xc,R2);
%                 case 'multi_K'
%                     if length(Search.constant) ~=1
%                         % the idea for wighted K with poly. harmonic spline
%                         [x y CS]=Constant_K_Search_FilterSQP_multi_K(x,inter_par_p, inter_par_g,xc,R2,Search);
%                     else
%                         error('Change the scheme')
%                     end
            end
            % check 
            if (CS<CSm)
                xm=x; ym=y; CSm = CS; indm=ind;
            end
            if (y<ym & CS==CSm)
                %          elseif (y<ym & CS==CSm)
                xm=x; ym=y; indm=ind;
            end
            % Adaptive K method
        elseif Search.method ==2
            
            %TODO only search on the optimum simplex
%              keyboard
            scheme = 'exp_fitting_direct';
            switch scheme
                case 'exp_fitting_direct'
                    [x y CS]=Constraint_Adaptive_K_Search(x,inter_par_p,inter_par_g,xc,R2, Search);
                case 'exp_fitting_inverse'
                    [x y CS]=Constraint_Adaptive_K_Search_inverse(x,inter_par_p,inter_par_g,xc,R2, Search);
            end
            %
            if CS ==0
                % there is a point that pk < y0
                xm=x; ym=y; CSm = CS; indm=ind;
                break
            end
            if y < ym & CS ==1
%             if  ym-y > Eps/10 & CS ==1
                xm=x; ym=y; CSm = CS; indm=ind;
            end
         
%             keyboard
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%
        
        
        end
    
    % go to the next simplex till you find a point that satisfies the
    % abovementioned criteria
       
end
end

function [ind, ym] = find_min_simplex(xi, tri, inter_par_p, inter_par_g)
% this function finds the index of the simplex over all the simplecies
% which the value of its cicumcenter max{(p-y0,g)/e} is minimum
global n ms MAXRAD Search
% solve the problem in each simplex
for ind=1:size(tri,1)
%      keyboard
    [xc,R2]=circhyp(xi(:,tri(ind,:)), n);
    if sqrt(R2)<MAXRAD % bounded circumradius criteria for convergance
        x=xi(:,tri(ind,:))*ones(n+1,1)/(n+1); % midpoint 
        [xmf,ymf(ind,1)]=inter_min(x,inter_par_p,1);
         for ii=1:ms
          [xmc,ymc(ind,ii)]=inter_min(x,inter_par_g{ii},1);
         end
      [ym,ind0] = min(max([ymc,ymf-Search.constant]'));
    end
end
if isempty(ind0),
    warning('no simplex was found ....!!!!')
    ind=1:size(tri,1);
end
           
end