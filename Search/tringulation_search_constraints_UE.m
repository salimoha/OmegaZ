function [xm ym CSm indm] = tringulation_search_constraints(inter_par_p,inter_par_g,xi,A,b)
% This function evalutes the
global n Search tri bnd1 bnd2
%  keyboard
ym=inf;
CSm=1;
tri=delaunayn(xi.');
for ind=1:size(tri,1)
    %          keyboard
    [xc,R2]=circhyp(xi(:,tri(ind,:)), n);
    if sqrt(R2)<2e3
        % body center of simplex as an initial point
        x=xi(:,tri(ind,:))*ones(n+1,1)/(n+1);
        if Search.method ==1
            %                keyboard
            [x y CS]=Constant_K_Search_FilterSQP(x,inter_par_p, inter_par_g,xc,R2,Search);
            %
            if (CS<CSm)
                xm=x; ym=y; CSm = CS; indm=ind;
            end
            if (y<ym & CS==CSm)
                xm=x; ym=y; indm=ind;
            end
            %%
        elseif Search.method ==2
            %       keyboard
            [x,A1,b1]= min_projection_control_feas(x,bnd1,bnd2,A,b);
            [x y CS]= Constraint_Adaptive_K_Search_UE(x,inter_par_p,inter_par_g,xc,R2, Search,A1,b1);
            
            if CS ==0
                % there is a point that pk < y0
                xm=x; ym=y; CSm = CS; indm=ind;
                break
            end
            if y < ym & CS ==1
                %           keyboard
                xm=x; ym=y; CSm = CS; indm=ind;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%
    end
    % go to the next simplex till you find a point that satisfies the
    % abovementioned criteria
end
end

