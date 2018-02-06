function [Output] = Score_all(data,ns)

[n,ncases] = size(data);

type = cell(1,n);
params = cell(1,n);
for i=1:n
type{i} = 'tabular';
params{i} = { 'prior_type', 'dirichlet', 'dirichlet_weight', 1 };
end
discrete = 1:n;
clamped = zeros(n, ncases);
scoring_fn='bayesian';


[n_nodes, thrash] = size(data);

Prior = zeros(n_nodes,1);
for i=1:(n_nodes-1)
	Prior(i+1,1) = Prior(i,1) + log(i) - log(n_nodes-i);
end

for child = 1:n_nodes 

% Score of the empty parent-set    
ps      = [];
Score   = score_local(child, ps, ns, data,discrete,params,clamped,scoring_fn,type);
Score   = Score + Prior(1);
Counter = 1;    
Structure.parents(Counter,1:3) = [0 0 0];
Structure.scores(Counter) = Score;
Counter = Counter +1;

% Local Bayesian scores of parent-sets of size 1:
    for parent_1 = 1:n_nodes  
    if(child~=parent_1)
    ps    = parent_1;
    Score = score_local(child, ps, ns, data,discrete,params,clamped,scoring_fn,type);
    Score = Score + Prior(2);
    Structure.parents(Counter,1:3) = [0, 0, parent_1];
    Structure.scores(Counter)      = Score;
    Counter = Counter +1;
    end
    end  

       
% Local Bayesian scores of parent-sets of size 2:
    for parent_1 = 1:(n_nodes-1)
        for parent_2 = (parent_1+1):n_nodes
            if(child~=parent_1 && child~=parent_2)
            ps    = [parent_1;parent_2]';
            Score =score_local(child, ps, ns, data,discrete,params,clamped,scoring_fn,type);
            Score = Score + Prior(3);
            Structure.parents(Counter,1:3) = [0, parent_1, parent_2];
            Structure.scores(Counter)      = Score;
            Counter = Counter +1;
            end
        end
    end

% Local Bayesian scores of parent-sets of size 3: 


for parent_1 = 1:(n_nodes-2)
     for parent_2 = (parent_1+1):(n_nodes-1)
       for parent_3 = (parent_2+1):n_nodes
            if(child~=parent_1 && child~=parent_2 &&child~=parent_3)
            ps    = [parent_1;parent_2;parent_3]';
            Score = score_local(child, ps, ns, data,discrete,params,clamped,scoring_fn,type)
            Score   = Score + Prior(4);
            Structure.parents(Counter,1:3) = [parent_1, parent_2, parent_3];
            Structure.scores(Counter)      = Score;
            Counter = Counter +1;
            end
       end
     end
end

% Delete all illegal parent-sets (containing the child-node itself):
indicis = find(Structure.parents(:,1)~=child & Structure.parents(:,2)~=child & Structure.parents(:,3)~=child);
Structure.parents = Structure.parents(indicis,:);
Structure.scores  = Structure.scores(indicis);


[thrash, indicis] = sort(Structure.scores,'descend');

Output{child}.parents       = Structure.parents(indicis,:);
Output{child}.log_scores    = Structure.scores(indicis);
clear Structure

end 

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function score = score_local(j, parents, ns, data,discrete,params,clamped,scoring_fn,type)

u = find(clamped(j,:)==0);
score=score_family(j, parents, type{j}, scoring_fn, ns, discrete, data(:,u), params{j});

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% 
% for node = 1:n_nodes 
%     for child= 1:n_nodes  
%         if(node~=child)   
%          temp_ps=output{child}.parents;
%          
%          temp_score=output{child}.log_scores;
%          
%          a=[0 0 0];
%          [a1,b1,c1]=intersect(output{node}.parents,a,'rows');
%          
%          single_score=temp_score(b1);
%          
%         
%          
%          [r,c,v]=find(temp_ps==node);
%          ps=temp_ps(r,:);
%          scores=temp_score(r);%-single_score;
%          
%          [thrash, indicis] = sort(scores,'descend');
% 
%          output{node}.Children{child}.parents=ps(indicis,:);
%          output{node}.Children{child}.scores=scores(indicis);
%         end
%     end
% end
% 







