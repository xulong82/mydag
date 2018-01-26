function [new_dag,R,n_child] = MBR_move(ns,output,dag)


%%%output is the pre_computed scores, ns is node size, dag is curretn dag
[new_dag,node]=get_new_MB_move_dag(output,ns,dag);

n_child=length(children);
GRAPH{1}=dag;
GRAPH{2}=new_dag;
R=Acceptance_Probability_computation_PMBR(GRAPH,output,node,children);




function [new_dag,node,children] = get_new_MB_move_dag(output,ns,dag)


new_dag=dag;
node=randi(length(ns)); %%%randomly select a node
new_dag(:,node)=0;%%orphan the selected node

pso = find(dag(:,node))';%%store its orgianl parents
children = find(dag(node,:))';%%%store its children nodes

log_scores=output{node}.log_scores;%%%restrive pre-computed scores for the node

%%%%find descend nodes
A_NEW = expm(new_dag') - eye(length(ns));
A_NEW = (A_NEW>0);

descendants_of_node = find(A_NEW(:,node));

         
%%%exclude descendants nodes and original parent nodes   
for k = [descendants_of_node',pso]                                   
    x1 = find(output{node}.parents(:,1)==k);
    x2 = find(output{node}.parents(:,2)==k);
    x3 = find(output{node}.parents(:,3)==k);
    x = [x1;x2;x3];
    log_scores(x)=-Inf;
end

%%%%sample new parent set 
log_max=max(log_scores)-100;
p=exp(log_scores-log_max);
cummulative1=cumsum(p)/sum(p);
x= rand;
indicis = find((cummulative1>=x));     
pick1  = indicis(1);
ps=output{node}.parents(pick1,:);
ps=ps(find(ps));
new_dag(ps,node)=1;%%%add new parent set
log_score=output{node}.log_scores(pick1);



%%%%%%%%%%% sample parents for children nodes%%%%%%%%%

new_dag(:,children')=0;%%orphan children node


if ~isempty(children)
% obj children's parents
 for child=children
     
      log_scores=output{child}.log_scores;
     
      A_NEW = expm(new_dag') - eye(length(ns));
      A_NEW = (A_NEW>0);
      descendants = find(A_NEW(:,child))';
    %exclude its descendants node
            for k = descendants                                  
                x1 = find(output{child}.parents(:,1)==k);
                x2 = find(output{child}.parents(:,2)==k);
                x3 = find(output{child}.parents(:,3)==k);
                x = [x1;x2;x3];
                log_scores(x)=-Inf;
            end
   
   %include its previous parent node
   Indicator_Vector = zeros(length(log_scores),1);         
   x_node = find(output{child}.parents(:,1)==node | output{child}.parents(:,2)==node | output{child}.parents(:,3)==node);
   Indicator_Vector(x_node,1)=1; 
   indicis_0 = find(Indicator_Vector==0);            
   log_scores(indicis_0) = -Inf;
   
   %%%sample new parent set        
    log_max=max(log_scores)-100;
    p=exp(log_scores-log_max);
    cummulative1=cumsum(p)/sum(p);
    x= rand;
    indicis = find((cummulative1>=x));     
    pick1  = indicis(1);


    childrens_parents=output{child}.parents(pick1,:);
    if ~isempty(childrens_parents)
        new_dag(childrens_parents(find(childrens_parents)),child)=1;


    end

end

end



