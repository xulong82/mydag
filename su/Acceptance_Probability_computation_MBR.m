function [ACCEPTANCE] = Acceptance_Probability_computation_PMBR(GRAPH,Output,node,children);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GRAPH_1  = GRAPH{1}; %old dag
GRAPH_O  = GRAPH{1}; %old dag
GRAPH_2  = GRAPH{2}; %new dag
GRAPH_N  = GRAPH{2}; %new dag
n        = length(GRAPH_1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%forward move%%%%%%%%%%%%%%%%%%%%%%%
children = find(GRAPH_1(node,:))'; %find old dag's childre
ps1 = find(GRAPH_1(:,node))'; %find old dag's parents
new_parents_node = find(GRAPH_2(:,node));  %find new dag's parent nodes

GRAPH_1(:,node)=0; %orphan the selected node

%%%find descendant nodes
A_1 = expm(GRAPH_1') - eye(n);
A_1 = (A_1>0);
descendants_of_node = find(A_1(:,node));

SCORES_COPY_1 = Output{node}.log_scores;


% exclude descendants_of_node and original parent nodes 
for k = [descendants_of_node',ps1]                                 
    x1 = find(Output{node}.parents(:,1)==k);
    x2 = find(Output{node}.parents(:,2)==k);
    x3 = find(Output{node}.parents(:,3)==k);
    x = [x1;x2;x3];
    SCORES_COPY_1(x)=-Inf;
end

%%add new parent node 
GRAPH_1(new_parents_node,node)=1;
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
GRAPH_1(:,children)=0;        
log_sum_node_2=0;
max_2=0;
        
        if ~isempty(children)
           for child = children        
      
           A_1 = expm(GRAPH_1') - eye(n);
           A_1 = (A_1>0);
           
           descendants_2 = find(A_1(:,child));% all descendants of node_a
               
           SCORES_COPY_2 = Output{child}.log_scores;
           
           %%exclude descendants node to avoid cycles 
                for k=descendants_2'
                    x1 = find(Output{child}.parents(:,1)==k);
                    x2 = find(Output{child}.parents(:,2)==k);
                    x3 = find(Output{child}.parents(:,3)==k);
                    x = [x1;x2;x3];
                    SCORES_COPY_2(x) = -Inf;
                end    
         
            % include the selected node 
            Indicator_Vector = zeros(length(SCORES_COPY_2),1);
            node_x = find(Output{child}.parents(:,1)==node | Output{child}.parents(:,2)==node | Output{child}.parents(:,3)==node);
            Indicator_Vector(node_x,1)=1;  % indicates relevant parent-sets      
            indicis_0 = find(Indicator_Vector==0);
            SCORES_COPY_2(indicis_0) = -Inf;
            
            
            max_2 = max_2+max(SCORES_COPY_2);
            log_sum_node_2 =log_sum_node_2+ log(sum(exp(SCORES_COPY_2-max(SCORES_COPY_2))));%%%sum of local scores of valid parent sets of children nodes
            SCORES_COPY_2=[];
            new_parents = find(GRAPH_N(:,child));
            GRAPH_1(new_parents,child)=1; %%%add new parent set to the graph
           end
           
        end        



        
 
%%%%%%%%%%%%%%%%   reverse move  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
        ps = find(GRAPH_2(:,node))'; %find parents of the selected node in the new dag
        GRAPH_2(:,node)= 0;%orphan the selected node

        %find descendant nodes
        A_2 = expm(GRAPH_2') - eye(n);
        A_2 = (A_2>0);           
        
   
           
        SCORES_COPY_b = Output{node}.log_scores;
        descendants_of_node = find(A_2(:,node));
        
        %exclude descendant nodes and its current original nodes
        for k = [descendants_of_node',ps]                                  
            x1 = find(Output{node}.parents(:,1)==k);
            x2 = find(Output{node}.parents(:,2)==k);
            x3 = find(Output{node}.parents(:,3)==k);
            x = [x1;x2;x3];
            SCORES_COPY_b(x)=-Inf;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        old_parents_node = find(GRAPH_1(:,node));
        GRAPH_2(old_parents_node,node)=1;%%%add old parent node
         
        GRAPH_2(:,children)= 0;
        log_sum_node_a=0;
        max_a=0;
        
        if ~isempty(children)
           for child = children       
      
           %find all descendants of a child
           A_2 = expm(GRAPH_2') - eye(n);
           A_2 = (A_2>0);
           
           descendants_a = find(A_2(:,child));% 
           SCORES_COPY_a = Output{child}.log_scores;
           
           %exclude all descendants of a child node
                for k=descendants_a'
                    x1 = find(Output{child}.parents(:,1)==k);
                    x2 = find(Output{child}.parents(:,2)==k);
                    x3 = find(Output{child}.parents(:,3)==k);
                    x = [x1;x2;x3];
                    SCORES_COPY_a(x) = -Inf;
                end    
         
            %%%incldue the selected node
            Indicator_Vector = zeros(length(SCORES_COPY_a),1);
            x = find(Output{child}.parents(:,1)==node | Output{child}.parents(:,2)==node | Output{child}.parents(:,3)==node );
            Indicator_Vector(x,1)=1;  % indicates relevant parent-sets
            x=[];       
            indicis_0 = find(Indicator_Vector==0);
            SCORES_COPY_a(indicis_0) = -Inf;
            
            
            max_a = max_a+max(SCORES_COPY_a);
            log_sum_node_a =log_sum_node_a+ log(sum(exp(SCORES_COPY_a-max(SCORES_COPY_a))));%%%sum of local scores of valid parent sets of children nodes
            SCORES_COPY_a=[];
            old_parents = find(GRAPH_O(:,child));
            GRAPH_2(old_parents,child)=1;
           end
           
        end        
                
                

 		 max_1 = max(SCORES_COPY_1);
		 max_b = max(SCORES_COPY_b);
         
		 
         log_sum_node_1 = log(sum(exp(SCORES_COPY_1-max_1)));%%sum of local scores of selected node in the old dag
		 log_sum_node_b = log(sum(exp(SCORES_COPY_b-max_b))); %%sum of local scores of selected node in the new dag

     
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
log_A_nominator   = log_sum_node_1 + log_sum_node_2 + max_1 + max_2;          
log_A_denominator = log_sum_node_a + log_sum_node_b + max_a + max_b;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
log_ACCEPTANCE = log_A_nominator - log_A_denominator;

ACCEPTANCE = exp(log_ACCEPTANCE);


return
