

IdentifiabilityGain<-function(input_g,cur_observ_Vec,nodeId,anscestor,desendent,immediate)
{
   #  i<-1
#     nodeId<-node_state_vec[i]
#     anscestor<-ance_Nodes[[node_state_vec[i]]]
#     desendent<-desc_Nodes[[node_state_vec[i]]]
#     immediate<-imme_Nodes[[node_state_vec[i]]]
#     input_g<-g
#     plot(input_g)
#
     #----------------------------------------
     #remove all the out-edges of descendent observed nodes and immediate observed nodes
     if(length(desendent)>0)
     {
        for(i in 1:length(desendent))
        {
            if(cur_observ_Vec[desendent[i]]==1)
            {
                out_neighbors<-neighbors(input_g,as.vector(desendent[i]),mode="out")
                if(length(out_neighbors)>0)
                {
                    for(j in 1:length(out_neighbors))
                    {
                        input_g<-input_g-edge(get.edge.ids(input_g,c(as.vector(desendent[i]),as.vector(out_neighbors[j]))))
                    } 
                } 
            }
        }
     }
     if(length(immediate)>0)
     {
        for(i in 1:length(immediate))
        {
            if(cur_observ_Vec[immediate[i]]==1)
            {
                out_neighbors<-neighbors(input_g,as.vector(immediate[i]),mode="out")
                if(length(out_neighbors)>0)
                {
                    for(j in 1:length(out_neighbors))
                    {
                          input_g<-input_g-edge(get.edge.ids(input_g,c(as.vector(immediate[i]),as.vector(out_neighbors[j]))))
                    } 
                } 
            }
        }
     }
     #remove all the in-edges of ancestor observed nodes  and  all the out-edge to immediate nodes of ancestor observed nodes
     if(length(anscestor)>0)
     {
        for(i in 1:length(anscestor))
        {
            if(cur_observ_Vec[anscestor[i]]==1)
            {
                in_neighbors<-neighbors(input_g,as.vector(anscestor[i]),mode="in")
                if(length(in_neighbors)>0)
                {
                    for(j in 1:length(in_neighbors))
                    {
                        input_g<-input_g-edge(get.edge.ids(input_g,c(as.vector(in_neighbors[j]),as.vector(anscestor[i]))))
                    } 
                }
                out_neighbors<-neighbors(input_g,as.vector(anscestor[i]),mode="out")
                boundary_Res<-out_neighbors %in%  immediate
                if(length(boundary_Res)>0)
                {
                    for(j in 1:length(boundary_Res))
                    {
                        if(boundary_Res[j])
                        {
                             input_g<-input_g-edge(get.edge.ids(input_g,c(as.vector(anscestor[i]),as.vector(out_neighbors[j]))))
                        }
                    } 
                } 
            }
        }
     }
     
     #get the number of observed nodes which are connected with nodeId by Wright's path
      
     cluster_Res<-components(input_g, mode = "weak")
     tmp_ance_Nodes<-ego(input_g,node_Num,as.vector(nodeId),mode="in",mindist = 1)
     added_number<-0
     for(i in 1:vcount(input_g))
     {
          if((cluster_Res[["membership"]][i]==cluster_Res[["membership"]][nodeId])&&(cur_observ_Vec[i]==1))
          {
               cur_ance_Nodes<-ego(input_g,node_Num,i,mode="in",mindist = 1)
               path_Res<-all_simple_paths(input_g,as.vector(nodeId),i,mode="in")
               path_Res2<-all_simple_paths(input_g,as.vector(nodeId),i,mode="out")
               path_Res<-c(path_Res,path_Res2)
               common_anstor_Flag<-any(cur_ance_Nodes %in% tmp_ance_Nodes)
               if(common_anstor_Flag || length(path_Res)>0)
               { 
                    added_number<-added_number+1
               }
          }
     }
     
     return(added_number)
}