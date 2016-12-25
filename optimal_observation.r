require("igraph")
source("identifiability_gain-fun.r")
 
#g<-read_graph("DAG-8-nodes.gml",format="gml")
#obs_file<-file("DAG-8-observation.txt")

#cat("0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n", file = "Influenza.txt")
#g<-read_graph("DAG-9-nodes.gml",format="gml")
#obs_file<-file("DAG-9-observation.txt")
#g<-read_graph("DCG-9-Nodes.gml",format="gml")
#obs_file<-file("DCG-9-Nodes.txt")
g<-read_graph("Influenza.gml",format="gml")
obs_file<-file("Influenza.txt")

observ<-read.table(obs_file,header = FALSE, sep =" ", nrows =1)
observ<-as.numeric(observ)

connected_Flag<-is_connected(g,mode="weak")
if(!connected_Flag)
{
    print("The graph is not connected, program will quit!")
    stopifnot(connected_Flag)
}

need_Fun_Num<-ecount(g)
node_Num<-vcount(g)
plot(g)

#find ancestor nodes, descendant nodes, immediate nodes and boundary nodes of each node
ance_Nodes<-list()
desc_Nodes<-list()
imme_Nodes<-list()
necesary_in<-vector()
necesary_out0<-vector()
necesary_out1<-vector()

for(i in 1:node_Num)
{
    in_neighbor<-ego(g,node_Num,V(g)[i],mode="in",mindist = 1)
    ance_Nodes<-c(ance_Nodes,list(in_neighbor[[1]]))
    out_neighbor<-ego(g,node_Num,V(g)[i],mode="out",mindist = 1)
    desc_Nodes<-c(desc_Nodes,list(out_neighbor[[1]]))
    immediate_nodes<-V(g)[-c(i,in_neighbor[[1]],out_neighbor[[1]])]
    imme_Nodes<-c(imme_Nodes,list(immediate_nodes))    
}

#the nodes with in-degree=0 and out-degree<3; out-degree=0; out-degree=1 must be observed
inD<-degree(g,v=V(g),mode="in",loops=FALSE,normalized = FALSE)
outD<-degree(g,v=V(g),mode="out",loops=FALSE,normalized = FALSE)
observ_Node_in<-which(inD %in% 0)

if(length(observ_Node_in)>0)
{
    necesary_in<-vector()
    for(i in 1:length(observ_Node_in))
    {
        if(outD[observ_Node_in[i]]<3)
        {
             necesary_in<-c(necesary_in,observ_Node_in[i])
        }
    }
    observ_Node_in<-necesary_in
}
necesary_node<-observ_Node_in
observ_Node_out<-which(outD %in% 0)
necesary_node<-c(necesary_node,observ_Node_out)
observ_Node_out_middle<-which(outD %in% 1)
necesary_node<-c(necesary_node,observ_Node_out_middle)
necesary_node<-unique(necesary_node)

valid_Fun_Num<-0
cur_observ_Vec<-rep(0,node_Num)

if(length(necesary_node)>0)
{
     cur_observ_Vec[necesary_node[1]]<-1
     if(length(necesary_node)>1)
     {
          for(i in 2:length(necesary_node))
          {
              added_equation_num<-IdentifiabilityGain(g,cur_observ_Vec,necesary_node[i],ance_Nodes[[necesary_node[i]]],desc_Nodes[[necesary_node[i]]],imme_Nodes[[necesary_node[i]]])
              valid_Fun_Num<-valid_Fun_Num+added_equation_num
              cur_observ_Vec[necesary_node[i]]<-1
          }
     }      
}
cur_observ_Vec<-as.numeric(cur_observ_Vec|observ)

#dynamic programming
found_flag<-FALSE
res_list<-list() 
if(valid_Fun_Num>=need_Fun_Num)
{
    res_list<-c(res_list,list(cur_observ_Vec))
    found_flag<-TRUE 
}else
{
    node_state_vec<-V(g)[-necesary_node]
    observed_list<-list()
    total_Fun_vec<-vector()
    if(length(node_state_vec)>0)
    {
        for(stage_id in 1:length(node_state_vec))
        {
            if(stage_id==1)
            {
                for(i in 1:length(node_state_vec))
                {
                     added_equation_num<-IdentifiabilityGain(g,cur_observ_Vec,node_state_vec[i],ance_Nodes[[node_state_vec[i]]],desc_Nodes[[node_state_vec[i]]],imme_Nodes[[node_state_vec[i]]])
                     total_Fun<-valid_Fun_Num+added_equation_num
                     total_Fun_vec<-c(total_Fun_vec,total_Fun)
                     observed_vec<-cur_observ_Vec
                     observed_vec[node_state_vec[i]]<-1
                     observed_list<-c(observed_list,list(observed_vec))
                }
            }
            else
            {
                next_total_Fun<-total_Fun_vec
                next_observed_list<-observed_list
                for(i in stage_id:length(node_state_vec))
                {
                     tmp_total_fun_vec<-vector()
                     for(j in (stage_id-1):(i-1))
                     {
                          added_equation_num<-IdentifiabilityGain(g,observed_list[[j]],node_state_vec[j],ance_Nodes[[node_state_vec[j]]],desc_Nodes[[node_state_vec[j]]],imme_Nodes[[node_state_vec[j]]])
                          total_Fun<-total_Fun_vec[j]+added_equation_num
                          tmp_total_fun_vec<-c(tmp_total_fun_vec,total_Fun)
                     }
                     next_total_Fun[i]<-max(tmp_total_fun_vec)
                     chosed_index<-match(max(tmp_total_fun_vec),tmp_total_fun_vec)
                     chosed_index<-chosed_index-1+stage_id-1 
                     new_observed_vec<-observed_list[[chosed_index]]
                     new_observed_vec[node_state_vec[i]]<-1
                     next_observed_list[[i]]<-new_observed_vec
                }
                total_Fun_vec<-next_total_Fun
                observed_list<-next_observed_list
            }
            if(max(total_Fun_vec)>=need_Fun_Num)
            {
                #output the optimal observation result
                for(i in 1:length(node_state_vec))
                {
                    if(total_Fun_vec[i]>=need_Fun_Num)
                    {
                        res_list<-c(res_list,list(observed_list[[i]]))
                    }
                }
                found_flag<-TRUE
                break
            }    
        }
    }
}

if(found_flag)
{
    print("The optimal observation strategy is:")
    print(res_list)
}else  #output failure to find the optimal observation strategy
{
    print("Failure to find the optimal observation strategy") 
}


