/*------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/GridBuilder.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/IGMeshIOService.h>
/*------------------------------------------------------------------------*/
#include <vector>
#include <set>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/
void reset_path(Mesh*AM, Variable<int>* APath){
   for(auto e:AM->edges())
       APath->set(e,0);
}
/*------------------------------------------------------------------------*/
int get_nb_adj_edge_in_path(const Node& AN,Variable<int>* APath){
    std::vector<TCellID > n_edge_ids = AN.getIDs<Edge>();
    int nb_edges_in_path=0;
    for(auto n_edge_i:n_edge_ids){
        if(APath->value(n_edge_i)==1){
            nb_edges_in_path++;
        }
    }
    return nb_edges_in_path;
}
/*------------------------------------------------------------------------*/
std::set<TCellID> init_bnd_path(Mesh*AM, Variable<int>* APath){
    std::set<TCellID> bnd_path;
    for(auto e:AM->edges()){
        Edge cur_edge = AM->get<Edge>(e);
        std::vector<Node> cur_nodes = cur_edge.get<Node>();
        for(auto& n:cur_nodes){
            auto nb_edges_in_path = get_nb_adj_edge_in_path(n, APath);
            if(nb_edges_in_path==0)
                throw GMDSException("error during path boundary init");
            if(nb_edges_in_path==1)
                bnd_path.insert(n.id());
        }
    }
}
/*------------------------------------------------------------------------*/
Edge select_next_edge_in_path(const Node& AN,
                Variable<double>* APheromone,
                Variable<int>* APath)
{
   std::vector<Edge> es = AN.get<Edge>(), candidates, in_path;
   for(auto e:es){
       if(APath->value(e.id())==0){
           candidates.push_back(e);
       }
       else{
           in_path.push_back(e);
       }
   }
   if(in_path.size()!=1)
        throw GMDSException("select_next_edge_in_path: Error in the path");

   //Now we pick one edge among the remainding one, considering the pheromon
}
/*------------------------------------------------------------------------*/
bool is_boundary_edge(const Edge& AE){
    return AE.get<Face>().size()==1;
}
/*------------------------------------------------------------------------*/
bool is_boundary_node(const Node& AN){
    return AN.get<Face>().size()!=4;
}
/*------------------------------------------------------------------------*/
void build_path(Mesh*AM,
                Variable<int>* AConstraint,
                Variable<double>* APheromone,
                Variable<int>* APath){
    // make the path empty
    reset_path(AM,APath);
    // add the constraint into the path
    for(auto e:AM->edges()){
        APath->set(e,AConstraint->value(e));
    }
    //We store the boundary nodes of the path.
    std::set<TCellID> bnd_path = init_bnd_path(AM, APath);
    while (!bnd_path.empty()){
        //get the first element
        TCellID bnd_node_i = *bnd_path.begin();
        //and remove it from the set
        bnd_path.erase(bnd_path.begin(),bnd_path.begin());

        //Now we pick the next edge
        Edge next_edge = select_next_edge_in_path(AM->get<Node>(bnd_node_i),APheromone,APath);
        //Add the edge in the path
        APath->set(next_edge.id(),1);
        //And the other extremity is added in the boundary if we didn't reach a dead end, that can be:
        // 1. we have a loop
        // 2. reach a boundary node from an inner edge
        std::vector<Node> next_edge_nodes = next_edge.get<Node>();
        Node next_node = next_edge_nodes[0];
        if(next_node.id()==bnd_node_i)
            next_node= next_edge_nodes[1];

        if(get_nb_adj_edge_in_path(next_node, APath)>=2){
            // 1. loop if 2, otherwise it is an invalid path (non-manifold)
            // so the node must already be in the bnd_path so
            bnd_path.erase(next_node.id());
        }
        else { // it is equal to 1
            //if the node is on the boundary and the incoming edge inside the domain, we stop too
            //case 2
            if(is_boundary_node(next_node) && !is_boundary_edge(next_edge)){
                //do nothing
            }
            else
                bnd_path.insert(next_node.id());
        }
    }
}
/*------------------------------------------------------------------------*/
double quality_path(){

}
/*------------------------------------------------------------------------*/
bool is_valid_path(){

}
/*------------------------------------------------------------------------*/
int main(){
    std::cout<<"Ant Colony"<<std::endl;
    double max_pheromone= 100;
    int nb_ants=10;

    Mesh m(MeshModel(DIM3|N|F|E|F2N|E2N|N2E|E2F|F2E));
    GridBuilder gb(&m,2);
    gb.execute(8,1,8,1);

    m.deleteFace(4);
    m.deleteFace(5);
    m.deleteFace(6);
    m.deleteFace(12);
    m.deleteFace(13);
    MeshDoctor doc(&m);
    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    Variable<int>* var_constraint = m.newVariable<int, GMDS_EDGE>("constraint");
    var_constraint->set(23,1);
    var_constraint->set(24,1);

    Variable<double>* var_pheromone = m.newVariable<double, GMDS_EDGE>("pheromone");
    Variable<int>* var_path = m.newVariable<int, GMDS_EDGE>("path");

    //================================================================
    //init the pheromone
    //================================================================
    for(auto e:m.edges())
        var_pheromone->set(e,max_pheromone);

    //================================================================
    // Loop over ants
    //================================================================
    for(int i=0;i<nb_ants;i++){
        //We build a path for an ant
        build_path(&m,var_constraint, var_pheromone, var_path);
    }
    IGMeshIOService ios(&m);
    VTKWriter wf(&ios);
    wf.setCellOptions(N|F);
    wf.setDataOptions(N|F);
    wf.write("ant_f.vtk");

    VTKWriter w(&ios);
    w.setCellOptions(N|E);
    w.setDataOptions(N|E);
    w.write("ant.vtk");
    std::cout<<"-- it's done! --"<<std::endl;
}