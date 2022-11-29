//
// Created by rochec on 18/11/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/Front_3D.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
Front_3D::Front_3D(int front_id, std::vector<TCellID> nodes_Id, std::vector<TCellID> faces_Id) {
	m_FrontID = front_id;
	m_nodesId = nodes_Id;
	m_facesId = faces_Id;
}
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
void Front_3D::setFrontID(int layer_id){
	m_FrontID = layer_id;
}
/*-------------------------------------------------------------------*/
int Front_3D::getFrontID(){
	return m_FrontID;
}
/*-------------------------------------------------------------------*/
std::vector<TCellID> Front_3D::getNodes(){
	return m_nodesId;
};
/*-------------------------------------------------------------------*/
std::vector<TCellID> Front_3D::getFaces(){
	return m_facesId;
};
/*-------------------------------------------------------------------*/
void Front_3D::addNodeId(TCellID n_id){
	m_nodesId.push_back(n_id);
}
/*-------------------------------------------------------------------*/
void Front_3D::addFaceId(TCellID f_id){
	m_facesId.push_back(f_id);
}
/*-------------------------------------------------------------------*/
std::vector<TCellID>
Front_3D::orderedFrontEdgesAroundNode(Mesh *m, TCellID n_id)
{
	std::vector<TCellID> n_ordered_edges_on_Front;

	Node n = m->get<Node>(n_id);
	std::vector<Edge> n_edges = n.get<Edge>();
	Variable<int>* var_node_couche_id = m->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");

	// Get the front edges connected to n
	std::vector<Edge> n_edges_on_Front;
	for (auto e:n_edges)
	{
		Node n_opp = e.getOppositeNode(n);
		if (var_node_couche_id->value(n_opp.id()) == m_FrontID)
		{
			n_edges_on_Front.push_back(e);
		}
	}

	int mark_isTreated = m->newMark<Face>();
	n_ordered_edges_on_Front.push_back(n_edges_on_Front[0].id());	// Choose a first edge, and a first face

	for (int i=1;i<=n_edges_on_Front.size()-1;i++)
	{
		TCellID e_id = n_ordered_edges_on_Front[n_ordered_edges_on_Front.size()-1];
		Edge e = m->get<Edge>(e_id);
		Node n_opp = e.getOppositeNode(n_id);

		Face f;
		std::vector<Face> e_faces = e.get<Face>() ;
		for (auto f_loc:e_faces)
		{
			std::vector<Node> f_loc_nodes = f_loc.get<Node>();

			if (!m->isMarked(f_loc, mark_isTreated)
			    && var_node_couche_id->value(f_loc_nodes[0].id()) == m_FrontID
			    && var_node_couche_id->value(f_loc_nodes[1].id()) == m_FrontID
			    && var_node_couche_id->value(f_loc_nodes[2].id()) == m_FrontID
			    && var_node_couche_id->value(f_loc_nodes[3].id()) == m_FrontID)
			{
				f = f_loc;
			}
		}
		m->mark(f, mark_isTreated);

		std::vector<Edge> f_edges = f.get<Edge>();
		for (auto e_loc:f_edges)
		{
			std::vector<Node> e_loc_nodes = e_loc.get<Node>();
			if ( (e_loc_nodes[0].id() == n_id && e_loc_nodes[1].id() != n_opp.id() )
			    || (e_loc_nodes[1].id() == n_id && e_loc_nodes[0].id() != n_opp.id() ) )
			{
				n_ordered_edges_on_Front.push_back(e_loc.id());
			}
		}

	}

	m->unmarkAll<Face>(mark_isTreated);
	m->freeMark<Face>(mark_isTreated);

	return n_ordered_edges_on_Front;
}
/*-------------------------------------------------------------------*/
std::vector<TCellID>
Front_3D::edgeFacesOnFront(Mesh *m, TCellID e_id)
{
	std::vector<TCellID> e_faces_on_Front;

	Edge e = m->get<Edge>(e_id);
	std::vector<Face> e_faces = e.get<Face>();
	for (auto f:e_faces)
	{
		for (auto f_front_id:m_facesId)
		{
			if (f.id()==f_front_id)
			{
				e_faces_on_Front.push_back(f.id());
			}
		}
	}

	// e_faces_on_Front is supposed to be sized 2

	return e_faces_on_Front;
}
/*-------------------------------------------------------------------*/