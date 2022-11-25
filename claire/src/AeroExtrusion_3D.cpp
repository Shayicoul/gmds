//
// Created by rochec on 18/11/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/claire/Front_3D.h>
#include <gmds/claire/Utils.h>
#include <gmds/claire/AeroException.h>
#include <gmds/claire/AeroExtrusion_3D.h>
#include <gmds/claire/AdvectedPointRK4_3D.h>
#include <gmds/claire/AeroMeshQuality.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

AeroExtrusion_3D::AeroExtrusion_3D(Mesh *AMeshT, Mesh *AMeshH, ParamsAero Aparams_aero, Variable<double>* A_DistanceField, Variable<math::Vector3d>* A_VectorField) {
	m_meshT = AMeshT;
	m_meshH = AMeshH;
	m_params_aero = Aparams_aero;
	m_DistanceField = A_DistanceField;
	m_VectorField = A_VectorField;
}


/*------------------------------------------------------------------------*/
AeroExtrusion_3D::STATUS
AeroExtrusion_3D::execute()
{
	double pas_couche = 1.0/m_params_aero.nbr_couches ;

	std::vector<TCellID> surface_block_corners_Id;
	std::vector<TCellID> surface_block_faces_Id;

	for (auto n_id:m_meshH->nodes())
	{
		surface_block_corners_Id.push_back(n_id);
	}
	for (auto f_id:m_meshH->faces())
	{
		surface_block_faces_Id.push_back(f_id);
	}

	Front_3D Front_Geom = Front_3D(0, surface_block_corners_Id, surface_block_faces_Id);
	//Front_3D Current_Front = Front_Geom;

	Front_3D Current_Front = Compute1stLayer(Front_Geom, m_meshT->getVariable<double,GMDS_NODE>("GMDS_Distance_Int"), m_VectorField);

	// Compute the successive layers
	for (int i=2; i <= m_params_aero.nbr_couches; i++) {
		Current_Front = ComputeLayer(Current_Front, m_meshT->getVariable<double,GMDS_NODE>("GMDS_Distance"), i*pas_couche,
		                             m_VectorField);
	}

	Variable<int>* var_node_layer = m_meshH->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	Variable<int>* var_face_layer = m_meshH->getOrCreateVariable<int, GMDS_FACE>("GMDS_FACE_Couche_Id");
	for (auto f_id:m_meshH->faces())
	{
		Face f = m_meshH->get<Face>(f_id);
		std::vector<Node> nodes = f.get<Node>();
		int max_layer_index = std::max(var_node_layer->value(nodes[0].id()), var_node_layer->value(nodes[1].id()));
		max_layer_index = std::max(max_layer_index, var_node_layer->value(nodes[2].id()));
		max_layer_index = std::max(max_layer_index, var_node_layer->value(nodes[3].id()));
		var_face_layer->set(f_id, max_layer_index);
	}

	return AeroExtrusion_3D::SUCCESS;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
std::map<TCellID, TCellID>
AeroExtrusion_3D::ComputeIdealPositions(Front_3D AFront, double dist_cible, Variable<double>* A_distance, Variable<math::Vector3d>* A_vectors)
{
	std::map<TCellID, TCellID> map_idealNextNodes;
	std::vector<TCellID> front_nodes = AFront.getNodes();

	for (auto n_id:front_nodes){
		Node n = m_meshH->get<Node>(n_id);
		math::Point M = n.point();
		AdvectedPointRK4_3D advpoint(m_meshT, M, dist_cible, A_distance, A_vectors);
		advpoint.execute();
		Node n_new = m_meshH->newNode(advpoint.getPend());
		map_idealNextNodes[n_id] = n_new.id() ;
	}

	return map_idealNextNodes;
}
/*------------------------------------------------------------------------*/
Front_3D
AeroExtrusion_3D::ComputeLayer(Front_3D Front_IN, Variable<double>* A_distance, double dist_cible, Variable<math::Vector3d>* A_vectors){

	std::cout << "---------> build layer: " << Front_IN.getFrontID()+1 << std::endl;

	std::map<TCellID, TCellID> map_new_nodes = ComputeIdealPositions(Front_IN, dist_cible, A_distance, A_vectors);

	// Mise à jour de l'indice de couche
	Variable<int>* couche_id = m_meshH->getVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	for (auto n_id:Front_IN.getNodes()){
		couche_id->set(map_new_nodes[n_id], Front_IN.getFrontID()+1);
	}

	std::vector<TCellID> front_nodes = Front_IN.getNodes();
	std::vector<TCellID> front_faces = Front_IN.getFaces();


	// Ajout des hex restants
	for (auto f_id:front_faces){
		Face f = m_meshH->get<Face>(f_id);
		std::vector<Node> nodes = f.get<Node>();
		CreateNormalHexa(f_id, Front_IN, map_new_nodes);
	}

	// Supression des noeuds non utilisés
	math::Utils::MeshCleaner(m_meshH);



	// Initialisation du front de sortie

	std::vector<TCellID> new_front_nodes_id;
	std::vector<TCellID> new_front_faces_id;

	Variable<int>* var_node_couche_id = m_meshH->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	Variable<int>* var_face_couche_id = m_meshH->getOrCreateVariable<int, GMDS_FACE>("GMDS_FACE_Couche_Id");

	for (auto n_id:m_meshH->nodes())
	{
		if (var_node_couche_id->value(n_id) == Front_IN.getFrontID()+1)
		{
			new_front_nodes_id.push_back(n_id);
		}
	}
	for (auto f_id:m_meshH->faces())
	{
		if (var_face_couche_id->value(f_id) == Front_IN.getFrontID()+1)
		{
			new_front_faces_id.push_back(f_id);
		}
	}
	Front_3D Front_OUT = Front_3D(Front_IN.getFrontID()+1, new_front_nodes_id, new_front_faces_id);

	return Front_OUT;
}
/*------------------------------------------------------------------------*/
Front_3D
AeroExtrusion_3D::Compute1stLayer(Front_3D A_Front_IN, Variable<double>* A_distance, Variable<math::Vector3d>* A_vectors)
{
	std::cout << "---------> build layer: " << A_Front_IN.getFrontID()+1 << std::endl;

	std::map<TCellID, TCellID> map_new_nodes = ComputeIdealPositions(A_Front_IN, m_params_aero.delta_cl, A_distance, A_vectors);

	// Mise à jour de l'indice de couche
	Variable<int>* couche_id = m_meshH->getVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	for (auto n_id:A_Front_IN.getNodes()){
		couche_id->set(map_new_nodes[n_id], A_Front_IN.getFrontID()+1);
	}

	std::vector<TCellID> front_nodes = A_Front_IN.getNodes();
	std::vector<TCellID> front_faces = A_Front_IN.getFaces();

	// Ajout des hex restants
	for (auto f_id:front_faces){
		Face f = m_meshH->get<Face>(f_id);
		std::vector<Node> nodes = f.get<Node>();
		CreateNormalHexa(f_id, A_Front_IN, map_new_nodes);
	}

	// Supression des noeuds non utilisés
	math::Utils::MeshCleaner(m_meshH);

	// Initialisation du front de sortie
	std::vector<TCellID> new_front_nodes_id;
	std::vector<TCellID> new_front_faces_id;

	Variable<int>* var_node_couche_id = m_meshH->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	Variable<int>* var_face_couche_id = m_meshH->getOrCreateVariable<int, GMDS_FACE>("GMDS_FACE_Couche_Id");

	for (auto n_id:m_meshH->nodes())
	{
		if (var_node_couche_id->value(n_id) == A_Front_IN.getFrontID()+1)
		{
			new_front_nodes_id.push_back(n_id);
		}
	}
	for (auto f_id:m_meshH->faces())
	{
		if (var_face_couche_id->value(f_id) == A_Front_IN.getFrontID()+1)
		{
			new_front_faces_id.push_back(f_id);
		}
	}
	Front_3D Front_OUT = Front_3D(A_Front_IN.getFrontID()+1, new_front_nodes_id, new_front_faces_id);

	return Front_OUT;

}
/*------------------------------------------------------------------------*/
void
AeroExtrusion_3D::CreateNormalHexa(TCellID f_id, Front_3D &Front_IN, std::map<TCellID, TCellID> map_new_nodes)
{
	Variable<int>* var_node_couche_id = m_meshH->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	Variable<int>* var_face_couche_id = m_meshH->getOrCreateVariable<int, GMDS_FACE>("GMDS_FACE_Couche_Id");

	Face f = m_meshH->get<Face>(f_id);
	std::vector<Node> nodes = f.get<Node>();

	TCellID n0_id = map_new_nodes[nodes[0].id()] ;
	TCellID n1_id = map_new_nodes[nodes[1].id()] ;
	TCellID n2_id = map_new_nodes[nodes[2].id()] ;
	TCellID n3_id = map_new_nodes[nodes[3].id()] ;

	Node n0 = m_meshH->get<Node>(n0_id);
	Node n1 = m_meshH->get<Node>(n1_id);
	Node n2 = m_meshH->get<Node>(n2_id);
	Node n3 = m_meshH->get<Node>(n3_id);
	var_node_couche_id->set(n0_id, Front_IN.getFrontID()+1);
	var_node_couche_id->set(n1_id, Front_IN.getFrontID()+1);
	var_node_couche_id->set(n2_id, Front_IN.getFrontID()+1);
	var_node_couche_id->set(n3_id, Front_IN.getFrontID()+1);

	TCellID r_id = math::Utils::CreateHexaNConnectivities(m_meshH, nodes[0], nodes[1], nodes[2], nodes[3], n0, n1, n2, n3);

	TCellID f_new_layer_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshH, n0.id(), n1.id(), n2.id(), n3.id());
	var_face_couche_id->set(f_new_layer_id, Front_IN.getFrontID()+1);


}
/*------------------------------------------------------------------------*/