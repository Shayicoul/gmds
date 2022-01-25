//
// Created by rochec on 19/01/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/LevelSetFromIntToOut.h>
//#include <limits>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

LevelSetFromIntToOut::LevelSetFromIntToOut(Mesh *AMesh, int AmarkFrontNodesInt, int AmarkFrontNodesOut) {
	m_mesh = AMesh;
	m_markFrontNodesInt = AmarkFrontNodesInt;
	m_markFrontNodesOut = AmarkFrontNodesOut;
	m_distance = m_mesh->newVariable<double,GMDS_NODE>("distance_combined");
	m_distance_Int = m_mesh->newVariable<double,GMDS_NODE>("distance_int");
	m_distance_Out = m_mesh->newVariable<double,GMDS_NODE>("distance_out");

}





/*------------------------------------------------------------------------*/
LevelSetFromIntToOut::STATUS
LevelSetFromIntToOut::execute()
{
	initialisationDistancesInt();
	initialisationDistancesOut();

	combineDistanceFields();

	return LevelSetFromIntToOut::SUCCESS;
}
/*------------------------------------------------------------------------*/



/*-------------------------------------------------------------------*/
void
LevelSetFromIntToOut::setValue(TCellID n_id, double v0){
	m_distance->value(n_id) = v0 ;
};
/*-------------------------------------------------------------------*/



/*-------------------------------------------------------------------*/
void
LevelSetFromIntToOut::combineDistanceFields() {
	double distInt;
	double distOut;
	for (auto id:m_mesh->nodes()){
		distInt = m_distance_Int->value(id);
		distOut = m_distance_Out->value(id);
		if (distInt+distOut == 0) {
			setValue(id, -1 );
		}
		else {
			setValue(id, distInt / (distInt + distOut));     // Attention, pas de barrière si distInt+dOut =0 (c'est à dire si les deux fronts ont un noeud commun)
		}
	}
};
/*-------------------------------------------------------------------*/



/*-------------------------------------------------------------------*/
void
LevelSetFromIntToOut::initialisationDistancesInt() {
	LevelSetNaif lsInt(m_mesh, m_markFrontNodesInt);
	lsInt.execute();
	double distInt;
	for (auto id:m_mesh->nodes()){
		lsInt.getValue(id, distInt);
		m_distance_Int->set(id, distInt);
	}
	m_mesh->deleteVariable(GMDS_NODE, "distance");
};
/*-------------------------------------------------------------------*/



/*-------------------------------------------------------------------*/
void
LevelSetFromIntToOut::initialisationDistancesOut() {
	LevelSetNaif lsOut(m_mesh, m_markFrontNodesOut);
	lsOut.execute();
	double distOut;
	for (auto id:m_mesh->nodes()){
		lsOut.getValue(id, distOut);
		m_distance_Out->set(id, distOut);
	}
	m_mesh->deleteVariable(GMDS_NODE, "distance");
};
/*-------------------------------------------------------------------*/