//
// Created by rochec on 13/01/2022.
//
/*------------------------------------------------------------------------*/
#include <gmds/claire/LevelSet2D.h>
#include <limits>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/


LevelSet2D::LevelSet2D(Mesh *AMesh, int AmarkFrontNodes) {
	m_mesh = AMesh;
	m_markFrontNodes = AmarkFrontNodes;
	m_distance = m_mesh->newVariable<double,GMDS_NODE>("distance");

}


/*------------------------------------------------------------------------*/
LevelSet2D::STATUS LevelSet2D::execute()
{
	initialisationDistances();

	// Tant qu'il y a des noeuds dans la map
	while(!m_DistanceMap.isEmpty()){
		TCellID n0_id;
		double v0;
		m_DistanceMap.getAndRemoveFirst(v0, n0_id);
	}

	return LevelSet2D::SUCCESS;
}
/*------------------------------------------------------------------------*/




/*-------------------------------------------------------------------*/
void LevelSet2D::initialisationDistances(){
	for (auto id:m_mesh->nodes()){
		if(m_mesh->isMarked<Node>(id, m_markFrontNodes)){
			m_DistanceMap.add(0, id);
			m_distance->set(id, 0);
		}
		else{
			m_DistanceMap.add(std::numeric_limits<double>::max(), id);
			m_distance->set(id, 1000);
		}
	}
};
/*-------------------------------------------------------------------*/





/*-------------------------------------------------------------------*/
void LevelSet2D::getValue(TCellID n_id, double &v0){
	v0 = m_distance->value(n_id);
};
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
void LevelSet2D::setValue(TCellID n_id, double v0){
	m_distance->value(n_id) = v0 ;
};
/*-------------------------------------------------------------------*/