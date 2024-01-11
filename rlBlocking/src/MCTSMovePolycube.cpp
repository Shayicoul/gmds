/*----------------------------------------------------------------------------*/
#include <gmds/rlBlocking/MCTSMovePolycube.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
MCTSMovePolycube::~MCTSMovePolycube()
{}
/*----------------------------------------------------------------------------*/
MCTSMovePolycube::MCTSMovePolycube(TCellID AIdEdge, TCellID AIdBlock, double AParamCut, bool ATypeMove)
		:m_AIdEdge(AIdEdge),m_AIdBlock(AIdBlock), m_AParamCut(AParamCut),m_typeMove(ATypeMove)
{}
/*----------------------------------------------------------------------------*/
bool
MCTSMovePolycube::operator==(const gmds::MCTSMove &AOther) const
{
	const MCTSMovePolycube &o = (const MCTSMovePolycube &) AOther;        // Note: Casting necessary
	return m_AIdEdge == o.m_AIdEdge && m_AIdBlock == o.m_AIdBlock
	       && m_AParamCut == o.m_AParamCut && m_typeMove== o.m_typeMove;
}
/*----------------------------------------------------------------------------*/
