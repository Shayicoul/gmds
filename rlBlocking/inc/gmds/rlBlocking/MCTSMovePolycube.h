//
// Created by bourmaudp on 02/12/22.
//
/*----------------------------------------------------------------------------------------*/
#ifndef GMDS_MCTSMOVE_POLYCUBE_H
#define GMDS_MCTSMOVE_POLYCUBE_H
/*----------------------------------------------------------------------------------------*/
#include "LIB_GMDS_RLBLOCKING_export.h"
#include <gmds/rlBlocking/MCTSMove.h>
/*----------------------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------------------*/
/** @class  MCTSMove
 *  @brief  Structure that provides ....
 */
struct LIB_GMDS_RLBLOCKING_API MCTSMovePolycube: public MCTSMove {
	/*------------------------------------------------------------------------*/
	/** @brief  Destructor
	 */
	virtual ~MCTSMovePolycube();
	/*------------------------------------------------------------------------*/
	/** @brief  Overloaded ==
	 */
	virtual bool operator==(const MCTSMove& AOther) const;
};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------------------*/
#endif     // GMDS_MCTSMOVE_POLYCUBE_H
/*----------------------------------------------------------------------------------------*/
