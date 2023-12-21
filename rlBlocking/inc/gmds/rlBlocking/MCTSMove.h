//
// Created by bourmaudp on 02/12/22.
//
/*----------------------------------------------------------------------------------------*/
#ifndef GMDS_MCTSMOVE_H
#define GMDS_MCTSMOVE_H
/*----------------------------------------------------------------------------------------*/
#include "LIB_GMDS_RLBLOCKING_export.h"
/*----------------------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------------------*/
/** @class  MCTSMove
 *  @brief  Structure that provides ....
 */
struct LIB_GMDS_RLBLOCKING_API MCTSMove {
	/*------------------------------------------------------------------------*/
	/** @brief  Destructor
	 */
	virtual ~MCTSMove() = default;
	/*------------------------------------------------------------------------*/
	/** @brief  Overloaded ==
	 */
	virtual bool operator==(const MCTSMove& AOther) const = 0;
};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------------------*/
#endif     // GMDS_MCTSMOVE_H
/*----------------------------------------------------------------------------------------*/
