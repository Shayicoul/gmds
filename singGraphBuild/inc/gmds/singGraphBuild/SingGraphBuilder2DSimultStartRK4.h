/*----------------------------------------------------------------------------*/
/*
 * SingGraphBuilder2DSimultStartRK4.h
 *
 *  Created on: April 10 2014
 *      Author: F. Ledoux
 */
/*----------------------------------------------------------------------------*/
#pragma once
/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_SINGGRAPHBUILD_export.h"
#include <gmds/singGraphBuild/SingularityGraph.h>
#include <gmds/singGraphBuild/SingularityGraphBuilder2D.h>
#include <gmds/singGraphBuild/Tools.h>
/*----------------------------------------------------------------------------*/

namespace gmds {

class LIB_GMDS_SINGGRAPHBUILD_API SingGraphBuilder2DSimultStartRK4 : public SingularityGraphBuilder2D
{
 public:
	SingGraphBuilder2DSimultStartRK4(Mesh *AMesh, Variable<math::Cross2D> *AField, const bool ABuildGeomSing = true);

	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Creation of singularity lines starting simultaneously from all singularities
	 *	(using RK4 computation)
	 *   formely SingularityGraphBuilder2D::createSingularityLinesSimultaneousStartRK4
	 */
	void createSingularityLines();
};
}     // namespace gmds
