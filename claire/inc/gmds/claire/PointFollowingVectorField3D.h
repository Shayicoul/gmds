//
// Created by rochec on 27/01/2022.
//

#ifndef GMDS_POINTFOLLOWINGVECTORFIELD3D_H
#define GMDS_POINTFOLLOWINGVECTORFIELD3D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/ig/Mesh.h>
#include <string>
#include <map>
#include <fstream>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_CLAIRE_API PointFollowingVectorField3D
{
 public:
	/*--------------------------------------------------------------------*/
	/** @enum  Status code for executing algorithms
	 */
	typedef enum {
		FAIL,
		SUCCESS
	} STATUS;

	/*-------------------------------------------------------------------*/
	/** @brief Constructor.
         *  @param AMesh the mesh where we work on
	 */
	PointFollowingVectorField3D(Mesh *AMesh, math::Point A_Pstart, double A_distance, Variable<math::Vector3d>* A_gradient3D);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();

 private:
	/*-------------------------------------------------------------------*/
	/** @brief Return true if the point M is in the triangle
	 */
	bool isInTriangle(TCellID face_id, math::Point M);
	/*-------------------------------------------------------------------*/
	/** @brief Return in which triangle M is
	 */
	TCellID inWhichTriangle(math::Point M);
	/*-------------------------------------------------------------------*/
	/** @brief Compute the minimal edge's lenght
	 */
	double minEdgeLenght();
	/*-------------------------------------------------------------------*/

 private:
	/** mesh we work on */
	Mesh *m_mesh;
	/** starting point */
	math::Point m_Pstart ;
	/** ending point */
	math::Point m_Pend ;
	/** the distance */
	double m_distance ;
	/** the vector field to follow */
	Variable<math::Vector3d>* m_gradient3D ;

};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_POINTFOLLOWINGVECTORFIELD3D_H
