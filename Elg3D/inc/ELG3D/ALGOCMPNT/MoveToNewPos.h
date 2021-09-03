/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux and N. Le Goff (2015)
 *
 * franck.ledoux@cea.fr
 * nicolas.le-goff@cea.fr
 *
 * This software is a computer program whose purpose is to provide a set of
 * functionnalities to represent and handle any type of meshes (2D, 3D,
 * triangles, tetrahedra, quad, hexa, polygons, polyhedra, etc.) and write
 * meshing algorithms. So it gathers many mathematical objects like points,
 * segment, quaternions, etc. and basic algorithms useful to build more evolved
 * ones.
 *
 * This software is governed by the CeCILL-C license under French law and
 * abiding by the rules of distribution of free software.  You can  use,
 * modify and/ or redistribute the software under the terms of the CeCILL-C
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability.
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and, more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 */
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    MoveToNewPos.h
 *  \author  legoff
 *  \date    05/18/2018
 */
/*----------------------------------------------------------------------------*/
#ifndef ELG3D_MOVETONEWPOS_H_
#define ELG3D_MOVETONEWPOS_H_
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <KM/DS/Mesh.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
#include "ELG3D/DATACMPNT/FracPres.h"
#include "ELG3D/DATACMPNT/MaterialAssignment.h"
/*----------------------------------------------------------------------------*/
namespace elg3d {
/*----------------------------------------------------------------------------*/

    const int MoveToNewPos_NBMAXCOLORS = 20;


    /*------------------------------------------------------------------------*/
    /** \brief  Move the nodes from the selection to their provided new position.
     *          Additionally projects unto the geometric association of the node, if any,
     *          and does not move if the new position provided is incorrect (ie NaN).
     *
     *  \param[in] ASelectionNodes the selection of the nodes that will be moved
     *  \param[in,out]  AMesh the mesh
     *  \param[in] AVarNodeDestination the new prospective positions
     *  \param[in] AVarNodeGeomAssociation the geometric association of the nodes, for those which have one
     *  \param[out] AVarNodeNbMoves the number of effective displacements
     */
    void moveToNewPos_basicMove(const kmds::GrowingView<kmds::TCellID>* ASelectionInterfaceNodes,
                                kmds::Mesh* AMesh,
                                const kmds::Variable<gmds::math::Point>* AVarNewPos,
                                const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation);


    /*------------------------------------------------------------------------*/
    /** \brief  Move the nodes from the selection to their provided new position, by increments.
     *          Additionally projects unto the geometric association of the node along the way, if any.
     *
     *  \param[in] ANbIter the displacement will be done by this number of increments
     *  \param[in] AQualThreshold a node will stay in place if the minimum scaled jacobian of
     *             the adjacent cells falls below this threshold
     *  \param[in] ASelectionNodes the selection of the nodes that will be moved
     *  \param[in,out]  AMesh the mesh
     *  \param[in] AVarNodeDestination the new prospective positions
     *  \param[in] AVarNodeGeomAssociation the geometric association of the nodes, for those which have one
     *  \param[out] AVarNodeNbMoves the number of effective displacements for each node (for debug purposes)
     *
     */
    void moveToNewPos_noBadMove_2D(const int ANbIter,
                                   const kmds::TCoord qualThreshold,
                                   const kmds::GrowingView<kmds::TCellID>* ASelectionNodes,
                                   kmds::Mesh* AMesh,
                                   const kmds::Connectivity* c_N2F,
                                   const kmds::Variable<gmds::math::Point>* AVarNodeDestination,
                                   const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation,
                                   kmds::Variable<int>* AVarNodeNbMoves);

    void moveToNewPos_noBadMove_3D(const int ANbIter,
                                   const kmds::TCoord qualThreshold,
                                   const kmds::GrowingView<kmds::TCellID>* ASelectionNodes,
                                   kmds::Mesh* AMesh,
                                   const kmds::Connectivity* c_N2R,
                                   const kmds::Variable<gmds::math::Point>* AVarNodeDestination,
                                   const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation,
                                   kmds::Variable<int>* AVarNodeNbMoves);


/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
#endif /* ELG3D_MOVETONEWPOS_H_ */
