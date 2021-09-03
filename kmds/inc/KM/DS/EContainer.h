/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux and N. Le Goff (2015)
 *
 * franck.ledoux@cea.fr
 * nicolas.le-goff@cea.fr
 *
 * The KMDS library is a computer program whose purpose is to provide a set of
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
 * data to be ensured and,  more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 */
/*----------------------------------------------------------------------------*/
/*
 * NodeContainer.h
 *
 *  Created on: 14 feb 2018
 *      Author: legoff
 */
/*----------------------------------------------------------------------------*/
#ifndef KMDS_ECONTAINER_H_
#define KMDS_ECONTAINER_H_
/*----------------------------------------------------------------------------*/
// STL headers
#include <iostream>
/*----------------------------------------------------------------------------*/
// Kokkos headers
#include <Kokkos_Core.hpp>
/*----------------------------------------------------------------------------*/
// KMDS headers
#include <KM/Utils/KTypes.h>
/*----------------------------------------------------------------------------*/
namespace kmds {
/*----------------------------------------------------------------------------*/
    class Mesh;
/*----------------------------------------------------------------------------*/
    class EXPORT_KMDS EContainer
    {
    public:
        /*------------------------------------------------------------------------*/
        /** \brief Constructor
                *
                * \param[in] AOwner the mesh that owns this edge container
                * \param[in] ACapacity initial capacity of the container
                                */
        EContainer(Mesh* AOwner, const TInt32 ACapacity = 16);

        /*------------------------------------------------------------------------*/
        /** \brief  Destructor
         */
        virtual ~EContainer();

        /*------------------------------------------------------------------------*/
        /** \brief Indicate if this container contains a node of id AID
         *
         *  \param AID a cell id
         */
        bool has(const TCellID AID) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Add/allocate \p ANb new cells
         *
         * \param[in] ANb the number of cells to be created
         * \return    The index of the first created cell. Others are contiguous righ behind.
         */
        TCellID add(const TInt32 ANb);

        /*------------------------------------------------------------------------*/
        /** \brief  Add/allocate one
         *
         * \return    The index of the created cell
         */
        TCellID add();
        TCellID add_unsafe();

        TCellID newEdge(const TCellID AN1, const TCellID AN2);
        TCellID newEdge_unsafe(const TCellID AN1, const TCellID AN2);

        /*------------------------------------------------------------------------*/
        /** \brief  Remove/Deallocate cell of id \p AId.
                 *
         * \param[in] AId The index of the cell to be removed.
         */
        void remove(const TCellID AId);

        /*------------------------------------------------------------------------*/
        /** \brief  Set the nodes of cell \p AId
         *
         *  \param[in] AI cell id we want to get the nodes
         *  \param[in] AN new nodes of the cell
         */
        void set(const TCellID AId, const TCellID AN1, const TCellID AN2);

        /*------------------------------------------------------------------------*/
        /** \brief  Get the nodes of cell \p AI
         *
         *  \param[in] AI cell id we want to get the nodes of
         *  \param[in] AN nodes of the cell
         */
        void get(const TCellID AI, Kokkos::View<TCellID*>& AV) const;
        void get(const TCellID AId, TCellID& AN1, TCellID& AN2) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Returns the capacity in terms of array size (including holes)
         */
        TSize capacity() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Returns the top id (including holes)
         */
        TSize top() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Returns the number of cells really stored
         */
        TSize nbCells() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Resize the container
         */
        void resize(const TInt32 ASize);

        /*------------------------------------------------------------------------*/
        /** \brief  Empty the container
         */
        void removeAll();

    private:
        /** mesh owner of this cell container */
        Mesh* m_mesh;

        /** store cell info*/
        Kokkos::View<bool*> m_cells;
        Kokkos::View<TCellID * [2]> m_E2N;
        //, Kokkos::MemoryTraits<Kokkos::Atomic> I remove the memory traits since I cannot resize otherwise

        /** top of the heap including holes so */
        TSize m_top;
    };
/*----------------------------------------------------------------------------*/
}  // namespace KMDS
/*----------------------------------------------------------------------------*/
#endif /* KMDS_ECONTAINER_H_ */
/*----------------------------------------------------------------------------*/
