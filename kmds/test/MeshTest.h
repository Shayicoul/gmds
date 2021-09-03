/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux and N. Le Goff (2015)
 *
 * franck.ledoux@cea.fr
 * nicolas.le-goff@cea.fr
 *
 * The GMDS library is a computer program whose purpose is to provide a set of
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
// Google Test headers
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
// Kokkos headers
#include <Kokkos_Core.hpp>
/*----------------------------------------------------------------------------*/
// GMDS headers
#include <KM/DS/Mesh.h>
/*----------------------------------------------------------------------------*/
// using namespace kmds;
/*----------------------------------------------------------------------------*/
class MeshTest : public ::testing::Test
{
 protected:
        MeshTest()
        {
                ;
        }
        virtual ~MeshTest()
        {
                ;
        }

    static void
    SetUpTestCase()
    {
            // Kokkos::Serial::initialize();
            // Kokkos::Threads::initialize();
            Kokkos::InitArguments kargs;
            kargs.num_threads = 3;
//            int num_threads = 4;
//            int use_numa = 1;
//            int use_core = 1;
//            Kokkos::OpenMP::initialize(num_threads, use_numa, use_core);
            Kokkos::initialize(kargs);
    }

    static void
    TearDownTestCase()
    {
            // Kokkos::Serial::finalize();
            // Kokkos::Threads::finalize();
//            Kokkos::OpenMP::finalize();
            Kokkos::finalize();
    }
};
/*----------------------------------------------------------------------------*/
TEST_F(MeshTest, init)
{
        kmds::Mesh m;

        EXPECT_EQ(0, m.getNbNodes());
        EXPECT_EQ(0, m.getNbFaces());
        EXPECT_EQ(0, m.getNbRegions());
}
/*----------------------------------------------------------------------------*/
TEST_F(MeshTest, nodesCreation)
{
        kmds::Mesh m;

        const int nx = 3;
        const int ny = 3;
        const int nz = 3;

        m.updateNodeCapacity(nx*ny*nz);

        for(int i=0; i<nx; i++) {
                for(int j=0; j<ny; j++) {
                        for(int k=0; k<nz; k++) {
                                int a = m.addNode();
                                m.setNodeLocation(a, i, j, k);
                        }
                }
        }

        EXPECT_EQ(nx*ny*nz, m.getNbNodes());
}
/*----------------------------------------------------------------------------*/
TEST_F(MeshTest, FaceCreation)
{
        kmds::Mesh m;

        const int nx = 3;
        const int ny = 3;

        const int ni = nx -1;
        const int nj = ny -1;

        m.updateNodeCapacity(nx*ny);
        m.updateFaceCapacity(ni*nj);

        for(int i=0; i<nx; i++) {
                for(int j=0; j<ny; j++) {
                        int a = m.addNode();
                        m.setNodeLocation(a, i, j, 0.);
                }
        }

        for(int i=0; i<ni; i++) {
                for(int j=0; j<nj; j++) {

                        m.newQuad(i*ny+j,
                                  (i+1)*ny+j,
                                  (i+1)*ny+(j+1),
                                  i*ny+(j+1)
                        );

                }
        }

        EXPECT_EQ(nx*ny, m.getNbNodes());
        EXPECT_EQ(ni*nj, m.getNbFaces());
        EXPECT_EQ(m.getNbQuads(), m.getNbQuads());
}
/*----------------------------------------------------------------------------*/
TEST_F(MeshTest, regionsCreation)
{
        kmds::Mesh m;

        const int nx = 3;
        const int ny = 3;
        const int nz = 3;

        const int ni = nx -1;
        const int nj = ny -1;
        const int nk = nz -1;

        m.updateNodeCapacity(nx*ny*nz);
        m.updateRegionCapacity(ni*nj*nk);

        for(int i=0; i<nx; i++) {
                for(int j=0; j<ny; j++) {
                        for(int k=0; k<nz; k++) {
                                int a = m.addNode();
                                m.setNodeLocation(a, i, j, k);
                        }
                }
        }

        for(int i=0; i<ni; i++) {
                for(int j=0; j<nj; j++) {
                        for(int k=0; k<nk; k++) {
                                m.newHexahedron(i*ny*nz+j*nz+k,
                                                (i+1)*ny*nz+j*nz+k,
                                                (i+1)*ny*nz+(j+1)*nz+k,
                                                i*ny*nz+(j+1)*nz+k,
                                                i*ny*nz+j*nz+k+1,
                                                (i+1)*ny*nz+j*nz+k+1,
                                                (i+1)*ny*nz+(j+1)*nz+k+1,
                                                i*ny*nz+(j+1)*nz+k+1

                                );
                        }
                }
        }

        EXPECT_EQ(nx*ny*nz, m.getNbNodes());
        EXPECT_EQ(ni*nj*nk, m.getNbRegions());
        EXPECT_EQ(m.getNbHexahedra(), m.getNbRegions());

        kmds::Region r = m.getRegion(0);
}
/*----------------------------------------------------------------------------*/
TEST_F(MeshTest, regionsDeletion)
{
        kmds::Mesh m;

        const int nx = 3;
        const int ny = 3;
        const int nz = 3;

        const int ni = nx -1;
        const int nj = ny -1;
        const int nk = nz -1;

        m.updateNodeCapacity(nx*ny*nz);
        m.updateRegionCapacity(ni*nj*nk);

        for(int i=0; i<nx; i++) {
                for(int j=0; j<ny; j++) {
                        for(int k=0; k<nz; k++) {
                                int a = m.addNode();
                                m.setNodeLocation(a, i, j, k);
                        }
                }
        }

        for(int i=0; i<ni; i++) {
                for(int j=0; j<nj; j++) {
                        for(int k=0; k<nk; k++) {
                                m.newHexahedron(i*ny*nk+j*nz+k,
                                                (i+1)*ny*nz+j*nz+k,
                                                (i+1)*ny*nz+(j+1)*nz+k,
                                                i*ny*nz+(j+1)*nz+k,
                                                i*ny*nz+j*nz+k+1,
                                                (i+1)*ny*nz+j*nz+k+1,
                                                (i+1)*ny*nz+(j+1)*nz+k+1,
                                                i*ny*nz+(j+1)*nz+k+1

                                );
                        }
                }
        }

        EXPECT_EQ(nx*ny*nz, m.getNbNodes());
        EXPECT_EQ(ni*nj*nk, m.getNbRegions());

        m.removeRegion(ni-1);
        EXPECT_EQ(ni*nj*nk-1, m.getNbRegions());
}
/*----------------------------------------------------------------------------*/
TEST_F(MeshTest, regionsAccessor)
{
        kmds::Mesh m;

        const int nx = 3;
        const int ny = 3;
        const int nz = 3;

        const int ni = nx -1;
        const int nj = ny -1;
        const int nk = nz -1;

        m.updateNodeCapacity(nx*ny*nz);
        m.updateRegionCapacity(ni*nj*nk);

        for(int i=0; i<nx; i++) {
                for(int j=0; j<ny; j++) {
                        for(int k=0; k<nz; k++) {
                                int a = m.addNode();
                                m.setNodeLocation(a, i, j, k);
                        }
                }
        }

        for(int i=0; i<ni; i++) {
                for(int j=0; j<nj; j++) {
                        for(int k=0; k<nk; k++) {
                                m.newHexahedron(i*ny*nk+j*nz+k,
                                                (i+1)*ny*nz+j*nz+k,
                                                (i+1)*ny*nz+(j+1)*nz+k,
                                                i*ny*nz+(j+1)*nz+k,
                                                i*ny*nz+j*nz+k+1,
                                                (i+1)*ny*nz+j*nz+k+1,
                                                (i+1)*ny*nz+(j+1)*nz+k+1,
                                                i*ny*nz+(j+1)*nz+k+1
                                );
                        }
                }
        }

        kmds::GrowingView<kmds::TCellID> cells("CELLS", m.getNbRegions());
        m.getRegionIDs(&cells);
        EXPECT_EQ(cells.getNbElems(), m.getNbRegions());

        m.removeRegion(ni-1);
        kmds::GrowingView<kmds::TCellID> cells_bis("CELLS", m.getNbRegions());
        m.getRegionIDs(&cells_bis);
        EXPECT_EQ(cells_bis.getNbElems(), m.getNbRegions());

        bool found = false;
        for(int i=0; i<cells_bis.getNbElems(); i++) {
                if(cells_bis.get(i) == ni-1) {
                        found = true;
                }
        }
        EXPECT_FALSE(found);
}
/*----------------------------------------------------------------------------*/