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
#include <KM/IO/VTKWriter.h>
/*----------------------------------------------------------------------------*/
// using namespace kmds;
/*----------------------------------------------------------------------------*/
class WriterTest : public ::testing::Test
{
protected:
    WriterTest()
    {
        ;
    }
    virtual ~WriterTest()
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
//        int num_threads = 4;
//        int use_numa = 1;
//        int use_core = 1;
//        Kokkos::OpenMP::initialize(num_threads, use_numa, use_core);
        Kokkos::initialize(kargs);
    }

    static void
    TearDownTestCase()
    {
        // Kokkos::Serial::finalize();
        // Kokkos::Threads::finalize();
//        Kokkos::OpenMP::finalize();
        Kokkos::finalize();
    }
};
/*----------------------------------------------------------------------------*/
TEST_F(WriterTest, N2F_quadtri)
{
    kmds::Mesh m;

    const int nx = 3;
    const int ny = 3;

    const int ni = nx - 1;
    const int nj = ny - 1;

    m.updateNodeCapacity(nx*ny+1);
    m.updateFaceCapacity(ni*nj+2);

    kmds::TCellID nodeIDs[nx*ny+1];
    int index = 0;

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            nodeIDs[index] = m.newNode(i, j, 0.);
            index++;
        }
    }
    int a = m.addNode();
    m.setNodeLocation(a, nx, 1., 0.);

    for (int i = 0; i < ni; i++) {
        for (int j = 0; j < nj; j++) {

            m.newQuad(i * ny + j,
                      (i + 1) * ny + j,
                      (i + 1) * ny + (j + 1),
                      i * ny + (j + 1)
            );

        }
    }
    m.newTriangle((nx-1)*ny,
                  a,
                  (nx-1)*ny+1
    );
    m.newTriangle((nx-1)*ny+1,
                  a,
                  (nx-1)*ny+2
    );

    kmds::VTKWriter<kmds::Mesh> w(m);
    w.write("vtkwritertest_2D", kmds::F);
}
/*----------------------------------------------------------------------------*/
TEST_F(WriterTest, N2R_hexpyrtetprism)
{
    kmds::Mesh m;

    m.updateNodeCapacity(13);
    m.updateRegionCapacity(4);

    kmds::TCellID nodeIDs[13];
    nodeIDs[0] = m.newNode(0., 0., 0.);
    nodeIDs[1] = m.newNode(1., 0., 0.);
    nodeIDs[2] = m.newNode(1., 1., 0.);
    nodeIDs[3] = m.newNode(0., 1., 0.);
    nodeIDs[4] = m.newNode(0., 0., 1.);
    nodeIDs[5] = m.newNode(1., 0., 1.);
    nodeIDs[6] = m.newNode(1., 1., 1.);
    nodeIDs[7] = m.newNode(0., 1., 1.);

    nodeIDs[8] = m.newNode(0.5, 0.5, 2.);
    nodeIDs[9] = m.newNode(1.5, 0.5, 1.5);

    nodeIDs[10] = m.newNode(0., 1., 1.);
    nodeIDs[11] = m.newNode(0., 1., 1.);
    nodeIDs[12] = m.newNode(0., 1., 1.);


    kmds::TCellID hexId = m.newHexahedron(
            nodeIDs[0],
            nodeIDs[1],
            nodeIDs[2],
            nodeIDs[3],
            nodeIDs[4],
            nodeIDs[5],
            nodeIDs[6],
            nodeIDs[7]
    );

    kmds::TCellID pyrId = m.newPyramid(
            nodeIDs[4],
            nodeIDs[5],
            nodeIDs[6],
            nodeIDs[7],
            nodeIDs[8]
    );

    kmds::TCellID tetId = m.newTetrahedron(
            nodeIDs[5],
            nodeIDs[6],
            nodeIDs[8],
            nodeIDs[9]
    );

    kmds::TCellID prismId = m.newPrism3(
            nodeIDs[7],
            nodeIDs[4],
            nodeIDs[8],
            nodeIDs[10],
            nodeIDs[11],
            nodeIDs[12]
    );

    kmds::Variable<int>* varint = m.createVariable<int>(-1, kmds::KMDS_REGION, "varint");
    (*varint)[hexId] = 17;
    (*varint)[pyrId] = 27;
    (*varint)[tetId] = 37;
    (*varint)[prismId] = 47;

    kmds::Variable<double>* vardouble = m.createVariable<double>(-1., kmds::KMDS_REGION, "vardouble");
    (*vardouble)[hexId] = 7.1;
    (*vardouble)[pyrId] = 7.2;
    (*vardouble)[tetId] = 7.3;
    (*vardouble)[prismId] = 7.4;

    kmds::Variable<gmds::math::Vector>* varvec = m.createVariable<gmds::math::Vector>(gmds::math::Vector(0., 0., 0.), kmds::KMDS_REGION, "varvec");
    (*varvec)[hexId] = gmds::math::Vector (1.,0.,0.);
    (*varvec)[pyrId] = gmds::math::Vector (1.,1.,0.);
    (*varvec)[tetId] = gmds::math::Vector (1.,1.,1.);
    (*varvec)[prismId] = gmds::math::Vector (2.,0.,0.);

    kmds::VTKWriter<kmds::Mesh> w(m);
    w.write("vtkwritertest_3D", kmds::R);
}
/*----------------------------------------------------------------------------*/