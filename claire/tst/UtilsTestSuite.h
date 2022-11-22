//
// Created by rochec on 22/03/2022.
//

#include <gmds/math/Line.h>
#include <gmds/claire/Utils.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/BoundaryOperator2D.h>
#include <gmds/igalgo/GridBuilder.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <gtest/gtest.h>
#include <iostream>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/

TEST(UtilsTestClass, Utils_Test1)
{
	// Test de la méthode math::Utils::distFromNodeIds

	Mesh m(MeshModel(DIM3 | R | F | N | F2N | N2F));

	GridBuilder gb(&m,2);
	gb.execute(3,1.0, 4, 1.0);

	double eps(pow(10,-5));

	ASSERT_NEAR(math::Utils::distFromNodeIds(&m, 0, 8), 2.0, eps);
	ASSERT_NEAR(math::Utils::distFromNodeIds(&m, 0, 3), 3.0, eps);
	ASSERT_NEAR(math::Utils::distFromNodeIds(&m, 0, 11), 3.60555, eps);
	ASSERT_NEAR(math::Utils::distFromNodeIds(&m, 3, 1), 2.0, eps);
	ASSERT_NEAR(math::Utils::distFromNodeIds(&m, 7, 9), 2.23607, eps);
	ASSERT_NEAR(math::Utils::distFromNodeIds(&m, 1, 7), 2.23607, eps);
	ASSERT_NEAR(math::Utils::distFromNodeIds(&m, 10, 3), 2.23607, eps);

	//std::cout << math::Utils::distFromNodeIds(&m, 10, 3) << std::endl;

	IGMeshIOService ioService_geom(&m);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("Utils_Test1.vtk");

}


TEST(UtilsTestClass, Utils_Test2)
{
	// Test de la méthode math::Utils::CommonEdge

	// WE READ
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::E | gmds::N2E | gmds::N2F | gmds::F2N | gmds::E2N | gmds::F2E | gmds::E2F));

	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir + "/Aero/Poubelle/Carre.vtk";

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N | gmds::F);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	ASSERT_EQ(math::Utils::CommonEdge(&m, 0, 4), 212);
	ASSERT_EQ(math::Utils::CommonEdge(&m, 51, 54), 141);
	ASSERT_EQ(math::Utils::CommonEdge(&m, 62, 85), 14);
	ASSERT_EQ(math::Utils::CommonEdge(&m, 16, 17), 115);
	ASSERT_EQ(math::Utils::CommonEdge(&m, 78, 67), 1);
	ASSERT_EQ(math::Utils::CommonEdge(&m, 17, 40), 116);
	ASSERT_EQ(math::Utils::CommonEdge(&m, 86, 51), 203);

	ASSERT_EQ(math::Utils::CommonEdge(&m, 44, 84), NullID);
	ASSERT_EQ(math::Utils::CommonEdge(&m, 56, 97), NullID);
	ASSERT_EQ(math::Utils::CommonEdge(&m, 45, 38), NullID);
	ASSERT_EQ(math::Utils::CommonEdge(&m, 87, 11), NullID);
	ASSERT_EQ(math::Utils::CommonEdge(&m, 20, 16), NullID);
	ASSERT_EQ(math::Utils::CommonEdge(&m, 24, 25), NullID);
	ASSERT_EQ(math::Utils::CommonEdge(&m, 14, 12), NullID);

	/*
	TCellID id1 = 86;
	TCellID id2 = 51;
	std::cout << math::Utils::CommonEdge(&m, id1, id2) << std::endl;
	Node n1 = m.get<Node>(id1);
	Node n2 = m.get<Node>(id2);
	std::vector<Edge> adj_edges_1 = n1.get<Edge>() ;
	std::vector<Edge> adj_edges_2 = n2.get<Edge>();
	std::cout << "---------" << std::endl;
	for (auto e:adj_edges_1){
		std::cout << e.id() << std::endl;
	}
	std::cout << "---------" << std::endl;
	for (auto e:adj_edges_2){
		std::cout << e.id() << std::endl;
	}
	 */

	IGMeshIOService ioService_geom(&m);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("Utils_Test2.vtk");

}


TEST(UtilsTestClass, Utils_Test3)
{
	// Test de la méthode math::Utils::MeshCleaner

	// WE READ
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::E | gmds::N2E | gmds::N2F | gmds::F2N | gmds::E2N | gmds::F2E | gmds::E2F));

	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir + "/Aero/Poubelle/Carre.vtk";

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N | gmds::F);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	math::Point p0(-1,-1,0);
	Node n0 = m.newNode(p0);

	math::Point p1(-60,-12,0);
	Node n1 = m.newNode(p1);

	math::Point p2(3,10,0);
	Node n2 = m.newNode(p2);

	math::Point p3(-13,4,0);
	Node n3 = m.newNode(p3);

	math::Utils::MeshCleaner(&m);

	for (auto n_id:m.nodes()){
		ASSERT_NE(n_id, n0.id());
		ASSERT_NE(n_id, n1.id());
		ASSERT_NE(n_id, n2.id());
		ASSERT_NE(n_id, n3.id());
	}

	IGMeshIOService ioService_geom(&m);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("Utils_Test3.vtk");

}


TEST(ClaireTestClass, Utils_AdjacentNodes)
{
	// Test
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::E | gmds::N2E | gmds::N2F | gmds::F2N | gmds::E2N | gmds::F2E | gmds::E2F));

	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir + "/Aero/Poubelle/Carre.vtk";

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N | gmds::F);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	TCellID n_id = 1;
	Node n = m.get<Node>(n_id);
	std::vector<Node> adj_nodes = math::Utils::AdjacentNodes(&m, n);
	ASSERT_EQ(adj_nodes.size(),3);
	ASSERT_EQ(adj_nodes[0].id(),11);
	ASSERT_EQ(adj_nodes[1].id(),91);
	ASSERT_EQ(adj_nodes[2].id(),10);

	n_id = 4;
	n = m.get<Node>(n_id);
	adj_nodes = math::Utils::AdjacentNodes(&m, n);
	ASSERT_EQ(adj_nodes.size(),4);
	ASSERT_EQ(adj_nodes[0].id(),5);
	ASSERT_EQ(adj_nodes[1].id(),79);
	ASSERT_EQ(adj_nodes[2].id(),0);
	ASSERT_EQ(adj_nodes[3].id(),90);

	n_id = 92;
	n = m.get<Node>(n_id);
	adj_nodes = math::Utils::AdjacentNodes(&m, n);
	ASSERT_EQ(adj_nodes.size(), 6);
	ASSERT_EQ(adj_nodes[0].id(), 69);
	ASSERT_EQ(adj_nodes[1].id(), 70);
	ASSERT_EQ(adj_nodes[2].id(), 71);
	ASSERT_EQ(adj_nodes[3].id(), 50);
	ASSERT_EQ(adj_nodes[4].id(), 86);
	ASSERT_EQ(adj_nodes[5].id(), 84);
	//std::cout << adj_nodes[5].id() << std::endl;

	IGMeshIOService ioService_geom(&m);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("Utils_Test4.vtk");

}


TEST(ClaireTestClass, Utils_BuildMesh2DFromBlocking2D)
{
	// Test
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::E | gmds::N2E | gmds::N2F | gmds::F2N | gmds::E2N | gmds::F2E | gmds::E2F));

	Blocking2D b;
	Node n1 = b.newBlockCorner(0,0);
	Node n2 = b.newBlockCorner(1,0);
	Node n3 = b.newBlockCorner(1,1);
	Node n4=  b.newBlockCorner(0,1);

	Blocking2D::Block b0 = b.newBlock(n1,n2,n3,n4);

	Node n5 = b.newBlockCorner(2,0,0);
	Node n6 = b.newBlockCorner(2,1.5,0);
	Blocking2D::Block b1 = b.newBlock(n2,n5,n6,n3);

	ASSERT_EQ(b0.id(), b.block(0).id());
	b0.setNbDiscretizationI(11);
	b0.setNbDiscretizationJ(11);
	b1.setNbDiscretizationI(11);
	b1.setNbDiscretizationJ(11);
	b.initializeGridPoints();

	int mark_block_nodes = m.newMark<Node>();
	int mark_first_layer = m.newMark<Node>();
	int mark_farfield_nodes = m.newMark<Node>();
	math::Utils::BuildMesh2DFromBlocking2D(&b, &m, mark_block_nodes, mark_first_layer, mark_farfield_nodes);
	m.unmarkAll<Node>(mark_block_nodes);
	m.freeMark<Node>(mark_block_nodes);
	m.unmarkAll<Node>(mark_first_layer);
	m.freeMark<Node>(mark_first_layer);
	m.unmarkAll<Node>(mark_farfield_nodes);
	m.freeMark<Node>(mark_farfield_nodes);


	ASSERT_EQ(m.getNbFaces(), 200);
	ASSERT_EQ(m.getNbNodes(), 231);

	IGMeshIOService ioService_geom(&b);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("Utils_BuildMesh2DFromBlocking2D_Blocking.vtk");

	IGMeshIOService ioService_geom_mesh(&m);
	VTKWriter writer_geom_mesh(&ioService_geom_mesh);
	writer_geom_mesh.setCellOptions(N|F);
	writer_geom_mesh.setDataOptions(N|F);
	writer_geom_mesh.write("Utils_BuildMesh2DFromBlocking2D_Mesh.vtk");

}


TEST(ClaireTestClass, Utils_WeightedPointOnBranch)
{
	{
		math::Point A({0.0, 0.0, 0.0});
		math::Point B({1.0, 0.0, 0.0});
		math::Point C({3.0, 0.0, 0.0});
		math::Point D = math::Utils::WeightedPointOnBranch(A, B, C, 0.5);
		ASSERT_FLOAT_EQ(D.X(), 1.5);
		ASSERT_FLOAT_EQ(D.Y(), 0.0);
	}

	{
		math::Point A({0.0, 0.0, 0.0});
		math::Point B({1.0, 0.0, 0.0});
		math::Point C({4.0, 0.0, 0.0});
		math::Point D = math::Utils::WeightedPointOnBranch(A, B, C, 0.75);
		ASSERT_FLOAT_EQ(D.X(), 3.0);
		ASSERT_FLOAT_EQ(D.Y(), 0.0);
	}

	{
		math::Point A({0.0, 0.0, 0.0});
		math::Point B({1.0, 0.0, 0.0});
		math::Point C({1.0, 3.0, 0.0});
		math::Point D = math::Utils::WeightedPointOnBranch(A, B, C, 0.75);
		ASSERT_FLOAT_EQ(D.X(), 1.0);
		ASSERT_FLOAT_EQ(D.Y(), 2.0);
	}

}


TEST(ClaireTestClass, Utils_CreateQuadAndConnectivitiesN2F)
{
	// Test
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::E | gmds::N2E | gmds::N2F | gmds::F2N | gmds::E2N | gmds::F2E | gmds::E2F));

	Node n0 = m.newNode({0,0,0});
	Node n1 = m.newNode({0,1,0});
	Node n2 = m.newNode({1,1,0});
	Node n3 = m.newNode({1,0,0});

	ASSERT_EQ(m.getNbFaces(), 0);

	TCellID f_id = math::Utils::CreateQuadAndConnectivitiesN2F(&m, n0.id(), n1.id(), n2.id(), n3.id());
	ASSERT_EQ(m.getNbFaces(), 1);

	std::vector<Face> n0_faces = n0.get<Face>();
	ASSERT_EQ(n0_faces.size(), 1);

}