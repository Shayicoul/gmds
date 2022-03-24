/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/ig/Blocking2D.h>
#include <iostream>

/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(Blocking2DTestSuite, test_blocking2D_1)
{
    Blocking2D m;
    Node n1 = m.newBlockCorner(0,0);
    Node n2 = m.newBlockCorner(1,0);
    Node n3 = m.newBlockCorner(1,1);
    Node n4=  m.newBlockCorner(0,1);

    Blocking2D::Block b0 = m.newBlock(n1,n2,n3,n4);

    Node n5 = m.newBlockCorner(2,0,0);
    Node n6 = m.newBlockCorner(2,1.5,0);
    Blocking2D::Block b1 = m.newBlock(n2,n5,n6,n3);

    ASSERT_EQ(b0.id(), m.block(0).id());
    b0.setNbDiscretizationI(11);
    b0.setNbDiscretizationJ(11);
    b1.setNbDiscretizationI(11);
    b1.setNbDiscretizationJ(11);

    ASSERT_EQ(n1.id(),b0.origin());
    ASSERT_FLOAT_EQ(b0.getUnitVectorI().X(),1.0);
    ASSERT_FLOAT_EQ(b0.getUnitVectorI().Y(),0.0);
    ASSERT_FLOAT_EQ(b0.getUnitVectorI().Z(),0.0);
    ASSERT_FLOAT_EQ(b0.getUnitVectorJ().X(),0.0);
    ASSERT_FLOAT_EQ(b0.getUnitVectorJ().Y(),1.0);
    ASSERT_FLOAT_EQ(b0.getUnitVectorJ().Z(),0.0);

    ASSERT_EQ(b0.getNbDiscretizationI(),11);
    ASSERT_EQ(b0.getNbDiscretizationJ(),11);

    m.initializeGridPoints();
    b0 = m.block(0);
    math::Point p00 = b0(0,0).point();
    ASSERT_FLOAT_EQ(p00.X(),0);
    ASSERT_FLOAT_EQ(p00.Y(),0);

    math::Point p04 = b0(0,4).point();
    ASSERT_FLOAT_EQ(p04.X(),0);
    ASSERT_FLOAT_EQ(p04.Y(),0.4);

    math::Point p31 = b0(3,1).point();
    ASSERT_FLOAT_EQ(p31.X(),0.3);
    ASSERT_FLOAT_EQ(p31.Y(),0.1);

    math::Point p76 = b0(7,6).point();
    ASSERT_FLOAT_EQ(p76.X(),0.7);
    ASSERT_FLOAT_EQ(p76.Y(),0.6);

}
/*----------------------------------------------------------------------------*/
TEST(Blocking2DTestSuite, test_blocking2D_2)
{
	Mesh mesh(MeshModel(DIM3|F|E|N|F2N|E2N|F2E|E2F|N2E));

	std::vector<TCellID> nodes_for_face;

	Node n1 = mesh.newNode(0,0);
	Node n2 = mesh.newNode(1,0);
	Node n3 = mesh.newNode(1,1);
	Node n4 = mesh.newNode(0,1);

	nodes_for_face.push_back(n1.id());
	nodes_for_face.push_back(n2.id());
	nodes_for_face.push_back(n3.id());
	nodes_for_face.push_back(n4.id());

	mesh.newFace(nodes_for_face);
	nodes_for_face.clear();

	Node n5 = mesh.newNode(2,0,0);
	Node n6 = mesh.newNode(2,1.5,0);

	nodes_for_face.push_back(n2.id());
	nodes_for_face.push_back(n5.id());
	nodes_for_face.push_back(n6.id());
	nodes_for_face.push_back(n3.id());

	mesh.newFace(nodes_for_face);

	Blocking2D blocking(mesh);

	int nbBlocks = blocking.allBlocks().size();

	ASSERT_EQ(nbBlocks, 2);

	Blocking2D::Block b0 = blocking.block(0);

	ASSERT_EQ(n1.id(),b0.origin());
	ASSERT_FLOAT_EQ(b0.getUnitVectorI().X(),1.0);
	ASSERT_FLOAT_EQ(b0.getUnitVectorI().Y(),0.0);
	ASSERT_FLOAT_EQ(b0.getUnitVectorI().Z(),0.0);
	ASSERT_FLOAT_EQ(b0.getUnitVectorJ().X(),0.0);
	ASSERT_FLOAT_EQ(b0.getUnitVectorJ().Y(),1.0);
	ASSERT_FLOAT_EQ(b0.getUnitVectorJ().Z(),0.0);

	ASSERT_EQ(b0.getNbDiscretizationI(),11);
	ASSERT_EQ(b0.getNbDiscretizationJ(),11);

	blocking.initializeGridPoints();

	b0 = blocking.block(0);

	math::Point p00 = b0(0,0).point();
	ASSERT_FLOAT_EQ(p00.X(),0);
	ASSERT_FLOAT_EQ(p00.Y(),0);

	math::Point p04 = b0(0,4).point();
	ASSERT_FLOAT_EQ(p04.X(),0);
	ASSERT_FLOAT_EQ(p04.Y(),0.4);

	math::Point p31 = b0(3,1).point();
	ASSERT_FLOAT_EQ(p31.X(),0.3);
	ASSERT_FLOAT_EQ(p31.Y(),0.1);

	math::Point p76 = b0(7,6).point();
	ASSERT_FLOAT_EQ(p76.X(),0.7);
	ASSERT_FLOAT_EQ(p76.Y(),0.6);

}