//
// Created by Paul Bourmaud on 23/03/2022.
//
#include <cstdlib>
#include <iostream>
#include <ctime>
#include "gmds/paul/Grid.h"
#include "gmds/ig/Mesh.h"
#include "gmds/igalgo/VolFracComputation.h"
#include "gmds/io/IGMeshIOService.h"
#include "gmds/io/VTKWriter.h"
#include "gmds/paul/Actions_Agent.h"
#include "gmds/paul/Tools.h"

using namespace gmds;

int main(){

	//Create the agent environment
	std::cout<<"=== Test Grid ==="<<std::endl;
	Grid g(2,4);
	std::cout<<g.getX()<<" - "<<g.getY()<<std::endl;

	//Create triangle reference for BuildGridAround
	Mesh m(MeshModel(DIM3|F|N|F2N|N2F));
	Node n0 = m.newNode(math::Point(0,0,0));
	Node n1 = m.newNode(math::Point(0,3,0));
	Node n2 = m.newNode(math::Point(3,3,0));
	/*Variable<double> *v = m.newVariable<double,GMDS_NODE>("val");
	v->set(n0.id(),10);
	v->set(n1.id(),40);*/
	m.newTriangle(n0,n1,n2);
	for(auto face_id:m.faces()){
		Face f = m.get<Face>(face_id);
		std::vector<Node> f_nodes = f.get<Node>();
		for(auto n:f_nodes){
			std::cout<<n<<std::endl;
		}
	}

	std::cout<<"Creation triangle target"<<std::endl;
	//Create triangle target
	Mesh mT(MeshModel(DIM3|F|N|F2N|N2F));
	Node n0b = mT.newNode(math::Point(0,0,0));
	Node n1b = mT.newNode(math::Point(0,3,0));
	Node n2b = mT.newNode(math::Point(3,3,0));
	Variable<double> *vb = mT.newVariable<double,GMDS_NODE>("valTarget");
	vb->set(n0b.id(),10);
	vb->set(n1b.id(),30);
	mT.newTriangle(n0b,n1b,n2b);
	Variable<int> *titi = mT.newVariable<int,gmds::GMDS_FACE>("titi");
	titi->set(0,1);

	for (auto face_id:mT.faces()){
		Face fT=mT.get<Face>(face_id);
		std::vector<Node> f_nodes = fT.get<Node>();
		for (auto n : f_nodes){
			std::cout<<n.id()<<std::endl;
			vb->set(n.id(),40); //attribution val aux noeuds
		}
	}

	//test VolFrac function create by nicolas legoff (pas possible pour le moment parce que pas bon type
	// de maillage en entrée, besoin maillage tétraédrique je crois)
	/*std::cout<<"=== Valeur Fraction ==="<<std::endl;
	Variable<double> *volFrac = mT.newVariable<double,GMDS_FACE>("Valeur Fraction");
	gmds::volfraccomputation_2d(&mT,&m,volFrac);
	std::cout<<volFrac<<std::endl;*/

	//Create the grid around the target (triangle in this case)
	std::cout<<"=== Test Grid Builder Around ==="<<std::endl;
	GridBuilderAround gba(&m,2);
	int nbNodes = 5;
	double sizeStep = 5;
	gba.executeGrid2D(nbNodes,sizeStep,nbNodes,sizeStep);
	Actions action(&gba);
	//Variable<int> *activate = m.newVariable<int,GMDS_FACE>("exist");

	for(auto face_id:gba.m_mesh.faces()){
		Face f =gba.m_mesh.get<Face>(face_id);
		std::vector<Node>f_nodes = f.get<Node>();
		//std::cout<<"Avant bool activate"<<std::endl;
		//std::cout<<activate<<std::endl;
		//std::cout<<"Apres bool activate"<<std::endl;
		for(auto n:f_nodes){
			//std::cout<<n<<std::endl;
			//std::cout<<"print de f.id()"<<std::endl;
			//std::cout<<f.id()<<std::endl;
			//std::cout<<"================================================="<<std::endl;
			//gba.activate->set(f.id(), 1);
		}
	}

	/*
	for(int id=0;id <9;id++) {
		action.executeDeleteFace(id);
	}*/

	Tools tool(&gba);

	//======================================== TEST CUT NODE GLIDE =============================================
	/*
	gba.m_mesh.get<Node>(4).X()=0.5;
	gba.m_mesh.get<Node>(4).Y()=4.5;

	action.executeCutEdge(gba.m_mesh.get<Node>(4),gba.m_mesh.get<Node>(9));
	*/
	/*
	action.executeCutEdge(tool.g_grid.m_mesh.get<Node>(18),tool.g_grid.m_mesh.get<Node>(13));
	action.executeCutEdge(tool.g_grid.m_mesh.get<Node>(1),tool.g_grid.m_mesh.get<Node>(2));
	action.executeCutEdge(tool.g_grid.m_mesh.get<Node>(3),tool.g_grid.m_mesh.get<Node>(2));
	action.executeCutEdge(tool.g_grid.m_mesh.get<Node>(14),tool.g_grid.m_mesh.get<Node>(25));
	action.executeCutEdge(tool.g_grid.m_mesh.get<Node>(19),tool.g_grid.m_mesh.get<Node>(25));
	action.executeCutEdge(tool.g_grid.m_mesh.get<Node>(4),tool.g_grid.m_mesh.get<Node>(3));
	*/

	//====================== TEST =====================
	std::srand(std::time(nullptr));
	for (int j=0; j<21;j++) {

		for (int i = 0; i < 1000; i++) {
			int i1 = rand() % 82;
			int i2 = rand() % 82;
			Node node1 = gba.m_mesh.get<Node>(i1);
			Node node2 = gba.m_mesh.get<Node>(i2);
			action.executeCutEdge(node1, node2);
		}

		for (int i; i < 18; i++) {
			int idFace = rand() % 50;
			action.executeDeleteFace(idFace);
		}
		IGMeshIOService ioService23(&m);
		VTKWriter vtkWriter23(&ioService23);
		vtkWriter23.setCellOptions(gmds::N|gmds::F);
		vtkWriter23.setDataOptions(gmds::N|gmds::F);
		std::string s = std::to_string(j);
		s.push_back('.');
		s.push_back('v');
		s.push_back('t');
		s.push_back('k');
		vtkWriter23.write(s);

	}





// Save Triangle Generation
	IGMeshIOService ioService(&mT);
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("lulu.vtk");


// Save GridBuilder Generation
	IGMeshIOService ioService2(&m);
	VTKWriter vtkWriter2(&ioService2);
	vtkWriter2.setCellOptions(gmds::N|gmds::F);
	vtkWriter2.setDataOptions(gmds::N|gmds::F);
	vtkWriter2.write("lili.vtk");

	exit(3);

}