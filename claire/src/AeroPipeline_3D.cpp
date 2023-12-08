//
// Created by rochec on 09/02/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/Utils.h>
#include <gmds/claire/AeroPipeline_3D.h>
#include <gmds/claire/AeroBoundaries_3D.h>
#include <gmds/claire/AeroExtrusion_3D.h>
#include <gmds/claire/LevelSetCombined.h>
#include <gmds/claire/LeastSquaresGradientComputation.h>
#include <gmds/claire/IntervalAssignment_3D.h>
#include <gmds/claire/MFEMMeshWriter.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <gmds/smoothy/LaplacianSmoother.h>
//#include <gmds/blocking/CurvedBlocking.h>
#include <gmds/math/BezierHex.h>
#include <gmds/claire/AdvectedPointRK4_3D.h>
#include <gmds/claire/TransfiniteInterpolation_3D.h>
#include <iostream>
#include <chrono>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

AeroPipeline_3D::AeroPipeline_3D(std::string Aparams, std::string &Aworking_dir) :
	AbstractAeroPipeline(Aparams, Aworking_dir),
  m_linker_BG(new cad::GeomMeshLinker())
{
	m_meshTet = new Mesh(gmds::MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R |
	                                  F2E | E2F | R2E | E2R | N2R | N2F | N2E));
	m_meshHex = new Mesh(gmds::MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R |
	                                     F2E | E2F | R2E | E2R | N2R | N2F | N2E));
	m_couche_id = m_meshHex->newVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	m_meshHex->getOrCreateVariable<int, GMDS_FACE>("GMDS_FACE_Couche_Id");
	m_Bnd = new AeroBoundaries_3D(m_meshTet) ;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
AbstractAeroPipeline::STATUS
AeroPipeline_3D::execute(){

	clock_t t_start, t_end;

	LectureMaillage();

	t_start = clock();
	m_Bnd->execute();
	t_end = clock();
	std::cout << "........................................ temps : " << 1.0*double(t_end-t_start)/CLOCKS_PER_SEC << "s" << std::endl;
	//m_manager.initAndLinkFrom3DMesh(&m_mTetra,&m_linker_TG);
	PreTraitementMeshTet();

	// Calcul du level set
	std::cout << "-> Calcul des Level Sets" << std::endl;
	t_start = clock();
	m_meshTet->newVariable<double,GMDS_NODE>("GMDS_Distance");
	m_meshTet->newVariable<double,GMDS_NODE>("GMDS_Distance_Int");
	m_meshTet->newVariable<double,GMDS_NODE>("GMDS_Distance_Out");
	LevelSetCombined lsCombined(m_meshTet,
	                            m_Bnd->getMarkParoi(),
	                            m_Bnd->getMarkAmont(),
	                            m_meshTet->getVariable<double,GMDS_NODE>("GMDS_Distance"),
	                            m_meshTet->getVariable<double,GMDS_NODE>("GMDS_Distance_Int"),
	                            m_meshTet->getVariable<double,GMDS_NODE>("GMDS_Distance_Out"));
	lsCombined.execute();
	t_end = clock();
	std::cout << "........................................ temps : " << 1.0*double(t_end-t_start)/CLOCKS_PER_SEC << "s" << std::endl;


	// Calcul du gradient du champ de Level Set
	std::cout << "-> Calcul du gradient du champ des Level Sets" << std::endl;
	t_start = clock();
	ComputeVectorFieldForExtrusion();
	t_end = clock();
	std::cout << "........................................ temps : " << 1.0*double(t_end-t_start)/CLOCKS_PER_SEC << "s" << std::endl;

	GeometrySurfaceBlockingGeneration();	// Generate the blocking of the geometry surface.
	SurfaceBlockingClassification();			// Link the surface blocking to the geometry.

	// Write the surface blocking to check the classification
	if (m_params.with_debug_files)
	{
		gmds::IGMeshIOService ioService(m_meshHex);
		gmds::VTKWriter vtkWriter(&ioService);
		vtkWriter.setCellOptions(gmds::N | gmds::F);
		vtkWriter.setDataOptions(gmds::N | gmds::F);
		vtkWriter.write("Surface_3D_GEOM_CLASSIFICATION.vtk");
	}

	// Extrusion
	std::cout << "-> Extrusion" << std::endl;
	Variable<math::Vector3d>* var_VectorsForExtrusion = m_meshTet->getOrCreateVariable<math::Vector3d, GMDS_NODE>("VectorField_Extrusion");
	t_start = clock();
	AeroExtrusion_3D aero_extrusion(m_meshTet,
	                                m_meshHex,
	                                m_params,
	                                m_meshTet->getVariable<double,GMDS_NODE>("GMDS_Distance"),
	                                var_VectorsForExtrusion);
	aero_extrusion.execute();
	t_end = clock();
	std::cout << "........................................ temps : " << 1.0*double(t_end-t_start)/CLOCKS_PER_SEC << "s" << std::endl;
	std::cout << " " << std::endl;

	BlockingGeometricClassification();	// Geometric classification of the final blocking, according to the geometric classification of the input tet mesh
	//smoothy::LaplacianSmoother smoother(m_linker_HG);
	//smoother.smoothCurves();
	//smoother.smoothSurfaces();
	//smoother.smoothVolumes(2);


	// Smoothy
	int nbr_iter_smoothing(3);
	for (int i=0;i<nbr_iter_smoothing;i++)
	{
		std::map<TCellID, math::Point> new_pos;
		// Smooth node position in volume
		for (auto n_id:m_meshHex->nodes())
		{
			if (m_couche_id->value(n_id) > 1
			    && m_couche_id->value(n_id) < m_params.nbr_couches)
			{
				double theta = 0.3*m_couche_id->value(n_id)/m_params.nbr_couches ;
				Node n = m_meshHex->get<Node>(n_id);
				std::vector<Edge> n_edges = n.get<Edge>();
				//math::Point p(0.0, 0.0, 0.0);
				math::Point p(n.point());
				for (auto const& e:n_edges)
				{
					Node n_opp = e.getOppositeNode(n_id);
					p = p + n_opp.point();
				}
				p.setX(p.X()/(double(n_edges.size())+1));
				p.setY(p.Y()/(double(n_edges.size())+1));
				p.setZ(p.Z()/(double(n_edges.size())+1));
				new_pos[n_id] = theta*p+(1.0-theta)*n.point();
			}
		}

		// Smooth node position on last front (could be improve by using Front3D class)
		for (auto n_id:m_meshHex->nodes())
		{
			if (m_couche_id->value(n_id) == m_params.nbr_couches)
			{
				double theta = 0.3 ;
				Node n = m_meshHex->get<Node>(n_id);
				std::vector<Edge> n_edges = n.get<Edge>();
				math::Point p(0.0, 0.0, 0.0);
				int nbr_adj(0);
				for (auto const& e:n_edges)
				{
					Node n_opp = e.getOppositeNode(n_id);
					if (m_couche_id->value(n_opp.id()) == m_params.nbr_couches)
					{
						p = p + e.getOppositeNode(n_id).point();
						nbr_adj++;
					}
				}
				p.setX(p.X()/nbr_adj);
				p.setY(p.Y()/nbr_adj);
				p.setZ(p.Z()/nbr_adj);
				new_pos[n_id] = theta*p+(1.0-theta)*n.point();
			}
		}

		// Update positions
		for (auto n_update:new_pos)
		{
			Node n = m_meshHex->get<Node>(n_update.first);
			if (m_couche_id->value(n.id()) == m_params.nbr_couches)
			{
				// Re-projection on exterior front
				int geom_id = m_linker_HG->getGeomId<Node>(n.id()) ;
				int geom_dim = m_linker_HG->getGeomDim<Node>(n.id()) ;
				cad::GeomSurface* surface = m_manager->getSurface(geom_id);
				surface->project(n_update.second);
			}
			n.setPoint(n_update.second);
		}

	}


	// Convert to CurvedBlocking structure
	/*
	gmds::blocking::CurvedBlocking blocking(m_manager) ;
	blocking.init_from_mesh(*m_meshHex);

	gmds::Mesh m_out(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E | gmds::F | gmds::R | gmds::E2N | gmds::F2N | gmds::R2N));
	blocking.convert_to_mesh(m_out);

	gmds::IGMeshIOService ioService(&m_out);
	gmds::VTKWriter writer(&ioService);
	writer.setCellOptions(gmds::N | gmds::R);
	writer.setDataOptions(gmds::N | gmds::R);
	writer.write("TEST_CURVED.vtk");
	 */

	// Init the Blocking3D from the hex mesh
	std::cout << "-> Init the Blocking3D structure from the Hex mesh" << std::endl;
	t_start = clock();
	initBlocking3DfromMesh();
	updateLayerValues();
	t_end = clock();
	std::cout << "........................................ temps : " << 1.0*double(t_end-t_start)/CLOCKS_PER_SEC << "s" << std::endl;
	std::cout << " " << std::endl;


	// Stat of the blocking
	std::cout << "=============================================" << std::endl;
	std::cout << "			BLOCKING STATISTICS:" << std::endl;
	std::cout << "=============================================" << std::endl;
	std::cout << "|| Number of Blocks: " << m_meshHex->getNbRegions() << std::endl;
	std::cout << "|| Number of Block Faces: " << m_meshHex->getNbFaces() << std::endl;
	std::cout << "|| Number of Block Edges: " << m_meshHex->getNbEdges() << std::endl;
	std::cout << "|| Number of Block Nodes: " << m_meshHex->getNbNodes() << std::endl;
	std::cout << "=============================================" << std::endl;

	// Write the final mesh.
	EcritureMaillage();

	return AbstractAeroPipeline::SUCCESS;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void
AeroPipeline_3D::LectureMaillage(){

	// Lecture du maillage
	std::cout << "-> Lecture du maillage ..." << std::endl;

	gmds::IGMeshIOService ioService(m_meshTet);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::R);
	vtkReader.read(m_params.input_file);

	gmds::MeshDoctor doctor(m_meshTet);
	doctor.buildFacesAndR2F();
	doctor.buildEdgesAndX2E();
	doctor.updateUpwardConnectivity();

	// Erase node connected to nothing
	math::Utils::MeshCleaner(m_meshTet);

}
/*------------------------------------------------------------------------*/
void
AeroPipeline_3D::EcritureMaillage(){

	std::cout << "-> Ecriture du maillage ..." << std::endl;

	// Ecriture du maillage généré
	gmds::IGMeshIOService ioService(m_meshHex);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write(m_params.output_file);

	// Write the edges of the blocking to vizualize the edge
	// discretization
	gmds::IGMeshIOService ioService_edges(&m_Blocking3D);
	gmds::VTKWriter vtkWriter_edges(&ioService_edges);
	vtkWriter_edges.setCellOptions(gmds::N|gmds::E);
	vtkWriter_edges.setDataOptions(gmds::N|gmds::E);
	vtkWriter_edges.write("AeroPipeline3D_EdgesDiscretization.vtk");

	// Write the initial tet mesh (with fields computed on it)
	ioService = IGMeshIOService(m_meshTet);
	gmds::VTKWriter vtkWriter2(&ioService);
	vtkWriter2.setCellOptions(gmds::N|gmds::R);
	vtkWriter2.setDataOptions(gmds::N|gmds::R);
	vtkWriter2.write("AeroPipeline3D_Tetra.vtk");

	// Write the Blocking
	gmds::IGMeshIOService ioService_blocking(&m_Blocking3D);
	gmds::VTKWriter vtkWriter_blocking(&ioService_blocking);
	vtkWriter_blocking.setCellOptions(gmds::N|gmds::R);
	vtkWriter_blocking.setDataOptions(gmds::N|gmds::R);
	vtkWriter_blocking.write("AeroPipeline3D_Blocking.vtk");

	// Write the Blocking
	gmds::IGMeshIOService ioService_ctrlpts(&m_CtrlPts);
	gmds::VTKWriter vtkWriter_ctrlpts(&ioService_ctrlpts);
	vtkWriter_ctrlpts.setCellOptions(gmds::N|gmds::R);
	vtkWriter_ctrlpts.setDataOptions(gmds::N|gmds::R);
	vtkWriter_ctrlpts.write("AeroPipeline3D_ControlPoints.vtk");


	// TEST FOR MFEM
	/*
	if (m_params.with_debug_files)
	{
		MeshDoctor doc(m_meshHex);
		doc.buildFacesAndR2F();
		doc.updateUpwardConnectivity();
		doc.orient2DFaces();

		std::cout << "Nbr elements avant: " << m_meshHex->getNbRegions() << std::endl;
		for (auto r_id:m_meshHex->regions())
		{
			math::Utils::orientRegion(m_meshHex, m_meshHex->get<Region>(r_id));
		}
		std::cout << "Nbr elements après: " << m_meshHex->getNbRegions() << std::endl;

		gmds::IGMeshIOService ioService_MFEM(m_meshHex);
		gmds::VTKWriter vtkWriter_MFEM(&ioService_MFEM);
		vtkWriter_MFEM.setCellOptions(gmds::N | gmds::R);
		vtkWriter_MFEM.setDataOptions(gmds::N | gmds::R);
		std::string dir(".");
		vtkWriter_MFEM.write("Apollo_3D_Blocks.vtk");
	}
	 */

	{
		std::cout << "MFEM writing..." << std::endl;
		FastLocalize fl = FastLocalize(m_meshTet);
		Variable<int>* var_couche = m_meshHex->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");
		Variable<math::Vector3d>* var_VectorsForExtrusion = m_meshTet->getOrCreateVariable<math::Vector3d, GMDS_NODE>("VectorField_Extrusion");
		std::map<TCellID,TCellID> new_nodes;
		for (auto n_id:m_meshHex->nodes())
		{
			if (var_couche->value(n_id)==0)
			{
				Cell::Data data = fl.find(m_meshHex->get<Node>(n_id).point());
				//Node n_closest = m_meshHex->get<Node>(data.id);
				math::Vector3d v = var_VectorsForExtrusion->value(data.id).normalize();
				math::Point p_new = m_meshHex->get<Node>(n_id).point() + -1*(m_params.delta_cl/2.0)*v;
				Node n_new = m_meshHex->newNode(p_new);
				new_nodes[n_id] = n_new.id();
				var_couche->set(n_new.id(), -1);
			}
		}
		for (auto f_id:m_meshHex->faces())
		{
			Face f = m_meshHex->get<Face>(f_id);
			std::vector<Node> f_nodes = f.get<Node>();
			if (var_couche->value(f_nodes[0].id())==0
			    && var_couche->value(f_nodes[1].id())==0
			    && var_couche->value(f_nodes[2].id())==0
			    && var_couche->value(f_nodes[3].id())==0)
			{
				Region r = m_meshHex->newHex(f_nodes[0].id(), f_nodes[1].id(), f_nodes[2].id(), f_nodes[3].id(),
				                  new_nodes[f_nodes[0].id()], new_nodes[f_nodes[1].id()], new_nodes[f_nodes[2].id()], new_nodes[f_nodes[3].id()]);
				var_couche->set(r.id(), -1);
			}
		}

		/*
		math::Point p;
		int compteur(0);
		for (auto n_id:m_meshHex->nodes())
		{
			if (var_couche->value(n_id)==0)
			{
				Node n = m_meshHex->get<Node>(n_id);
				p = p + n.point();
				compteur++;
			}
		}
		if (compteur != 0)
		{
			p.setX(p.X() / compteur);
			p.setY(p.Y() / compteur);
		}
		Node n_new = m_meshHex->newNode(p);
		for (auto f_id:m_meshHex->faces())
		{
			Face f = m_meshHex->get<Face>(f_id);
			std::vector<Node> f_nodes = f.get<Node>();
			if (var_couche->value(f_nodes[0].id())==0
			    && var_couche->value(f_nodes[1].id())==0
			    && var_couche->value(f_nodes[2].id())==0
			    && var_couche->value(f_nodes[3].id())==0)
			{
				m_meshHex->newPyramid(f_nodes[0], f_nodes[1], f_nodes[2], f_nodes[3], n_new);
			}
		}
		*/

		if (m_params.with_debug_files)
		{
			MeshDoctor doc(m_meshHex);
			doc.buildFacesAndR2F();
			doc.updateUpwardConnectivity();
			doc.orient2DFaces();

			std::cout << "Nbr elements avant: " << m_meshHex->getNbRegions() << std::endl;
			for (auto r_id:m_meshHex->regions())
			{
				math::Utils::orientRegion(m_meshHex, m_meshHex->get<Region>(r_id));
			}
			std::cout << "Nbr elements après: " << m_meshHex->getNbRegions() << std::endl;

			gmds::IGMeshIOService ioService_MFEM(m_meshHex);
			gmds::VTKWriter vtkWriter_MFEM(&ioService_MFEM);
			vtkWriter_MFEM.setCellOptions(gmds::N | gmds::R);
			vtkWriter_MFEM.setDataOptions(gmds::N | gmds::R);
			std::string dir(".");
			vtkWriter_MFEM.write("AeroPipeline3D_Blocking_withPyramids.vtk");
		}

		MFEMMeshWriter mfemwriter = MFEMMeshWriter(m_meshHex, "AeroPipeline3D_Blocking_toFit");
		mfemwriter.execute();
	}

	// Write the Blocking3D as a mesh
	m_meshHex->clear();
	math::Utils::BuildMesh3DFromBlocking3D(&m_Blocking3D, m_meshHex);
	gmds::IGMeshIOService ioService_mesh3D(m_meshHex);
	gmds::VTKWriter vtkWriter_mesh3D(&ioService_mesh3D);
	vtkWriter_mesh3D.setCellOptions(gmds::N|gmds::R);
	vtkWriter_mesh3D.setDataOptions(gmds::N|gmds::R);
	vtkWriter_mesh3D.write("AeroPipeline3D_HexMesh.vtk");

	// Stat of the blocking
	std::cout << "=============================================" << std::endl;
	std::cout << "	FINAL HEX MESH STATISTICS:" << std::endl;
	std::cout << "=============================================" << std::endl;
	std::cout << "|| Number of Hex: " << m_meshHex->getNbRegions() << std::endl;
	//std::cout << "|| Number of Block Faces: " << m_meshHex->getNbFaces() << std::endl;
	//std::cout << "|| Number of Block Edges: " << m_meshHex->getNbEdges() << std::endl;
	std::cout << "|| Number of Nodes: " << m_meshHex->getNbNodes() << std::endl;
	std::cout << "=============================================" << std::endl;

	// Write a trick file to visualize the curved block edges
	m_meshHex->clear();
	math::Utils::CurveBlockEdgesReveal3D(&m_CtrlPts, m_meshHex, 10);
	gmds::IGMeshIOService ioService_visuCurved(m_meshHex);
	gmds::VTKWriter vtkWriter_visuCurved(&ioService_visuCurved);
	vtkWriter_visuCurved.setCellOptions(gmds::N|gmds::F);
	vtkWriter_visuCurved.setDataOptions(gmds::N|gmds::F);
	vtkWriter_visuCurved.write("AeroPipeline3D_CurvedBlockEdges.vtk");

}
/*------------------------------------------------------------------------*/
void
AeroPipeline_3D::GeometrySurfaceBlockingGeneration()
{
	if (m_params.block_surface_3D==0)
	{
		// Lecture du maillage
		std::cout << "-> Read the input quad surface blocking ..." << std::endl;

		gmds::IGMeshIOService ioService(m_meshHex);
		gmds::VTKReader vtkReader(&ioService);
		vtkReader.setCellOptions(gmds::N|gmds::F);
		vtkReader.read(m_params.input_file_3D_surface);

		for (auto n_id : m_meshHex->nodes()) {
			m_couche_id->set(n_id, 0);
		}

		// Build the edges and the connectivities
		math::Utils::buildEfromFandConnectivies(m_meshHex);

	}

	//-------------------------------------//
	// Special case of the C2_3D geometry  //
	//	Surface Blocking 1						//
	//-------------------------------------//
	else if (m_params.block_surface_3D==1)
	{
		Node n0 = m_meshHex->newNode({-0.5, -0.5, -0.5});
		Node n1 = m_meshHex->newNode({-0.5, 0.5, -0.5});
		Node n2 = m_meshHex->newNode({0.5, 0.5, -0.5});
		Node n3 = m_meshHex->newNode({0.5, -0.5, -0.5});

		Node n4 = m_meshHex->newNode({-0.5, -0.5, 0.5});
		Node n5 = m_meshHex->newNode({-0.5, 0.5, 0.5});
		Node n6 = m_meshHex->newNode({0.5, 0.5, 0.5});
		Node n7 = m_meshHex->newNode({0.5, -0.5, 0.5});

		// Creates the faces
		TCellID f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n0.id(), n1.id(), n2.id(), n3.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n4.id(), n5.id(), n6.id(), n7.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n0.id(), n1.id(), n5.id(), n4.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n0.id(), n3.id(), n7.id(), n4.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n3.id(), n2.id(), n6.id(), n7.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n2.id(), n1.id(), n5.id(), n6.id());

		for (auto n_id : m_meshHex->nodes()) {
			m_couche_id->set(n_id, 0);
		}
	}

	//-------------------------------------//
	// Special case of the C2_3D geometry  //
	//	Surface Blocking 2						//
	//-------------------------------------//

	else if (m_params.block_surface_3D==2)
	{
		Node n1 = m_meshHex->newNode({-0.5, -0.5, -0.5});
		Node n2 = m_meshHex->newNode({0.0, -0.5, -0.5});
		Node n3 = m_meshHex->newNode({0.5, -0.5, -0.5});

		Node n4 = m_meshHex->newNode({-0.5, 0, -0.5});
		Node n5 = m_meshHex->newNode({0.0, 0, -0.5});
		Node n6 = m_meshHex->newNode({0.5, 0, -0.5});

		Node n7 = m_meshHex->newNode({-0.5, 0.5, -0.5});
		Node n8 = m_meshHex->newNode({0.0, 0.5, -0.5});
		Node n9 = m_meshHex->newNode({0.5, 0.5, -0.5});

		Node n10 = m_meshHex->newNode({-0.5, -0.5, 0});
		Node n11 = m_meshHex->newNode({0.0, -0.5, 0});
		Node n12 = m_meshHex->newNode({0.5, -0.5, 0});

		Node n13 = m_meshHex->newNode({-0.5, 0, 0});
		Node n14 = m_meshHex->newNode({0.5, 0, 0});

		Node n15 = m_meshHex->newNode({-0.5, 0.5, 0});
		Node n16 = m_meshHex->newNode({0.0, 0.5, 0});
		Node n17 = m_meshHex->newNode({0.5, 0.5, 0});

		Node n18 = m_meshHex->newNode({-0.5, -0.5, 0.5});
		Node n19 = m_meshHex->newNode({0.0, -0.5, 0.5});
		Node n20 = m_meshHex->newNode({0.5, -0.5, 0.5});

		Node n21 = m_meshHex->newNode({-0.5, 0, 0.5});
		Node n22 = m_meshHex->newNode({0.0, 0, 0.5});
		Node n23 = m_meshHex->newNode({0.5, 0, 0.5});

		Node n24 = m_meshHex->newNode({-0.5, 0.5, 0.5});
		Node n25 = m_meshHex->newNode({0.0, 0.5, 0.5});
		Node n26 = m_meshHex->newNode({0.5, 0.5, 0.5});

		// Create the faces
		TCellID f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n1.id(), n2.id(), n5.id(), n4.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n2.id(), n3.id(), n6.id(), n5.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n4.id(), n5.id(), n8.id(), n7.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n5.id(), n6.id(), n9.id(), n8.id());

		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n18.id(), n19.id(), n22.id(), n21.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n19.id(), n20.id(), n23.id(), n22.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n21.id(), n22.id(), n25.id(), n24.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n22.id(), n23.id(), n26.id(), n25.id());

		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n1.id(), n2.id(), n11.id(), n10.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n2.id(), n3.id(), n12.id(), n11.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n10.id(), n11.id(), n19.id(), n18.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n11.id(), n12.id(), n20.id(), n19.id());

		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n3.id(), n6.id(), n14.id(), n12.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n6.id(), n9.id(), n17.id(), n14.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n12.id(), n14.id(), n23.id(), n20.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n14.id(), n17.id(), n26.id(), n23.id());

		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n9.id(), n8.id(), n16.id(), n17.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n8.id(), n7.id(), n15.id(), n16.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n17.id(), n16.id(), n25.id(), n26.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n16.id(), n15.id(), n24.id(), n25.id());

		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n4.id(), n1.id(), n10.id(), n13.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n7.id(), n4.id(), n13.id(), n15.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n13.id(), n10.id(), n18.id(), n21.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n15.id(), n13.id(), n21.id(), n24.id());

		for (auto n_id : m_meshHex->nodes()) {
			m_couche_id->set(n_id, 0);
		}
	}


	//-------------------------------------//
	// Special case of the C3_3D geometry  //
	//	Surface Blocking 1						//
	//-------------------------------------//

	else if (m_params.block_surface_3D==3)
	{
		Node n0 = m_meshHex->newNode({0, 0, 0});

		Node n1 = m_meshHex->newNode({-0.5, -0.5, 0.5});
		Node n2 = m_meshHex->newNode({0.5, -0.5, 0.5});
		Node n3 = m_meshHex->newNode({0.0, 0.0, 0.5});
		Node n4 = m_meshHex->newNode({0.5, 0.0, 0.5});
		Node n5 = m_meshHex->newNode({-0.5, 0.5, 0.5});
		Node n6 = m_meshHex->newNode({0.0, 0.5, 0.5});

		Node n7 = m_meshHex->newNode({0.5, 0.0, 0.0});
		Node n8 = m_meshHex->newNode({0.0, 0.5, 0.0});
		Node n9 = m_meshHex->newNode({0.5, 0.5, 0.0});

		Node n10 = m_meshHex->newNode({-0.5, -0.5, -0.5});
		Node n11 = m_meshHex->newNode({0.5, -0.5, -0.5});
		Node n12 = m_meshHex->newNode({-0.5, 0.5, -0.5});
		Node n13 = m_meshHex->newNode({0.5, 0.5, -0.5});

		// Create the faces
		TCellID f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n1.id(), n2.id(), n11.id(), n10.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n1.id(), n10.id(), n12.id(), n5.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n10.id(), n11.id(), n13.id(), n12.id());

		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n3.id(), n4.id(), n7.id(), n0.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n0.id(), n7.id(), n9.id(), n8.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n3.id(), n0.id(), n8.id(), n6.id());

		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n1.id(), n3.id(), n6.id(), n5.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n1.id(), n2.id(), n4.id(), n3.id());

		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n5.id(), n6.id(), n8.id(), n12.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n8.id(), n9.id(), n13.id(), n12.id());

		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n2.id(), n4.id(), n7.id(), n11.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n7.id(), n9.id(), n13.id(), n11.id());

		for (auto n_id : m_meshHex->nodes()) {
			m_couche_id->set(n_id, 0);
		}

	}


	//-------------------------------------//
	// Special case of the C3_3D geometry  //
	//	Surface Blocking 2						//
	//-------------------------------------//

	else if (m_params.block_surface_3D==4)
	{
		Node n0 = m_meshHex->newNode({0, 0, 0});

		Node n1 = m_meshHex->newNode({-0.5, -0.5, 0.5});
		Node n2 = m_meshHex->newNode({0.5, -0.5, 0.5});
		Node n3 = m_meshHex->newNode({0.0, 0.0, 0.5});
		Node n4 = m_meshHex->newNode({0.5, 0.0, 0.5});
		Node n5 = m_meshHex->newNode({-0.5, 0.5, 0.5});
		Node n6 = m_meshHex->newNode({0.0, 0.5, 0.5});

		Node n7 = m_meshHex->newNode({0.5, 0.0, 0.0});
		Node n8 = m_meshHex->newNode({0.0, 0.5, 0.0});
		Node n9 = m_meshHex->newNode({0.5, 0.5, 0.0});

		Node n10 = m_meshHex->newNode({-0.5, -0.5, -0.5});
		Node n11 = m_meshHex->newNode({0.5, -0.5, -0.5});
		Node n12 = m_meshHex->newNode({-0.5, 0.5, -0.5});
		Node n13 = m_meshHex->newNode({0.5, 0.5, -0.5});

		Node n20 = m_meshHex->newNode({0, -0.5, 0.5});
		Node n21 = m_meshHex->newNode({-0.5, -0.5, 0});
		Node n22 = m_meshHex->newNode({0.0, -0.5, 0});
		Node n23 = m_meshHex->newNode({0.5, -0.5, 0});
		Node n24 = m_meshHex->newNode({0.0, -0.5, -0.5});
		Node n25 = m_meshHex->newNode({0.25, 0.0, 0.5});
		Node n26 = m_meshHex->newNode({0.0, 0.0, 0.25});
		Node n27 = m_meshHex->newNode({0.25, 0.0, 0.25});
		Node n28 = m_meshHex->newNode({0.5, 0.0, 0.25});
		Node n29 = m_meshHex->newNode({0.25, 0.0, 0.0});
		Node n30 = m_meshHex->newNode({-0.5, 0.0, 0.5});
		Node n31 = m_meshHex->newNode({0.0, 0.25, 0.5});
		Node n32 = m_meshHex->newNode({-0.5, 0.0, 0.0});
		Node n33 = m_meshHex->newNode({0.0, 0.25, 0.25});
		Node n34 = m_meshHex->newNode({0.0, 0.25, 0.0});
		Node n35 = m_meshHex->newNode({0.25, 0.25, 0.0});
		Node n36 = m_meshHex->newNode({0.5, 0.25, 0.0});
		Node n37 = m_meshHex->newNode({0.5, 0.0, -0.5});
		Node n38 = m_meshHex->newNode({-0.5, 0.0, -0.5});
		Node n39 = m_meshHex->newNode({0.0, 0.0, -0.5});
		Node n40 = m_meshHex->newNode({0.0, 0.5, 0.25});
		Node n41 = m_meshHex->newNode({-0.5, 0.5, 0.0});
		Node n42 = m_meshHex->newNode({0.25, 0.5, 0.0});
		Node n43 = m_meshHex->newNode({0.0, 0.5, -0.5});

		// Create the faces
		TCellID f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n1.id(), n3.id(), n31.id(), n30.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n30.id(), n31.id(), n6.id(), n5.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n1.id(), n20.id(), n25.id(), n3.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n20.id(), n2.id(), n4.id(), n25.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n0.id(), n29.id(), n35.id(), n34.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n29.id(), n7.id(), n36.id(), n35.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n34.id(), n35.id(), n42.id(), n8.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n35.id(), n36.id(), n9.id(), n42.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n10.id(), n24.id(), n39.id(), n38.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n24.id(), n11.id(), n37.id(), n39.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n12.id(), n43.id(), n39.id(), n38.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n39.id(), n37.id(), n13.id(), n43.id());

		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n2.id(), n4.id(), n28.id(), n23.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n11.id(), n7.id(), n28.id(), n23.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n11.id(), n7.id(), n36.id(), n37.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n13.id(), n9.id(), n36.id(), n37.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n3.id(), n31.id(), n33.id(), n26.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n0.id(), n34.id(), n33.id(), n26.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n40.id(), n33.id(), n34.id(), n8.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n40.id(), n33.id(), n31.id(), n6.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n1.id(), n30.id(), n32.id(), n21.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n10.id(), n38.id(), n32.id(), n21.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n12.id(), n38.id(), n32.id(), n41.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n5.id(), n30.id(), n32.id(), n41.id());

		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n5.id(), n6.id(), n40.id(), n41.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n12.id(), n8.id(), n40.id(), n41.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n12.id(), n8.id(), n42.id(), n43.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n13.id(), n9.id(), n42.id(), n43.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n0.id(), n29.id(), n27.id(), n26.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n3.id(), n25.id(), n27.id(), n26.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n28.id(), n27.id(), n25.id(), n4.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n28.id(), n27.id(), n29.id(), n7.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n24.id(), n22.id(), n21.id(), n10.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n24.id(), n22.id(), n23.id(), n11.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n20.id(), n22.id(), n23.id(), n2.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n20.id(), n22.id(), n21.id(), n1.id());

		for (auto n_id : m_meshHex->nodes()) {
			m_couche_id->set(n_id, 0);
		}

	}

	else if (m_params.block_surface_3D==5)
	{
		Node n0 = m_meshHex->newNode({-1, -1, 1});
		Node n1 = m_meshHex->newNode({0, -1, 1});
		Node n2 = m_meshHex->newNode({1, -1, 1});
		Node n3 = m_meshHex->newNode({-1, 0, 1});
		Node n4 = m_meshHex->newNode({0, 0, 1});
		Node n5 = m_meshHex->newNode({1, 0, 1});
		Node n6 = m_meshHex->newNode({-1, -1, 0});
		Node n7 = m_meshHex->newNode({0, -1, 0});
		Node n8 = m_meshHex->newNode({1, -1, 0});
		Node n9 = m_meshHex->newNode({-1, 0, 0});
		Node n10 = m_meshHex->newNode({0, 0, 0});
		Node n11 = m_meshHex->newNode({1, 0, 0});
		Node n12 = m_meshHex->newNode({-1, 1, 0});
		Node n13 = m_meshHex->newNode({0, 1, 0});
		Node n14 = m_meshHex->newNode({-1, -1, -1});
		Node n15 = m_meshHex->newNode({0, -1, -1});
		Node n16 = m_meshHex->newNode({1, -1, -1});
		Node n17 = m_meshHex->newNode({-1, 0, -1});
		Node n18 = m_meshHex->newNode({0, 0, -1});
		Node n19 = m_meshHex->newNode({1, 0, -1});
		Node n20 = m_meshHex->newNode({-1, 1, -1});
		Node n21 = m_meshHex->newNode({0, 1, -1});
		Node n30 = m_meshHex->newNode({-0.5, -1, 1});
		Node n31 = m_meshHex->newNode({0.5, -1, 1});
		Node n32 = m_meshHex->newNode({-1, -1, 0.5});
		Node n33 = m_meshHex->newNode({-0.5, -1, 0.5});
		Node n34 = m_meshHex->newNode({0, -1, 0.5});
		Node n35 = m_meshHex->newNode({0.5, -1, 0.5});
		Node n36 = m_meshHex->newNode({1, -1, 0.5});
		Node n37 = m_meshHex->newNode({-0.5, -1, 0});
		Node n38 = m_meshHex->newNode({0.5, -1, 0});
		Node n39 = m_meshHex->newNode({-1, -1, -0.5});
		Node n40 = m_meshHex->newNode({-0.5, -1, -0.5});
		Node n41 = m_meshHex->newNode({0, -1, -0.5});
		Node n42 = m_meshHex->newNode({0.5, -1, -0.5});
		Node n43 = m_meshHex->newNode({1, -1, -0.5});
		Node n44 = m_meshHex->newNode({-0.5, -1, -1});
		Node n45 = m_meshHex->newNode({0.5, -1, -1});
		Node n46 = m_meshHex->newNode({-0.5, 0, 1});
		Node n47 = m_meshHex->newNode({0.5, 0, 1});
		Node n48 = m_meshHex->newNode({-1, 0, 0.5});
		Node n49 = m_meshHex->newNode({-0.5, 0, 0.5});
		Node n50 = m_meshHex->newNode({0, 0, 0.5});
		Node n51 = m_meshHex->newNode({0.5, 0, 0.5});
		Node n52 = m_meshHex->newNode({1, 0, 0.5});
		Node n53 = m_meshHex->newNode({-0.5, 0, 0});
		Node n54 = m_meshHex->newNode({0.5, 0, 0});
		Node n55 = m_meshHex->newNode({-1, 0, -0.5});
		Node n56 = m_meshHex->newNode({0, 0, -0.5});
		Node n57 = m_meshHex->newNode({0.5, 0, -0.5});
		Node n58 = m_meshHex->newNode({1, 0, -0.5});
		Node n59 = m_meshHex->newNode({-0.5, 0, -1});
		Node n60 = m_meshHex->newNode({0.5, 0, -1});
		Node n61 = m_meshHex->newNode({-1, 0.5, 0});
		Node n62 = m_meshHex->newNode({-0.5, 0.5, 0});
		Node n63 = m_meshHex->newNode({0, 0.5, 0});
		Node n64 = m_meshHex->newNode({-1, 0.5, -0.5});
		Node n65 = m_meshHex->newNode({0, 0.5, -0.5});
		Node n66 = m_meshHex->newNode({-1, 0.5, -1});
		Node n67 = m_meshHex->newNode({-0.5, 0.5, -1});
		Node n68 = m_meshHex->newNode({0, 0.5, -1});
		Node n69 = m_meshHex->newNode({-0.5, 1, 0});
		Node n70 = m_meshHex->newNode({-1, 1, -0.5});
		Node n71 = m_meshHex->newNode({-0.5, 1, -0.5});
		Node n72 = m_meshHex->newNode({0, 1, -0.5});
		Node n73 = m_meshHex->newNode({-0.5, 1, -1});

		// Creates the faces
		// x=-1
		TCellID f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n0.id(), n32.id(), n48.id(), n3.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n48.id(), n32.id(), n6.id(), n9.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n6.id(), n9.id(), n55.id(), n39.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n14.id(), n17.id(), n55.id(), n39.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n55.id(), n17.id(), n66.id(), n64.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n61.id(), n9.id(), n55.id(), n64.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n61.id(), n12.id(), n70.id(), n64.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n66.id(), n20.id(), n70.id(), n64.id());
		// x=1
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n2.id(), n5.id(), n52.id(), n36.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n8.id(), n11.id(), n52.id(), n36.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n8.id(), n11.id(), n58.id(), n43.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n16.id(), n19.id(), n58.id(), n43.id());
		// x=0
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n10.id(), n63.id(), n65.id(), n56.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n18.id(), n68.id(), n65.id(), n56.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n21.id(), n68.id(), n65.id(), n72.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n13.id(), n63.id(), n65.id(), n72.id());

		// z=-1
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n17.id(), n14.id(), n44.id(), n59.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n18.id(), n15.id(), n44.id(), n59.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n18.id(), n15.id(), n45.id(), n60.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n19.id(), n16.id(), n45.id(), n60.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n17.id(), n59.id(), n67.id(), n66.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n20.id(), n73.id(), n67.id(), n66.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n21.id(), n73.id(), n67.id(), n68.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n18.id(), n59.id(), n67.id(), n68.id());
		//z=1
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n0.id(), n3.id(), n46.id(), n30.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n1.id(), n4.id(), n46.id(), n30.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n1.id(), n4.id(), n47.id(), n31.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n2.id(), n5.id(), n47.id(), n31.id());
		// z=0
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n10.id(), n63.id(), n62.id(), n53.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n9.id(), n61.id(), n62.id(), n53.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n12.id(), n61.id(), n62.id(), n69.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n13.id(), n63.id(), n62.id(), n69.id());

		// y=-1
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n14.id(), n44.id(), n40.id(), n39.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n6.id(), n37.id(), n40.id(), n39.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n6.id(), n37.id(), n33.id(), n32.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n0.id(), n30.id(), n33.id(), n32.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n1.id(), n30.id(), n33.id(), n34.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n7.id(), n37.id(), n33.id(), n34.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n7.id(), n37.id(), n40.id(), n41.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n15.id(), n44.id(), n40.id(), n41.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n42.id(), n45.id(), n15.id(), n41.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n42.id(), n41.id(), n7.id(), n38.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n35.id(), n34.id(), n7.id(), n38.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n35.id(), n34.id(), n1.id(), n31.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n35.id(), n36.id(), n2.id(), n31.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n35.id(), n36.id(), n8.id(), n38.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n42.id(), n43.id(), n8.id(), n38.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n42.id(), n43.id(), n16.id(), n45.id());
		// y = 0
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n3.id(), n46.id(), n49.id(), n48.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n9.id(), n53.id(), n49.id(), n48.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n10.id(), n53.id(), n49.id(), n50.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n4.id(), n46.id(), n49.id(), n50.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n4.id(), n47.id(), n51.id(), n50.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n5.id(), n47.id(), n51.id(), n52.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n11.id(), n54.id(), n51.id(), n52.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n11.id(), n54.id(), n57.id(), n58.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n19.id(), n60.id(), n57.id(), n58.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n18.id(), n60.id(), n57.id(), n56.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n10.id(), n54.id(), n57.id(), n56.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n10.id(), n54.id(), n51.id(), n50.id());
		// y = 1
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n12.id(), n69.id(), n71.id(), n70.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n20.id(), n73.id(), n71.id(), n70.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n21.id(), n73.id(), n71.id(), n72.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n13.id(), n69.id(), n71.id(), n72.id());


		for (auto n_id : m_meshHex->nodes()) {
			m_couche_id->set(n_id, 0);
		}

	}

	// Write the surface block structure
	gmds::IGMeshIOService ioService(m_meshHex);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("Surface_3D.vtk");

}
/*------------------------------------------------------------------------*/
void
AeroPipeline_3D::PreTraitementMeshTet()
{
	for (auto r_id:m_meshTet->regions())
	{
		Region r = m_meshTet->get<Region>(r_id);
		std::vector<Node> r_nodes = r.get<Node>();
		if ( m_Bnd->getNodeColor(r_nodes[0].id()) != 0
		    && (m_Bnd->getNodeColor(r_nodes[0].id()) == m_Bnd->getNodeColor(r_nodes[1].id()))
		    && (m_Bnd->getNodeColor(r_nodes[0].id()) == m_Bnd->getNodeColor(r_nodes[2].id()))
		    && (m_Bnd->getNodeColor(r_nodes[0].id()) == m_Bnd->getNodeColor(r_nodes[3].id())) )
		{
			//std::cout << "Tetra " << r_id << " pré traité" << std::endl;

			math::Point P_bar = r.center();
			m_meshTet->deleteRegion(r);
			Node n_new = m_meshTet->newNode(P_bar);

			r_nodes[0].remove<Region>(r);
			r_nodes[1].remove<Region>(r);
			r_nodes[2].remove<Region>(r);
			r_nodes[3].remove<Region>(r);

			// Get the 4 faces of the tet
			TCellID f1_id = math::Utils::CommonFace3Nodes(m_meshTet, r_nodes[0].id(), r_nodes[1].id(), r_nodes[2].id());
			TCellID f2_id = math::Utils::CommonFace3Nodes(m_meshTet, r_nodes[1].id(), r_nodes[2].id(), r_nodes[3].id());
			TCellID f3_id = math::Utils::CommonFace3Nodes(m_meshTet, r_nodes[0].id(), r_nodes[2].id(), r_nodes[3].id());
			TCellID f4_id = math::Utils::CommonFace3Nodes(m_meshTet, r_nodes[0].id(), r_nodes[1].id(), r_nodes[2].id());

			// New tets
			Region r_1 = m_meshTet->newTet(n_new, r_nodes[0], r_nodes[1], r_nodes[2]);
			Region r_2 = m_meshTet->newTet(n_new, r_nodes[1], r_nodes[2], r_nodes[3]);
			Region r_3 = m_meshTet->newTet(n_new, r_nodes[0], r_nodes[2], r_nodes[3]);
			Region r_4 = m_meshTet->newTet(n_new, r_nodes[0], r_nodes[1], r_nodes[3]);

			// The 3 new faces
			Face f_new_01 = m_meshTet->newTriangle(n_new, r_nodes[0], r_nodes[1]);
			Face f_new_02 = m_meshTet->newTriangle(n_new, r_nodes[0], r_nodes[2]);
			Face f_new_03 = m_meshTet->newTriangle(n_new, r_nodes[0], r_nodes[3]);
			Face f_new_12 = m_meshTet->newTriangle(n_new, r_nodes[1], r_nodes[2]);
			Face f_new_13 = m_meshTet->newTriangle(n_new, r_nodes[1], r_nodes[3]);
			Face f_new_23 = m_meshTet->newTriangle(n_new, r_nodes[2], r_nodes[3]);

			// The 4 new edges
			Edge e_new_0 = m_meshTet->newEdge(n_new, r_nodes[0]);
			Edge e_new_1 = m_meshTet->newEdge(n_new, r_nodes[1]);
			Edge e_new_2 = m_meshTet->newEdge(n_new, r_nodes[2]);
			Edge e_new_3 = m_meshTet->newEdge(n_new, r_nodes[3]);

			// Get the old tet edges
			TCellID e_01_id = math::Utils::CommonEdge(m_meshTet, r_nodes[0].id(), r_nodes[1].id()) ;
			TCellID e_02_id = math::Utils::CommonEdge(m_meshTet, r_nodes[0].id(), r_nodes[2].id()) ;
			TCellID e_03_id = math::Utils::CommonEdge(m_meshTet, r_nodes[0].id(), r_nodes[2].id()) ;
			TCellID e_12_id = math::Utils::CommonEdge(m_meshTet, r_nodes[1].id(), r_nodes[2].id()) ;
			TCellID e_13_id = math::Utils::CommonEdge(m_meshTet, r_nodes[1].id(), r_nodes[3].id()) ;
			TCellID e_23_id = math::Utils::CommonEdge(m_meshTet, r_nodes[2].id(), r_nodes[3].id()) ;

			// Update the connectivities N->R
			r_nodes[0].add<Region>(r_1);
			r_nodes[0].add<Region>(r_3);
			r_nodes[0].add<Region>(r_4);

			r_nodes[1].add<Region>(r_1);
			r_nodes[1].add<Region>(r_2);
			r_nodes[1].add<Region>(r_4);

			r_nodes[2].add<Region>(r_1);
			r_nodes[2].add<Region>(r_2);
			r_nodes[2].add<Region>(r_3);

			r_nodes[3].add<Region>(r_2);
			r_nodes[3].add<Region>(r_3);
			r_nodes[3].add<Region>(r_4);

			// Update the connectivities N->R (x4) for the new node
			n_new.add<Region>(r_1);
			n_new.add<Region>(r_2);
			n_new.add<Region>(r_3);
			n_new.add<Region>(r_4);

			// Update the connectivities, N->F (x6) and N->E (x6)
			n_new.add<Face>(f_new_01);
			n_new.add<Face>(f_new_02);
			n_new.add<Face>(f_new_03);
			n_new.add<Face>(f_new_12);
			n_new.add<Face>(f_new_13);
			n_new.add<Face>(f_new_23);

			n_new.add<Edge>(e_01_id);
			n_new.add<Edge>(e_02_id);
			n_new.add<Edge>(e_03_id);
			n_new.add<Edge>(e_12_id);
			n_new.add<Edge>(e_13_id);
			n_new.add<Edge>(e_23_id);

			// Update the connectivities for the nodes of the old tet
			r_nodes[0].add<Edge>(e_new_0);
			r_nodes[0].add<Face>(f_new_01);
			r_nodes[0].add<Face>(f_new_02);
			r_nodes[0].add<Face>(f_new_03);

			r_nodes[1].add<Edge>(e_new_1);
			r_nodes[1].add<Face>(f_new_01);
			r_nodes[1].add<Face>(f_new_12);
			r_nodes[1].add<Face>(f_new_13);

			r_nodes[2].add<Edge>(e_new_2);
			r_nodes[2].add<Face>(f_new_02);
			r_nodes[2].add<Face>(f_new_12);
			r_nodes[2].add<Face>(f_new_23);

			r_nodes[3].add<Edge>(e_new_3);
			r_nodes[3].add<Face>(f_new_03);
			r_nodes[3].add<Face>(f_new_13);
			r_nodes[3].add<Face>(f_new_23);

			// Update the connectivities R->F (x3) and R->E (x3)
			r_1.add<Face>(f1_id);
			r_1.add<Face>(f_new_01.id());
			r_1.add<Face>(f_new_12.id());
			r_1.add<Edge>(e_01_id);
			r_1.add<Edge>(e_12_id);
			r_1.add<Edge>(e_02_id);

			r_2.add<Face>(f2_id);
			r_2.add<Face>(f_new_12.id());
			r_2.add<Face>(f_new_23.id());
			r_2.add<Edge>(e_12_id);
			r_2.add<Edge>(e_02_id);
			r_2.add<Edge>(e_01_id);

			r_3.add<Face>(f3_id);
			r_3.add<Face>(f_new_02.id());
			r_3.add<Face>(f_new_23.id());
			r_3.add<Edge>(e_02_id);
			r_3.add<Edge>(e_23_id);
			r_3.add<Edge>(e_03_id);

			r_4.add<Face>(f4_id);
			r_4.add<Face>(f_new_01.id());
			r_4.add<Face>(f_new_13.id());
			r_4.add<Edge>(e_01_id);
			r_4.add<Edge>(e_13_id);
			r_4.add<Edge>(e_03_id);

		}
	}

	// Ecriture du maillage initial (tetra)
	gmds::IGMeshIOService ioService(m_meshTet);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write("AeroPipeline3D_Tetra_PreTraite.vtk");
}
/*------------------------------------------------------------------------*/
void
AeroPipeline_3D::SurfaceBlockingClassification()
{
	std::cout << "-> Classify the input quad surface blocking on the geometry ..." << std::endl;

	// Init the geometry manager and the linker
	m_manager->initAndLinkFrom3DMesh(m_meshTet, m_linker_TG);

	// Write the initial tet mesh (with fields computed on it)
	gmds::IGMeshIOService ioService = IGMeshIOService(m_meshTet);
	gmds::VTKWriter vtkWriter2(&ioService);
	vtkWriter2.setCellOptions(gmds::N|gmds::F);
	vtkWriter2.setDataOptions(gmds::N|gmds::F);
	vtkWriter2.write("AeroPipeline3D_Tetra_CLASSIFICATION.vtk");

	// Init the linker for the Blocking
	m_linker_HG->setGeometry(m_manager);
	m_linker_HG->setMesh(m_meshHex);

	//------------------------------//
	//	BLOCK CORNERS CLASSIFICATION //
	//------------------------------//
	for (auto n_id:m_meshHex->nodes())
	{
		std::vector<cad::GeomPoint*> points;
		m_manager->getPoints(points);

		std::vector<cad::GeomCurve*> curves;
		m_manager->getCurves(curves);

		std::vector<cad::GeomSurface*> surfaces;
		m_manager->getSurfaces(surfaces);

		Node n = m_meshHex->get<Node>(n_id);
		math::Point p = n.point();
		math::Point p_proj = n.point();

		double min_dist(std::numeric_limits<double>::max());
		int geom_dim = surfaces[0]->dim();
		int geom_id = 3;

		// Check on the geometrical points first
		for (auto point:points)
		{
			p_proj = n.point();
			point->project(p_proj);
			if ( (n.point()-p_proj).norm() < min_dist )
			{
				min_dist = (n.point()-p_proj).norm();
				geom_dim = 1;
				geom_id = point->id();
			}
		}

		// Check on the geometrical curves
		for (auto curve:curves)
		{
			p_proj = n.point();
			curve->project(p_proj);
			if ( (n.point()-p_proj).norm() < min_dist )
			{
				min_dist = (n.point()-p_proj).norm();
				geom_dim = 2;
				geom_id = curve->id();
			}
		}

		// Check on the geometrical surfaces
		for (auto surface:surfaces)
		{
			p_proj = n.point();
			surface->project(p_proj);
			if ( (n.point()-p_proj).norm() < min_dist
			    && min_dist > pow(10,-20))		// Need to add a tolerence here.
			{
				min_dist = (n.point()-p_proj).norm();
				geom_dim = 3;
				geom_id = surface->id();
			}
		}

		// Link the node of the Blocking to the geometrical entity that
		// minimize the distance
		if(geom_dim == 1){
			m_linker_HG->linkNodeToPoint(n_id, geom_id);
		}
		else if(geom_dim==2){
			m_linker_HG->linkNodeToCurve(n_id, geom_id);
		}
		else if(geom_dim==3){
			m_linker_HG->linkNodeToSurface(n_id, geom_id);
		}

	}

	//------------------------------//
	//	 BLOCK EDGES CLASSIFICATION  //
	//------------------------------//
	for (auto e_id:m_meshHex->edges())
	{
		std::vector<Node> e_nodes = m_meshHex->get<Edge>(e_id).get<Node>();
		int geom_id_n0 = m_linker_HG->getGeomId<Node>(e_nodes[0].id()) ;
		int geom_id_n1 = m_linker_HG->getGeomId<Node>(e_nodes[1].id()) ;
		int geom_dim_n0 = m_linker_HG->getGeomDim<Node>(e_nodes[0].id()) ;
		int geom_dim_n1 = m_linker_HG->getGeomDim<Node>(e_nodes[1].id());

		int geom_dim;
		int geom_id;

		if (geom_dim_n0 > geom_dim_n1)
		{
			geom_dim = geom_dim_n0;
			geom_id = geom_id_n0;
		}
		else if (geom_dim_n1 > geom_dim_n0)
		{
			geom_dim = geom_dim_n1;
			geom_id = geom_id_n1;
		}
		else if (geom_dim_n0 == geom_dim_n1
		         && geom_id_n0 == geom_id_n1)
		{
			geom_dim = geom_dim_n0;
			geom_id = geom_id_n0;
		}
		else if (geom_dim_n0 == geom_dim_n1
		         && geom_dim_n1 == 1)	// Find the common curve between the two geometrical points
		{
			cad::GeomPoint* geom_p0 = m_manager->getPoint(geom_id_n0);
			cad::GeomPoint* geom_p1 = m_manager->getPoint(geom_id_n1);
			geom_id = m_manager->getCommonCurve(geom_p0, geom_p1 );
			geom_dim = 2;
		}

		// Link the edge of the Blocking to the geometrical entity
		if(geom_dim==2){
			m_linker_HG->linkEdgeToCurve(e_id, geom_id);
		}
		else if(geom_dim==3){
			m_linker_HG->linkEdgeToSurface(e_id, geom_id);
		}

	}

	//------------------------------//
	//	 BLOCK FACES CLASSIFICATION  //
	//------------------------------//
	for (auto f_id:m_meshHex->faces())
	{
		std::vector<Edge> f_edges = m_meshHex->get<Face>(f_id).get<Edge>();
		int geom_id_e0 = m_linker_HG->getGeomId<Edge>(f_edges[0].id()) ;
		int geom_id_e1 = m_linker_HG->getGeomId<Edge>(f_edges[1].id()) ;
		int geom_id_e2 = m_linker_HG->getGeomId<Edge>(f_edges[2].id()) ;
		int geom_id_e3 = m_linker_HG->getGeomId<Edge>(f_edges[3].id()) ;
		int geom_dim_e0 = m_linker_HG->getGeomDim<Edge>(f_edges[0].id()) ;
		int geom_dim_e1 = m_linker_HG->getGeomDim<Edge>(f_edges[1].id()) ;
		int geom_dim_e2 = m_linker_HG->getGeomDim<Edge>(f_edges[2].id()) ;
		int geom_dim_e3 = m_linker_HG->getGeomDim<Edge>(f_edges[3].id()) ;

		int geom_dim(3);
		int geom_id;

		/*
		if (geom_dim_e0 == 3)
		{
			geom_id = geom_id_e0;
		}
		else if (geom_dim_e1 == 3)
		{
			geom_id = geom_id_e1;
		}
		else if (geom_dim_e2 == 3)
		{
			geom_id = geom_id_e2;
		}
		else if (geom_dim_e3 == 3)
		{
			geom_id = geom_id_e3;
		}
		else
		{
		 */
			std::vector<cad::GeomSurface*> surfs = m_manager->getSurfaces();
			Face f = m_meshHex->get<Face>(f_id);
			math::Point p = f.center() ;
			surfs[0]->project(p);
			double dist = (f.center()-p).norm();
			geom_id = surfs[0]->id();
			for (int i=1;i<surfs.size();i++)
			{
				p = f.center();
				surfs[i]->project(p);
				if ((f.center()-p).norm() < dist)
				{
					dist = (f.center()-p).norm();
					geom_id = surfs[i]->id();
				}
			}
			/*
			cad::GeomCurve* curve_0 = m_manager->getCurve(geom_id_e0);
			cad::GeomCurve* curve_1 = m_manager->getCurve(geom_id_e1);
			geom_id = m_manager->getCommonSurface(curve_0, curve_1);
			 */
		//}
		// Link the face to the surface
		m_linker_HG->linkFaceToSurface(f_id, geom_id);
	}

}
/*------------------------------------------------------------------------*/
void
AeroPipeline_3D::ComputeVectorFieldForExtrusion()
{
	Variable<math::Vector3d>* var_VectorsForExtrusion = m_meshTet->getOrCreateVariable<math::Vector3d, GMDS_NODE>("VectorField_Extrusion");

	if (m_params.vectors_field <= 0 || m_params.vectors_field > 4)
	{
		//m_meshTet->newVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient");
		LeastSquaresGradientComputation grad3D(m_meshTet,
		                                       m_meshTet->getVariable<double, GMDS_NODE>("GMDS_Distance"),
		   												var_VectorsForExtrusion);
		grad3D.execute();
	}

	else if (m_params.vectors_field == 1)
	{
		//m_meshTet->newVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient");
		LeastSquaresGradientComputation grad3D(m_meshTet,
		                                       m_meshTet->getVariable<double, GMDS_NODE>("GMDS_Distance_Int"),
		   									var_VectorsForExtrusion);
		grad3D.execute();
	}

	else if (m_params.vectors_field == 2)
	{
		Variable<math::Vector3d>* var_VectorField_1 = m_meshTet->getOrCreateVariable<math::Vector3d, GMDS_NODE>("VectorField_1");
		Variable<math::Vector3d>* var_VectorField_2 = m_meshTet->getOrCreateVariable<math::Vector3d, GMDS_NODE>("VectorField_2");

		LeastSquaresGradientComputation grad2D_1(m_meshTet, m_meshTet->getVariable<double,GMDS_NODE>("GMDS_Distance_Int"),
		                                         var_VectorField_1);
		grad2D_1.execute();
		LeastSquaresGradientComputation grad2D_2(m_meshTet, m_meshTet->getVariable<double,GMDS_NODE>("GMDS_Distance"),
		                                         var_VectorField_2);
		grad2D_2.execute();

		for (auto n_id:m_meshTet->nodes())
		{
			math::Vector3d vec_1 = var_VectorField_1->value(n_id);
			math::Vector3d vec_2 = var_VectorField_2->value(n_id);

			vec_1.normalize();
			vec_2.normalize();

			Node n = m_meshTet->get<Node>(n_id);
			math::Point p = n.point();

			if (p.X() <= m_params.x_VectorField_Z1)
			{
				var_VectorsForExtrusion->set(n_id, vec_1);
			}
			else if (p.X() >= m_params.x_VectorField_Z2)
			{
				var_VectorsForExtrusion->set(n_id, vec_2);
			}
			else
			{
				// Compute the transition field
				double alpha = (p.X() - m_params.x_VectorField_Z1)/(m_params.x_VectorField_Z2-m_params.x_VectorField_Z1) ;
				math::Vector3d v_transit = alpha*vec_2 + (1.0-alpha)*vec_1 ;
				v_transit.normalize();
				var_VectorsForExtrusion->set(n_id, v_transit);
			}

		}

		m_meshTet->deleteVariable(GMDS_NODE, var_VectorField_1);
		m_meshTet->deleteVariable(GMDS_NODE, var_VectorField_2);

	}

}
/*------------------------------------------------------------------------*/
void
AeroPipeline_3D::BlockingGeometricClassification()
{
	std::vector<cad::GeomSurface*> surfaces;
	m_manager->getSurfaces(surfaces);

	// Get the exterior surface
	double xyz_min[3];
	double xyz_max[3];
	cad::GeomSurface* surf_ext = surfaces[0];
	surf_ext->computeBoundingBox(xyz_min, xyz_max);
	for (auto const surf:surfaces)
	{
		double xyz_min_surf[3];
		double xyz_max_surf[3];
		surf->computeBoundingBox(xyz_min_surf, xyz_max_surf);
		if (xyz_min_surf[0] <= xyz_min[0]
		    && xyz_max[0] <= xyz_max_surf[0]
		    && xyz_min_surf[1] <= xyz_min[1]
		    && xyz_max[1] <= xyz_max_surf[1]
		    && xyz_min_surf[2] <= xyz_min[2]
		    && xyz_max[2] <= xyz_max_surf[2])
		{
			surf_ext = surf;
			xyz_min[0] = xyz_min_surf[0];
			xyz_min[1] = xyz_min_surf[1];
			xyz_min[2] = xyz_min_surf[2];
			xyz_max[0] = xyz_max_surf[0];
			xyz_max[1] = xyz_max_surf[1];
			xyz_max[2] = xyz_max_surf[2];
		}
	}

	//--------------------------------//
	//	  BLOCK NODES CLASSIFICATION   //
	//--------------------------------//
	for (auto const n_id:m_meshHex->nodes())
	{
		if (m_couche_id->value(n_id) == m_params.nbr_couches)
		{
			m_linker_HG->linkNodeToSurface(n_id, surf_ext->id());
		}
	}

	//--------------------------------//
	//	  BLOCK EDGES CLASSIFICATION   //
	//--------------------------------//
	for (auto const e_id:m_meshHex->edges())
	{
		std::vector<Node> e_nodes = m_meshHex->get<Edge>(e_id).get<Node>();
		if (m_couche_id->value(e_nodes[0].id()) == m_params.nbr_couches
		    && m_couche_id->value(e_nodes[1].id()) == m_params.nbr_couches)
		{
			m_linker_HG->linkEdgeToSurface(e_id, surf_ext->id());
		}
	}

	//--------------------------------//
	//	  BLOCK FACES CLASSIFICATION   //
	//--------------------------------//
	for (auto const f_id:m_meshHex->faces())
	{
		std::vector<Node> f_nodes = m_meshHex->get<Face>(f_id).get<Node>();
		if (m_couche_id->value(f_nodes[0].id()) == m_params.nbr_couches
		    && m_couche_id->value(f_nodes[1].id()) == m_params.nbr_couches
		    && m_couche_id->value(f_nodes[2].id()) == m_params.nbr_couches
		    && m_couche_id->value(f_nodes[3].id()) == m_params.nbr_couches)
		{
			m_linker_HG->linkFaceToSurface(f_id, surf_ext->id());
		}
	}

}
/*------------------------------------------------------------------------*/
void
AeroPipeline_3D::initBlocking3DfromMesh()
{
	std::cout << "-> init Blocking3D structure from the mesh..." << std::endl;

	m_linker_BG->setGeometry(m_manager);
	m_linker_BG->setMesh(&m_Blocking3D);

	//std::cout << "Nbr curves " << m_manager->getCurves().size() << std::endl;
	//std::cout << "Nbr surfaces " << m_manager->getSurfaces().size() << std::endl;

	Variable<int>* var_couche_blocking = m_Blocking3D.newVariable<int, GMDS_NODE>("GMDS_Couche");
	Variable<int>* var_couche_ctrlpts = m_CtrlPts.newVariable<int, GMDS_NODE>("GMDS_Couche");
	std::map<TCellID,TCellID> map_new_nodes_IDS_blocking;
	std::map<TCellID,TCellID> map_new_nodes_IDS_ctrlpts;

	//=========================================
	// Create the Block Corners and Blocks in
	// m_Blocking3D
	//=========================================
	// Copy the block corners
	for (auto n_id:m_meshHex->nodes())
	{
		Node n = m_meshHex->get<Node>(n_id);
		Node n_blocking = m_Blocking3D.newBlockCorner(n.point());
		Node n_ctrlpts = m_CtrlPts.newBlockCorner(n.point());

		map_new_nodes_IDS_blocking[n.id()] = n_blocking.id() ;
		map_new_nodes_IDS_ctrlpts[n.id()] = n_ctrlpts.id() ;

		var_couche_blocking->set(n_blocking.id(), m_couche_id->value(n_id));
		var_couche_ctrlpts->set(n_ctrlpts.id(), m_couche_id->value(n_id));

		math::Utils::UpdateLinker3D(m_linker_HG, n, m_linker_BG, n_blocking); // Init linker_BG
	}

	// Create the Blocks from the Block Corners
	for (auto r_id:m_meshHex->regions())
	{
		Region r = m_meshHex->get<Region>(r_id);
		std::vector<Node> r_nodes = r.get<Node>();
		Node n0 = m_meshHex->get<Node>(map_new_nodes_IDS_blocking[r_nodes[0].id()]);
		Node n1 = m_meshHex->get<Node>(map_new_nodes_IDS_blocking[r_nodes[1].id()]);
		Node n2 = m_meshHex->get<Node>(map_new_nodes_IDS_blocking[r_nodes[2].id()]);
		Node n3 = m_meshHex->get<Node>(map_new_nodes_IDS_blocking[r_nodes[3].id()]);
		Node n4 = m_meshHex->get<Node>(map_new_nodes_IDS_blocking[r_nodes[4].id()]);
		Node n5 = m_meshHex->get<Node>(map_new_nodes_IDS_blocking[r_nodes[5].id()]);
		Node n6 = m_meshHex->get<Node>(map_new_nodes_IDS_blocking[r_nodes[6].id()]);
		Node n7 = m_meshHex->get<Node>(map_new_nodes_IDS_blocking[r_nodes[7].id()]);
		Blocking3D::Block b_blocking = m_Blocking3D.newBlock(n0, n1, n2, n3, n4, n5, n6, n7);

		n0 = m_meshHex->get<Node>(map_new_nodes_IDS_ctrlpts[r_nodes[0].id()]);
		n1 = m_meshHex->get<Node>(map_new_nodes_IDS_ctrlpts[r_nodes[1].id()]);
		n2 = m_meshHex->get<Node>(map_new_nodes_IDS_ctrlpts[r_nodes[2].id()]);
		n3 = m_meshHex->get<Node>(map_new_nodes_IDS_ctrlpts[r_nodes[3].id()]);
		n4 = m_meshHex->get<Node>(map_new_nodes_IDS_ctrlpts[r_nodes[4].id()]);
		n5 = m_meshHex->get<Node>(map_new_nodes_IDS_ctrlpts[r_nodes[5].id()]);
		n6 = m_meshHex->get<Node>(map_new_nodes_IDS_ctrlpts[r_nodes[6].id()]);
		n7 = m_meshHex->get<Node>(map_new_nodes_IDS_ctrlpts[r_nodes[7].id()]);
		Blocking3D::Block b_ctrlpts = m_CtrlPts.newBlock(n0, n1, n2, n3, n4, n5, n6, n7);

	}


	//=========================================
	// Update classification of Block edges
	// and Block faces of m_Blocking3D
	//=========================================
	// Update the classification of the block edges
	for (auto e_id:m_meshHex->edges())
	{
		Edge e = m_meshHex->get<Edge>(e_id);
		std::vector<Node> e_nodes = e.get<Node>();
		if ( (m_couche_id->value(e_nodes[0].id()) == 0
		    && m_couche_id->value(e_nodes[1].id()) == 0)
		    || (m_couche_id->value(e_nodes[0].id()) == m_params.nbr_couches
		    && m_couche_id->value(e_nodes[1].id()) == m_params.nbr_couches) )
		{
			TCellID e_blocking_id = math::Utils::CommonEdge(&m_Blocking3D, map_new_nodes_IDS_blocking[e_nodes[0].id()], map_new_nodes_IDS_blocking[e_nodes[1].id()]);
			int geom_dim = m_linker_HG->getGeomDim<Edge>(e_id);
			if (geom_dim==2)
			{
				m_linker_BG->linkEdgeToCurve(e_blocking_id, m_linker_HG->getGeomId<Edge>(e_id));
			}
			else if (geom_dim==3)
			{
				m_linker_BG->linkEdgeToSurface(e_blocking_id, m_linker_HG->getGeomId<Edge>(e_id));
			}
		}
	}

	// Update the classification of the block faces
	for (auto f_id:m_meshHex->faces())
	{
		Face f = m_meshHex->get<Face>(f_id);
		std::vector<Node> f_nodes = f.get<Node>();
		if ( (m_couche_id->value(f_nodes[0].id()) == 0
		    && m_couche_id->value(f_nodes[1].id()) == 0
		    && m_couche_id->value(f_nodes[2].id()) == 0
		    && m_couche_id->value(f_nodes[3].id()) == 0)
		    || (m_couche_id->value(f_nodes[0].id()) == m_params.nbr_couches
		    && m_couche_id->value(f_nodes[1].id()) == m_params.nbr_couches
		    && m_couche_id->value(f_nodes[2].id()) == m_params.nbr_couches
		    && m_couche_id->value(f_nodes[3].id()) == m_params.nbr_couches))
		{
			TCellID f_blocking_id = math::Utils::CommonFace(&m_Blocking3D,
			                                                map_new_nodes_IDS_blocking[f_nodes[0].id()],
			                                                map_new_nodes_IDS_blocking[f_nodes[1].id()],
			                                                map_new_nodes_IDS_blocking[f_nodes[2].id()],
			                                                map_new_nodes_IDS_blocking[f_nodes[3].id()]);
			m_linker_BG->linkFaceToSurface(f_blocking_id, m_linker_HG->getGeomId<Face>(f_id));
		}
	}


	// Init the control points of each Block
	int degree_Bezier(3);
	for (auto bloc:m_CtrlPts.allBlocks())
	{
		bloc.setNbDiscretizationI(degree_Bezier+1);
		bloc.setNbDiscretizationJ(degree_Bezier+1);
		bloc.setNbDiscretizationK(degree_Bezier+1);
	}
	m_CtrlPts.initializeGridPoints();

	// Compute the positions of the control points around boundaries
	// to interpolate the boundary geometry
	computeControlPointstoInterpolateBoundaries();


	//===========================================
	// Interval Assignment:
	// Compute the discretization of each Block
	//===========================================

	std::cout << "-> Interval Assignment" << std::endl;
	IntervalAssignment_3D intAss(&m_Blocking3D,
	                             &m_CtrlPts,
	                             m_params,
	                             m_Blocking3D.newVariable<int,GMDS_EDGE>("GMDS_EdgeDiscretization"));
	intAss.execute();


	// Temporary: set the discretization of each block in an uniform way
	/*
	for (auto b:m_Blocking3D.allBlocks())
	{
		b.setNbDiscretizationI(11);
		b.setNbDiscretizationJ(11);
		b.setNbDiscretizationK(11);
	}
	 */

	// Init the grid points (the inner nodes of each block edge, face and hex)
	m_Blocking3D.initializeGridPoints();

	//=========================================
	// Update classification of Block edges
	// and Block faces inner nodes of
	// m_Blocking3D
	//=========================================
	for (auto b:m_Blocking3D.allBlocks())
	{
		int Nb_i = b.getNbDiscretizationI();
		int Nb_j = b.getNbDiscretizationJ();
		int Nb_k = b.getNbDiscretizationK();
		//--------------------------------
		// Classification of the node on
		// block faces
		//--------------------------------
		// Face k=0
		if ( (var_couche_blocking->value(b.getNode(0).id()) == 0
		    && var_couche_blocking->value(b.getNode(1).id()) == 0
		    && var_couche_blocking->value(b.getNode(2).id()) == 0
		    && var_couche_blocking->value(b.getNode(3).id()) == 0)
		    || (var_couche_blocking->value(b.getNode(0).id()) == m_params.nbr_couches
		    && var_couche_blocking->value(b.getNode(1).id()) == m_params.nbr_couches
		    && var_couche_blocking->value(b.getNode(2).id()) == m_params.nbr_couches
		    && var_couche_blocking->value(b.getNode(3).id()) == m_params.nbr_couches))
		{
			TCellID f0123_id = math::Utils::CommonFace(&m_Blocking3D,
			                                           b.getNode(0).id(),
			                                           b.getNode(1).id(),
			                                           b.getNode(2).id(),
			                                           b.getNode(3).id()) ;
			for (auto i=1;i<Nb_i-1;i++)
			{
				for (auto j=1;j<Nb_j-1;j++)
				{
					m_linker_BG->linkNodeToSurface(b(i,j,0).id(), m_linker_BG->getGeomId<Face>(f0123_id));
				}
			}
		}
		// Face k=Nb_k-1
		if ( (var_couche_blocking->value(b.getNode(4).id()) == 0
		    && var_couche_blocking->value(b.getNode(5).id()) == 0
		    && var_couche_blocking->value(b.getNode(6).id()) == 0
		    && var_couche_blocking->value(b.getNode(7).id()) == 0)
		    || (var_couche_blocking->value(b.getNode(4).id()) == m_params.nbr_couches
		    && var_couche_blocking->value(b.getNode(5).id()) == m_params.nbr_couches
		    && var_couche_blocking->value(b.getNode(6).id()) == m_params.nbr_couches
		    && var_couche_blocking->value(b.getNode(7).id()) == m_params.nbr_couches))
		{
			TCellID f4567_id = math::Utils::CommonFace(&m_Blocking3D,
			                                           b.getNode(4).id(),
			                                           b.getNode(5).id(),
			                                           b.getNode(6).id(),
			                                           b.getNode(7).id()) ;
			for (auto i=1;i<Nb_i-1;i++)
			{
				for (auto j=1;j<Nb_j-1;j++)
				{
					m_linker_BG->linkNodeToSurface(b(i,j,Nb_k-1).id(), m_linker_BG->getGeomId<Face>(f4567_id));
				}
			}
		}
		// Face j=0
		if ( (var_couche_blocking->value(b.getNode(0).id()) == 0
		     && var_couche_blocking->value(b.getNode(1).id()) == 0
		     && var_couche_blocking->value(b.getNode(5).id()) == 0
		     && var_couche_blocking->value(b.getNode(4).id()) == 0)
		    || (var_couche_blocking->value(b.getNode(0).id()) == m_params.nbr_couches
		        && var_couche_blocking->value(b.getNode(1).id()) == m_params.nbr_couches
		        && var_couche_blocking->value(b.getNode(5).id()) == m_params.nbr_couches
		        && var_couche_blocking->value(b.getNode(4).id()) == m_params.nbr_couches))
		{
			TCellID f0154_id = math::Utils::CommonFace(&m_Blocking3D,
			                                           b.getNode(0).id(),
			                                           b.getNode(1).id(),
			                                           b.getNode(5).id(),
			                                           b.getNode(4).id()) ;
			for (auto i=1;i<Nb_i-1;i++)
			{
				for (auto k=1;k<Nb_k-1;k++)
				{
					m_linker_BG->linkNodeToSurface(b(i,0,k).id(), m_linker_BG->getGeomId<Face>(f0154_id));
				}
			}
		}
		// Face j=Nb_j-1
		if ( (var_couche_blocking->value(b.getNode(3).id()) == 0
		     && var_couche_blocking->value(b.getNode(2).id()) == 0
		     && var_couche_blocking->value(b.getNode(6).id()) == 0
		     && var_couche_blocking->value(b.getNode(7).id()) == 0)
		    || (var_couche_blocking->value(b.getNode(3).id()) == m_params.nbr_couches
		        && var_couche_blocking->value(b.getNode(2).id()) == m_params.nbr_couches
		        && var_couche_blocking->value(b.getNode(6).id()) == m_params.nbr_couches
		        && var_couche_blocking->value(b.getNode(7).id()) == m_params.nbr_couches))
		{
			TCellID f3267_id = math::Utils::CommonFace(&m_Blocking3D,
			                                           b.getNode(3).id(),
			                                           b.getNode(2).id(),
			                                           b.getNode(6).id(),
			                                           b.getNode(7).id()) ;
			for (auto i=1;i<Nb_i-1;i++)
			{
				for (auto k=1;k<Nb_k-1;k++)
				{
					m_linker_BG->linkNodeToSurface(b(i,Nb_j-1,k).id(), m_linker_BG->getGeomId<Face>(f3267_id));
				}
			}
		}
		// Face i=0
		if ( (var_couche_blocking->value(b.getNode(0).id()) == 0
		     && var_couche_blocking->value(b.getNode(3).id()) == 0
		     && var_couche_blocking->value(b.getNode(7).id()) == 0
		     && var_couche_blocking->value(b.getNode(4).id()) == 0)
		    || (var_couche_blocking->value(b.getNode(0).id()) == m_params.nbr_couches
		        && var_couche_blocking->value(b.getNode(3).id()) == m_params.nbr_couches
		        && var_couche_blocking->value(b.getNode(7).id()) == m_params.nbr_couches
		        && var_couche_blocking->value(b.getNode(4).id()) == m_params.nbr_couches))
		{
			TCellID f0374_id = math::Utils::CommonFace(&m_Blocking3D,
			                                           b.getNode(0).id(),
			                                           b.getNode(3).id(),
			                                           b.getNode(7).id(),
			                                           b.getNode(4).id()) ;
			for (auto j=1;j<Nb_j-1;j++)
			{
				for (auto k=1;k<Nb_k-1;k++)
				{
					m_linker_BG->linkNodeToSurface(b(0,j,k).id(), m_linker_BG->getGeomId<Face>(f0374_id));
				}
			}
		}
		// Face i=Nb_i-1
		if ( (var_couche_blocking->value(b.getNode(1).id()) == 0
		     && var_couche_blocking->value(b.getNode(2).id()) == 0
		     && var_couche_blocking->value(b.getNode(6).id()) == 0
		     && var_couche_blocking->value(b.getNode(5).id()) == 0)
		    || (var_couche_blocking->value(b.getNode(1).id()) == m_params.nbr_couches
		        && var_couche_blocking->value(b.getNode(2).id()) == m_params.nbr_couches
		        && var_couche_blocking->value(b.getNode(6).id()) == m_params.nbr_couches
		        && var_couche_blocking->value(b.getNode(5).id()) == m_params.nbr_couches))
		{
			TCellID f1265_id = math::Utils::CommonFace(&m_Blocking3D,
			                                           b.getNode(1).id(),
			                                           b.getNode(2).id(),
			                                           b.getNode(6).id(),
			                                           b.getNode(5).id()) ;
			for (auto j=1;j<Nb_j-1;j++)
			{
				for (auto k=1;k<Nb_k-1;k++)
				{
					m_linker_BG->linkNodeToSurface(b(Nb_i-1,j,k).id(), m_linker_BG->getGeomId<Face>(f1265_id));
				}
			}
		}
		//--------------------------------
		// Classification of the node on
		// block edges
		//--------------------------------
		// Edge I: 0,1
		if ( (var_couche_blocking->value(b.getNode(0).id()) == 0
		     && var_couche_blocking->value(b.getNode(1).id()) == 0)
		    || (var_couche_blocking->value(b.getNode(0).id()) == m_params.nbr_couches
		        && var_couche_blocking->value(b.getNode(1).id()) == m_params.nbr_couches))
		{
			TCellID e01_id = math::Utils::CommonEdge(&m_Blocking3D,
			                                         b.getNode(0).id(),
			                                         b.getNode(1).id());
			int geom_dim = m_linker_BG->getGeomDim<Edge>(e01_id);
			if (geom_dim == 2)
			{
				for (auto i=1;i<Nb_i-1;i++)
				{
					m_linker_BG->linkNodeToCurve(b(i,0,0).id(), m_linker_BG->getGeomId<Edge>(e01_id));
				}
			}
			else if (geom_dim == 3)
			{
				for (auto i=1;i<Nb_i-1;i++)
				{
					m_linker_BG->linkNodeToSurface(b(i,0,0).id(), m_linker_BG->getGeomId<Edge>(e01_id));
				}
			}
		}
		// Edge I: 2,3
		if ( (var_couche_blocking->value(b.getNode(2).id()) == 0
		     && var_couche_blocking->value(b.getNode(3).id()) == 0)
		    || (var_couche_blocking->value(b.getNode(2).id()) == m_params.nbr_couches
		        && var_couche_blocking->value(b.getNode(3).id()) == m_params.nbr_couches))
		{
			TCellID e23_id = math::Utils::CommonEdge(&m_Blocking3D,
			                                         b.getNode(2).id(),
			                                         b.getNode(3).id());
			int geom_dim = m_linker_BG->getGeomDim<Edge>(e23_id);
			if (geom_dim == 2)
			{
				for (auto i=1;i<Nb_i-1;i++)
				{
					m_linker_BG->linkNodeToCurve(b(i,Nb_j-1,0).id(), m_linker_BG->getGeomId<Edge>(e23_id));
				}
			}
			else if (geom_dim == 3)
			{
				for (auto i=1;i<Nb_i-1;i++)
				{
					m_linker_BG->linkNodeToSurface(b(i,Nb_j-1,0).id(), m_linker_BG->getGeomId<Edge>(e23_id));
				}
			}
		}
		// Edge I: 4,5
		if ( (var_couche_blocking->value(b.getNode(4).id()) == 0
		     && var_couche_blocking->value(b.getNode(5).id()) == 0)
		    || (var_couche_blocking->value(b.getNode(4).id()) == m_params.nbr_couches
		        && var_couche_blocking->value(b.getNode(5).id()) == m_params.nbr_couches))
		{
			TCellID e45_id = math::Utils::CommonEdge(&m_Blocking3D,
			                                         b.getNode(4).id(),
			                                         b.getNode(5).id());
			int geom_dim = m_linker_BG->getGeomDim<Edge>(e45_id);
			if (geom_dim == 2)
			{
				for (auto i=1;i<Nb_i-1;i++)
				{
					m_linker_BG->linkNodeToCurve(b(i,0,Nb_k-1).id(), m_linker_BG->getGeomId<Edge>(e45_id));
				}
			}
			else if (geom_dim == 3)
			{
				for (auto i=1;i<Nb_i-1;i++)
				{
					m_linker_BG->linkNodeToSurface(b(i,0,Nb_k-1).id(), m_linker_BG->getGeomId<Edge>(e45_id));
				}
			}
		}
		// Edge I: 7,6
		if ( (var_couche_blocking->value(b.getNode(7).id()) == 0
		     && var_couche_blocking->value(b.getNode(6).id()) == 0)
		    || (var_couche_blocking->value(b.getNode(7).id()) == m_params.nbr_couches
		        && var_couche_blocking->value(b.getNode(6).id()) == m_params.nbr_couches))
		{
			TCellID e76_id = math::Utils::CommonEdge(&m_Blocking3D,
			                                         b.getNode(7).id(),
			                                         b.getNode(6).id());
			int geom_dim = m_linker_BG->getGeomDim<Edge>(e76_id);
			if (geom_dim == 2)
			{
				for (auto i=1;i<Nb_i-1;i++)
				{
					m_linker_BG->linkNodeToCurve(b(i,Nb_j-1,Nb_k-1).id(), m_linker_BG->getGeomId<Edge>(e76_id));
				}
			}
			else if (geom_dim == 3)
			{
				for (auto i=1;i<Nb_i-1;i++)
				{
					m_linker_BG->linkNodeToSurface(b(i,Nb_j-1,Nb_k-1).id(), m_linker_BG->getGeomId<Edge>(e76_id));
				}
			}
		}
		// Edge J: 0,3
		if ( (var_couche_blocking->value(b.getNode(0).id()) == 0
		     && var_couche_blocking->value(b.getNode(3).id()) == 0)
		    || (var_couche_blocking->value(b.getNode(0).id()) == m_params.nbr_couches
		        && var_couche_blocking->value(b.getNode(3).id()) == m_params.nbr_couches))
		{
			TCellID e03_id = math::Utils::CommonEdge(&m_Blocking3D,
			                                         b.getNode(0).id(),
			                                         b.getNode(3).id());
			int geom_dim = m_linker_BG->getGeomDim<Edge>(e03_id);
			if (geom_dim == 2)
			{
				for (auto j=1;j<Nb_j-1;j++)
				{
					m_linker_BG->linkNodeToCurve(b(0,j,0).id(), m_linker_BG->getGeomId<Edge>(e03_id));
				}
			}
			else if (geom_dim == 3)
			{
				for (auto j=1;j<Nb_j-1;j++)
				{
					m_linker_BG->linkNodeToSurface(b(0,j,0).id(), m_linker_BG->getGeomId<Edge>(e03_id));
				}
			}
		}
		// Edge J: 1,2
		if ( (var_couche_blocking->value(b.getNode(1).id()) == 0
		     && var_couche_blocking->value(b.getNode(2).id()) == 0)
		    || (var_couche_blocking->value(b.getNode(1).id()) == m_params.nbr_couches
		        && var_couche_blocking->value(b.getNode(2).id()) == m_params.nbr_couches))
		{
			TCellID e12_id = math::Utils::CommonEdge(&m_Blocking3D,
			                                         b.getNode(1).id(),
			                                         b.getNode(2).id());
			int geom_dim = m_linker_BG->getGeomDim<Edge>(e12_id);
			if (geom_dim == 2)
			{
				for (auto j=1;j<Nb_j-1;j++)
				{
					m_linker_BG->linkNodeToCurve(b(Nb_i-1,j,0).id(), m_linker_BG->getGeomId<Edge>(e12_id));
				}
			}
			else if (geom_dim == 3)
			{
				for (auto j=1;j<Nb_j-1;j++)
				{
					m_linker_BG->linkNodeToSurface(b(Nb_i-1,j,0).id(), m_linker_BG->getGeomId<Edge>(e12_id));
				}
			}
		}
		// Edge J: 5,6
		if ( (var_couche_blocking->value(b.getNode(5).id()) == 0
		     && var_couche_blocking->value(b.getNode(6).id()) == 0)
		    || (var_couche_blocking->value(b.getNode(5).id()) == m_params.nbr_couches
		        && var_couche_blocking->value(b.getNode(6).id()) == m_params.nbr_couches))
		{
			TCellID e56_id = math::Utils::CommonEdge(&m_Blocking3D,
			                                         b.getNode(5).id(),
			                                         b.getNode(6).id());
			int geom_dim = m_linker_BG->getGeomDim<Edge>(e56_id);
			if (geom_dim == 2)
			{
				for (auto j=1;j<Nb_j-1;j++)
				{
					m_linker_BG->linkNodeToCurve(b(Nb_i-1,j,Nb_k-1).id(), m_linker_BG->getGeomId<Edge>(e56_id));
				}
			}
			else if (geom_dim == 3)
			{
				for (auto j=1;j<Nb_j-1;j++)
				{
					m_linker_BG->linkNodeToSurface(b(Nb_i-1,j,Nb_k-1).id(), m_linker_BG->getGeomId<Edge>(e56_id));
				}
			}
		}
		// Edge J: 4,7
		if ( (var_couche_blocking->value(b.getNode(4).id()) == 0
		     && var_couche_blocking->value(b.getNode(7).id()) == 0)
		    || (var_couche_blocking->value(b.getNode(4).id()) == m_params.nbr_couches
		        && var_couche_blocking->value(b.getNode(7).id()) == m_params.nbr_couches))
		{
			TCellID e47_id = math::Utils::CommonEdge(&m_Blocking3D,
			                                         b.getNode(4).id(),
			                                         b.getNode(7).id());
			int geom_dim = m_linker_BG->getGeomDim<Edge>(e47_id);
			if (geom_dim == 2)
			{
				for (auto j=1;j<Nb_j-1;j++)
				{
					m_linker_BG->linkNodeToCurve(b(0,j,Nb_k-1).id(), m_linker_BG->getGeomId<Edge>(e47_id));
				}
			}
			else if (geom_dim == 3)
			{
				for (auto j=1;j<Nb_j-1;j++)
				{
					m_linker_BG->linkNodeToSurface(b(0,j,Nb_k-1).id(), m_linker_BG->getGeomId<Edge>(e47_id));
				}
			}
		}
		// Edge K: 0,4
		if ( (var_couche_blocking->value(b.getNode(0).id()) == 0
		     && var_couche_blocking->value(b.getNode(4).id()) == 0)
		    || (var_couche_blocking->value(b.getNode(0).id()) == m_params.nbr_couches
		        && var_couche_blocking->value(b.getNode(4).id()) == m_params.nbr_couches))
		{
			TCellID e04_id = math::Utils::CommonEdge(&m_Blocking3D,
			                                         b.getNode(0).id(),
			                                         b.getNode(4).id());
			int geom_dim = m_linker_BG->getGeomDim<Edge>(e04_id);
			if (geom_dim == 2)
			{
				for (auto k=1;k<Nb_k-1;k++)
				{
					m_linker_BG->linkNodeToCurve(b(0,0,k).id(), m_linker_BG->getGeomId<Edge>(e04_id));
				}
			}
			else if (geom_dim == 3)
			{
				for (auto k=1;k<Nb_k-1;k++)
				{
					m_linker_BG->linkNodeToSurface(b(0,0,k).id(), m_linker_BG->getGeomId<Edge>(e04_id));
				}
			}
		}
		// Edge K: 1,5
		if ( (var_couche_blocking->value(b.getNode(1).id()) == 0
		     && var_couche_blocking->value(b.getNode(5).id()) == 0)
		    || (var_couche_blocking->value(b.getNode(1).id()) == m_params.nbr_couches
		        && var_couche_blocking->value(b.getNode(5).id()) == m_params.nbr_couches))
		{
			TCellID e15_id = math::Utils::CommonEdge(&m_Blocking3D,
			                                         b.getNode(1).id(),
			                                         b.getNode(5).id());
			int geom_dim = m_linker_BG->getGeomDim<Edge>(e15_id);
			if (geom_dim == 2)
			{
				for (auto k=1;k<Nb_k-1;k++)
				{
					m_linker_BG->linkNodeToCurve(b(Nb_i-1,0,k).id(), m_linker_BG->getGeomId<Edge>(e15_id));
				}
			}
			else if (geom_dim == 3)
			{
				for (auto k=1;k<Nb_k-1;k++)
				{
					m_linker_BG->linkNodeToSurface(b(Nb_i-1,0,k).id(), m_linker_BG->getGeomId<Edge>(e15_id));
				}
			}
		}
		// Edge K: 2,6
		if ( (var_couche_blocking->value(b.getNode(2).id()) == 0
		     && var_couche_blocking->value(b.getNode(6).id()) == 0)
		    || (var_couche_blocking->value(b.getNode(2).id()) == m_params.nbr_couches
		        && var_couche_blocking->value(b.getNode(6).id()) == m_params.nbr_couches))
		{
			TCellID e26_id = math::Utils::CommonEdge(&m_Blocking3D,
			                                         b.getNode(2).id(),
			                                         b.getNode(6).id());
			int geom_dim = m_linker_BG->getGeomDim<Edge>(e26_id);
			if (geom_dim == 2)
			{
				for (auto k=1;k<Nb_k-1;k++)
				{
					m_linker_BG->linkNodeToCurve(b(Nb_i-1,Nb_j-1,k).id(), m_linker_BG->getGeomId<Edge>(e26_id));
				}
			}
			else if (geom_dim == 3)
			{
				for (auto k=1;k<Nb_k-1;k++)
				{
					m_linker_BG->linkNodeToSurface(b(Nb_i-1,Nb_j-1,k).id(), m_linker_BG->getGeomId<Edge>(e26_id));
				}
			}
		}
		// Edge K: 3,7
		if ( (var_couche_blocking->value(b.getNode(3).id()) == 0
		     && var_couche_blocking->value(b.getNode(7).id()) == 0)
		    || (var_couche_blocking->value(b.getNode(3).id()) == m_params.nbr_couches
		        && var_couche_blocking->value(b.getNode(7).id()) == m_params.nbr_couches))
		{
			TCellID e37_id = math::Utils::CommonEdge(&m_Blocking3D,
			                                         b.getNode(3).id(),
			                                         b.getNode(7).id());
			int geom_dim = m_linker_BG->getGeomDim<Edge>(e37_id);
			if (geom_dim == 2)
			{
				for (auto k=1;k<Nb_k-1;k++)
				{
					m_linker_BG->linkNodeToCurve(b(0,Nb_j-1,k).id(), m_linker_BG->getGeomId<Edge>(e37_id));
				}
			}
			else if (geom_dim == 3)
			{
				for (auto k=1;k<Nb_k-1;k++)
				{
					m_linker_BG->linkNodeToSurface(b(0,Nb_j-1,k).id(), m_linker_BG->getGeomId<Edge>(e37_id));
				}
			}
		}

	}



	// Project boundary nodes onto the geometry
	for (auto n_id:m_Blocking3D.nodes())
	{
		if (var_couche_blocking->value(n_id) == 0
		    || var_couche_blocking->value(n_id) == m_params.nbr_couches )
		{
			Node n =  m_Blocking3D.get<Node>(n_id);
			int geom_dim = m_linker_BG->getGeomDim<Node>(n_id);
			if (geom_dim == 2)
			{
				cad::GeomCurve* curve = m_manager->getCurve(m_linker_BG->getGeomId<Node>(n_id));
				math::Point p = n.point();
				curve->project(p);
				n.setPoint(p);
			}
			else if (geom_dim == 3)
			{
				cad::GeomSurface* surface = m_manager->getSurface(m_linker_BG->getGeomId<Node>(n_id));
				math::Point p = n.point();
				surface->project(p);
				n.setPoint(p);
			}
		}
	}

	// Test new method to compute the block nodes positions from ctrl points positions
	computeBlockNodesPositionsFromCtrlPoints();

}
/*------------------------------------------------------------------------*/
void
AeroPipeline_3D::updateLayerValues()
{
	Variable<int>* var_couche_blocking = m_Blocking3D.getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche");

	// Update the "Couche" variable on nodes
	for (auto b:m_Blocking3D.allBlocks())
	{
		int layer(0);
		layer = std::max(var_couche_blocking->value(b.getNode(0).id()), var_couche_blocking->value(b.getNode(1).id()));
		layer = std::max(layer, var_couche_blocking->value(b.getNode(2).id()));
		layer = std::max(layer, var_couche_blocking->value(b.getNode(3).id()));
		layer = std::max(layer, var_couche_blocking->value(b.getNode(4).id()));
		layer = std::max(layer, var_couche_blocking->value(b.getNode(5).id()));
		layer = std::max(layer, var_couche_blocking->value(b.getNode(6).id()));
		layer = std::max(layer, var_couche_blocking->value(b.getNode(7).id()));
		int Nb_i = b.getNbDiscretizationI();
		int Nb_j = b.getNbDiscretizationJ();
		int Nb_k = b.getNbDiscretizationK();
		for (auto i=1; i < Nb_i-1; i++)
		{
			for (auto j=1; j<Nb_j-1; j++)
			{
				for (auto k=1; k<Nb_k-1; k++)
				{
					var_couche_blocking->set(b(i,j,k).id(), layer);
				}
			}
		}

		// -------------------
		// Update face nodes
		// -------------------
		// Faces k=0 and k=Nb_k-1
		int layer_k0 = std::max(var_couche_blocking->value(b.getNode(0).id()), var_couche_blocking->value(b.getNode(1).id()));
		layer_k0 = std::max(layer_k0, var_couche_blocking->value(b.getNode(2).id()));
		layer_k0 = std::max(layer_k0, var_couche_blocking->value(b.getNode(3).id()));
		int layer_kmax = std::max(var_couche_blocking->value(b.getNode(4).id()), var_couche_blocking->value(b.getNode(5).id()));
		layer_kmax = std::max(layer_kmax, var_couche_blocking->value(b.getNode(6).id()));
		layer_kmax = std::max(layer_kmax, var_couche_blocking->value(b.getNode(7).id()));
		for (auto i=1; i < Nb_i-1; i++)
		{
			for (auto j=1; j<Nb_j-1; j++)
			{
				var_couche_blocking->set(b(i,j,0).id(), layer_k0);
				var_couche_blocking->set(b(i,j,Nb_k-1).id(), layer_kmax);
			}
		}

		// Faces j=0 and j=Nb_j-1
		int layer_j0 = std::max(var_couche_blocking->value(b.getNode(0).id()), var_couche_blocking->value(b.getNode(1).id()));
		layer_j0 = std::max(layer_j0, var_couche_blocking->value(b.getNode(5).id()));
		layer_j0 = std::max(layer_j0, var_couche_blocking->value(b.getNode(4).id()));
		int layer_jmax = std::max(var_couche_blocking->value(b.getNode(3).id()), var_couche_blocking->value(b.getNode(2).id()));
		layer_jmax = std::max(layer_jmax, var_couche_blocking->value(b.getNode(6).id()));
		layer_jmax = std::max(layer_jmax, var_couche_blocking->value(b.getNode(7).id()));
		for (auto i=1; i < Nb_i-1; i++)
		{
			for (auto k=1; k<Nb_k-1; k++)
			{
				var_couche_blocking->set(b(i,0,k).id(), layer_j0);
				var_couche_blocking->set(b(i,Nb_j-1,k).id(), layer_jmax);
			}
		}

		// Faces i=0 and i=Nb_i-1
		int layer_i0 = std::max(var_couche_blocking->value(b.getNode(0).id()), var_couche_blocking->value(b.getNode(3).id()));
		layer_i0 = std::max(layer_i0, var_couche_blocking->value(b.getNode(7).id()));
		layer_i0 = std::max(layer_i0, var_couche_blocking->value(b.getNode(4).id()));
		int layer_imax = std::max(var_couche_blocking->value(b.getNode(1).id()), var_couche_blocking->value(b.getNode(2).id()));
		layer_imax = std::max(layer_imax, var_couche_blocking->value(b.getNode(6).id()));
		layer_imax = std::max(layer_imax, var_couche_blocking->value(b.getNode(5).id()));
		for (auto j=1; j < Nb_j-1; j++)
		{
			for (auto k=1; k<Nb_k-1; k++)
			{
				var_couche_blocking->set(b(0,j,k).id(), layer_i0);
				var_couche_blocking->set(b(Nb_i-1,j,k).id(), layer_imax);
			}
		}


		// -------------------
		// Update edge nodes
		// -------------------
		// Edges i
		int layer_e01 =  std::max(var_couche_blocking->value(b.getNode(0).id()), var_couche_blocking->value(b.getNode(1).id()));
		int layer_e23 =  std::max(var_couche_blocking->value(b.getNode(2).id()), var_couche_blocking->value(b.getNode(3).id()));
		int layer_e67 =  std::max(var_couche_blocking->value(b.getNode(6).id()), var_couche_blocking->value(b.getNode(7).id()));
		int layer_e45 =  std::max(var_couche_blocking->value(b.getNode(4).id()), var_couche_blocking->value(b.getNode(5).id()));
		for (auto i=1; i<Nb_i-1; i++)
		{
			var_couche_blocking->set(b(i,0,0).id(), layer_e01);
			var_couche_blocking->set(b(i,Nb_j-1,0).id(), layer_e23);
			var_couche_blocking->set(b(i,Nb_j-1,Nb_k-1).id(), layer_e67);
			var_couche_blocking->set(b(i,0,Nb_k-1).id(), layer_e45);

		}

		// Edges j
		int layer_e03 =  std::max(var_couche_blocking->value(b.getNode(0).id()), var_couche_blocking->value(b.getNode(3).id()));
		int layer_e12 =  std::max(var_couche_blocking->value(b.getNode(1).id()), var_couche_blocking->value(b.getNode(2).id()));
		int layer_e56 =  std::max(var_couche_blocking->value(b.getNode(5).id()), var_couche_blocking->value(b.getNode(6).id()));
		int layer_e47 =  std::max(var_couche_blocking->value(b.getNode(4).id()), var_couche_blocking->value(b.getNode(7).id()));
		for (auto j=1; j<Nb_j-1; j++)
		{
			var_couche_blocking->set(b(0,j,0).id(), layer_e03);
			var_couche_blocking->set(b(Nb_i-1,j,0).id(), layer_e12);
			var_couche_blocking->set(b(Nb_i-1,j,Nb_k-1).id(), layer_e56);
			var_couche_blocking->set(b(0,j,Nb_k-1).id(), layer_e47);
		}

		// Edges k
		int layer_e04 =  std::max(var_couche_blocking->value(b.getNode(0).id()), var_couche_blocking->value(b.getNode(4).id()));
		int layer_e15 =  std::max(var_couche_blocking->value(b.getNode(1).id()), var_couche_blocking->value(b.getNode(5).id()));
		int layer_e62 =  std::max(var_couche_blocking->value(b.getNode(6).id()), var_couche_blocking->value(b.getNode(2).id()));
		int layer_e37 =  std::max(var_couche_blocking->value(b.getNode(3).id()), var_couche_blocking->value(b.getNode(7).id()));
		for (auto k=1; k<Nb_k-1; k++)
		{
			var_couche_blocking->set(b(0,0,k).id(), layer_e04);
			var_couche_blocking->set(b(Nb_i-1,0,k).id(), layer_e15);
			var_couche_blocking->set(b(Nb_i-1,Nb_j-1,k).id(), layer_e62);
			var_couche_blocking->set(b(0,Nb_j-1,k).id(), layer_e37);
		}

	}
}
/*------------------------------------------------------------------------*/
void
AeroPipeline_3D::computeBlockNodesPositionsFromCtrlPoints()
{
	for (auto b:m_Blocking3D.allBlocks())
	{
		Blocking3D::Block b_ctrlpoints = m_CtrlPts.block(b.id());	// Because we suppose the numbering of blocks is the same between m_Blocking3D and m_CtrlPts, but we could use a variable to store the corresponding numbering
		Array3D<math::Point> ctrl_points(b_ctrlpoints.getNbDiscretizationI(), b_ctrlpoints.getNbDiscretizationJ(), b_ctrlpoints.getNbDiscretizationK());
		for (auto i=0;i<b_ctrlpoints.getNbDiscretizationI();i++)
		{
			for (auto j=0;j<b_ctrlpoints.getNbDiscretizationJ();j++)
			{
				for (auto k=0;k<b_ctrlpoints.getNbDiscretizationK();k++)
				{
					ctrl_points(i,j,k) = b_ctrlpoints(i,j,k).point();
				}
			}
		}

		math::BezierHex bezier_hex(ctrl_points);
		Array3D<math::Point> b_points = bezier_hex.getDiscretization(b.getNbDiscretizationI(), b.getNbDiscretizationJ(), b.getNbDiscretizationK());
		for (auto i=0;i<b.getNbDiscretizationI();i++)
		{
			for (auto j=0;j<b.getNbDiscretizationJ();j++)
			{
				for (auto k=0;k<b.getNbDiscretizationK();k++)
				{
					//b(i,j,k).setPoint(b_points(i,j,k));
					double u = 1.0*i/(b.getNbDiscretizationI()-1.0) ;
					double v = 1.0*j/(b.getNbDiscretizationJ()-1.0) ;
					double w = 1.0*k/(b.getNbDiscretizationK()-1.0) ;
					b(i,j,k).setPoint(bezier_hex(u,v,w));
				}
			}
		}

	}
}
/*------------------------------------------------------------------------*/
void
AeroPipeline_3D::computeControlPointstoInterpolateBoundaries()
{
	Variable<int>* var_couche_b = m_Blocking3D.getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche");

	int degree = m_CtrlPts.block(0).getNbDiscretizationI()-1;

	// On each face of the boundary, there are (degree+1)*(degree+1) points to interpolate
	Eigen::MatrixXd mat_B((degree+1)*(degree+1), (degree+1)*(degree+1));
	Eigen::VectorXd ctrl_points_x((degree+1)*(degree+1));
	Eigen::VectorXd ctrl_points_y((degree+1)*(degree+1));
	Eigen::VectorXd ctrl_points_z((degree+1)*(degree+1));
	Eigen::VectorXd interp_points_x((degree+1)*(degree+1));
	Eigen::VectorXd interp_points_y((degree+1)*(degree+1));
	Eigen::VectorXd interp_points_z((degree+1)*(degree+1));

	for (auto b:m_Blocking3D.allBlocks())
	{
		Blocking3D::Block b_ctrlpts = m_CtrlPts.block(b.id()) ;

		// Face K=0
		if ( (var_couche_b->value(b.getNode(0).id()) == 0
		    && var_couche_b->value(b.getNode(1).id()) == 0
		    && var_couche_b->value(b.getNode(2).id()) == 0
		    && var_couche_b->value(b.getNode(3).id()) == 0)
		    || (var_couche_b->value(b.getNode(0).id()) == m_params.nbr_couches
		    && var_couche_b->value(b.getNode(1).id()) == m_params.nbr_couches
		    && var_couche_b->value(b.getNode(2).id()) == m_params.nbr_couches
		    && var_couche_b->value(b.getNode(3).id()) == m_params.nbr_couches) )
		{
			TCellID f0123_id = math::Utils::CommonFace(&m_Blocking3D, b.getNode(0).id(), b.getNode(1).id(), b.getNode(2).id(), b.getNode(3).id());
			int geom_id = m_linker_BG->getGeomId<Face>(f0123_id) ;
			cad::GeomSurface* surf = m_manager->getSurface(geom_id);
			// Compute positions of points to interpolate and fill interp_points_* vectors
			for (int i=0;i<b_ctrlpts.getNbDiscretizationI();i++)
			{
				for (int j=0;j<b_ctrlpts.getNbDiscretizationJ();j++)
				{
					math::Point p = b_ctrlpts(i, j, 0).point();
					if ( !(i==0 && j==0)
					    && !(i==0 && j==b_ctrlpts.getNbDiscretizationJ()-1)
					    && !(i==b_ctrlpts.getNbDiscretizationI()-1 && j==0 )
					    && !(i==b_ctrlpts.getNbDiscretizationI()-1 && j==b_ctrlpts.getNbDiscretizationJ()-1) )
					{
						surf->project(p);
					}
					interp_points_x[i+j*(b_ctrlpts.getNbDiscretizationJ())] = p.X();
					interp_points_y[i+j*(b_ctrlpts.getNbDiscretizationJ())] = p.Y();
					interp_points_z[i+j*(b_ctrlpts.getNbDiscretizationJ())] = p.Z();
				}
			}
			// Matrix Assembly
			for (int i=0;i<b_ctrlpts.getNbDiscretizationI();i++)
			{
				for (int j=0;j<b_ctrlpts.getNbDiscretizationJ();j++)
				{
					double u = 1.0*i/(b_ctrlpts.getNbDiscretizationI()-1) ;
					double v = 1.0*j/(b_ctrlpts.getNbDiscretizationJ()-1) ;
					for (int i2=0;i2<b_ctrlpts.getNbDiscretizationI();i2++)
					{
						for (int j2=0;j2<b_ctrlpts.getNbDiscretizationJ();j2++)
						{
							double Bi2n = math::Utils::BernsteinPolynomial(b_ctrlpts.getNbDiscretizationI()-1, i2, u) ;
							double Bj2m = math::Utils::BernsteinPolynomial(b_ctrlpts.getNbDiscretizationJ()-1, j2, v) ;
							mat_B(i+j*b_ctrlpts.getNbDiscretizationJ(), i2+j2*b_ctrlpts.getNbDiscretizationJ()) = Bi2n*Bj2m ;
						}
					}
				}
			}
			Eigen::MatrixXd mat_B_inv = mat_B.inverse();
			ctrl_points_x = mat_B_inv*interp_points_x;
			ctrl_points_y = mat_B_inv*interp_points_y;
			ctrl_points_z = mat_B_inv*interp_points_z;
			for (int i=0;i<b_ctrlpts.getNbDiscretizationI();i++)
			{
				for (int j=0;j<b_ctrlpts.getNbDiscretizationJ();j++)
				{
					// Set the same offset on the opposite face control points
					math::Point p_opp = b_ctrlpts(i,j,b_ctrlpts.getNbDiscretizationK()-1).point();
					p_opp.setX(p_opp.X() + (ctrl_points_x[i+j*b_ctrlpts.getNbDiscretizationJ()] - b_ctrlpts(i,j,0).X()) );
					p_opp.setY(p_opp.Y() + (ctrl_points_y[i+j*b_ctrlpts.getNbDiscretizationJ()] - b_ctrlpts(i,j,0).Y()) );
					p_opp.setZ(p_opp.Z() + (ctrl_points_z[i+j*b_ctrlpts.getNbDiscretizationJ()] - b_ctrlpts(i,j,0).Z()) );
					b_ctrlpts(i,j,b_ctrlpts.getNbDiscretizationK()-1).setPoint(p_opp);
					// Set the value of the new control point on the boundary
					b_ctrlpts(i, j, 0).setX(ctrl_points_x[i+j*b_ctrlpts.getNbDiscretizationJ()]);
					b_ctrlpts(i, j, 0).setY(ctrl_points_y[i+j*b_ctrlpts.getNbDiscretizationJ()]);
					b_ctrlpts(i, j, 0).setZ(ctrl_points_z[i+j*b_ctrlpts.getNbDiscretizationJ()]);
				}
			}

			b_ctrlpts.computeFaceNodesPoints(1,2,5,6);
			b_ctrlpts.computeFaceNodesPoints(2,3,6,7);
			b_ctrlpts.computeFaceNodesPoints(0,3,4,7);
			b_ctrlpts.computeFaceNodesPoints(0,1,4,5);
		}
		// Face K=max
		if ( (var_couche_b->value(b.getNode(4).id()) == 0
		     && var_couche_b->value(b.getNode(5).id()) == 0
		     && var_couche_b->value(b.getNode(6).id()) == 0
		     && var_couche_b->value(b.getNode(7).id()) == 0)
		    || (var_couche_b->value(b.getNode(4).id()) == m_params.nbr_couches
		        && var_couche_b->value(b.getNode(5).id()) == m_params.nbr_couches
		        && var_couche_b->value(b.getNode(6).id()) == m_params.nbr_couches
		        && var_couche_b->value(b.getNode(7).id()) == m_params.nbr_couches) )
		{
			TCellID f4567_id = math::Utils::CommonFace(&m_Blocking3D, b.getNode(4).id(), b.getNode(5).id(),
			                                           b.getNode(6).id(), b.getNode(7).id());
			int geom_id = m_linker_BG->getGeomId<Face>(f4567_id) ;
			cad::GeomSurface* surf = m_manager->getSurface(geom_id);
			// Compute positions of points to interpolate and fill interp_points_* vectors
			for (int i=0;i<b_ctrlpts.getNbDiscretizationI();i++)
			{
				for (int j=0;j<b_ctrlpts.getNbDiscretizationJ();j++)
				{
					math::Point p = b_ctrlpts(i, j, b_ctrlpts.getNbDiscretizationK()-1).point();
					if ( !(i==0 && j==0)
					    && !(i==0 && j==b_ctrlpts.getNbDiscretizationJ()-1)
					    && !(i==b_ctrlpts.getNbDiscretizationI()-1 && j==0 )
					    && !(i==b_ctrlpts.getNbDiscretizationI()-1 && j==b_ctrlpts.getNbDiscretizationJ()-1) )
					{
						surf->project(p);
					}
					interp_points_x[i+j*(b_ctrlpts.getNbDiscretizationJ())] = p.X();
					interp_points_y[i+j*(b_ctrlpts.getNbDiscretizationJ())] = p.Y();
					interp_points_z[i+j*(b_ctrlpts.getNbDiscretizationJ())] = p.Z();
				}
			}
			// Matrix Assembly
			for (int i=0;i<b_ctrlpts.getNbDiscretizationI();i++)
			{
				for (int j=0;j<b_ctrlpts.getNbDiscretizationJ();j++)
				{
					double u = 1.0*i/(b_ctrlpts.getNbDiscretizationI()-1) ;
					double v = 1.0*j/(b_ctrlpts.getNbDiscretizationJ()-1) ;
					for (int i2=0;i2<b_ctrlpts.getNbDiscretizationI();i2++)
					{
						for (int j2=0;j2<b_ctrlpts.getNbDiscretizationJ();j2++)
						{
							double Bi2n = math::Utils::BernsteinPolynomial(b_ctrlpts.getNbDiscretizationI()-1, i2, u) ;
							double Bj2m = math::Utils::BernsteinPolynomial(b_ctrlpts.getNbDiscretizationJ()-1, j2, v) ;
							mat_B(i+j*b_ctrlpts.getNbDiscretizationJ(), i2+j2*b_ctrlpts.getNbDiscretizationJ()) = Bi2n*Bj2m ;
						}
					}
				}
			}
			Eigen::MatrixXd mat_B_inv = mat_B.inverse();
			ctrl_points_x = mat_B_inv*interp_points_x;
			ctrl_points_y = mat_B_inv*interp_points_y;
			ctrl_points_z = mat_B_inv*interp_points_z;
			for (int i=0;i<b_ctrlpts.getNbDiscretizationI();i++)
			{
				for (int j=0;j<b_ctrlpts.getNbDiscretizationJ();j++)
				{
					// Set the same offset on the opposite face control points
					math::Point p_opp = b_ctrlpts(i,j,0).point();
					p_opp.setX(p_opp.X() + (ctrl_points_x[i+j*b_ctrlpts.getNbDiscretizationJ()] - b_ctrlpts(i,j,b_ctrlpts.getNbDiscretizationK()-1).X()) );
					p_opp.setY(p_opp.Y() + (ctrl_points_y[i+j*b_ctrlpts.getNbDiscretizationJ()] - b_ctrlpts(i,j,b_ctrlpts.getNbDiscretizationK()-1).Y()) );
					p_opp.setZ(p_opp.Z() + (ctrl_points_z[i+j*b_ctrlpts.getNbDiscretizationJ()] - b_ctrlpts(i,j,b_ctrlpts.getNbDiscretizationK()-1).Z()) );
					b_ctrlpts(i,j,0).setPoint(p_opp);
					// Set the value of the new control point on the boundary
					b_ctrlpts(i, j, b_ctrlpts.getNbDiscretizationK()-1).setX(ctrl_points_x[i+j*b_ctrlpts.getNbDiscretizationJ()]);
					b_ctrlpts(i, j, b_ctrlpts.getNbDiscretizationK()-1).setY(ctrl_points_y[i+j*b_ctrlpts.getNbDiscretizationJ()]);
					b_ctrlpts(i, j, b_ctrlpts.getNbDiscretizationK()-1).setZ(ctrl_points_z[i+j*b_ctrlpts.getNbDiscretizationJ()]);
				}
			}

			b_ctrlpts.computeFaceNodesPoints(1,2,5,6);
			b_ctrlpts.computeFaceNodesPoints(2,3,6,7);
			b_ctrlpts.computeFaceNodesPoints(0,3,4,7);
			b_ctrlpts.computeFaceNodesPoints(0,1,4,5);
		}
		// Face J=0
		if ( (var_couche_b->value(b.getNode(0).id()) == 0
		     && var_couche_b->value(b.getNode(1).id()) == 0
		     && var_couche_b->value(b.getNode(5).id()) == 0
		     && var_couche_b->value(b.getNode(4).id()) == 0)
		    || (var_couche_b->value(b.getNode(0).id()) == m_params.nbr_couches
		        && var_couche_b->value(b.getNode(1).id()) == m_params.nbr_couches
		        && var_couche_b->value(b.getNode(5).id()) == m_params.nbr_couches
		        && var_couche_b->value(b.getNode(4).id()) == m_params.nbr_couches) )
		{
			TCellID f0154_id = math::Utils::CommonFace(&m_Blocking3D, b.getNode(0).id(), b.getNode(1).id(),
			                                           b.getNode(5).id(), b.getNode(4).id());
			int geom_id = m_linker_BG->getGeomId<Face>(f0154_id) ;
			cad::GeomSurface* surf = m_manager->getSurface(geom_id);
			// Compute positions of points to interpolate and fill interp_points_* vectors
			for (int i=0;i<b_ctrlpts.getNbDiscretizationI();i++)
			{
				for (int k=0;k<b_ctrlpts.getNbDiscretizationK();k++)
				{
					math::Point p = b_ctrlpts(i, 0, k).point();
					if ( !(i==0 && k==0)
					    && !(i==0 && k==b_ctrlpts.getNbDiscretizationK()-1)
					    && !(i==b_ctrlpts.getNbDiscretizationI()-1 && k==0 )
					    && !(i==b_ctrlpts.getNbDiscretizationI()-1 && k==b_ctrlpts.getNbDiscretizationK()-1) )
					{
						surf->project(p);
					}
					interp_points_x[i+k*(b_ctrlpts.getNbDiscretizationK())] = p.X();
					interp_points_y[i+k*(b_ctrlpts.getNbDiscretizationK())] = p.Y();
					interp_points_z[i+k*(b_ctrlpts.getNbDiscretizationK())] = p.Z();
				}
			}
			// Matrix Assembly
			for (int i=0;i<b_ctrlpts.getNbDiscretizationI();i++)
			{
				for (int k=0;k<b_ctrlpts.getNbDiscretizationK();k++)
				{
					double u = 1.0*i/(b_ctrlpts.getNbDiscretizationI()-1) ;
					double v = 1.0*k/(b_ctrlpts.getNbDiscretizationK()-1) ;
					for (int i2=0;i2<b_ctrlpts.getNbDiscretizationI();i2++)
					{
						for (int k2=0;k2<b_ctrlpts.getNbDiscretizationK();k2++)
						{
							double Bi2n = math::Utils::BernsteinPolynomial(b_ctrlpts.getNbDiscretizationI()-1, i2, u) ;
							double Bk2m = math::Utils::BernsteinPolynomial(b_ctrlpts.getNbDiscretizationK()-1, k2, v) ;
							mat_B(i+k*b_ctrlpts.getNbDiscretizationK(), i2+k2*b_ctrlpts.getNbDiscretizationK()) = Bi2n*Bk2m ;
						}
					}
				}
			}
			Eigen::MatrixXd mat_B_inv = mat_B.inverse();
			ctrl_points_x = mat_B_inv*interp_points_x;
			ctrl_points_y = mat_B_inv*interp_points_y;
			ctrl_points_z = mat_B_inv*interp_points_z;
			for (int i=0;i<b_ctrlpts.getNbDiscretizationI();i++)
			{
				for (int k=0;k<b_ctrlpts.getNbDiscretizationK();k++)
				{
					// Set the same offset on the opposite face control points
					math::Point p_opp = b_ctrlpts(i,b_ctrlpts.getNbDiscretizationJ()-1,k).point();
					p_opp.setX(p_opp.X() + (ctrl_points_x[i+k*b_ctrlpts.getNbDiscretizationK()] - b_ctrlpts(i,0,k).X()) );
					p_opp.setY(p_opp.Y() + (ctrl_points_y[i+k*b_ctrlpts.getNbDiscretizationK()] - b_ctrlpts(i,0,k).Y()) );
					p_opp.setZ(p_opp.Z() + (ctrl_points_z[i+k*b_ctrlpts.getNbDiscretizationK()] - b_ctrlpts(i,0,k).Z()) );
					b_ctrlpts(i,b_ctrlpts.getNbDiscretizationJ()-1,k).setPoint(p_opp);
					// Set the value of the new control point on the boundary
					b_ctrlpts(i, 0, k).setX(ctrl_points_x[i+k*b_ctrlpts.getNbDiscretizationK()]);
					b_ctrlpts(i, 0, k).setY(ctrl_points_y[i+k*b_ctrlpts.getNbDiscretizationK()]);
					b_ctrlpts(i, 0, k).setZ(ctrl_points_z[i+k*b_ctrlpts.getNbDiscretizationK()]);
				}
			}

			b_ctrlpts.computeFaceNodesPoints(0,1,2,3);
			b_ctrlpts.computeFaceNodesPoints(1,2,5,6);
			b_ctrlpts.computeFaceNodesPoints(4,5,6,7);
			b_ctrlpts.computeFaceNodesPoints(0,3,4,7);
		}
		// Face J=max
		if ( (var_couche_b->value(b.getNode(3).id()) == 0
		     && var_couche_b->value(b.getNode(2).id()) == 0
		     && var_couche_b->value(b.getNode(6).id()) == 0
		     && var_couche_b->value(b.getNode(7).id()) == 0)
		    || (var_couche_b->value(b.getNode(3).id()) == m_params.nbr_couches
		        && var_couche_b->value(b.getNode(2).id()) == m_params.nbr_couches
		        && var_couche_b->value(b.getNode(6).id()) == m_params.nbr_couches
		        && var_couche_b->value(b.getNode(7).id()) == m_params.nbr_couches) )
		{
			TCellID f3267_id = math::Utils::CommonFace(&m_Blocking3D, b.getNode(3).id(), b.getNode(2).id(),
			                                           b.getNode(6).id(), b.getNode(7).id());
			int geom_id = m_linker_BG->getGeomId<Face>(f3267_id) ;
			cad::GeomSurface* surf = m_manager->getSurface(geom_id);
			// Compute positions of points to interpolate and fill interp_points_* vectors
			for (int i=0;i<b_ctrlpts.getNbDiscretizationI();i++)
			{
				for (int k=0;k<b_ctrlpts.getNbDiscretizationK();k++)
				{
					math::Point p = b_ctrlpts(i, b_ctrlpts.getNbDiscretizationJ()-1, k).point();
					if ( !(i==0 && k==0)
					    && !(i==0 && k==b_ctrlpts.getNbDiscretizationK()-1)
					    && !(i==b_ctrlpts.getNbDiscretizationI()-1 && k==0 )
					    && !(i==b_ctrlpts.getNbDiscretizationI()-1 && k==b_ctrlpts.getNbDiscretizationK()-1) )
					{
						surf->project(p);
					}
					interp_points_x[i+k*(b_ctrlpts.getNbDiscretizationK())] = p.X();
					interp_points_y[i+k*(b_ctrlpts.getNbDiscretizationK())] = p.Y();
					interp_points_z[i+k*(b_ctrlpts.getNbDiscretizationK())] = p.Z();
				}
			}
			// Matrix Assembly
			for (int i=0;i<b_ctrlpts.getNbDiscretizationI();i++)
			{
				for (int k=0;k<b_ctrlpts.getNbDiscretizationK();k++)
				{
					double u = 1.0*i/(b_ctrlpts.getNbDiscretizationI()-1) ;
					double v = 1.0*k/(b_ctrlpts.getNbDiscretizationK()-1) ;
					for (int i2=0;i2<b_ctrlpts.getNbDiscretizationI();i2++)
					{
						for (int k2=0;k2<b_ctrlpts.getNbDiscretizationK();k2++)
						{
							double Bi2n = math::Utils::BernsteinPolynomial(b_ctrlpts.getNbDiscretizationI()-1, i2, u) ;
							double Bk2m = math::Utils::BernsteinPolynomial(b_ctrlpts.getNbDiscretizationK()-1, k2, v) ;
							mat_B(i+k*b_ctrlpts.getNbDiscretizationK(), i2+k2*b_ctrlpts.getNbDiscretizationK()) = Bi2n*Bk2m ;
						}
					}
				}
			}
			Eigen::MatrixXd mat_B_inv = mat_B.inverse();
			ctrl_points_x = mat_B_inv*interp_points_x;
			ctrl_points_y = mat_B_inv*interp_points_y;
			ctrl_points_z = mat_B_inv*interp_points_z;
			for (int i=0;i<b_ctrlpts.getNbDiscretizationI();i++)
			{
				for (int k=0;k<b_ctrlpts.getNbDiscretizationK();k++)
				{
					// Set the same offset on the opposite face control points
					math::Point p_opp = b_ctrlpts(i,0,k).point();
					p_opp.setX(p_opp.X() + (ctrl_points_x[i+k*b_ctrlpts.getNbDiscretizationK()] - b_ctrlpts(i,b_ctrlpts.getNbDiscretizationJ()-1,k).X()) );
					p_opp.setY(p_opp.Y() + (ctrl_points_y[i+k*b_ctrlpts.getNbDiscretizationK()] - b_ctrlpts(i,b_ctrlpts.getNbDiscretizationJ()-1,k).Y()) );
					p_opp.setZ(p_opp.Z() + (ctrl_points_z[i+k*b_ctrlpts.getNbDiscretizationK()] - b_ctrlpts(i,b_ctrlpts.getNbDiscretizationJ()-1,k).Z()) );
					b_ctrlpts(i,0,k).setPoint(p_opp);
					// Set the value of the new control point on the boundary
					b_ctrlpts(i, b_ctrlpts.getNbDiscretizationJ()-1, k).setX(ctrl_points_x[i+k*b_ctrlpts.getNbDiscretizationK()]);
					b_ctrlpts(i, b_ctrlpts.getNbDiscretizationJ()-1, k).setY(ctrl_points_y[i+k*b_ctrlpts.getNbDiscretizationK()]);
					b_ctrlpts(i, b_ctrlpts.getNbDiscretizationJ()-1, k).setZ(ctrl_points_z[i+k*b_ctrlpts.getNbDiscretizationK()]);
				}
			}

			b_ctrlpts.computeFaceNodesPoints(0,1,2,3);
			b_ctrlpts.computeFaceNodesPoints(1,2,5,6);
			b_ctrlpts.computeFaceNodesPoints(4,5,6,7);
			b_ctrlpts.computeFaceNodesPoints(0,3,4,7);
		}
		// Face I=0
		if ( (var_couche_b->value(b.getNode(0).id()) == 0
		     && var_couche_b->value(b.getNode(3).id()) == 0
		     && var_couche_b->value(b.getNode(7).id()) == 0
		     && var_couche_b->value(b.getNode(4).id()) == 0)
		    || (var_couche_b->value(b.getNode(0).id()) == m_params.nbr_couches
		        && var_couche_b->value(b.getNode(3).id()) == m_params.nbr_couches
		        && var_couche_b->value(b.getNode(7).id()) == m_params.nbr_couches
		        && var_couche_b->value(b.getNode(4).id()) == m_params.nbr_couches) )
		{
			TCellID f0374_id = math::Utils::CommonFace(&m_Blocking3D, b.getNode(0).id(), b.getNode(3).id(),
			                                           b.getNode(7).id(), b.getNode(4).id());
			int geom_id = m_linker_BG->getGeomId<Face>(f0374_id) ;
			cad::GeomSurface* surf = m_manager->getSurface(geom_id);
			// Compute positions of points to interpolate and fill interp_points_* vectors
			for (int j=0;j<b_ctrlpts.getNbDiscretizationJ();j++)
			{
				for (int k=0;k<b_ctrlpts.getNbDiscretizationK();k++)
				{
					math::Point p = b_ctrlpts(0, j, k).point();
					if ( !(k==0 && j==0)
					    && !(k==0 && j==b_ctrlpts.getNbDiscretizationJ()-1)
					    && !(k==b_ctrlpts.getNbDiscretizationK()-1 && j==0 )
					    && !(k==b_ctrlpts.getNbDiscretizationK()-1 && j==b_ctrlpts.getNbDiscretizationJ()-1) )
					{
						surf->project(p);
					}
					interp_points_x[j+k*(b_ctrlpts.getNbDiscretizationK())] = p.X();
					interp_points_y[j+k*(b_ctrlpts.getNbDiscretizationK())] = p.Y();
					interp_points_z[j+k*(b_ctrlpts.getNbDiscretizationK())] = p.Z();
				}
			}
			// Matrix Assembly
			for (int j=0;j<b_ctrlpts.getNbDiscretizationJ();j++)
			{
				for (int k=0;k<b_ctrlpts.getNbDiscretizationK();k++)
				{
					double u = 1.0*j/(b_ctrlpts.getNbDiscretizationJ()-1) ;
					double v = 1.0*k/(b_ctrlpts.getNbDiscretizationK()-1) ;
					for (int j2=0;j2<b_ctrlpts.getNbDiscretizationJ();j2++)
					{
						for (int k2=0;k2<b_ctrlpts.getNbDiscretizationK();k2++)
						{
							double Bj2n = math::Utils::BernsteinPolynomial(b_ctrlpts.getNbDiscretizationJ()-1, j2, u) ;
							double Bk2m = math::Utils::BernsteinPolynomial(b_ctrlpts.getNbDiscretizationK()-1, k2, v) ;
							mat_B(j+k*b_ctrlpts.getNbDiscretizationK(), j2+k2*b_ctrlpts.getNbDiscretizationK()) = Bj2n*Bk2m ;
						}
					}
				}
			}
			Eigen::MatrixXd mat_B_inv = mat_B.inverse();
			ctrl_points_x = mat_B_inv*interp_points_x;
			ctrl_points_y = mat_B_inv*interp_points_y;
			ctrl_points_z = mat_B_inv*interp_points_z;
			for (int j=0;j<b_ctrlpts.getNbDiscretizationJ();j++)
			{
				for (int k=0;k<b_ctrlpts.getNbDiscretizationK();k++)
				{
					// Set the same offset on the opposite face control points
					math::Point p_opp = b_ctrlpts(b_ctrlpts.getNbDiscretizationI()-1,j,k).point();
					p_opp.setX(p_opp.X() + (ctrl_points_x[j+k*b_ctrlpts.getNbDiscretizationK()] - b_ctrlpts(0,j,k).X()) );
					p_opp.setY(p_opp.Y() + (ctrl_points_y[j+k*b_ctrlpts.getNbDiscretizationK()] - b_ctrlpts(0,j,k).Y()) );
					p_opp.setZ(p_opp.Z() + (ctrl_points_z[j+k*b_ctrlpts.getNbDiscretizationK()] - b_ctrlpts(0,j,k).Z()) );
					b_ctrlpts(b_ctrlpts.getNbDiscretizationI()-1,j,k).setPoint(p_opp);
					// Set the value of the new control point on the boundary
					b_ctrlpts(0, j, k).setX(ctrl_points_x[j+k*b_ctrlpts.getNbDiscretizationK()]);
					b_ctrlpts(0, j, k).setY(ctrl_points_y[j+k*b_ctrlpts.getNbDiscretizationK()]);
					b_ctrlpts(0, j, k).setZ(ctrl_points_z[j+k*b_ctrlpts.getNbDiscretizationK()]);
				}
			}

			b_ctrlpts.computeFaceNodesPoints(0,1,4,5);
			b_ctrlpts.computeFaceNodesPoints(0,1,2,3);
			b_ctrlpts.computeFaceNodesPoints(4,5,6,7);
			b_ctrlpts.computeFaceNodesPoints(2,3,6,7);
		}
		// Face I=max
		if ( (var_couche_b->value(b.getNode(1).id()) == 0
		     && var_couche_b->value(b.getNode(2).id()) == 0
		     && var_couche_b->value(b.getNode(6).id()) == 0
		     && var_couche_b->value(b.getNode(5).id()) == 0)
		    || (var_couche_b->value(b.getNode(1).id()) == m_params.nbr_couches
		        && var_couche_b->value(b.getNode(2).id()) == m_params.nbr_couches
		        && var_couche_b->value(b.getNode(6).id()) == m_params.nbr_couches
		        && var_couche_b->value(b.getNode(5).id()) == m_params.nbr_couches) )
		{
			TCellID f1265_id = math::Utils::CommonFace(&m_Blocking3D, b.getNode(1).id(), b.getNode(2).id(),
			                                           b.getNode(6).id(), b.getNode(5).id());
			int geom_id = m_linker_BG->getGeomId<Face>(f1265_id) ;
			cad::GeomSurface* surf = m_manager->getSurface(geom_id);
			// Compute positions of points to interpolate and fill interp_points_* vectors
			for (int j=0;j<b_ctrlpts.getNbDiscretizationJ();j++)
			{
				for (int k=0;k<b_ctrlpts.getNbDiscretizationK();k++)
				{
					math::Point p = b_ctrlpts(b_ctrlpts.getNbDiscretizationI()-1, j, k).point();
					if ( !(k==0 && j==0)
					    && !(k==0 && j==b_ctrlpts.getNbDiscretizationJ()-1)
					    && !(k==b_ctrlpts.getNbDiscretizationK()-1 && j==0 )
					    && !(k==b_ctrlpts.getNbDiscretizationK()-1 && j==b_ctrlpts.getNbDiscretizationJ()-1) )
					{
						surf->project(p);
					}
					interp_points_x[j+k*(b_ctrlpts.getNbDiscretizationK())] = p.X();
					interp_points_y[j+k*(b_ctrlpts.getNbDiscretizationK())] = p.Y();
					interp_points_z[j+k*(b_ctrlpts.getNbDiscretizationK())] = p.Z();
				}
			}
			// Matrix Assembly
			for (int j=0;j<b_ctrlpts.getNbDiscretizationJ();j++)
			{
				for (int k=0;k<b_ctrlpts.getNbDiscretizationK();k++)
				{
					double u = 1.0*j/(b_ctrlpts.getNbDiscretizationJ()-1) ;
					double v = 1.0*k/(b_ctrlpts.getNbDiscretizationK()-1) ;
					for (int j2=0;j2<b_ctrlpts.getNbDiscretizationJ();j2++)
					{
						for (int k2=0;k2<b_ctrlpts.getNbDiscretizationK();k2++)
						{
							double Bj2n = math::Utils::BernsteinPolynomial(b_ctrlpts.getNbDiscretizationJ()-1, j2, u) ;
							double Bk2m = math::Utils::BernsteinPolynomial(b_ctrlpts.getNbDiscretizationK()-1, k2, v) ;
							mat_B(j+k*b_ctrlpts.getNbDiscretizationK(), j2+k2*b_ctrlpts.getNbDiscretizationK()) = Bj2n*Bk2m ;
						}
					}
				}
			}
			Eigen::MatrixXd mat_B_inv = mat_B.inverse();
			ctrl_points_x = mat_B_inv*interp_points_x;
			ctrl_points_y = mat_B_inv*interp_points_y;
			ctrl_points_z = mat_B_inv*interp_points_z;
			for (int j=0;j<b_ctrlpts.getNbDiscretizationJ();j++)
			{
				for (int k=0;k<b_ctrlpts.getNbDiscretizationK();k++)
				{
					// Set the same offset on the opposite face control points
					math::Point p_opp = b_ctrlpts(0,j,k).point();
					p_opp.setX(p_opp.X() + (ctrl_points_x[j+k*b_ctrlpts.getNbDiscretizationK()] - b_ctrlpts(b_ctrlpts.getNbDiscretizationI()-1,j,k).X()) );
					p_opp.setY(p_opp.Y() + (ctrl_points_y[j+k*b_ctrlpts.getNbDiscretizationK()] - b_ctrlpts(b_ctrlpts.getNbDiscretizationI()-1,j,k).Y()) );
					p_opp.setZ(p_opp.Z() + (ctrl_points_z[j+k*b_ctrlpts.getNbDiscretizationK()] - b_ctrlpts(b_ctrlpts.getNbDiscretizationI()-1,j,k).Z()) );
					b_ctrlpts(0,j,k).setPoint(p_opp);
					// Set the value of the new control point on the boundary
					b_ctrlpts(b_ctrlpts.getNbDiscretizationI()-1, j, k).setX(ctrl_points_x[j+k*b_ctrlpts.getNbDiscretizationK()]);
					b_ctrlpts(b_ctrlpts.getNbDiscretizationI()-1, j, k).setY(ctrl_points_y[j+k*b_ctrlpts.getNbDiscretizationK()]);
					b_ctrlpts(b_ctrlpts.getNbDiscretizationI()-1, j, k).setZ(ctrl_points_z[j+k*b_ctrlpts.getNbDiscretizationK()]);
				}
			}

			b_ctrlpts.computeFaceNodesPoints(0,1,4,5);
			b_ctrlpts.computeFaceNodesPoints(0,1,2,3);
			b_ctrlpts.computeFaceNodesPoints(4,5,6,7);
			b_ctrlpts.computeFaceNodesPoints(2,3,6,7);
		}
	}

	// Recompute the positions of the inner control points
	for (auto b:m_CtrlPts.allBlocks())
	{
		b.computeInnerBlockNodesPoints();
		/*
		auto nb_I = b.getNbDiscretizationI();
		auto nb_J = b.getNbDiscretizationJ();
		auto nb_K = b.getNbDiscretizationK();
		Array3D<math::Point> pnts(nb_I, nb_J, nb_K);
		for (auto i=0;i<nb_I;i++)
		{
			for (auto j=0;j<nb_J;j++)
			{
				pnts(i,j,0) = b(i,j,0).point();
				pnts(i,j,nb_K-1) = b(i,j,nb_K-1).point();
			}
		}
		for (auto i=0;i<nb_I;i++)
		{
			for (auto k=0;k<nb_K;k++)
			{
				pnts(i,0,k) = b(i,0,k).point();
				pnts(i,nb_J-1,k) = b(i,nb_J-1,k).point();
			}
		}
		for (auto j=0;j<nb_J;j++)
		{
			for (auto k=0;k<nb_K;k++)
			{
				pnts(0,j,k) = b(0,j,k).point();
				pnts(nb_I-1,j,k) = b(nb_I-1,j,k).point();
			}
		}
		TransfiniteInterpolation_3D::computeHex(pnts);
		for(auto i=1; i<nb_I-1;i++)
		{
			for(auto j=1; j<nb_J-1;j++)
			{
				for (auto k=1; k<nb_K-1;k++)
				{
					b(i,j,k).setPoint(pnts(i,j,k));
				}
			}
		}
		 */
	}


}
/*------------------------------------------------------------------------*/