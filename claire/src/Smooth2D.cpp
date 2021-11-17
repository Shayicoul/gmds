/*------------------------------------------------------------------------*/
//
// Created by Claire Roche on 21/10/2021.
//
/*------------------------------------------------------------------------*/
#include <gmds/claire/Smooth2D.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/
Smooth2D::Smooth2D(Mesh *AMesh,
                   const Variable<int>* AVarBnd,
                   const int ANbIterations) {
	m_mesh = AMesh;
	m_node_constrained = AVarBnd;
	m_nb_max_iterations = ANbIterations;

}
/*------------------------------------------------------------------------*/
void Smooth2D::setNbIterations(const int ANbIterations)
{
	m_nb_max_iterations=ANbIterations;

}
/*------------------------------------------------------------------------*/
Smooth2D::STATUS Smooth2D::execute()
{
    //we traverse all the nodes of the whole mesh to get the nodes that
    //are not constrained, i.e are so free
	for(auto n_id:m_mesh->nodes()){
		if(m_node_constrained->value(n_id)!=1){
			m_free_nodes.push_back(n_id);
		}
	}

	// Check if the mesh is structured or not
	bool CheckMesh = CheckStructuredMesh();

	if (CheckMesh) {
		// for all free nodes, we build and store stencils locally
		buildStencils();

		// Un noeud a un point, une position géométrique
		// Variable est un grand tableau qui permet de stocker une valeur à chaque noeud du maillage
		// Ici, la valeur stockée est de type math::Point
		// On définit un objet math::Point pour chaque sommet du maillage
		// On accède aux valeurs en utilisant les numéros des id
		Variable<math::Point> *old_coords = NULL;
		// Ici, le tableau est alloué directement avec le nombre de sommets du maillage.
		// Attention : si un élément est supprimé, sa case n'est pas supprimée
		old_coords = m_mesh->newVariable<math::Point,GMDS_NODE>("old_coords") ;
		// Exemple : How to access to elements
		//math::Point test_point = old_coords->value(0) ;

		for (int iteration = 1; iteration <= m_nb_max_iterations; iteration++) {
			// ATTENTION : a local copy of the mesh is needed to compute the position of each node at the iteration n+1
			// on the positions of the nodes at the iteration n

			// Idée 1 :
			// Mesh old_mesh ; // NE FONCTIONNE PAS
			// old_mesh = *m_mesh ;

			// Idée 2 :
			// Mesh *old_mesh ; // old_mesh est un pointeur vers un maillage
			// *old_mesh = *m_mesh ; // Ne fonctionne pas
			// |-----> Je souhaite que le contenu de la case pointée par old_mesh prenne la même
			// valeur que le contenu de la case pointée par m_mesh

			// Idée 3 :
			// Mesh old_mesh(MeshModel(DIM3 | R | F | N | F2N | N2F)) ; // On créé un maillage nommé old_mesh
			// old_mesh = *m_mesh ; // NE FONCTIONNE PAS !
			// |----> Je souhaite que old_mesh prenne la même valeur que la case pointée par m_mesh

			// Boucle sur l'ensemble du maillage pour stocker les coordonnées de chaque noeud
			for(auto n_id:m_mesh->nodes()) {
				old_coords->set(n_id, m_mesh->get<Node>(n_id).point());
			}

			Mesh *new_mesh = m_mesh;
			Node n_new;
			// PENSER A REMPLACER DANS CETTE BOUCLE POUR CALCULER AVEC LES COORDS DU TABLEAU OLD COORDS
			for (auto n_id : m_free_nodes) {
				math::Point H1, H2, H3;
				math::Point V1, V2, V3;
				math::Point A, B, C;
				Node noeud_traite;
				Node noeud_voisin;
				// std::cout << "Noeud :" << n_id << std::endl;
				noeud_traite = m_mesh->get<Node>(n_id);

				// Noeuds de la première branche (verticale)
				noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[0][0]);
				std::cout << "NOEUD NOEUD :" << noeud_voisin << std::endl ;
				//A = noeud_voisin.point();
				A = old_coords->value(noeud_voisin.id()) ;
				noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[1][0]);
				//B = noeud_voisin.point();
				B = old_coords->value(noeud_voisin.id()) ;
				noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[2][0]);
				//C = noeud_voisin.point();
				C = old_coords->value(noeud_voisin.id()) ;
				V1 = FindMidBranche(A, B, C);

				noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[0][1]);
				//A = noeud_voisin.point();
				A = old_coords->value(noeud_voisin.id()) ;
				noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[1][1]);
				//B = noeud_voisin.point();
				B = old_coords->value(noeud_voisin.id()) ;
				noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[2][1]);
				//C = noeud_voisin.point();
				C = old_coords->value(noeud_voisin.id()) ;
				V2 = FindMidBranche(A, B, C);

				noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[0][2]);
				//A = noeud_voisin.point();
				A = old_coords->value(noeud_voisin.id()) ;
				noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[1][2]);
				//B = noeud_voisin.point();
				B = old_coords->value(noeud_voisin.id()) ;
				noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[2][2]);
				//C = noeud_voisin.point();
				C = old_coords->value(noeud_voisin.id()) ;
				V3 = FindMidBranche(A, B, C);

				// Noeuds de la seconde branche (horizontale)
				noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[0][0]);
				//A = noeud_voisin.point();
				A = old_coords->value(noeud_voisin.id()) ;
				noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[0][1]);
				//B = noeud_voisin.point();
				B = old_coords->value(noeud_voisin.id()) ;
				noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[0][2]);
				//C = noeud_voisin.point();
				C = old_coords->value(noeud_voisin.id()) ;
				H1 = FindMidBranche(A, B, C);

				noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[1][0]);
				//A = noeud_voisin.point();
				A = old_coords->value(noeud_voisin.id()) ;
				noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[1][1]);
				//B = noeud_voisin.point();
				B = old_coords->value(noeud_voisin.id()) ;
				noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[1][2]);
				//C = noeud_voisin.point();
				C = old_coords->value(noeud_voisin.id()) ;
				H2 = FindMidBranche(A, B, C);

				noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[2][0]);
				//A = noeud_voisin.point();
				A = old_coords->value(noeud_voisin.id()) ;
				noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[2][1]);
				//B = noeud_voisin.point();
				B = old_coords->value(noeud_voisin.id()) ;
				noeud_voisin = m_mesh->get<Node>(m_stencil[n_id].val[2][2]);
				//C = noeud_voisin.point();
				C = old_coords->value(noeud_voisin.id()) ;
				H3 = FindMidBranche(A, B, C);

				// Recherche de l'intersection entre les 4 segments
				bool intersection_trouvee(false);
				math::Point M(0, 0, 0);
				math::Segment Seg_Vert_1(V1, V2);
				math::Segment Seg_Vert_2(V2, V3);
				math::Segment Seg_Hori_1(H1, H2);
				math::Segment Seg_Hori_2(H2, H3);

				intersection_trouvee = Seg_Vert_1.intersect2D(Seg_Hori_1, M);
				if (!intersection_trouvee) {
					intersection_trouvee = Seg_Vert_1.intersect2D(Seg_Hori_2, M);
				}
				if (!intersection_trouvee) {
					intersection_trouvee = Seg_Vert_2.intersect2D(Seg_Hori_1, M);
				}
				if (!intersection_trouvee) {
					intersection_trouvee = Seg_Vert_2.intersect2D(Seg_Hori_2, M);
				}

				/*
				// PRINT POUR DEBUG
				if (iteration==m_nb_max_iterations) {
				   std::cout << "--------------------------------" << std::endl;
				   std::cout << "Iterations :" << iteration << std::endl;
				   std::cout << "Noeud :" << n_id << std::endl;
				   std::cout << "Intersection trouvee ?" << intersection_trouvee << std::endl;
				   std::cout << "Pos actuelle :" << noeud_traite.point() << std::endl;
				   std::cout << "Nouvelle position :" << M << std::endl;
				}
				// FIN DES PRINT POUR DEBUG
				*/

				n_new = new_mesh->get<Node>(n_id);
				// Stocker l'ancienne position dans la variable old_coords avant de changer sa position
				// old_coords->set(n_id, n_new.point());
				if (intersection_trouvee) {
					n_new.setPoint(M);
				}

				if ( (n_id == 5) && (iteration == 1) ) {
					// std::cout << "HELLOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOW" << std::endl ;
					write_debug_txt(n_id, old_coords, H1, H2, H3, V1, V2, V3, n_new.point(), "test.txt");
				}

			}

			m_mesh->deleteVariable(GMDS_NODE, old_coords);
			m_mesh = new_mesh;
		}
	}


	    return Smooth2D::SUCCESS;
}
/*------------------------------------------------------------------------*/
void Smooth2D::buildStencils() {
    for(auto n_id:m_free_nodes){
        Node current_node = m_mesh->get<Node>(n_id);
        std::vector<Face> current_faces = current_node.get<Face>();
        // FIRST FACE
        Face f0 = current_faces[0];
        Node f0_n0;
        Node f0_n1;
        f0.getAdjacentNodes(current_node,f0_n0,f0_n1);
        stencil current_stencil;
        current_stencil.val[1][1] = current_node.id();
        current_stencil.val[1][0] = f0_n0.id();
        current_stencil.val[0][1] = f0_n1.id();
        std::vector<TCellID> f_nodes=f0.getIDs<Node>();
        for(auto f_n_id:f_nodes){
            if(f_n_id!=current_node.id() &&
               f_n_id!=f0_n0.id() &&
               f_n_id!=f0_n1.id()) {
                current_stencil.val[0][0] = f_n_id;
            }
        }
        //SECOND FACE
        std::vector<TCellID> faces_11_10 = m_mesh->getCommonFaces(current_node,f0_n0);
        TCellID f1_id = NullID;
        if(faces_11_10[0]==f0.id()){
            f1_id=faces_11_10[1];
        }
        else{
            f1_id=faces_11_10[0];
        }
        Face f1 = m_mesh->get<Face>(f1_id);
        Node f1_n0;
        Node f1_n1;
        f1.getAdjacentNodes(current_node,f1_n0,f1_n1);
        if(f1_n0.id()==current_stencil.val[1][0])
            current_stencil.val[2][1]=f1_n1.id();
        else
            current_stencil.val[2][1]=f1_n0.id();

        f_nodes=f1.getIDs<Node>();
        for(auto f_n_id:f_nodes){
            if(f_n_id!=current_node.id() &&
               f_n_id!=f1_n0.id() &&
               f_n_id!=f1_n1.id()) {
                current_stencil.val[2][0] = f_n_id;
            }
        }
        //THIRD FACE
        std::vector<TCellID> faces_11_21 =
                m_mesh->getCommonFaces(current_node,
                                       m_mesh->get<Node>(current_stencil.val[2][1]));
        TCellID f2_id = NullID;
        if(faces_11_21[0]==f1.id()){
            f2_id=faces_11_21[1];
        }
        else{
            f2_id=faces_11_21[0];
        }
        Face f2 = m_mesh->get<Face>(f2_id);
        Node f2_n0;
        Node f2_n1;
        f2.getAdjacentNodes(current_node,f2_n0,f2_n1);

        if(f2_n0.id()==current_stencil.val[2][1])
            current_stencil.val[1][2]=f2_n1.id();
        else
            current_stencil.val[1][2]=f2_n0.id();

        f_nodes=f2.getIDs<Node>();
        for(auto f_n_id:f_nodes){
            if(f_n_id!=current_node.id() &&
               f_n_id!=f2_n0.id() &&
               f_n_id!=f2_n1.id()) {
                current_stencil.val[2][2] = f_n_id;
            }
        }

        //FOURTH FACE
        std::vector<TCellID> faces_11_12 =
                m_mesh->getCommonFaces(current_node,
                                       m_mesh->get<Node>(current_stencil.val[1][2]));
        TCellID f3_id = NullID;
        if(faces_11_12[0]==f2.id()){
            f3_id=faces_11_12[1];
        }
        else{
            f3_id=faces_11_12[0];
        }
        Face f3 = m_mesh->get<Face>(f3_id);

        f_nodes=f3.getIDs<Node>();
        for(auto f_n_id:f_nodes){
            if(f_n_id!=current_node.id() &&
               f_n_id!=current_stencil.val[1][2] &&
               f_n_id!=current_stencil.val[0][1]) {
                current_stencil.val[0][2] = f_n_id;
            }
        }

        m_stencil.insert(std::make_pair(n_id, current_stencil));
		  // Affichage du stencil pour chaque noeud intérieur
		  /*
        std::cout<<"Node "<<n_id<<std::endl;
        for(auto i=0;i<3;i++){
            for(auto j=0;j<3;j++){
                std::cout << m_stencil[n_id].val[i][j] << " ";
            }
            std::cout<<std::endl;
        }
        */
    }
}
/*------------------------------------------------------------------------*/


/*
// ------------ PERTURBATION ALEATOIRE DU MAILLAGE ------------
void Smooth2D::PerturbationMaillage(const Variable<int>* var_bnd, const double dx, const double dy) {
	// Initialisation de la variable n qui est un noeud et n_coords qui récupère les coordonnées du noeud
	Node n = m_mesh->get<Node>(1);
	math::Point n_coords = n.point();
	double c1 = 0.0;
	double c2 = 0.0;
	for(auto id:m_mesh->nodes()){
		if(var_bnd->value(id)==0){
			n = m_mesh->get<Node>(id);
			n_coords = n.point();
			c1 = rand()/ double(RAND_MAX) ;
			c2 = rand()/ double(RAND_MAX) ;
			//std::cout << "Numéro node :" << n << std::endl ;
			//std::cout << "Nombre aléatoire c1 =" << c1 << std::endl ;
			n.setX( n_coords.X() + (2.0*c1-1.0) * dx);
			n.setY( n_coords.Y() + (2.0*c2-1.0) * dy);
		}
	}
	// ---------------------------------------------------------------
}
*/




/*------------------------------------------------------------------------*/
// Fonction FindMidBranche : Si on considère une branche composée de 3 points, cette fonction retourne le point
//                            positionné au milieu de cette branche
// En entrée : A, B, C -> 3 points
// En sortie : le point milieu
math::Point Smooth2D::FindMidBranche(const math::Point A, const math::Point B, const math::Point C) {
	math::Point Point_Milieu ;
	math::Vector3d Vec_AB = B-A ;
	math::Vector3d Vec_BC = C-B ;
	double norme_1 = Vec_AB.norm() ;
	double norme_2 = Vec_BC.norm() ;
	double norme_branche = norme_1 + norme_2 ;
	double norme_milieu = norme_branche / 2.0 ;

	if (norme_milieu <= norme_1){
		Vec_AB.normalize();
		Point_Milieu = A + norme_milieu*Vec_AB ;
	}
	else if (norme_milieu > norme_1){
		math::Vector3d Vec_CB = - Vec_BC ;
		Vec_CB.normalize();
		Point_Milieu = C + norme_milieu*Vec_CB ;
	}

	return Point_Milieu;
}
/*------------------------------------------------------------------------*/






/*------------------------------------------------------------------------*/
// Fonction CheckStructuredMesh : Prend un maillage en entrée et vérifie si il est structuré ou non
// En entrée :
// En sortie :
bool Smooth2D::CheckStructuredMesh() {
	bool checkmesh(true);
	// Check if each inner node have 4 adjacent faces
	for(auto n_id:m_free_nodes) {
		if (checkmesh) {
			Node current_node = m_mesh->get<Node>(n_id);
			std::vector<Face> current_faces = current_node.get<Face>();
			int Nbr_Faces = current_faces.size();
			if (Nbr_Faces != 4) {
				checkmesh = false;
			}
		}
	}

	// Check if each face is a quad
	for(auto n_id:m_mesh->faces()){
		if (checkmesh) {
			Face current_face = m_mesh->get<Face>(n_id);
			std::vector<Node> current_nodes = current_face.get<Node>();
			int Nbr_Noeuds = current_nodes.size();
			if (Nbr_Noeuds != 4) {
				checkmesh = false;
			}
		}
	}

	//std::cout << "Maillage structuré ? " << checkmesh << std::endl ;
	return checkmesh;
}
/*------------------------------------------------------------------------*/







/*------------------------------------------------------------------------*/
void Smooth2D::write_debug_txt(int n_id, const Variable<math::Point> *old_coords,
                          math::Point H1, math::Point H2, math::Point H3,
                          math::Point V1, math::Point V2, math::Point V3,
                          math::Point Point_Intersection,
                          std::string AFileName){
	//===============================================================
	// First, we create the file where we are going to store the info
	std::ofstream stream= std::ofstream(AFileName, std::ios::out);
	//set the numerical precision (number of digits)
	stream.precision(15);
	//Header indicating which type of file it is
	stream << "# Claire PHD Debug Version 1.0\n\n";

	//Write the id of the node
	stream << "# NODE " << n_id << "\n";

	// For each node, I write the coordinates of the stencil
	stream << "# COORDINATES OF THE STENCIL " << "\n";
	for(auto i=0;i<3;i++){
		for(auto j=0;j<3;j++){
			Node n_local = m_mesh->get<Node>(m_stencil[n_id].val[i][j]);
			//math::Point point_local = n_local.point();
			math::Point point_local = old_coords->value(n_local.id()) ;
			stream << point_local.X() << " " << point_local.Y() << " " << point_local.Z() << "\n";
		}
		stream << "\n" ;
	}

	for(auto j=0;j<3;j++){
		for(auto i=0;i<3;i++){
			Node n_local = m_mesh->get<Node>(m_stencil[n_id].val[i][j]);
			//math::Point point_local = n_local.point();
			math::Point point_local = old_coords->value(n_local.id()) ;
			stream << point_local.X() << " " << point_local.Y() << " " << point_local.Z() << "\n";
		}
		stream << "\n" ;
	}
	stream << "\n" ;

	// Write the 6 nodes for the Line-Sweeping algorithm
	stream << "# FIRST BRANCH OF MID POINTS " << "\n";
	stream << H1.X() << " " << H1.Y() << " " << H1.Z() << "\n";
	stream << H2.X() << " " << H2.Y() << " " << H2.Z() << "\n";
	stream << H3.X() << " " << H3.Y() << " " << H3.Z() << "\n";
	stream << "\n";
	stream << "\n";

	stream << "# SECOND BRANCH OF MID POINTS " << "\n";
	stream << V1.X() << " " << V1.Y() << " " << V1.Z() << "\n";
	stream << V2.X() << " " << V2.Y() << " " << V2.Z() << "\n";
	stream << V3.X() << " " << V3.Y() << " " << V3.Z() << "\n";
	stream << "\n";
	stream << "\n";

	// Write the intersection point
	stream << "# INTERSECTION POINT " << "\n";
	stream << Point_Intersection.X() << " " << Point_Intersection.Y() << " " << Point_Intersection.Z() << "\n";
	stream << "\n";
	stream << "\n";

	stream << "# FOR GNUPLOT :" << "\n";
	stream << "plot 'test.txt' index 0 w lp linecolor rgb 'black' title 'stencil' " << "\n";
	stream << "replot 'test.txt' index 1 w lp linecolor rgb 'red' title 'first mid branch' " << "\n";
	stream << "replot 'test.txt' index 2 w lp linecolor rgb 'green' title 'second mid branch' " << "\n";
	stream << "replot 'test.txt' index 3 w lp linecolor rgb 'pink' title 'intersection point' " << "\n";
	stream << "\n";


	stream.close();
}
/*------------------------------------------------------------------------*/
