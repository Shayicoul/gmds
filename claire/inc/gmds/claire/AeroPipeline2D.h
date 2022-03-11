//
// Created by rochec on 09/02/2022.
//

#ifndef GMDS_AEROPIPELINE2D_H
#define GMDS_AEROPIPELINE2D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/claire/AbstractAeroPipeline.h>
#include <gmds/claire/Params.h>
/*----------------------------------------------------------------------------*/
namespace  gmds {
/*----------------------------------------------------------------------------*/
/** \class  AeroPipeline2D
 *  \brief  Pipeline de génération de maillages 2D pour l'aéro.
 */
class LIB_GMDS_CLAIRE_API AeroPipeline2D: public AbstractAeroPipeline {
 public:

	/*------------------------------------------------------------------------*/
	/** \brief Default constructor
	 */
	AeroPipeline2D(ParamsAero Aparams);

	/*------------------------------------------------------------------------*/
	/** \brief Function to be called for mesh generation
	 */
	virtual void execute();
	/*------------------------------------------------------------------------*/

 private:
	/*----------------------------------------------------------------------------*/
	/** @brief Lecture du fichier de maillage au format .vtk
	 */
	void LectureMaillage();
	/*----------------------------------------------------------------------------*/
	/** @brief Retire les noeuds qui ne sont connectés à rien
	 */
	void MeshCleaner();
	/*----------------------------------------------------------------------------*/
	/** @brief Ecriture du fichier de maillage au format .vtk
	 */
	void EcritureMaillage(Mesh* p_mesh);
	/*----------------------------------------------------------------------------*/
	/** @brief Initialisation des marques sur les fronts
	 */
	void InitialisationFronts();
	/*----------------------------------------------------------------------------*/
	/** @brief Initialisation du maillage quad généré
	 */
	void InitialisationMeshGen();
	/*----------------------------------------------------------------------------*/
	/** @brief Initialisation du maillage quad généré, test 2
	 */
	void InitialisationMeshParoi();
	/*----------------------------------------------------------------------------*/
	/** @brief Génère une couche de noeuds du maillage
	 */
	void GenerationCouche(int couche_id, double dist);
	/*----------------------------------------------------------------------------*/
	/** @brief Donne le noeud de la couche i construit à partir du noeud n0 de la
	 * couche i-1. On suppose ici, pour l'instant, qu'il est unique.
	 */
	Node SuccessorNode(Node n0);
	/*----------------------------------------------------------------------------*/
	/** @brief Donne le noeud de la couche i construit à partir du noeud n0 de la
	 * couche i+1. On suppose ici, pour l'instant, qu'il est unique.
	 */
	Node AnteriorNode(Node n0);
	/*----------------------------------------------------------------------------*/
	/** @brief Donne le vecteur des noeuds adjacents à n0 qui sont dans la couche i.
	 */
	std::vector<Node> AdjNodesInLayer(Node n0, int couche_i);
	/*----------------------------------------------------------------------------*/
	/** @brief Créé une face de type quad dans la couche avec les noeuds n0 de couche i,
	 * n1 son antécédant dans la couche i-1, n2 un noeud adjacent à n1 dans la couche i-1.
	 * Le dernier noeud n3 de la couche i est obtenu à l'aide de n2.
	 */
	void CreateQuadAndConnectivities(Node n0, Node n1, Node n2);
	/*----------------------------------------------------------------------------*/
	/** @brief Vérifie si une face de type quad est créée. Face correspondant
	 * aux noeuds n0 de couche i, n1 son antécédant dans la couche i-1, n2 un noeud
	 *  adjacent à n1 dans la couche i-1. Le dernier noeud n3 de la couche i est
	 *  obtenu à l'aide de n2.
	 */
	bool isQuadCreated(Node n0, Node n1, Node n2);
	/*----------------------------------------------------------------------------*/
	/** @brief Retourne la longueur d'un bord (périmètre) à partir d'un noeud au
	 * hasard
	 */
	double computeBoundaryLength(Node n0);
	/*----------------------------------------------------------------------------*/
 protected:
	/** mesh we work on */
	Mesh m_m;
	/** mesh quad generated */
	Mesh m_mGen;

};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif     // GMDS_AEROPIPELINE2D_H
