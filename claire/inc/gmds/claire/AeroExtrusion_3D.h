//
// Created by rochec on 18/11/2022.
//

#ifndef GMDS_AEROEXTRUSION_3D_H
#define GMDS_AEROEXTRUSION_3D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include "gmds/ig/Mesh.h"
#include <gmds/claire/AeroException.h>
#include <gmds/claire/Front_3D.h>
#include <gmds/claire/Params.h>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_CLAIRE_API AeroExtrusion_3D
{
 public:
	/*--------------------------------------------------------------------*/
	/** @enum  Status code for executing algorithms
	 */
	typedef enum {
		FAIL,
		SUCCESS
	} STATUS;

	/*-------------------------------------------------------------------*/
	/** @brief Constructor.
         *  @param[in] AMeshT the triangular mesh where we work on
         *  @param[in] AMeshH the quad mesh to generate
         *  @param[in] Aparams_aero parameters for aero algorithm
         *  @param[in] A_DistanceField distance field for extrusion
         *  @param[in] A_VectorField vector field for extrusion
         *
	 */
	AeroExtrusion_3D(Mesh *AMeshT, Mesh *AMeshH, ParamsAero Aparams_aero, Variable<double>* A_DistanceField, Variable<math::Vector3d>* A_VectorField);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();
	/*-------------------------------------------------------------------*/

 private:
	/*-------------------------------------------------------------------*/
	/** @brief Calcule la position idéale du prochain noeud pour chaque
	 	* noeud du front.
	 	* \param[in] AFront the front
   	* \param[in] A_distance the distance variable on the first mesh
		* \param[in] dist_cible the distance wanted for the first layer
		*
		* \return  a map with (TCellID, TCellID) for ideal positions
	 */
	std::map<TCellID, TCellID> ComputeIdealPositions(Front_3D AFront, double dist_cible, Variable<double>* A_distance, Variable<math::Vector3d>* A_vectors);
	/*-------------------------------------------------------------------*/
	/** @brief Construit une couche de mailles à partir d'un front. Ici,
	 	* des mailles peuvent être fusionnées ou insérées.
	 	* \param[in] Front_IN the front before
   	* \param[in] A_distance the distance variable on the first mesh
		* \param[in] dist_cible the distance wanted for the first layer
	 	* \param[in] A_vectors le champ de vecteurs à utiliser
		*
		* \return  the front computed
	 */
	Front_3D ComputeLayer(Front_3D Front_IN, Variable<double>* A_distance, double dist_cible, Variable<math::Vector3d>* A_vectors);
	/*-------------------------------------------------------------------*/
	/** @brief Construit la première couche de blocs. Pour cette couche,
	 	* les conditions sont particulières.
	 	* \param[in] A_Front_IN front faces and nodes
   	* \param[in] A_distance the distance variable on the first mesh
	 	* \param[in] A_vectors le champ de vecteurs à utiliser
		*
		* \return  the first front computed
	 */
	Front_3D Compute1stLayer(Front_3D A_Front_IN, Variable<double>* A_distance, Variable<math::Vector3d>* A_vectors);
	/*-------------------------------------------------------------------*/
	/** @brief Créé un hax normal sur la couche à partir d'une face
	 	* \param[in] f_id la face concernée
	 	* \param[in] Front_IN front en entrée
		*
		* \return
	 */
	void CreateNormalHexa(TCellID f_id, Front_3D &Front_IN, std::map<TCellID, TCellID> map_new_nodes);
	/*-------------------------------------------------------------------*/
	/** @brief Classification of one edge of the front
	 	* \param[in] e_id l'arête concernée
	 	* \param[in] Front front en entrée
		*
		* \return the edge classification
	 */
	int SingleEdgeClassification(TCellID e_id, Front_3D &Front);
	/*-------------------------------------------------------------------*/
	/** @brief Classification of all the edges on the front
	 	* \param[in] Front front en entrée
		*
		* \return the edge classification
	 */
	Variable<int>* FrontEdgesClassification(Front_3D &Front);
	/*-------------------------------------------------------------------*/
 private:
	/** triangular mesh we work on */
	Mesh *m_meshT;
	/** quad mesh to generate */
	Mesh *m_meshH;
	/** Params pour l'aéro */
	ParamsAero m_params_aero;
	/** Distance Field for extrusion */
	Variable<double>* m_DistanceField;
	/** Vector Field for extrusion */
	Variable<math::Vector3d>* m_VectorField;

};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_AEROEXTRUSION_3D_H
