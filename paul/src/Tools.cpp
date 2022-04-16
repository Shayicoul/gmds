//
// Created by Paul Bourmaud on 12/04/2022.
//

#include "gmds/paul/Tools.h"

using namespace gmds;
Tools::Tools(GridBuilderAround *AGrid)
	:g_grid(*AGrid){;}

int Tools::getValueActivateFace(const int faceIDChecked)
{
	int valueActivate;
	valueActivate = g_grid.activate->value(faceIDChecked);
	return valueActivate;

}

bool Tools::checkExistEdge(const int i1, const int i2, const int faceID)
{
	if(Tools::checkFollowIdNode(i1,i2,faceID) && Tools::checkCommonFace(i1,i2)){
		std::cout<<"LES 2 NOEUDS ONT UNE ARETE !"<<std::endl;
		return true;
	}
	else{
		std::cout<<"LES 2 NOEUDS N'ONT PAS UNE ARETE !"<<std::endl;
		return false;
	}

}

bool Tools::checkCommonFace(const int i1, const int i2)
{
	std::vector<Face> list_f= Tools::getListFacesOfNode(i1);
	for (auto f : list_f){
		std::vector<Node> list_n = Tools::getListNodesOfFace(f.id());
		for (auto n : list_n){
			if (n.id() == i2 ){
				//std::cout<<"LES 2 NOEUDS ONT UNE FACE COMMUNE"<<std::endl;
				return true;
			}
		}
	}
	//std::cout<<"LES 2 NOEUDS N'ONT PAS UNE FACE COMMUNE"<<std::endl;
	return false;
}

bool Tools::checkFollowIdNode(const int i1, const int i2,const int faceID)
{
	if (i1 == Tools::getIdNextNode(i2,faceID) || i1 == Tools::getIdPreviousNode(i2,faceID)){
		//std::cout<<"LES 2 NOEUDS SE SUIVENT !"<<std::endl;
		return true;
	}
	else{
		//std::cout<<"LES 2 NOEUDS NE SE SUIVENT PAS !"<<std::endl;
		return false;
	}
}

std::vector<Node> Tools::getListNodesOfFace(const int faceID)
{
	Face f = g_grid.m_mesh.get<Face>(faceID);
	std::vector<Node>list_nodes = f.get<Node>();
	/*for (auto n : list_nodes){
		std::cout<<"Face :"<< faceID <<" avec les noeuds :"<< n<<std::endl;
	}*/
	return list_nodes;
}

std::vector<Face> Tools::getListFacesOfNode(const int nodeID)
{
	Node n =g_grid.m_mesh.get<Node>(nodeID);
	std::vector<Face> list_faces = n.get<Face>();
	/*for (auto n : list_faces){
		std::cout<<"Le noeud " << nodeID << " avec les faces"<< n <<std::endl;
	}*/
	return list_faces;
}

std::vector<Face> Tools::getFacesCommon(const int i1, const int i2)
{
	std::vector<Face> list_Face_Common;
	std::vector<Face> list_f_i1 = Tools::getListFacesOfNode(i1);
	std::vector<Face> list_f_i2 = Tools::getListFacesOfNode(i2);
	for (auto f1 : list_f_i1){
		for (auto f2  : list_f_i2){
			if (f2 == f1 ){
				list_Face_Common.insert(list_Face_Common.end(),f1);
			}
		}
	}
	if(list_Face_Common.empty()){
		std::cout<<"Erreur, noeuds pas de face en commun "<<std::endl;
		return {};
	}
	/*std::cout<<"List des faces en commun entre "<<i1<<" et "<<i2<<" :"<<std::endl;
	for (auto f : list_Face_Common){
		std::cout<< f <<std::endl;
	}*/
	return list_Face_Common;
}

int Tools::getIdOneCommonFace(const int i1, const int i2)
{
	std::vector<Face> list_f= Tools::getListFacesOfNode(i1);
	for (auto f : list_f){
		std::vector<Node> list_n = Tools::getListNodesOfFace(f.id());
		for (auto n : list_n){
			if (n.id() == i2 ){
				//std::cout<<"LES 2 NOEUDS ONT UNE FACE COMMUNE"<<std::endl;
				return f.id();
			}
		}
	}
	//std::cout<<"LES 2 NOEUDS N'ONT PAS UNE FACE COMMUNE"<<std::endl;
	return {};

}

int Tools::getIdPreviousNode(const int idNode,const int idFaceNode)
{
	std::vector<Node> listNodeFace = Tools::getListNodesOfFace(idFaceNode);

	for(int i = 0; i < listNodeFace.size();i++){
		//std::cout<<"L'element i : "<< i <<" du vecteur est :"<<listNodeFace[i].id()<<std::endl;
		//std::cout<<"Le dernier element du vec est "<<listNodeFace.back()<<std::endl;
		if (listNodeFace[i].id() == idNode){
			if(i == 0){
				std::cout << "Le noeud precedent est " << listNodeFace.back().id() << std::endl;
				return listNodeFace.back().id();
			}
			else {
				std::cout << "Le noeud precedent est " << listNodeFace[i - 1].id() << std::endl;
				return listNodeFace[i - 1].id();
			}
		}
	}
}

int Tools::getIdNextNode(const int idNode, const int idFaceNode)
{
	std::vector<Node> listNodeFace = Tools::getListNodesOfFace(idFaceNode);

	for(int i = 0; i < listNodeFace.size();i++){
		//std::cout<<"L'element i : "<< i <<" du vecteur est :"<<listNodeFace[i].id()<<std::endl;
		//std::cout<<"Le dernier element du vec est "<<listNodeFace.back()<<std::endl;
		if (listNodeFace[i].id() == idNode){
			if(i == listNodeFace.size()-1){
				std::cout << "Le noeud suivant est " << listNodeFace.front().id() << std::endl;
				return listNodeFace.front().id();
			}
			else {
				std::cout << "Le noeud suivant est " << listNodeFace[i + 1].id() << std::endl;
				return listNodeFace[i + 1].id();
			}
		}
	}
}

std::vector<std::vector<int>> Tools::getOtherNodes(const int i1, const int i2)
{
	std::vector<Face> listFace = Tools::getFacesCommon(i1,i2);
	std::vector<Node> listOtherNodes;
	std::vector<std::vector<int>> listPairNodes;
	std::vector<int> listElementAdd;

	for(auto f : listFace){
		std::vector<Node> nodesFace = Tools::getListNodesOfFace(f.id());
		for (auto n : nodesFace){
			if(n.id() != i1 && n.id() != i2){
				listOtherNodes.insert(listOtherNodes.end(),n);
			}
		}
	}

	for (auto n : listOtherNodes){
		for(auto n1 : listOtherNodes){
			if (n != n1 && Tools::checkCommonFace(n.id(),n1.id())){
				if (std::any_of(listElementAdd.begin(),listElementAdd.end(),[&](const int& elem) { return elem == n.id(); })
					 && std::any_of(listElementAdd.begin(),listElementAdd.end(),[&](const int& elem) { return elem == n1.id(); })){
					//std::cout<<"DEJA DANS LA LISTE :"<<n<< " & "<<n1<<std::endl;
				}
				else{
					std::vector<int> pairNodes;
					listElementAdd.push_back(n.id());
					listElementAdd.push_back(n1.id());
					pairNodes.insert(pairNodes.end(),n.id());
					pairNodes.insert(pairNodes.end(),n1.id());
					listPairNodes.push_back(pairNodes);
				}
			}
		}
	}
	/*for (auto n : listPairNodes){
		std::cout<<"UNE PAIRE : "<<std::endl;
		for(auto i : n) {
			std::cout << "============"<< "Elements list pair : " << i << "===========" << std::endl;
		}
	}*/
	return listPairNodes;
}

std::vector<Face> Tools::getAllFacesChain(const int i1, const int i2)
{
	std::vector<Face> listAllFaces = Tools::getFacesCommon(i1,i2);
	std::vector<std::vector<int>> listAllPairNodes = Tools::getOtherNodes(i1,i2);
	std::vector<int> nodeCheck = {i1,i2};
	auto listPairToCheck = listAllPairNodes;

	while(listPairToCheck.empty() != 1){
		auto listParcour = listPairToCheck;
		for (auto p : listParcour){
			std::cout<<"LES ID QUI SONT CHECK"<<std::endl;
			for (auto n : p){
				std::cout<<n<<std::endl;
			}
			if (*find(nodeCheck.begin(), nodeCheck.end(), p.front()) == p.front() ||
				 *find(nodeCheck.begin(), nodeCheck.end(), p.back()) == p.back()) {
				std::cout<<"Noeud DEJA Check"<<std::endl;

			}
			else{
				nodeCheck.push_back(p.front());
				nodeCheck.push_back(p.back());

			}
			std::cout<<"Noeuds Check :"<<std::endl;
			for (auto j : nodeCheck){
				std::cout<<j<<std::endl;
			}
			std::vector<std::vector<int>> otherNodes = Tools::getOtherNodes(p.front(),p.back());


			/*for(auto f : Tools::getFacesCommon(p.front(),p.back())){
				listAllFaces.push_back(f);
				std::cout<<"Face ADD : "<<f<<std::endl;
			}*/
			for(auto n : otherNodes){
				if (*find(nodeCheck.begin(), nodeCheck.end(), n.front()) == n.front() ||
					 *find(nodeCheck.begin(), nodeCheck.end(), n.back()) == n.back()) {
					std::cout<<"Other nodes deja check"<<std::endl;

				}
				else{
					listPairToCheck.push_back(n);
				}
			}
			auto newFaces = Tools::getFacesCommon(p.front(),p.back());
			for (auto f : newFaces){

				if (*find(listAllFaces.begin(), listAllFaces.end(),f) == f){
					std::cout<<"Other nodes deja check"<<std::endl;
				}
				else{
					listAllFaces.push_back(f);
				}
			}
			listPairToCheck.erase(std::remove(listPairToCheck.begin(), listPairToCheck.end(), p), listPairToCheck.end());
			for (auto n : listPairToCheck){
				for (auto j : n){
					//std::cout<<"LIST EN COURS DE SUPP : "<<j<<std::endl;
				}
			}

		}
	}
	std::cout<<"LA LISTE DES FACES "<<std::endl;
	for(auto f : listAllFaces){
		std::cout<<f<<std::endl;
	}
	return listAllFaces;
}

std::vector<std::vector<int>> Tools::getAllNodesChain(const int i1, const int i2)
{
	//auto listAllFaces = Tools::getAllFacesChain(i1,i2);
	std::vector<int> pairNodes={i1,i2};
	auto nodeCheck = pairNodes;
	std::vector<std::vector<int>> listAllPairNodesChain= {pairNodes};
	std::vector<std::vector<int>> listPairToCheck = Tools::getOtherNodes(i1,i2);

	while(listPairToCheck.empty()!=1) {
		auto listParcour = listPairToCheck;
		for (auto p : listParcour) {
			listAllPairNodesChain.push_back(p);
			pairNodes.push_back(p.front());
			pairNodes.push_back(p.back());
			auto newPairToAdd = Tools::getOtherNodes(p.front(), p.back());
			for (auto n : newPairToAdd) {
				if ((*find(pairNodes.begin(), pairNodes.end(), n.front()) == n.front()) || (*find(pairNodes.begin(), pairNodes.end(), n.back()) == n.back())) {
					std::cout << "DEJA DANS LA LISTE DES NOEUDS " << n.front() << "" << n.back() << std::endl;
				}
				else {
					listPairToCheck.push_back(n);
				}
				listPairToCheck.erase(std::remove(listPairToCheck.begin(), listPairToCheck.end(), p), listPairToCheck.end());
			}
		}
	}

		for (auto n : listAllPairNodesChain){
			std::cout<<"\nUne Pair : "<<std::endl;
			for (auto i : n){
				std::cout<<i<<std::endl;
			}
		}
		return listAllPairNodesChain;
	}



