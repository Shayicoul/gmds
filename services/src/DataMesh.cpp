/*----------------------------------------------------------------------------*/
#include "gmds/services/DataMesh.h"
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
DataMesh::DataMesh(gmds::Mesh *AMesh, const std::string AName)
        : AbstractData(AName), m_mesh(AMesh){}
/*----------------------------------------------------------------------------*/
Mesh* DataMesh::mesh() const { return m_mesh;}
/*----------------------------------------------------------------------------*/
void DataMesh::set(gmds::Mesh *AM) {m_mesh=AM;}
/*----------------------------------------------------------------------------*/
