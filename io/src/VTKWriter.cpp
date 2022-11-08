/*----------------------------------------------------------------------------*/
//
// Created by ledouxf on 1/22/19.
//
/*----------------------------------------------------------------------------*/
#include <gmds/io/VTKWriter.h>
#include <fstream>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
VTKWriter::VTKWriter(IMeshIOService *AMeshService)
        :IWriter(AMeshService)
{}
/*----------------------------------------------------------------------------*/
VTKWriter::~VTKWriter()
{}
/*----------------------------------------------------------------------------*/
void VTKWriter::initialize(const std::string &AFileName) {

    m_stream= new std::ofstream(AFileName, std::ios::out);
    if (!m_stream){
        std::string s ="Impossible to create a VTK File (Legacy format): "+AFileName;
        throw GMDSException(s);
    }

    *m_stream << "# vtk DataFile Version 3.0\n";
    *m_stream << "Generated by GMDS VTK Writer\n\n";
    *m_stream << "ASCII\n";
    *m_stream << "DATASET UNSTRUCTURED_GRID\n";

}
/*----------------------------------------------------------------------------*/
void VTKWriter::writeNodes() {
    std::vector<IMeshIOService::NodeInfo> nodes_info;
    m_mesh_service->getNodes(nodes_info);
    *m_stream  << "POINTS " << nodes_info.size() << " float\n";

	(*m_stream).precision(15);
    auto vtk_node_id=0;
    for (auto info : nodes_info) {
        math::Point p = info.point;
        *m_stream << p.X() <<  " " << p.Y() << " " << p.Z() << "\n";
        m_node_ids_mapping[info.id] = vtk_node_id++;
    }
	(*m_stream).precision(6);
}
/*----------------------------------------------------------------------------*/
void VTKWriter::writeEdges() {
    m_mesh_service->getEdges(m_edges_info);
}
/*----------------------------------------------------------------------------*/
void VTKWriter::writeFaces()
{
    m_mesh_service->getFaces(m_faces_info);
}
/*----------------------------------------------------------------------------*/
void VTKWriter::writeRegions()
{
    m_mesh_service->getRegions(m_regions_info);
}
/*----------------------------------------------------------------------------*/
void VTKWriter::writeDataNodes() {

    IMeshIOService::DataID data_id;
    std::vector<IMeshIOService::DataInt> data_int;
    std::vector<IMeshIOService::DataReal> data_real;
    std::vector<IMeshIOService::DataVector> data_vec;
    m_mesh_service->getDataNodes(data_id,data_int, data_real, data_vec);


    int nb_elems = data_id.values.size();

    if(nb_elems==0)
        return;

    m_data_node_stream << "POINT_DATA " << nb_elems << "\n";

    /* FIRST WE PUT THE GMDS IDS*/
    m_data_node_stream<< "SCALARS GMDS_ID int 1\n";
    m_data_node_stream<< "LOOKUP_TABLE default\n";

    for (auto i: data_id.values){
        m_data_node_stream<< i.second<<" ";
    }
    m_data_node_stream<<"\n";

    for(auto d: data_int){

        m_data_node_stream<< "SCALARS " << d.name.c_str() << " int 1\n";
        m_data_node_stream<< "LOOKUP_TABLE default\n";
        for (auto id : d.values)
            m_data_node_stream << id.second<<" ";

        m_data_node_stream << "\n";
    }
    for(auto d: data_real){

        m_data_node_stream<< "SCALARS " << d.name.c_str() << " float 1\n";
        m_data_node_stream << "LOOKUP_TABLE default\n";
        for (auto id : d.values)
            m_data_node_stream << id.second<<" ";

        m_data_node_stream << "\n";
    }
    for(auto d: data_vec){

        m_data_node_stream<< "VECTORS "<< d.name.c_str() << " float\n";
        for (auto v : d.values)
            m_data_node_stream << v.second.X()<<" "<< v.second.Y()<<" "<< v.second.Z()<<" ";

        m_data_node_stream << "\n";
    }
}
/*----------------------------------------------------------------------------*/
void VTKWriter::writeDataEdges()  {

    m_mesh_service->getDataEdges(m_cell_id[0],
                                 m_cell_var_int[0],
                                 m_cell_var_real[0],
                                 m_cell_var_vec[0]);

}
/*----------------------------------------------------------------------------*/
void VTKWriter::writeDataFaces()
{
    m_mesh_service->getDataFaces(m_cell_id[1],
                                 m_cell_var_int[1],
                                 m_cell_var_real[1],
                                 m_cell_var_vec[1]);
}
/*----------------------------------------------------------------------------*/
void VTKWriter::writeDataRegions()
{
    m_mesh_service->getDataRegions(m_cell_id[2],
                                   m_cell_var_int[2],
                                   m_cell_var_real[2],
                                   m_cell_var_vec[2]);
}
/*----------------------------------------------------------------------------*/
void VTKWriter::finalize() {
    int nb_cells = m_edges_info.size()+m_faces_info.size()+m_regions_info.size();

    if(nb_cells==0)
        return;

    *m_stream << "\n";
    *m_stream << "CELLS ";
    *m_stream << nb_cells << " ";


    int vtk_cell_size = 3*m_edges_info.size();
    // 3 because 1 for saying that an edge cell has 2 nodes + the two node ids
    // so something like "2 id1 id2"

    for(auto info : m_faces_info){
        vtk_cell_size += 1 + info.node_ids.size();
    }
    for(auto info : m_regions_info){
        vtk_cell_size += 1 + info.node_ids.size();
    }

    *m_stream << vtk_cell_size << "\n";

    int vtk_cell_id = 0;
    for(auto info:m_edges_info){
        *m_stream << "2 ";
        *m_stream << m_node_ids_mapping[info.node_ids[0]] << " ";
        *m_stream << m_node_ids_mapping[info.node_ids[1]] << "\n";
        m_cell_ids_mapping[info.id] = vtk_cell_id++;
    }

    for(auto info:m_faces_info){
        *m_stream << info.node_ids.size()<<" ";
        for(auto nid : info.node_ids){
            *m_stream<< m_node_ids_mapping[nid] << " ";
        }
        *m_stream << "\n";
        m_cell_ids_mapping[info.id] = vtk_cell_id++;
    }

    for(auto info:m_regions_info){
        *m_stream << info.node_ids.size()<<" ";
        for(auto nid : info.node_ids){
            *m_stream << m_node_ids_mapping[nid] << " ";
        }
        *m_stream << "\n";
        m_cell_ids_mapping[info.id] = vtk_cell_id++;
    }
    *m_stream << "\n";
    *m_stream << "CELL_TYPES " << nb_cells << "\n";

    for(auto info: m_edges_info) {
        *m_stream << "2 \n";
    }
    for(auto info:m_faces_info){
        if (info.type == GMDS_TRIANGLE)
            *m_stream << "5 \n";
        else if (info.type == GMDS_QUAD)
            *m_stream << "9 \n";
        else {
            *m_stream << "7 \n";
        }
    }
    for(auto info:m_regions_info){
        if (info.type == GMDS_TETRA)
            *m_stream << "10 \n";
        else if (info.type == GMDS_HEX)
            *m_stream << "12 \n";
        else if (info.type == GMDS_PRISM3)
            *m_stream << "13 \n";
        else if (info.type == GMDS_PYRAMID)
            *m_stream << "14 \n";
        else
            throw GMDSException("VTKWriter::writeCellRegions cell type not handled.");
    }


    *m_stream << m_data_node_stream.str()<<"\n\n";

    //NOW WE WRITE EDGE, FACE, REGION INFO
    auto nb_edges  = m_cell_id[0].values.size();
    auto nb_faces  = m_cell_id[1].values.size();
    auto nb_regions= m_cell_id[2].values.size();

    if(nb_edges+nb_faces+nb_regions==0)
        return;

    *m_stream << "CELL_DATA " << nb_edges+nb_faces+nb_regions<< "\n";

    /* FIRST WE PUT THE GMDS IDS*/
    *m_stream << "SCALARS GMDS_ID int 1\n";
    *m_stream << "LOOKUP_TABLE default\n";
    for(auto dim=0; dim<3; dim++)
        for (auto i:m_cell_id[dim].values){
            *m_stream << i.second << " ";
        }
     *m_stream<<"\n\n";

    //======================================================
    //DATA FIELD FOR EDGES
    for(auto d: m_cell_var_int[0]){

        *m_stream<< "SCALARS " << d.name.c_str() << " int 1\n";
        *m_stream<< "LOOKUP_TABLE default\n";

        for (auto id : d.values)
            *m_stream << id.second<<" ";

        for(unsigned int idle=0; idle<nb_faces+nb_regions;idle++)
            *m_stream<<"-1 ";

        *m_stream<<"\n\n";
    }
    for(auto d:  m_cell_var_real[0]){

        *m_stream<< "SCALARS " << d.name.c_str() << " float 1\n";
        *m_stream << "LOOKUP_TABLE default\n";

        for (auto id : d.values)
            *m_stream << id.second<<" ";

        for(unsigned int idle=0; idle<nb_faces+nb_regions;idle++)
            *m_stream<<"-1.0 ";

        *m_stream<<"\n\n";
    }

    for(auto d:  m_cell_var_vec[0]){

        *m_stream<< "VECTORS " << d.name.c_str() << " float\n";
        for (auto v : d.values)
            *m_stream << v.second.X()<<" "<< v.second.Y()<<" "<< v.second.Z()<<"\n";

        for(unsigned int idle=0; idle<nb_faces+nb_regions;idle++)
            *m_stream<<"0 0 0\n";

        *m_stream<<"\n\n";
    }
    //======================================================
    //DATA FIELD FOR FACES
    for(auto d: m_cell_var_int[1]){

        *m_stream<< "SCALARS " << d.name.c_str() << " int 1\n";
        *m_stream<< "LOOKUP_TABLE default\n";
        for(unsigned int idle=0; idle<nb_edges;idle++)
            *m_stream<<"-1 ";

        for (auto id : d.values)
            *m_stream << id.second<<" ";

        for(unsigned int idle=0; idle<nb_regions;idle++)
            *m_stream<<"-1 ";

        *m_stream<<"\n\n";
    }
    for(auto d:  m_cell_var_real[1]){
        *m_stream<< "SCALARS " << d.name.c_str() << " float 1\n";
        *m_stream << "LOOKUP_TABLE default\n";


        for(unsigned int idle=0; idle<nb_edges;idle++)
            *m_stream<<"-1.0 ";

        for (auto id : d.values)
            *m_stream << id.second<<" ";

        for(unsigned int idle=0; idle<nb_regions;idle++)
            *m_stream<<"-1.0 ";

        *m_stream<<"\n\n";
    }
    for(auto d:  m_cell_var_vec[1]){

        *m_stream<< "VECTORS " << d.name.c_str() << " float\n";
        for(unsigned int idle=0; idle<nb_edges;idle++)
            *m_stream<<"0 0 0\n";

        for (auto v : d.values)
            *m_stream << v.second.X()<<" "<< v.second.Y()<<" "<< v.second.Z()<<"\n";

        for(unsigned int idle=0; idle<nb_regions;idle++)
            *m_stream<<"0 0 0\n";

        *m_stream<<"\n\n";
    }

    //======================================================
    //DATA FIELD FOR REGIONS
    for(auto d: m_cell_var_int[2]){

        *m_stream<< "SCALARS " << d.name.c_str() << " int 1\n";
        *m_stream<< "LOOKUP_TABLE default\n";

        for(unsigned int idle=0; idle<nb_faces+nb_edges;idle++)
            *m_stream<<"-1 ";

        for (auto id : d.values)
            *m_stream << id.second<<" ";


        *m_stream<<"\n\n";
    }
    for(auto d:  m_cell_var_real[2]){

        *m_stream<< "SCALARS " << d.name.c_str() << " float 1\n";
        *m_stream << "LOOKUP_TABLE default\n";

        for(unsigned int idle=0; idle<nb_faces+nb_edges;idle++)
            *m_stream<<"-1.0 ";

        for (auto id : d.values)
            *m_stream << id.second<<" ";

        *m_stream<<"\n\n";
    }
    for(auto d:  m_cell_var_vec[2]){

        *m_stream<< "VECTORS " << d.name.c_str() << " float\n";
        for(unsigned int idle=0; idle<nb_faces+nb_edges;idle++)
            *m_stream<<"0 0 0\n";

        for (auto v : d.values)
            *m_stream << v.second.X()<<" "<< v.second.Y()<<" "<< v.second.Z()<<"\n";

        *m_stream<<"\n\n";
    }
    m_stream->close();

}
