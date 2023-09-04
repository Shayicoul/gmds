/*----------------------------------------------------------------------------*/
#include "gmds/blocking/CurvedBlocking.h"
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::blocking;
/*----------------------------------------------------------------------------*/
int CellInfo::m_counter_global_id = 0;

/*----------------------------------------------------------------------------*/
CurvedBlocking::CurvedBlocking(cad::GeomManager *AGeomModel, bool AInitAsBoundingBox) : m_geom_model(AGeomModel) {
    if (AInitAsBoundingBox) {
        TCoord min[3] = {MAXFLOAT, MAXFLOAT, MAXFLOAT};
        TCoord max[3] = {-MAXFLOAT, -MAXFLOAT, -MAXFLOAT};
        std::vector<cad::GeomVolume *> vols;
        m_geom_model->getVolumes(vols);
        for (auto v: vols) {
            TCoord v_min[3], v_max[3];
            v->computeBoundingBox(v_min, v_max);
            for (auto i = 0; i < 3; i++)
                if (v_min[i] < min[i]) min[i] = v_min[i];
            for (auto i = 0; i < 3; i++)
                if (v_max[i] > max[i]) max[i] = v_max[i];
        }
        math::Point p1(min[0], min[1], min[2]);
        math::Point p2(min[0], max[1], min[2]);
        math::Point p3(max[0], max[1], min[2]);
        math::Point p4(max[0], min[1], min[2]);
        math::Point p5(min[0], min[1], max[2]);
        math::Point p6(min[0], max[1], max[2]);
        math::Point p7(max[0], max[1], max[2]);
        math::Point p8(max[0], min[1], max[2]);
        create_block(p1, p2, p3, p4, p5, p6, p7, p8);
    }
}

/*----------------------------------------------------------------------------*/
CurvedBlocking::~CurvedBlocking() {}

/*----------------------------------------------------------------------------*/
GMap3 *
CurvedBlocking::gmap() {
    return &m_gmap;
}

/*----------------------------------------------------------------------------*/
cad::GeomManager *
CurvedBlocking::geom_model() {
    return m_geom_model;
}

/*----------------------------------------------------------------------------*/
CurvedBlocking::Node
CurvedBlocking::create_node(const int AGeomDim, const int AGeomId, const math::Point &APoint) {
    return m_gmap.create_attribute<0>(NodeInfo(AGeomDim, AGeomId, APoint));
}

/*----------------------------------------------------------------------------*/
CurvedBlocking::Edge
CurvedBlocking::create_edge(const int AGeomDim, const int AGeomId) {
    return m_gmap.create_attribute<1>(CellInfo(1, AGeomDim, AGeomId));
}

/*----------------------------------------------------------------------------*/
CurvedBlocking::Face
CurvedBlocking::create_face(const int AGeomDim, const int AGeomId) {
    return m_gmap.create_attribute<2>(CellInfo(2, AGeomDim, AGeomId));
}

/*----------------------------------------------------------------------------*/
CurvedBlocking::Block
CurvedBlocking::create_block(const int AGeomDim, const int AGeomId) {
    return m_gmap.create_attribute<3>(CellInfo(3, AGeomDim, AGeomId));
}

/*----------------------------------------------------------------------------*/
std::vector<CurvedBlocking::Face>
CurvedBlocking::get_faces_of_block(const CurvedBlocking::Block AB) {
    Dart3 d1 = AB->dart();
    std::vector<CurvedBlocking::Face> faces;
    faces.resize(6);

    faces[0] = m_gmap.attribute<2>(d1);
    faces[1] = m_gmap.attribute<2>(m_gmap.alpha<2, 1, 0, 1, 2>(d1));
    faces[2] = m_gmap.attribute<2>(m_gmap.alpha<2>(d1));
    faces[3] = m_gmap.attribute<2>(m_gmap.alpha<1, 0, 1, 2>(d1));
    faces[4] = m_gmap.attribute<2>(m_gmap.alpha<1, 2>(d1));
    faces[5] = m_gmap.attribute<2>(m_gmap.alpha<0, 1, 2>(d1));
    return faces;
}

/*----------------------------------------------------------------------------*/
std::vector<CurvedBlocking::Block>
CurvedBlocking::get_all_blocks() {
    std::vector<CurvedBlocking::Block> blocks;
    for (auto it = m_gmap.attributes<3>().begin(), itend = m_gmap.attributes<3>().end(); it != itend; ++it) {
        blocks.push_back(it);
    }
    return blocks;
}

/*----------------------------------------------------------------------------*/
std::vector<CurvedBlocking::Face>
CurvedBlocking::get_all_faces() {
    std::vector<CurvedBlocking::Face> faces;
    for (auto it = m_gmap.attributes<2>().begin(), itend = m_gmap.attributes<2>().end(); it != itend; ++it) {
        faces.push_back(it);
    }
    return faces;
}

/*----------------------------------------------------------------------------*/
std::vector<TCellID> CurvedBlocking::get_all_id_nodes() {
    auto nodes = get_all_nodes();
    std::vector<TCellID> ids;
    ids.reserve(nodes.size());
    for (auto n: nodes)
        ids.push_back(n->info().topo_id);
    return ids;
}

/*----------------------------------------------------------------------------*/
std::vector<TCellID> CurvedBlocking::get_all_id_edges() {
    auto edges = get_all_edges();
    std::vector<TCellID> ids;
    ids.reserve(edges.size());
    for (auto e: edges)
        ids.push_back(e->info().topo_id);
    return ids;
}

/*----------------------------------------------------------------------------*/
std::vector<TCellID> CurvedBlocking::get_all_id_faces() {
    auto faces = get_all_faces();
    std::vector<TCellID> ids;
    ids.reserve(faces.size());
    for (auto f: faces)
        ids.push_back(f->info().topo_id);
    return ids;
}

/*----------------------------------------------------------------------------*/
std::vector<TCellID> CurvedBlocking::get_all_id_blocks() {
    auto blocks = get_all_blocks();
    std::vector<TCellID> ids;
    ids.reserve(blocks.size());
    for (auto b: blocks)
        ids.push_back(b->info().topo_id);
    return ids;
}

/*----------------------------------------------------------------------------*/
std::vector<CurvedBlocking::Edge>
CurvedBlocking::get_edges_of_node(const Node AN) {
    std::vector<Edge> edges;
    Dart3 d = AN->dart();
    for (auto it = m_gmap.one_dart_per_incident_cell<1, 0>(d).begin(),
                 itend = m_gmap.one_dart_per_incident_cell<1, 0>(d).end(); it != itend; ++it) {
        edges.push_back(m_gmap.attribute<1>(it));
    }
    return edges;
}

/*----------------------------------------------------------------------------*/
std::vector<CurvedBlocking::Face>
CurvedBlocking::get_faces_of_node(const Node AN) {
    std::vector<Face> faces;
    Dart3 d = AN->dart();
    for (auto it = m_gmap.one_dart_per_incident_cell<2, 0>(d).begin(),
                 itend = m_gmap.one_dart_per_incident_cell<2, 0>(d).end(); it != itend; ++it) {
        faces.push_back(m_gmap.attribute<2>(it));
    }
    return faces;
}

/*----------------------------------------------------------------------------*/
std::vector<CurvedBlocking::Edge>
CurvedBlocking::get_edges_of_face(const Face AF) {
    std::vector<Edge> edges;
    Dart3 d = AF->dart();
    for (auto it = m_gmap.one_dart_per_incident_cell<1, 2>(d).begin(),
                 itend = m_gmap.one_dart_per_incident_cell<1, 2>(d).end(); it != itend; ++it) {
        edges.push_back(m_gmap.attribute<1>(it));
    }
    return edges;
}

/*----------------------------------------------------------------------------*/
std::vector<CurvedBlocking::Block>
CurvedBlocking::get_blocks_of_node(const Node AN) {
    std::vector<Block> blocks;
    Dart3 d = AN->dart();
    for (auto it = m_gmap.one_dart_per_incident_cell<3, 0>(d).begin(),
                 itend = m_gmap.one_dart_per_incident_cell<3, 0>(d).end(); it != itend; ++it) {
        blocks.push_back(m_gmap.attribute<3>(it));
    }
    return blocks;
}

/*----------------------------------------------------------------------------*/
std::vector<CurvedBlocking::Face>
CurvedBlocking::get_faces_of_edge(const Edge AE) {
    std::vector<Face> faces;
    Dart3 d = AE->dart();
    for (auto it = m_gmap.one_dart_per_incident_cell<2, 1>(d).begin(),
                 itend = m_gmap.one_dart_per_incident_cell<2, 1>(d).end(); it != itend; ++it) {
        faces.push_back(m_gmap.attribute<2>(it));
    }
    return faces;
}

/*----------------------------------------------------------------------------*/
std::vector<CurvedBlocking::Block>
CurvedBlocking::get_blocks_of_edge(const Edge AE) {
    std::vector<Block> blocks;
    Dart3 d = AE->dart();
    for (auto it = m_gmap.one_dart_per_incident_cell<3, 1>(d).begin(),
                 itend = m_gmap.one_dart_per_incident_cell<3, 1>(d).end(); it != itend; ++it) {
        blocks.push_back(m_gmap.attribute<3>(it));
    }
    return blocks;
}

/*----------------------------------------------------------------------------*/
std::vector<CurvedBlocking::Block>
CurvedBlocking::get_blocks_of_face(const Face AF) {
    std::vector<Block> blocks;
    Dart3 d = AF->dart();
    for (auto it = m_gmap.one_dart_per_incident_cell<3, 2>(d).begin(),
                 itend = m_gmap.one_dart_per_incident_cell<3, 2>(d).end(); it != itend; ++it) {
        blocks.push_back(m_gmap.attribute<3>(it));
    }
    return blocks;
}

/*----------------------------------------------------------------------------*/
std::vector<CurvedBlocking::Edge>
CurvedBlocking::get_all_edges() {
    std::vector<CurvedBlocking::Edge> edges;
    for (auto it = m_gmap.attributes<1>().begin(), itend = m_gmap.attributes<1>().end(); it != itend; ++it) {
        edges.push_back(it);
    }
    return edges;
}

/*----------------------------------------------------------------------------*/
std::vector<CurvedBlocking::Node>
CurvedBlocking::get_all_nodes() {
    std::vector<CurvedBlocking::Node> nodes;
    for (auto it = m_gmap.attributes<0>().begin(), itend = m_gmap.attributes<0>().end(); it != itend; ++it) {
        nodes.push_back(it);
    }
    return nodes;
}

/*----------------------------------------------------------------------------*/
std::vector<CurvedBlocking::Edge>
CurvedBlocking::get_edges_of_block(const CurvedBlocking::Block AB) {
    Dart3 d = AB->dart();
    std::vector<CurvedBlocking::Edge> edges;
    edges.resize(12);

    edges[0] = m_gmap.attribute<1>(d);
    edges[1] = m_gmap.attribute<1>(m_gmap.alpha<1, 0, 1>(d));
    edges[2] = m_gmap.attribute<1>(m_gmap.alpha<2, 1, 0, 1>(d));
    edges[3] = m_gmap.attribute<1>(m_gmap.alpha<1, 0, 1, 2, 1, 0, 1>(d));
    d = m_gmap.alpha<1>(d);
    edges[4] = m_gmap.attribute<1>(d);
    edges[5] = m_gmap.attribute<1>(m_gmap.alpha<1, 0, 1>(d));
    edges[6] = m_gmap.attribute<1>(m_gmap.alpha<2, 1, 0, 1>(d));
    edges[7] = m_gmap.attribute<1>(m_gmap.alpha<1, 0, 1, 2, 1, 0, 1>(d));
    d = m_gmap.alpha<1>(d);
    d = m_gmap.alpha<2, 1>(d);
    edges[8] = m_gmap.attribute<1>(d);
    edges[9] = m_gmap.attribute<1>(m_gmap.alpha<1, 0, 1>(d));
    edges[10] = m_gmap.attribute<1>(m_gmap.alpha<2, 1, 0, 1>(d));
    edges[11] = m_gmap.attribute<1>(m_gmap.alpha<1, 0, 1, 2, 1, 0, 1>(d));

    return edges;
}

/*----------------------------------------------------------------------------*/
std::vector<CurvedBlocking::Node>
CurvedBlocking::get_nodes_of_block(const CurvedBlocking::Block AB) {
    Dart3 d1 = AB->dart();
    std::vector<CurvedBlocking::Node> nodes;
    nodes.resize(8);
    nodes[0] = m_gmap.attribute<0>(d1);
    nodes[1] = m_gmap.attribute<0>(m_gmap.alpha<0, 1>(d1));
    nodes[2] = m_gmap.attribute<0>(m_gmap.alpha<0, 1, 0, 1>(d1));
    nodes[3] = m_gmap.attribute<0>(m_gmap.alpha<1, 0>(d1));

    Dart3 d2 = m_gmap.alpha<2, 1, 0, 1, 2>(d1);
    nodes[4] = m_gmap.attribute<0>(d2);
    nodes[5] = m_gmap.attribute<0>(m_gmap.alpha<0, 1>(d2));
    nodes[6] = m_gmap.attribute<0>(m_gmap.alpha<0, 1, 0, 1>(d2));
    nodes[7] = m_gmap.attribute<0>(m_gmap.alpha<1, 0>(d2));


    return nodes;
}

/*----------------------------------------------------------------------------*/
std::vector<CurvedBlocking::Node>
CurvedBlocking::get_nodes_of_face(const CurvedBlocking::Face AF) {
    Dart3 d1 = AF->dart();
    std::vector<CurvedBlocking::Node> nodes;
    Dart3 d = d1;
    do {
        nodes.push_back(m_gmap.attribute<0>(d));
        d = m_gmap.alpha<0, 1>(d);
    } while (d != d1);
    return nodes;
}

/*----------------------------------------------------------------------------*/
std::vector<CurvedBlocking::Node>
CurvedBlocking::get_nodes_of_edge(const CurvedBlocking::Edge AE) {
    std::vector<CurvedBlocking::Node> nodes;
    nodes.resize(2);
    Dart3 d = AE->dart();
    nodes[0] = m_gmap.attribute<0>(d);
    nodes[1] = m_gmap.attribute<0>(m_gmap.alpha<0>(d));
    return nodes;
}

/*----------------------------------------------------------------------------*/
math::Point
CurvedBlocking::get_center_of_edge(const Edge AE) {
    std::vector<Node> n = get_nodes_of_edge(AE);
    return 0.5 * (n[0]->info().point + n[1]->info().point);
}

/*----------------------------------------------------------------------------*/
math::Point
CurvedBlocking::get_center_of_face(const Face AF) {
    std::vector<Node> nodes = get_nodes_of_face(AF);
    math::Point center(0, 0, 0);
    for (auto n: nodes) {
        center = center + n->info().point;
    }
    return (1.0 / nodes.size()) * center;
}

/*----------------------------------------------------------------------------*/
math::Vector3d
CurvedBlocking::get_normal_of_face(const Face AF) {
    std::vector<Node> nodes = get_nodes_of_face(AF);
    math::Point center = get_center_of_face(AF);
    math::Vector3d n({0, 0, 0});
    for (auto i = 0; i < nodes.size(); i++) {
        Node prev_node = (i == 0) ? nodes[nodes.size() - 1] : nodes[i - 1];
        Node curr_node = nodes[i];
        math::Vector3d vp = prev_node->info().point - center;
        math::Vector3d vc = curr_node->info().point - center;
        n = n + vp.cross(vc);
    }
    n.normalize();
    return n;
}

/*----------------------------------------------------------------------------*/
math::Point
CurvedBlocking::get_center_of_block(const Block AB) {
    std::vector<Node> nodes = get_nodes_of_block(AB);
    math::Point center(0, 0, 0);
    for (auto n: nodes) {
        center = center + n->info().point;
    }
    return (1.0 / nodes.size()) * center;
}

/*----------------------------------------------------------------------------*/
void
CurvedBlocking::move_node(Node AN, math::Point &ALoc) {
    AN->info().point = ALoc;
    if (AN->info().geom_dim == 0) {
        // classified onto a point
        cad::GeomPoint *geom_pnt = m_geom_model->getPoint(AN->info().geom_id);
        AN->info().point = geom_pnt->point();
    } else if (AN->info().geom_dim == 1) {
        // classified onto a curve
        cad::GeomCurve *geom_curve = m_geom_model->getCurve(AN->info().geom_id);
        geom_curve->project(AN->info().point);
    } else if (AN->info().geom_dim == 2) {
        // classified onto a surface
        cad::GeomSurface *geom_surf = m_geom_model->getSurface(AN->info().geom_id);
        geom_surf->project(AN->info().point);
    }
    // otherwise nothing else to do
}

/*----------------------------------------------------------------------------*/
CurvedBlocking::Block
CurvedBlocking::create_block(
        const math::Point &AP1, const math::Point &AP2,
        const math::Point &AP3, const math::Point &AP4,
        const math::Point &AP5, const math::Point &AP6,
        const math::Point &AP7, const math::Point &AP8) {
    // Initialize attribute for the created hexahedron
    Dart3 d1 = m_gmap.make_combinatorial_hexahedron();

    // 0D attribute
    m_gmap.set_attribute<0>(d1, create_node(1, 4, AP1));
    m_gmap.set_attribute<0>(m_gmap.alpha<0>(d1), create_node(1, 4, AP2));
    m_gmap.set_attribute<0>(m_gmap.alpha<0, 1, 0>(d1), create_node(1, 4, AP3));
    m_gmap.set_attribute<0>(m_gmap.alpha<1, 0>(d1), create_node(1, 4, AP4));
    m_gmap.set_attribute<0>(m_gmap.alpha<2, 1, 0, 1, 2>(d1), create_node(1, 4, AP5));
    m_gmap.set_attribute<0>(m_gmap.alpha<2, 1, 0, 1, 2, 0, 1>(d1), create_node(1, 4, AP6));
    m_gmap.set_attribute<0>(m_gmap.alpha<2, 1, 0, 1, 2, 0, 1, 0, 1>(d1), create_node(1, 4, AP7));
    m_gmap.set_attribute<0>(m_gmap.alpha<2, 1, 0, 1, 2, 1, 0>(d1), create_node(1, 4, AP8));

    // go through all the edges
    for (auto it = m_gmap.one_dart_per_incident_cell<1, 3>(d1).begin(), itend = m_gmap.one_dart_per_incident_cell<1, 3>(
            d1).end(); it != itend; ++it) {
        m_gmap.set_attribute<1>(it, create_edge(4, NullID));
    }

    // go through all the faces
    for (auto it = m_gmap.one_dart_per_incident_cell<2, 3>(d1).begin(), itend = m_gmap.one_dart_per_incident_cell<2, 3>(
            d1).end(); it != itend; ++it) {
        m_gmap.set_attribute<2>(it, create_face(4, NullID));
    }
    Block b = create_block(4, NullID);
    m_gmap.set_attribute<3>(d1, b);
    return b;
}

/*----------------------------------------------------------------------------*/
void
CurvedBlocking::remove_block(CurvedBlocking::Block AB) {
    Dart3 d = AB->dart();
    m_gmap.remove_cell<3>(d);
}

/*----------------------------------------------------------------------------*/
std::tuple<int, int, math::Point>
CurvedBlocking::get_node_info(const int ANodeId) {
    for (auto it = m_gmap.attributes<0>().begin(), itend = m_gmap.attributes<0>().end(); it != itend; ++it) {
        if (it->info().topo_id == ANodeId) {
            return std::make_tuple(it->info().geom_dim, it->info().geom_id, it->info().point);
        }
    }
    return std::make_tuple(-1, -1, math::Point(0, 0, 0));
}

/*----------------------------------------------------------------------------*/
std::tuple<int, int>
CurvedBlocking::get_edge_info(const int AEdgeId) {
    for (auto it = m_gmap.attributes<1>().begin(), itend = m_gmap.attributes<1>().end(); it != itend; ++it) {
        if (it->info().topo_id == AEdgeId) {
            return std::make_tuple(it->info().geom_dim, it->info().geom_id);
        }
    }
    return std::make_tuple(-1, -1);
}

/*----------------------------------------------------------------------------*/
std::tuple<int, int>
CurvedBlocking::get_face_info(const int AFaceId) {
    for (auto it = m_gmap.attributes<2>().begin(), itend = m_gmap.attributes<2>().end(); it != itend; ++it) {
        if (it->info().topo_id == AFaceId) {
            return std::make_tuple(it->info().geom_dim, it->info().geom_id);
        }
    }
    return std::make_tuple(-1, -1);
}

/*----------------------------------------------------------------------------*/
std::tuple<int, int>
CurvedBlocking::get_block_info(const int ABlockId) {
    for (auto it = m_gmap.attributes<3>().begin(), itend = m_gmap.attributes<3>().end(); it != itend; ++it) {
        if (it->info().topo_id == ABlockId) {
            return std::make_tuple(it->info().geom_dim, it->info().geom_id);
        }
    }
    return std::make_tuple(-1, -1);
}

/*----------------------------------------------------------------------------*/
bool
CurvedBlocking::is_valid_topology() const {
    return m_gmap.is_valid();
}

/*----------------------------------------------------------------------------*/
std::string
CurvedBlocking::info() const {
    std::ostringstream mess;
    mess << "Blocking Info: " << std::endl;
    m_gmap.display_characteristics(mess);
    mess << ", validity=" << (is_valid_topology() ? "true" : "false");
    mess << "\nOd-attributes: ";
    for (auto at: m_gmap.attributes<0>()) {
        mess << "[" << at.info().point << ",(dim:" << at.info().geom_dim << ",id:" << at.info().geom_id << ")]; ";
    }
    mess << "\n1d-attributes: ";
    for (auto at: m_gmap.attributes<1>()) {
        mess << "(dim:" << at.info().geom_dim << ",id:" << at.info().geom_id << "); ";
    }
    mess << "\n2d-attributes: ";
    for (auto at: m_gmap.attributes<2>()) {
        mess << "(dim:" << at.info().geom_dim << ",id:" << at.info().geom_id << "); ";
    }
    mess << "\n3d-attributes: ";
    for (auto at: m_gmap.attributes<3>()) {
        mess << "(dim:" << at.info().geom_dim << ",id:" << at.info().geom_id << "); ";
    }
    return mess.str();
}

/*----------------------------------------------------------------------------*/
void
CurvedBlocking::convert_to_mesh(Mesh &AMesh) {
    MeshModel model = AMesh.getModel();
    if (!model.has(N) || !model.has(E) || !model.has(F) || !model.has(R) || !model.has(E2N) || !model.has(F2N) ||
        !model.has(R2N))
        throw GMDSException("Wrong mesh model for block->mesh conversion");

    AMesh.clear();

    Variable<int> *var_node_topo_id = AMesh.newVariable<int, GMDS_NODE>("blocking_topo_id");
    Variable<int> *var_node_topo_dim = AMesh.newVariable<int, GMDS_NODE>("blocking_topo_dim");
    Variable<int> *var_node_geom_id = AMesh.newVariable<int, GMDS_NODE>("blocking_geom_id");
    Variable<int> *var_node_geom_dim = AMesh.newVariable<int, GMDS_NODE>("blocking_geom_dim");

    Variable<int> *var_edge_topo_id = AMesh.newVariable<int, GMDS_EDGE>("blocking_topo_id");
    Variable<int> *var_edge_topo_dim = AMesh.newVariable<int, GMDS_EDGE>("blocking_topo_dim");
    Variable<int> *var_edge_geom_id = AMesh.newVariable<int, GMDS_EDGE>("blocking_geom_id");
    Variable<int> *var_edge_geom_dim = AMesh.newVariable<int, GMDS_EDGE>("blocking_geom_dim");

    Variable<int> *var_face_topo_id = AMesh.newVariable<int, GMDS_FACE>("blocking_topo_id");
    Variable<int> *var_face_topo_dim = AMesh.newVariable<int, GMDS_FACE>("blocking_topo_dim");
    Variable<int> *var_face_geom_id = AMesh.newVariable<int, GMDS_FACE>("blocking_geom_id");
    Variable<int> *var_face_geom_dim = AMesh.newVariable<int, GMDS_FACE>("blocking_geom_dim");

    Variable<int> *var_region_topo_id = AMesh.newVariable<int, GMDS_REGION>("blocking_topo_id");
    Variable<int> *var_region_topo_dim = AMesh.newVariable<int, GMDS_REGION>("blocking_topo_dim");
    Variable<int> *var_region_geom_id = AMesh.newVariable<int, GMDS_REGION>("blocking_geom_id");
    Variable<int> *var_region_geom_dim = AMesh.newVariable<int, GMDS_REGION>("blocking_geom_dim");

    // mapping from blocking node ids to mesh node ids
    std::map<int, TCellID> n2n;

    // nodes
    for (auto at: m_gmap.attributes<0>()) {
        gmds::Node n = AMesh.newNode(at.info().point);
        var_node_topo_id->set(n.id(), at.info().topo_id);
        var_node_topo_dim->set(n.id(), at.info().topo_dim);
        var_node_geom_id->set(n.id(), at.info().geom_id);
        var_node_geom_dim->set(n.id(), at.info().geom_dim);
        n2n[at.info().topo_id] = n.id();
    }
    // edges
    for (auto it = m_gmap.attributes<1>().begin(), itend = m_gmap.attributes<1>().end(); it != itend; ++it) {
        auto att = m_gmap.info_of_attribute<1>(it);
        std::vector<Node> cell_nodes = get_nodes_of_edge(it);
        gmds::Edge e = AMesh.newEdge(n2n[cell_nodes[0]->info().topo_id], n2n[cell_nodes[1]->info().topo_id]);
        var_edge_topo_id->set(e.id(), att.topo_id);
        var_edge_topo_dim->set(e.id(), att.topo_dim);
        var_edge_geom_id->set(e.id(), att.geom_id);
        var_edge_geom_dim->set(e.id(), att.geom_dim);
    }
    // faces
    for (auto it = m_gmap.attributes<2>().begin(), itend = m_gmap.attributes<2>().end(); it != itend; ++it) {
        auto att = m_gmap.info_of_attribute<2>(it);
        std::vector<Node> cell_nodes = get_nodes_of_face(it);
        if (cell_nodes.size() != 4) throw GMDSException("Only quad blocking faces can be converted into mesh");

        gmds::Face f = AMesh.newQuad(n2n[cell_nodes[0]->info().topo_id], n2n[cell_nodes[1]->info().topo_id],
                                     n2n[cell_nodes[2]->info().topo_id],
                                     n2n[cell_nodes[3]->info().topo_id]);

        var_face_topo_id->set(f.id(), att.topo_id);
        var_face_topo_dim->set(f.id(), att.topo_dim);
        var_face_geom_id->set(f.id(), att.geom_id);
        var_face_geom_dim->set(f.id(), att.geom_dim);
    }
    // blocks
    for (auto it = m_gmap.attributes<3>().begin(), itend = m_gmap.attributes<3>().end(); it != itend; ++it) {
        auto att = m_gmap.info_of_attribute<3>(it);
        std::vector<Node> cell_nodes = get_nodes_of_block(it);
        if (cell_nodes.size() != 8) throw GMDSException("Only hex blocks can be converted into mesh");

        gmds::Region r = AMesh.newHex(n2n[cell_nodes[0]->info().topo_id], n2n[cell_nodes[1]->info().topo_id],
                                      n2n[cell_nodes[2]->info().topo_id],
                                      n2n[cell_nodes[3]->info().topo_id], n2n[cell_nodes[4]->info().topo_id],
                                      n2n[cell_nodes[5]->info().topo_id],
                                      n2n[cell_nodes[6]->info().topo_id], n2n[cell_nodes[7]->info().topo_id]);

        var_region_topo_id->set(r.id(), att.topo_id);
        var_region_topo_dim->set(r.id(), att.topo_dim);
        var_region_geom_id->set(r.id(), att.geom_id);
        var_region_geom_dim->set(r.id(), att.geom_dim);
    }
}

/*----------------------------------------------------------------------------*/
void
CurvedBlocking::init_from_mesh(Mesh &AMesh) {
    MeshModel model = AMesh.getModel();
    if (!model.has(N) || !model.has(R) || !model.has(R2N))
        throw GMDSException("Wrong mesh model for mesh->block->mesh conversion");

    //map from gmap node ids to init mesh id
    std::map<TCellID, TCellID> n2n;
    for (auto r_id: AMesh.regions()) {
        Region r = AMesh.get<Region>(r_id);
        //only hex cells are converted
        if (r.type() == GMDS_HEX) {
            std::vector<gmds::Node> ns = r.get<gmds::Node>();
            Block b = create_block(ns[0].point(), ns[1].point(), ns[2].point(), ns[3].point(),
                                   ns[4].point(), ns[5].point(), ns[6].point(), ns[7].point());
            std::vector<Node> b_nodes = get_nodes_of_block(b);
            for (auto i = 0; i < 8; i++) {
                math::Point pi =ns[i].point();
                auto min_dist = pi.distance2(b_nodes[0]->info().point);
                auto min_id = b_nodes[0]->info().topo_id;
                for(auto j=1; j<8;j++){
                    auto dist_j= pi.distance2(b_nodes[j]->info().point);
                    if(dist_j<min_dist){
                        min_dist=dist_j;
                        min_id =b_nodes[j]->info().topo_id;
                    }
                }
                n2n[min_id] = ns[i].id();
            }
        }
    }

    // now we glue blocks;
    for (auto it_r = m_gmap.attributes<3>().begin(), it_rend = m_gmap.attributes<3>().end(); it_r != it_rend; ++it_r)
    {
        std::vector<Face> cell_faces = get_faces_of_block(it_r);
        for (auto f: cell_faces) {
            Dart3 d = f->dart();
            if (m_gmap.is_free<3>(d)) {
                //means the face is not connected
                std::vector<Node> f_nodes = get_nodes_of_face(f);
                //We go through all the 3-free faces and try to connect to f
                bool found_and_glue = false;
                for (auto it_f = m_gmap.attributes<2>().begin(), it_fend = m_gmap.attributes<2>().end();
                     it_f != it_fend; ++it_f) {
                    auto att = m_gmap.info_of_attribute<2>(it_f);
                    Dart3 d2 = it_f->dart();
                    if (d2 != d && m_gmap.is_free<3>(d2)) {
                        //free dart and different face, we check if this face can be connected
                        std::vector<Node> f2_nodes = get_nodes_of_face(it_f);
                        //same nodes?
                        bool one_missing = false;
                        for (auto i = 0; i < f_nodes.size() && !one_missing; i++) {
                            bool found_same = false;
                            for (auto j = 0; j < f2_nodes.size() && !found_same; j++) {
                                if (n2n[f_nodes[i]->info().topo_id] == n2n[f2_nodes[j]->info().topo_id]) {
                                    found_same = true;
                                }
                            }
                            if (!found_same)
                                one_missing = true;
                        }

                        if (!one_missing) {
                            //we have the same ids, we need to glue in the right order
                            found_and_glue = true;
                            //Dart of the initial face d, and the face to glue d2
                            TCellID ref_id = n2n[m_gmap.attribute<0>(d)->info().topo_id];
                            TCellID ref_id_0 = n2n[m_gmap.attribute<0>(m_gmap.alpha<0>(d))->info().topo_id];

                            //starting from d2, we move to find the one that correspond to ref_id
                            Dart3 d_glue = d2;
                            while (n2n[m_gmap.attribute<0>(d_glue)->info().topo_id] != ref_id) {
                                d_glue = m_gmap.alpha<1, 0>(d_glue);
                            }
                            if (n2n[m_gmap.attribute<0>(m_gmap.alpha<0>(d_glue))->info().topo_id] == ref_id_0) {
                                m_gmap.sew<3>(d_glue, d);
                            } else {
                                m_gmap.sew<3>(m_gmap.alpha<1>(d_glue), d);
                            }
                        }
                    }
                }
            }
        }
    }
}

/*----------------------------------------------------------------------------*/
std::vector<std::vector<CurvedBlocking::Edge> >
CurvedBlocking::get_all_sheet_edge_sets() {
    std::vector<std::vector<Edge> > edges;
    std::vector<Edge> all_edges = get_all_edges();
    std::map<TCellID, bool> edge_done;
    for (auto e: all_edges) {
        edge_done[e->info().topo_id] = false;
    }

    bool remain_edge_to_do = true;
    while (remain_edge_to_do) {
        remain_edge_to_do = false;
        //we try and find the first that is not already put into a sheet set
        bool found_edge = false;
        auto edge_index = 0;
        for (auto i = 0; i < all_edges.size() && !found_edge; i++) {
            //we found an edge to treat
            if (edge_done[all_edges[i]->info().topo_id] == false) {
                found_edge = true;
                edge_index = i;
            }
        }
        if (found_edge) {
            //work to do, we will do another iteration
            remain_edge_to_do = true;
            std::vector<Edge> sh_edges;
            get_all_sheet_edges(all_edges[edge_index], sh_edges);
            //we store the sheet edges
            edges.push_back(sh_edges);
            //now we mark them as treated
            for (auto e: sh_edges) {
                edge_done[e->info().topo_id] = true;
            }
        }
    }
    return edges;
}

/*----------------------------------------------------------------------------*/
void
CurvedBlocking::get_all_sheet_edges(const Edge AE, std::vector<Edge> &AEdges) {
    AEdges.clear();
    std::vector<Dart3> sheet_darts;
    get_all_sheet_darts(AE, sheet_darts);
    AEdges.resize(sheet_darts.size());
    for (auto i = 0; i < sheet_darts.size(); i++) {
        AEdges[i] = m_gmap.attribute<1>(sheet_darts[i]);
    }
}

/*----------------------------------------------------------------------------*/
void
CurvedBlocking::get_all_sheet_darts(const Edge AE, std::vector<Dart3> &ADarts) {
    ADarts.clear();
    // we allocate a mark to know all the edges we go through
    auto edge_mark = m_gmap.get_new_mark();

    std::vector<Dart3> front;
    front.push_back(AE->dart());
    // the current dart belongs to the final set of darts
    ADarts.push_back(AE->dart());
    // we mark all the dart of the inital edge to avoid to traverse it twice
    m_gmap.mark_cell<1>(AE->dart(), edge_mark);

    // Now we propagate along topological parallel edges in each adjacent hexahedral cell
    while (!front.empty()) {
        // we pick the last dart of the front
        Dart3 d = front.back();
        front.pop_back();

        // we traverse all the darts of the orbit<2,3> starting from d

        for (GMap3::Dart_of_orbit_range<2, 3>::iterator it(m_gmap.darts_of_orbit<2, 3>(d).begin()), itend(
                m_gmap.darts_of_orbit<2, 3>(d).end()); it != itend;
             ++it) {
            auto d_next_edge = m_gmap.alpha<1, 0, 1>(it);
            if (!m_gmap.is_marked(d_next_edge, edge_mark)) {
                // it means that the edge containing the dart d_next_edge has not been traversed already.
                // We mark the dart of the corresponding edge, and we add it to the front
                front.push_back(d_next_edge);
                m_gmap.mark_cell<1>(d_next_edge, edge_mark);
                // We also add it to the set of darts to return
                ADarts.push_back(d_next_edge);
            }
        }
    }
    // We must unmark all the marked edges. As we stored one dart per edge, it is straightforward
    for (auto d: ADarts) {
        m_gmap.unmark_cell<1>(d, edge_mark);
    }
    m_gmap.free_mark(edge_mark);
}
/*----------------------------------------------------------------------------*/
std::vector<CurvedBlocking::Block> CurvedBlocking::get_all_chord_blocks(const Face AF)
{
  std::vector<CurvedBlocking::Block> bls;
  auto block_mark = m_gmap.get_new_mark();

  std::vector<Dart3> face_darts;
  get_all_chord_darts(AF, face_darts);
  for(auto d: face_darts){
      if(!m_gmap.is_marked(d,block_mark)) {
          bls.push_back(m_gmap.attribute<3>(d));
          m_gmap.mark_cell<3>(d,block_mark);
      }
      if(!m_gmap.is_marked(m_gmap.alpha<3>(d),block_mark)) {
          bls.push_back(m_gmap.attribute<3>(m_gmap.alpha<3>(d)));
          m_gmap.mark_cell<3>(m_gmap.alpha<3>(d),block_mark);
      }
  }

    // We must unmark all the marked edges. As we stored one dart per edge, it is straightforward
    for (auto b: bls) {
        m_gmap.unmark_cell<3>(b->dart(), block_mark);
    }
    m_gmap.free_mark(block_mark);
  return bls;
}
/*----------------------------------------------------------------------------*/
void
CurvedBlocking::get_all_chord_darts(const Face AF, std::vector<Dart3> &ADarts)
{
    ADarts.clear();
    // we allocate a mark to know all the faces we go through
    auto face_mark = m_gmap.get_new_mark();

    std::vector<Dart3> front;
    front.push_back(AF->dart());
    if (!m_gmap.is_free<3>(AF->dart()))
        front.push_back(m_gmap.alpha<3>(AF->dart()));

    // the current dart belongs to the final set of darts
    ADarts.push_back(AF->dart());
    // we mark all the dart of the initial face to avoid to traverse it twice
    m_gmap.mark_cell<2>(AF->dart(), face_mark);

    // Now we propagate along topological parallel faces in each adjacent block
    while (!front.empty()) {
        // we pick the last dart of the front
        Dart3 d = front.back();
        front.pop_back();

        auto d_next_face = m_gmap.alpha<2, 1, 0, 1, 2, 3>(d);
        if (!m_gmap.is_marked(d_next_face, face_mark)) {
            // it means that the edge containing the dart d_next_edge has not been traversed already.
            // We mark the dart of the corresponding edge, and we add it to the front
            front.push_back(d_next_face);
            m_gmap.mark_cell<2>(d_next_face, face_mark);
            // We also add it to the set of darts to return
            ADarts.push_back(d_next_face);
            if (!m_gmap.is_free<3>(d_next_face))
                front.push_back(m_gmap.alpha<3>(d_next_face));

        }
    }
    // We must unmark all the marked edges. As we stored one dart per edge, it is straightforward
    for (auto d: ADarts) {
        m_gmap.unmark_cell<2>(d, face_mark);
    }
    m_gmap.free_mark(face_mark);
}

/*----------------------------------------------------------------------------*/
void
CurvedBlocking::cut_sheet(const Edge AE) {
    cut_sheet(AE, 0.5);
}

/*----------------------------------------------------------------------------*/
void
CurvedBlocking::cut_sheet(const Edge AE, const math::Point &AP) {
    std::vector<Node> n = get_nodes_of_edge(AE);
    math::Point p0 = n[0]->info().point;
    math::Point p1 = n[1]->info().point;
    math::Segment s01(p0, p1);
    math::Point p = s01.project(AP);
    double param = math::Segment(p0, p).computeLength() / s01.computeLength();
    cut_sheet(AE, param);
}

/*----------------------------------------------------------------------------*/
void
CurvedBlocking::cut_sheet(const Edge AE, const double AParam) {
    assert(AParam > 0 && AParam < 1);
    //Note: the parameterization starts from the first node of AE, which is the one containing AE->dart().
    //As a consequence, our orientation is consistent in the algorithm.

    // We get a dart per sheet edge
    std::vector<Dart3> sheet_darts;
    get_all_sheet_darts(AE, sheet_darts);

    // then we insert a node in the middle of each edge, and we keep one of its darts
    std::vector<Dart3> in_edge_darts;
    auto mark_edge_darts = m_gmap.get_new_mark();

    for (auto d: sheet_darts) {
        Dart3 d0 = m_gmap.alpha<0>(d);
        math::Point pa = m_gmap.attribute<0>(d)->info().point;
        math::Point pb = m_gmap.attribute<0>(d0)->info().point;
        auto edge_att = m_gmap.attribute<1>(d)->info();

        Dart3 d_edge = m_gmap.insert_cell_0_in_cell_1(d);     // d_edge = alpha<0,1>(d)
        auto darts_23 = m_gmap.darts_of_orbit<2, 3>(d_edge);
        for (auto d_edge_23 = darts_23.begin(); d_edge_23 != darts_23.end(); d_edge_23++) {
            in_edge_darts.push_back(d_edge_23);
            m_gmap.mark(d_edge_23, mark_edge_darts);
        }

        math::Point pc = (1 - AParam) * pa + AParam * pb;
        if (edge_att.geom_dim == 1) {
            auto curve = m_geom_model->getCurve(edge_att.geom_id);
            curve->project(pc);
        } else if (edge_att.geom_dim == 2) {
            auto surf = m_geom_model->getSurface(edge_att.geom_id);
            surf->project(pc);
        }
        m_gmap.set_attribute<0>(d_edge, create_node(edge_att.geom_dim, edge_att.geom_id, pc));
    }

    // then we cut faces by inserting an edge
    auto mark_done = m_gmap.get_new_mark();
    std::vector<Dart3> in_face_darts;
    for (auto d: in_edge_darts) {
        if (!m_gmap.is_marked(d, mark_done)) {
            // as the node is not marked, it means we didn't insert an edge in this face
            // WARNING wrong pattern for self-intersecting sheet
            Dart3 d2 = m_gmap.alpha<0, 1, 0, 1, 0>(d);
            if (!m_gmap.is_marked(d2, mark_edge_darts)) {
                d2 = m_gmap.alpha<1, 0>(d2);
            }
            auto face_att = m_gmap.attribute<2>(d)->info();

            Dart3 d_face = m_gmap.insert_cell_1_in_cell_2(d, d2);

            m_gmap.set_attribute<1>(d_face, create_edge(face_att.geom_dim, face_att.geom_id));

            in_face_darts.push_back(d_face);
            m_gmap.mark(d, mark_done);
            m_gmap.mark(d2, mark_done);
            //the face can be shared by two regions and so, we must mark the right darts
            if (!m_gmap.is_free<3>(d)) {
                m_gmap.mark(m_gmap.alpha<3>(d), mark_done);
                m_gmap.mark(m_gmap.alpha<3>(d2), mark_done);

            }
        }
    }
    // We unmark all the darts, and we reuse the mark mark_done for faces now
    m_gmap.unmark_all(mark_done);

    // and finally, we cut blocks
    for (auto d: in_face_darts) {
        if (!m_gmap.is_marked(d, mark_done)) {
            // we get the loop of darts we are going to use to cut the block
            // one dart per edge must be given (see the implementation of insert_cell_2_in_cell_3 in CGAL)
            std::vector<Dart3> loop_darts;
            Dart3 current_dart = d;
            do {
                loop_darts.push_back(current_dart);
                m_gmap.mark(current_dart, mark_done);
                m_gmap.mark(m_gmap.alpha<0>(current_dart), mark_done);
                current_dart = m_gmap.alpha<0, 1, 2, 1>(current_dart);
            } while (current_dart != d);
            auto face_att = m_gmap.attribute<2>(d)->info();
            Dart3 d_block = m_gmap.insert_cell_2_in_cell_3(loop_darts.begin(), loop_darts.end());
            m_gmap.set_attribute<2>(d_block, create_face(face_att.geom_dim, face_att.geom_id));
        }
    }
    // we free the mark after unmarking all the darts. This is not optimal.
    m_gmap.unmark_all(mark_edge_darts);
    m_gmap.free_mark(mark_edge_darts);
    m_gmap.unmark_all(mark_done);
    m_gmap.free_mark(mark_done);

}

/*----------------------------------------------------------------------------*/
std::vector<std::pair<double, double> >
CurvedBlocking::get_projection_info(math::Point &AP, std::vector<CurvedBlocking::Edge> &AEdges) {
    std::vector<std::pair<double, double> > dist_coord;
    for (auto e: AEdges) {
        std::vector<Node> end_points = get_nodes_of_edge(e);
        math::Point end0 = end_points[0]->info().point;
        math::Point end1 = end_points[1]->info().point;
        math::Vector3d v1 = end1 - end0;
        math::Vector3d v2 = AP - end0;
        double coord = 0.0;
        double distance = 0.0;
        auto a = v1.dot(v2);
        if (a <= 0.0) {
            coord = 0.0;
            distance = AP.distance(end0);
        } else {
            auto b = v1.dot(v1);
            if (a >= b) {
                coord = 1.0;
                distance = AP.distance(end1);
            } else {
                coord = a / b;
                distance = AP.distance(end0 + coord * v1);
            }
        }
        dist_coord.push_back(std::make_pair(distance, coord));
    }
    return dist_coord;

}

/*----------------------------------------------------------------------------*/
bool
CurvedBlocking::collapse_chord(const Face AF, const Node AN1, const Node AN2) {
    //first we check that AN1 and AN2 belongs to AF and are opposite in AF
    Dart3 dart_f = AF->dart();
    Dart3 dart_n1 = AN1->dart();
    Dart3 dart_n2 = AN2->dart();
    // Are n1 and n2 in the face?
    bool found_n1 = false, found_n2 = false;
    std::vector<Dart3> all_darts_n1, all_darts_n2;
    for (GMap3::Dart_of_orbit_range<1, 2, 3>::iterator it(m_gmap.darts_of_orbit<1, 2, 3>(dart_n1).begin()), itend(
            m_gmap.darts_of_orbit<1, 2, 3>(dart_n1).end()); it != itend; ++it) {
        all_darts_n1.push_back(it);
    }
    for (GMap3::Dart_of_orbit_range<1, 2, 3>::iterator it(m_gmap.darts_of_orbit<1, 2, 3>(dart_n2).begin()), itend(
            m_gmap.darts_of_orbit<1, 2, 3>(dart_n2).end()); it != itend; ++it) {
        all_darts_n2.push_back(it);
    }
        //We traverse all the darts of the face
    for (GMap3::Dart_of_orbit_range<0, 1, 3>::iterator it(m_gmap.darts_of_orbit<0, 1, 3>(dart_f).begin()), itend(
            m_gmap.darts_of_orbit<0, 1, 3>(dart_f).end()); it != itend;
         ++it) {
        //Now we go through the dart of the current node
        for(auto d:all_darts_n1) {
            if (it == d)
                found_n1 = true;
        }
        for(auto d:all_darts_n2) {
            if (it == d)
                found_n2 = true;
        }
    }
    if (!found_n1 || !found_n2)
        return false;

    // Are n1 and n2 opposite in the face ?
    bool n1_is_n2_opposite = false;
    for(auto d2:all_darts_n2) {
        for(auto d1:all_darts_n1) {
            if (m_gmap.alpha<0, 1, 0>(d1) == d2)
                n1_is_n2_opposite = true;
        }
    }
    if (!n1_is_n2_opposite)
        return false;

    // We need the pair of nodes to collapse in each face. To do that, we iteratively work in each
    // face traversed by the sheet
    std::vector<Dart3> chord_darts, chord_darts_opp;
    get_all_chord_darts(AF,chord_darts);
    for(auto d:chord_darts){
        chord_darts_opp.push_back(m_gmap.alpha<0,1,0>(d));
    }

    //Before going further, we check that all the pair of opposite nodes can be merged (geometric conditions)
    for(auto i=0; i<chord_darts.size();i++) {
        Node ni = m_gmap.attribute<0>(chord_darts[i]);
        Node nj = m_gmap.attribute<0>(chord_darts_opp[i]);
        if(ni->info().geom_dim==0 && nj->info().geom_dim==0){
            //two nodes can not be merged together
            return false;
        }
        else if(ni->info().geom_dim==1 && nj->info().geom_dim==1 && ni->info().geom_id!=nj->info().geom_id){
            //two nodes on different curves cannot be merged
            return false;
        }
        else if(ni->info().geom_dim==2 && nj->info().geom_dim==2 && ni->info().geom_id!=nj->info().geom_id){
            //two nodes on different surfaces cannot be merged
            return false;
        }
    }
    //Now we have pairs of darts that belongs to opposite nodes that must be collapsed. We get nodes of adjacent faces
    //to 3-sew
    std::vector<std::pair<Dart3,Dart3> > to_sew;
    for(auto i=0; i<chord_darts.size();i++){
        Dart3 di = chord_darts[i];
        Dart3 dj = chord_darts_opp[i];
        to_sew.push_back(std::make_pair(m_gmap.alpha<2,3>(di),m_gmap.alpha<2,3>(dj)));
        to_sew.push_back(std::make_pair(m_gmap.alpha<1,2,3>(di),m_gmap.alpha<1,2,3>(dj)));
        if(!m_gmap.is_free<3>(di)){
            to_sew.push_back(std::make_pair(m_gmap.alpha<3,2,3>(di),m_gmap.alpha<3,2,3>(dj)));
            to_sew.push_back(std::make_pair(m_gmap.alpha<3,1,2,3>(di),m_gmap.alpha<3,1,2,3>(dj)));
        }
    }
    //We remove chord blocks
    std::vector<Block> chord_blocks = get_all_chord_blocks(AF);
    for(auto b:chord_blocks) {
        remove_block(b);
    }
    //and now, we fill the void by sewing darts
    for(auto dd : to_sew){
        if (m_gmap.is_free<3>(dd.first)){
            m_gmap.sew<3>(dd.first,dd.second);
        }
    }
    return true;
}
