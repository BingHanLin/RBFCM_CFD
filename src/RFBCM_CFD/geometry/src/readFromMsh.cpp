#include "readFromMsh.hpp"
#include "messages.hpp"
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>
bool readFromMsh(const std::string& filePath, std::vector<vec3d<double>>& nodes)
{
    const auto invalidFormat = [](std::string messages) -> bool {
        ASSERT("Invalid format: " << messages);
        return false;
    };
    const auto notImplemented = [](std::string messages) -> bool {
        ASSERT("Not implemented:" << messages);
        return false;
    };

    const auto passWhiteSpace = [](std::ifstream& fin) {
        char next = fin.peek();
        while (next == '\n' || next == ' ' || next == '\t' || next == '\r')
        {
            fin.get();
            next = fin.peek();
        }
    };

    // =========================================================================
    // open file and check
    // =========================================================================
    std::ifstream fileStream(filePath.c_str(), std::ios::in | std::ios::binary);

    if (!fileStream.is_open())
    {
        ASSERT("msh file is not found in path: " << filePath);
        return false;
    }

    // =========================================================================
    // parse header and check
    // =========================================================================
    std::string buf;
    fileStream >> buf;
    if (buf != "$MeshFormat")
    {
        return invalidFormat("msh MeshFormat tag not found");
    }

    double versionNumber;
    int fileType;  // equal to 0 indicates the ASCII file format.
    size_t dataSize;
    fileStream >> versionNumber >> fileType >> dataSize;

    if (fileType != 0)
    {
        return notImplemented("only accept ASCII file format for .msh");
    }

    if (dataSize != 8)
    {
        return invalidFormat("data size must be 8 bytes");
    }

    fileStream >> buf;
    if (buf != "$EndMeshFormat")
    {
        return invalidFormat("msh EndMeshFormat tag not found");
    }

    // =========================================================================
    // parse other tag
    // =========================================================================
    while (!fileStream.eof())
    {
        buf.clear();
        fileStream >> buf;
        if (buf == "$Nodes")
        {
            parseNodes(fileStream, nodes);
            fileStream >> buf;
            if (buf != "$EndNodes")
            {
                return invalidFormat("msh EndNodes tag not found");
            }
        }
        // else if (buf == "$Elements")
        // {
        //     parse_elements(fileStream);
        //     fileStream >> buf;
        //     if (buf != "$EndElements")
        //     {
        //         return invalid_format();
        //     }
        // }
        // else if (buf == "$NodeData")
        // {
        //     parse_node_field(fileStream);
        //     fileStream >> buf;
        //     if (buf != "$EndNodeData")
        //     {
        //         return invalid_format();
        //     }
        // }
        // else if (buf == "$ElementData")
        // {
        //     parse_element_field(fileStream);
        //     fileStream >> buf;
        //     if (buf != "$EndElementData")
        //     {
        //         return invalid_format();
        //     }
        // }
        // else if (fileStream.eof())
        // {
        //     break;
        // }
        // else
        // {
        //     parse_unknown_field(fileStream, buf);
        // }
    }
    fileStream.close();

    // typedef typename DerivedV::Scalar Float;
    // typedef Eigen::Matrix<Float, Eigen::Dynamic, 1> VectorF;
    // typedef Eigen::Matrix<int, Eigen::Dynamic, 1> VectorI;
    // typedef std::map<std::string, VectorF> FieldMap;
    // typedef std::vector<std::string> FieldNames;
    // VectorF m_nodes;
    // VectorI m_elements;
    // FieldMap m_node_fields;
    // FieldMap m_element_fields;

    // bool m_binary;
    // size_t m_data_size;
    // size_t m_nodes_per_element;
    // size_t m_element_type;

    // const auto parse_elements = [&](std::ifstream& fin) {
    //     size_t num_elements;
    //     fin >> num_elements;

    //     // Tmp storage of elements;
    //     std::vector<int> triangle_element_idx;
    //     std::vector<int> triangle_elements;
    //     std::vector<int> quad_element_idx;
    //     std::vector<int> quad_elements;
    //     std::vector<int> tet_element_idx;
    //     std::vector<int> tet_elements;
    //     std::vector<int> hex_element_idx;
    //     std::vector<int> hex_elements;

    //     auto get_element_storage = [&](int elem_type) -> std::vector<int>* {
    //         switch (elem_type)
    //         {
    //             default:
    //                 assert(false && "Unsupported element type encountered");
    //             case 2:
    //                 return &triangle_elements;
    //             case 3:
    //                 return &quad_elements;
    //             case 4:
    //                 return &tet_elements;
    //             case 5:
    //                 return &hex_elements;
    //         };
    //     };

    //     auto get_element_idx_storage = [&](int elem_type) ->
    //     std::vector<int>* {
    //         switch (elem_type)
    //         {
    //             default:
    //                 assert(false && "Unsupported element type encountered");
    //             case 2:
    //                 return &triangle_element_idx;
    //             case 3:
    //                 return &quad_element_idx;
    //             case 4:
    //                 return &tet_element_idx;
    //             case 5:
    //                 return &hex_element_idx;
    //         };
    //     };

    //     size_t nodes_per_element;
    //     int glob_elem_type = -1;

    //     if (m_binary)
    //     {
    //         eat_white_space(fin);
    //         int elem_read = 0;
    //         while (elem_read < num_elements)
    //         {
    //             // Parse element header.
    //             int elem_type, num_elems, num_tags;
    //             fin.read((char*)&elem_type, sizeof(int));
    //             fin.read((char*)&num_elems, sizeof(int));
    //             fin.read((char*)&num_tags, sizeof(int));
    //             nodes_per_element = num_nodes_per_elem_type(elem_type);
    //             std::vector<int>& elements = *get_element_storage(elem_type);
    //             std::vector<int>& element_idx =
    //                 *get_element_idx_storage(elem_type);

    //             for (size_t i = 0; i < num_elems; i++)
    //             {
    //                 int elem_idx;
    //                 fin.read((char*)&elem_idx, sizeof(int));
    //                 elem_idx -= 1;
    //                 element_idx.push_back(elem_idx);

    //                 // Eat up tags.
    //                 for (size_t j = 0; j < num_tags; j++)
    //                 {
    //                     int tag;
    //                     fin.read((char*)&tag, sizeof(int));
    //                 }

    //                 // Element values.
    //                 for (size_t j = 0; j < nodes_per_element; j++)
    //                 {
    //                     int idx;
    //                     fin.read((char*)&idx, sizeof(int));
    //                     elements.push_back(idx - 1);
    //                 }
    //             }

    //             elem_read += num_elems;
    //         }
    //     }
    //     else
    //     {
    //         for (size_t i = 0; i < num_elements; i++)
    //         {
    //             // Parse per element header
    //             int elem_num, elem_type, num_tags;
    //             fin >> elem_num >> elem_type >> num_tags;
    //             for (size_t j = 0; j < num_tags; j++)
    //             {
    //                 int tag;
    //                 fin >> tag;
    //             }
    //             nodes_per_element = num_nodes_per_elem_type(elem_type);
    //             std::vector<int>& elements = *get_element_storage(elem_type);
    //             std::vector<int>& element_idx =
    //                 *get_element_idx_storage(elem_type);

    //             elem_num -= 1;
    //             element_idx.push_back(elem_num);

    //             // Parse node idx.
    //             for (size_t j = 0; j < nodes_per_element; j++)
    //             {
    //                 int idx;
    //                 fin >> idx;
    //                 elements.push_back(idx - 1);  // msh index starts from 1.
    //             }
    //         }
    //     }

    //     auto copy_to_array = [&](const std::vector<int>& elements,
    //                              const int nodes_per_element) {
    //         const size_t num_elements = elements.size() / nodes_per_element;
    //         if (elements.size() % nodes_per_element != 0)
    //         {
    //             assert(false && "parsing element failed");
    //             return;
    //         }
    //         m_elements.resize(elements.size());
    //         std::copy(elements.begin(), elements.end(), m_elements.data());
    //         m_nodes_per_element = nodes_per_element;
    //     };

    //     if (!tet_elements.empty())
    //     {
    //         copy_to_array(tet_elements, 4);
    //         m_element_type = 4;
    //     }
    //     else if (!hex_elements.empty())
    //     {
    //         copy_to_array(hex_elements, 8);
    //         m_element_type = 5;
    //     }
    //     else if (!triangle_elements.empty())
    //     {
    //         copy_to_array(triangle_elements, 3);
    //         m_element_type = 2;
    //     }
    //     else if (!quad_elements.empty())
    //     {
    //         copy_to_array(quad_elements, 4);
    //         m_element_type = 3;
    //     }
    //     else
    //     {
    //         // 0 elements, use triangle by default.
    //         m_element_type = 2;
    //     }
    // };

    // const auto parse_element_field = [&](std::ifstream& fin) {
    //     size_t num_string_tags;
    //     size_t num_real_tags;
    //     size_t num_int_tags;

    //     fin >> num_string_tags;
    //     std::string* str_tags = new std::string[num_string_tags];
    //     for (size_t i = 0; i < num_string_tags; i++)
    //     {
    //         eat_white_space(fin);
    //         if (fin.peek() == '\"')
    //         {
    //             // Handle field name between quoates.
    //             char buf[128];
    //             fin.get();  // remove the quote at the beginning.
    //             fin.getline(buf, 128, '\"');
    //             str_tags[i] = std::string(buf);
    //         }
    //         else
    //         {
    //             fin >> str_tags[i];
    //         }
    //     }

    //     fin >> num_real_tags;
    //     Float* real_tags = new Float[num_real_tags];
    //     for (size_t i = 0; i < num_real_tags; i++)
    //         fin >> real_tags[i];

    //     fin >> num_int_tags;
    //     int* int_tags = new int[num_int_tags];
    //     for (size_t i = 0; i < num_int_tags; i++)
    //         fin >> int_tags[i];

    //     if (num_string_tags <= 0 || num_int_tags <= 2)
    //     {
    //         delete[] str_tags;
    //         delete[] real_tags;
    //         delete[] int_tags;
    //         assert(false && "invalid format");
    //         return;
    //     }
    //     std::string fieldname = str_tags[0];
    //     int num_components = int_tags[1];
    //     int num_entries = int_tags[2];
    //     VectorF field(num_entries * num_components);

    //     delete[] str_tags;
    //     delete[] real_tags;
    //     delete[] int_tags;

    //     if (m_binary)
    //     {
    //         size_t num_bytes = (num_components * m_data_size + 4) *
    //         num_entries; char* data = new char[num_bytes];
    //         eat_white_space(fin);
    //         fin.read(data, num_bytes);
    //         for (size_t i = 0; i < num_entries; i++)
    //         {
    //             int elem_idx = *reinterpret_cast<int*>(
    //                 &data[i * (4 + num_components * m_data_size)]);
    //             elem_idx -= 1;
    //             size_t base_idx = i * (4 + num_components * m_data_size) + 4;
    //             for (size_t j = 0; j < num_components; j++)
    //             {
    //                 field[elem_idx * num_components + j] =
    //                     *reinterpret_cast<Float*>(
    //                         &data[base_idx + j * m_data_size]);
    //             }
    //         }
    //         delete[] data;
    //     }
    //     else
    //     {
    //         int elem_idx;
    //         for (size_t i = 0; i < num_entries; i++)
    //         {
    //             fin >> elem_idx;
    //             elem_idx -= 1;
    //             for (size_t j = 0; j < num_components; j++)
    //             {
    //                 fin >> field[elem_idx * num_components + j];
    //             }
    //         }
    //     }

    //     m_element_fields[fieldname] = field;
    // };

    // const auto parse_node_field = [&](std::ifstream& fin) {
    //     size_t num_string_tags;
    //     size_t num_real_tags;
    //     size_t num_int_tags;

    //     fin >> num_string_tags;
    //     std::string* str_tags = new std::string[num_string_tags];
    //     for (size_t i = 0; i < num_string_tags; i++)
    //     {
    //         eat_white_space(fin);
    //         if (fin.peek() == '\"')
    //         {
    //             // Handle field name between quoates.
    //             char buf[128];
    //             fin.get();  // remove the quote at the beginning.
    //             fin.getline(buf, 128, '\"');
    //             str_tags[i] = std::string(buf);
    //         }
    //         else
    //         {
    //             fin >> str_tags[i];
    //         }
    //     }

    //     fin >> num_real_tags;
    //     Float* real_tags = new Float[num_real_tags];
    //     for (size_t i = 0; i < num_real_tags; i++)
    //         fin >> real_tags[i];

    //     fin >> num_int_tags;
    //     int* int_tags = new int[num_int_tags];
    //     for (size_t i = 0; i < num_int_tags; i++)
    //         fin >> int_tags[i];

    //     if (num_string_tags <= 0 || num_int_tags <= 2)
    //     {
    //         delete[] str_tags;
    //         delete[] real_tags;
    //         delete[] int_tags;
    //         assert(false && "invalid format");
    //         return;
    //     }
    //     std::string fieldname = str_tags[0];
    //     int num_components = int_tags[1];
    //     int num_entries = int_tags[2];
    //     VectorF field(num_entries * num_components);

    //     delete[] str_tags;
    //     delete[] real_tags;
    //     delete[] int_tags;

    //     if (m_binary)
    //     {
    //         size_t num_bytes = (num_components * m_data_size + 4) *
    //         num_entries; char* data = new char[num_bytes];
    //         eat_white_space(fin);
    //         fin.read(data, num_bytes);
    //         for (size_t i = 0; i < num_entries; i++)
    //         {
    //             int node_idx = *reinterpret_cast<int*>(
    //                 &data[i * (4 + num_components * m_data_size)]);
    //             node_idx -= 1;
    //             size_t base_idx = i * (4 + num_components * m_data_size) + 4;
    //             for (size_t j = 0; j < num_components; j++)
    //             {
    //                 field[node_idx * num_components + j] =
    //                     *reinterpret_cast<Float*>(
    //                         &data[base_idx + j * m_data_size]);
    //             }
    //         }
    //         delete[] data;
    //     }
    //     else
    //     {
    //         int node_idx;
    //         for (size_t i = 0; i < num_entries; i++)
    //         {
    //             fin >> node_idx;
    //             node_idx -= 1;
    //             for (size_t j = 0; j < num_components; j++)
    //             {
    //                 fin >> field[node_idx * num_components + j];
    //             }
    //         }
    //     }

    //     m_node_fields[fieldname] = field;
    // };
    // const auto parse_unknown_field = [](std::ifstream& fin,
    //                                     const std::string& fieldname) {
    //     std::cerr << "Warning: \"" << fieldname
    //               << "\" not supported yet.  Ignored." << std::endl;
    //     std::string endmark = fieldname.substr(0, 1) + "End" +
    //                           fieldname.substr(1, fieldname.size() - 1);

    //     std::string buf("");
    //     while (buf != endmark && !fin.eof())
    //     {
    //         fin >> buf;
    //     }
    // };

    // V.resize(m_nodes.rows() / 3, 3);
    // for (int i = 0; i < m_nodes.rows() / 3; i++)
    // {
    //     for (int j = 0; j < 3; j++)
    //     {
    //         V(i, j) = m_nodes(i * 3 + j);
    //     }
    // }
    // int ss = num_nodes_per_elem_type(m_element_type);
    // T.resize(m_elements.rows() / ss, ss);
    // for (int i = 0; i < m_elements.rows() / ss; i++)
    // {
    //     for (int j = 0; j < ss; j++)
    //     {
    //         T(i, j) = m_elements(i * ss + j);
    //     }
    // }
    // return true;
}

void parseNodes(std::ifstream& fileStream, std::vector<vec3d<double>>& nodes)
{
    size_t numNodes;
    fileStream >> numNodes;
    nodes.resize(numNodes);

    int node_idx;
    for (size_t i = 0; i < numNodes; i++)
    {
        fileStream >> node_idx;
        node_idx -= 1;
        fileStream >> nodes[node_idx](0) >> nodes[node_idx](1) >>
            nodes[node_idx](2);
    }
};