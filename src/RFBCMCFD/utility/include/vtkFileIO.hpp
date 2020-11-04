#ifndef VTKFILEIO_HPP
#define VTKFILEIO_HPP

#include "pugixml.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <fstream>
template <typename T, int N>
void appendArrayToVTKNode(const std::vector<Eigen::Matrix<T, N, 1>> &values,
                          const std::string name, pugi::xml_node &node)
{
    pugi::xml_node dataArray = node.append_child("DataArray");
    std::stringstream positions;
    for (size_t i = 0; i < values.size(); ++i)
    {
        for (size_t j = 0; j < N; ++j)
        {
            positions << values[i][j] << " ";
        }
    }
    dataArray.append_child(pugi::node_pcdata)
        .set_value(positions.str().c_str());
    dataArray.append_attribute("Name") = name.c_str();

    if (std::is_same<T, double>::value)
    {
        dataArray.append_attribute("type") = "Float64";
    }
    else if (std::is_same<T, size_t>::value)
    {
        dataArray.append_attribute("type") = "Int64";
    }

    dataArray.append_attribute("NumberOfComponents") =
        std::to_string(N).c_str();
    dataArray.append_attribute("format") = "ascii";
}

template <typename T>
inline void appendScalarsToVTKNode(
    const Eigen::Matrix<T, Eigen::Dynamic, 1> &values, const std::string name,
    pugi::xml_node &node)
{
    pugi::xml_node dataArray = node.append_child("DataArray");
    std::stringstream valuesStream;
    for (size_t i = 0; i < values.size(); ++i)
    {
        valuesStream << values[i] << " ";
    }
    dataArray.append_child(pugi::node_pcdata)
        .set_value(valuesStream.str().c_str());
    dataArray.append_attribute("Name") = name.c_str();

    if (std::is_same<T, double>::value)
    {
        dataArray.append_attribute("type") = "Float64";
    }
    else if (std::is_same<T, size_t>::value)
    {
        dataArray.append_attribute("type") = "Int64";
    }

    dataArray.append_attribute("format") = "ascii";
}

void addCells(const size_t cellNum, pugi::xml_node &Node)
{
    // set connectivity
    pugi::xml_node connectivityArray = Node.append_child("DataArray");
    std::stringstream valuesStream;
    for (size_t i = 0; i < cellNum; ++i)
    {
        valuesStream << i << " ";
    }
    connectivityArray.append_child(pugi::node_pcdata)
        .set_value(valuesStream.str().c_str());
    connectivityArray.append_attribute("Name") = "connectivity";
    connectivityArray.append_attribute("type") = "Int32";
    connectivityArray.append_attribute("format") = "ascii";

    // set offsets
    pugi::xml_node offsetsArray = Node.append_child("DataArray");
    offsetsArray.append_child(pugi::node_pcdata)
        .set_value(valuesStream.str().c_str());
    offsetsArray.append_attribute("Name") = "offsets";
    offsetsArray.append_attribute("type") = "Int32";
    offsetsArray.append_attribute("format") = "ascii";

    // set types
    pugi::xml_node typesArray = Node.append_child("DataArray");
    std::stringstream valuesStreamTypes;
    for (size_t i = 0; i < cellNum; ++i)
    {
        valuesStreamTypes << 1 << " ";
    }
    typesArray.append_child(pugi::node_pcdata)
        .set_value(valuesStreamTypes.str().c_str());
    typesArray.append_attribute("Name") = "types";
    typesArray.append_attribute("type") = "Int32";
    typesArray.append_attribute("format") = "ascii";
}

void writeVTKGroupFile(const std::string gourpFileName,
                       const std::string childFileNmae, const double timeStep)
{
    std::ifstream vtkGroupFile(gourpFileName);
    pugi::xml_document doc;
    if (!vtkGroupFile)
    {
        pugi::xml_node VTKFile = doc.append_child("VTKFile");
        VTKFile.append_attribute("type") = "Collection";
        VTKFile.append_attribute("version") = "0.1";
        VTKFile.append_attribute("byte_order") = "LittleEndian";
        pugi::xml_node Collection = VTKFile.append_child("Collection");
    }
    else
    {
        pugi::xml_parse_result result = doc.load_file(gourpFileName.c_str());
    }

    pugi::xml_node DataSet =
        doc.child("VTKFile").child("Collection").append_child("DataSet");
    DataSet.append_attribute("timestep") = std::to_string(timeStep).c_str();
    DataSet.append_attribute("part") = 0;
    DataSet.append_attribute("file") = childFileNmae.c_str();

    doc.save_file(gourpFileName.c_str());
}
#endif
