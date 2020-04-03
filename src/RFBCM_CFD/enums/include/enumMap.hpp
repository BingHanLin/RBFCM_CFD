#ifndef ENUMMAP_HPP
#define ENUMMAP_HPP

#include "json.h"

enum meshTypeEnum
{
    RECTNAGLE,
    DEFAULT
};

NLOHMANN_JSON_SERIALIZE_ENUM(meshTypeEnum, {
                                               {RECTNAGLE, "rectangle"},
                                               {DEFAULT, "mshFile"},
                                           })

enum dimensionTypeEnum
{
    threeD,
    twoD
};

NLOHMANN_JSON_SERIALIZE_ENUM(dimensionTypeEnum, {
                                                    {threeD, "threeD"},
                                                    {twoD, "twoD"},
                                                })

enum boundaryCondition
{
    NONESLIP
};

enum elementType
{
    LINE,
    TRIANGLE,
    QUADRANGLE,
    TETRAHEDRON,
    HEXAHEDRON,
    NONE
};

#endif