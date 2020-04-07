#ifndef ENUMMAP_HPP
#define ENUMMAP_HPP

#include "json.h"

enum solverTypeEnum
{
    NAVIERSTOKES,
    POISSON
};

NLOHMANN_JSON_SERIALIZE_ENUM(solverTypeEnum, {
                                                 {NAVIERSTOKES, "navierStokes"},
                                                 {POISSON, "poisson"},
                                             })

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

enum boundaryConditionType
{
    constantValue
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