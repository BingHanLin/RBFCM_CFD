#ifndef ENUMMAP_HPP
#define ENUMMAP_HPP

#include "json.h"

enum solverTypeEnum
{
    NAVIERSTOKES,
    TRANSFEREQ
};

NLOHMANN_JSON_SERIALIZE_ENUM(solverTypeEnum, {
                                                 {NAVIERSTOKES, "navierStokes"},
                                                 {TRANSFEREQ, "transferEq"},
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

enum rbfOperatorType
{
    CONSTANT,
    LAPLACE,
    PARTIAL_D1,
    PARTIAL_D2,
    PARTIAL_D3,
    NEUMANN
};

#endif