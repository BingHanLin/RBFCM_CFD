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

enum initTypeEnum
{
    UNIFORM
};

NLOHMANN_JSON_SERIALIZE_ENUM(initTypeEnum, {{UNIFORM, "uniform"}})

enum meshTypeEnum
{
    RECTNAGLE,
    DEFAULT
};

NLOHMANN_JSON_SERIALIZE_ENUM(meshTypeEnum, {
                                               {RECTNAGLE, "rectangle"},
                                               {DEFAULT, "mshFile"},
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
    DIVERGENCE,
    PARTIAL_D1,
    PARTIAL_D2,
    PARTIAL_D3,
    NEUMANN
};

#endif