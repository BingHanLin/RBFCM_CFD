#ifndef ENUMMAP_HPP
#define ENUMMAP_HPP

#include "json.h"

enum solverType
{
    NAVIERSTOKES,
    TRANSFEREQ
};

NLOHMANN_JSON_SERIALIZE_ENUM(solverType, {
                                             {NAVIERSTOKES, "navierStokes"},
                                             {TRANSFEREQ, "transferEq"},
                                         })

enum systemSateType
{
    TRANSIENT,
    STEADY
};

NLOHMANN_JSON_SERIALIZE_ENUM(systemSateType, {
                                                 {TRANSIENT, "transient"},
                                                 {STEADY, "steady"},
                                             })

enum initConditionType
{
    UNIFORM
};

NLOHMANN_JSON_SERIALIZE_ENUM(initConditionType, {{UNIFORM, "uniform"}})

enum meshType
{
    RECTNAGLE,
    DEFAULT
};

NLOHMANN_JSON_SERIALIZE_ENUM(meshType, {
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