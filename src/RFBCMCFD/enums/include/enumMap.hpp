#ifndef ENUMMAP_HPP
#define ENUMMAP_HPP

#include "json.h"

enum solverType
{
    INCOMPRESSIBLE,
    SCALARTRANSPORT
};

NLOHMANN_JSON_SERIALIZE_ENUM(solverType,
                             {
                                 {INCOMPRESSIBLE, "incompressible"},
                                 {SCALARTRANSPORT, "scalarTransport"},
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

enum meshType
{
    RECTNAGLE,
    DEFAULT
};

NLOHMANN_JSON_SERIALIZE_ENUM(meshType, {
                                           {RECTNAGLE, "rectangle"},
                                           {DEFAULT, "mshFile"},
                                       })

enum class boundaryConditionType
{
    constantValue
};

enum class initConditionType
{
    UNIFORM
};

enum class elementType
{
    LINE,
    TRIANGLE,
    QUADRANGLE,
    TETRAHEDRON,
    HEXAHEDRON,
    NONE
};

enum class rbfOperatorType
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