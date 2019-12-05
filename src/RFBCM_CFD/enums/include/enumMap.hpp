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

enum boundaryCondition
{
    Inner,
    NoNSLIP
};

#endif