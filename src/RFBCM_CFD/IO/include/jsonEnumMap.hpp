#ifndef JSONENUMMAP_HPP
#define JSONENUMMAP_HPP

#include "json.h"

enum meshTypeEnum { RECTNAGLE, DEAFULT };

NLOHMANN_JSON_SERIALIZE_ENUM( meshTypeEnum, {
                                                { RECTNAGLE, "rectangle" },
                                                { DEAFULT, "inpFile" },
                                            } )

#endif