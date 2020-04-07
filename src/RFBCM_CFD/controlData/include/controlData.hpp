#ifndef CONTROLDATA_HPP
#define CONTROLDATA_HPP

#include <iostream>

class controlData
{
   public:
    controlData(const std::string workingDir);

   private:
    std::string workingDir_;
};

#endif
