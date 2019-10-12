#ifndef RECTANGLE_HPP
#define RECTANGLE_HPP

#include <vector>
/*************************************************************************


*************************************************************************/

class Rectangle
{
   private:
    int numNodeX_, numNodeY_;
    double sizeX_, sizeY_;
    double centerX_, centerY_;

    std::vector<double> getInnerNode(int i) const;
    std::vector<double> getBoundaryNode(int i) const;
    std::vector<double> getAllNode(int i) const;
    std::vector<double> getOutNormal(int i) const;

   public:
    int numInnNodes_, numBouNodes_, numAllNodes_;
    std::vector<std::vector<double> > nodes_;
    std::vector<std::vector<double> > normals_;

    Rectangle(int numNodeX, int numNodeY, double sizeX, double sizeY);

    ~Rectangle(){};
};
#endif
