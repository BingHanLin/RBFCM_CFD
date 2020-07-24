#ifndef RECTANGLE_HPP
#define RECTANGLE_HPP

#include <vector>
/*************************************************************************


*************************************************************************/

class Rectangle
{
   public:
    size_t numInnNodes_, numBouNodes_, numAllNodes_;

    Rectangle(size_t numNodeX, size_t numNodeY, double sizeX, double sizeY);

    ~Rectangle(){};

    const std::vector<std::vector<double> >& getNodes() const;
    const std::vector<std::vector<double> >& getNormals() const;

   private:
    size_t numNodeX_, numNodeY_;
    double sizeX_, sizeY_;
    double centerX_, centerY_;
    std::vector<std::vector<double> > nodes_;
    std::vector<std::vector<double> > normals_;

    std::vector<double> computeAllNodes(size_t i) const;
    std::vector<double> computeInnerNodes(size_t i) const;
    std::vector<double> computeBoundaryNodes(size_t i) const;
    std::vector<double> computeNormals(size_t i) const;
};
#endif
