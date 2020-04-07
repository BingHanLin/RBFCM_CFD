#ifndef RECTANGLE_HPP
#define RECTANGLE_HPP

#include <vector>
/*************************************************************************


*************************************************************************/

class Rectangle
{
   public:
    int numInnNodes_, numBouNodes_, numAllNodes_;

    Rectangle(int numNodeX, int numNodeY, double sizeX, double sizeY);

    ~Rectangle(){};

    const std::vector<std::vector<double> >& getNodes() const;
    const std::vector<std::vector<double> >& getNormals() const;

   private:
    int numNodeX_, numNodeY_;
    double sizeX_, sizeY_;
    double centerX_, centerY_;
    std::vector<std::vector<double> > nodes_;
    std::vector<std::vector<double> > normals_;

    std::vector<double> computeAllNodes(int i) const;
    std::vector<double> computeInnerNodes(int i) const;
    std::vector<double> computeBoundaryNodes(int i) const;
    std::vector<double> computeNormals(int i) const;
};
#endif
