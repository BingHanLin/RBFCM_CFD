#ifndef RECTANGLE_HPP
#define RECTANGLE_HPP

#include <vector>
/*************************************************************************


*************************************************************************/

class Rectangle {
private:
    int    numNodeX_, numNodeY_;
    double sizeX_, sizeY_;
    double centerX_, centerY_;

public:
    int numInn_, numBou_, numAll_;

    Rectangle( int numNodeX, int numNodeY, double sizeX, double sizeY );

    ~Rectangle(){};

    std::vector< double > getInnerNode( int i ) const;
    std::vector< double > getBoundaryNode( int i ) const;
    std::vector< double > getAllNode( int i ) const;
    std::vector< double > getOutNormal( int i ) const;
};
#endif
