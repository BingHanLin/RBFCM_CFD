# Radial basis function collocation method (RBFCM) in computational fluid dynamics

## Introduction
**RBFCM-IN-CFD** is a simple CFD code that solves some classical fluid prolbems. It is designed to be modified easily by students  interested in learning meshless method. More features will be added in the future.
***
## Brief theoretical description

***
## Program Structure

* **GMMMQBasis2D.h**
A class define the shape parameter when it is initialize. Compute the linear operation of the radial basis function and return it to **Collocation2D.h**.

* **Collocation2D.h**
A fucntion collocate the nodes near the target node and return the local vector of this node cloud.
note: the first index indicate the target node.

* **GMMRectangular.h**
A class generate a retangular domain with orthogonal node distribution. The nodes on the corner are neglected.

***
## Examples

***
## References
>> [1] D