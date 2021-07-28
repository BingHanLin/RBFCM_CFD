# Radial basis function collocation method (RBFCM) for computational fluid dynamics

## Introduction
**RBFCM_CFD** is a simple CFD project that solves some classical fluid prolbems. It is designed to be modified easily by students  size_terested in learning meshless method. More features will be added in the future.
***
## Brief theoretical description

The local radial-basis-function collocation method
(LRBFCM) that based on the multiquadric type radial basis function is used to discretize the spatial
derivitives of the governing equations. As a result, a meshless numerical method is developed for
solving problems. The meshless numerical method is much simpler compared with the traditional numerical methods, such as the finite difference method and the finite element method.

***
## Program Structure

* **controlData** \
Reads the configuration data and offers method to access the data.

* **mesh** \
Reads the mesh file and store the information of nodes.

* **rbfBasis** \
Defines the shape parameter and searching radius when it is initialize. Compute the linear operation of the radial basis function and return it.

* **conditionPool** \
    Collects the **initial conditions** and **boundary conditions** from configuration.
    note: the first index indicate the target node.

    * **boundaryConditions** \
        Defines the boundary conditions, which fill the coefficient matrix and source vector.

    * **initialConditions** \
        Defines the initial conditions, which fill solution vector at begining.

* **simulationDomain** \
Execute the simulation which includes matrix assemble and solving.

***
## Examples

1. **LAPLACE equation** \
LAPLACE equation is solved in a rectangular domain with dirichlet boundary conditions on all sides.
![](asset\laplace.png)

2. **Navier Stokes equation** \
The projection method is used to solve the Navier Stokes equation. A lid-driven square cavity flow case is presented.

![](asset\lidCavity.png) 
***

## Dependencies
The following open source libraries or third party functions are used by this project:
- [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) library for linear algebra.
- [nanoflann](https://github.com/jlblancoc/nanoflann) is a C++ header-only library for building KD-Trees.
***
## References
[1] [Tsai, C.-C., Lin, Z.-H., & Hsu, T.-W. (2015). Using a local radial basis function collocation method to approximate radiation boundary conditions. Ocean Engineering, 105, 231–241.](https://doi.org/10.1016/j.oceaneng.2015.06.030) 

[2] [Chen, W., Zhuo. & Chen, C. (2014). Recent advances on radial basis function collocation methods. Berlin: Springer](https://www.springer.com/gp/book/9783642395710).

[3] [Catch2](http://blog.guorongfei.com/2016/08/22/cpp-unit-test-catch/) 

[4] [Eigen Cheat sheet](https://gist.github.com/gocarlos/c91237b02c120c6319612e42fa196d77) 

[5] [Eigen Macros](https://eigen.tuxfamily.org/dox/TopicPreprocessorDirectives.html) 