#ifndef LAPLACE_HH
#define LAPLACE_HH

#include "hdnum.hh"
using namespace hdnum;

/** @brief Laplace Problem

    \tparam T a type representing time

    \tparam N a type representing states and f-values

    \tparam DF a type representing the domain function

    \tparam BF a type representing the boundary function

*/
template<class T, class N, class DF, class BF, int dimension>
class LaplaceCentralDifferences
{
public:
  /** \brief export size_type */
  typedef std::size_t size_type;

  /** \brief export time_type */
  typedef T time_type;

  /** \brief export number_type */
  typedef N number_type;

  /** \brief export domain function type */
  typedef DF DomainFunction;

  /** \brief export boundary function type */
  typedef BF BoundaryFunction;

  enum{ dim = dimension };

private:
  typedef ::SGrid<N,DF,dimension> SGrid;
  SGrid grid;
  const BoundaryFunction & bf;
  
public:

  // make the model 
  LaplaceCentralDifferences(const Vector<number_type> extent_, 
                            const Vector<size_type> size_, 
                            const DomainFunction & df_,
                            const BoundaryFunction & bf_) 
    : grid(extent_,size_,df_), bf(bf_)
  {}

  //! return number of componentes for the model
  std::size_t size () const
  {
    return grid.getNumberOfNodes();
  }

  //! set initial state including time value
  void initialize (T& , Vector<N>& ) const
  {
    // As this is a stationary model, we do not do anything
  }


  //! Returns the grid
  const SGrid & getGrid()
  {
    return grid;
  }

  //! Iterate the grid nodes and evaluate the function assuming U=0!
  //! As this is a linear problem, the choice of U is irelevant, as
  //! the Newton solver will solve in one step.
  void f (const T& t, const Vector<N>& /* U = 0 */, Vector<N>& result) const
  {
    // Reset rhs
    result = 0.;

    const Vector<number_type> h = grid.getCellWidth();

    // Iterate grid nodes
    size_type n_nodes = grid.getNumberOfNodes();
    for(size_type n = 0; n<n_nodes; ++n){

      // Check if node is on grid boundary
       if(grid.isBoundaryNode(n)){
         
         // Get node world coordinates
         const Vector<number_type> x = grid.getCoordinates(n);

         // Determine whether this is a Dirichlet boundary condition
         // (otherwise it is assumed to be a Neumann condition).
         if(bf.isDirichlet(x)){
           result[n] = - bf.getDirichletValue(t,x);
         }
         else{
           for(int d=0; d<dim; ++d){
             for(int s=0; s<2; ++s){
               int side = s*2 -1;
               const int bneighbor = grid.getNeighborIndex(n,d,side);
               if(size_type(bneighbor) == grid.invalid_node){
                 result[n] += number_type(side) * bf.getNeumannValue(t,x)[d] / h[d];
               }
             }// s
           }// d

         }
       }
    }
  }

  //! Iterate the grid nodes and evaluate the function derivatives
  //! assuming U=0!  As this is a linear problem, the choice of U is
  //! irelevant, as the Newton solver will solve in one step.
  void f_x (const T& , const Vector<N>& /* U = 0 */, Matrix<N>& result) const
  {
    // Reset matrix
    result = 0.;

    const Vector<number_type> h = grid.getCellWidth();

    // Iterate grid nodes
    size_type n_nodes = grid.getNumberOfNodes();
    for(size_type n = 0; n<n_nodes; ++n){

      // Check if node is NOT on domain boundary
      if(!grid.isBoundaryNode(n)){
        for(int d=0; d<dim; ++d){
          result[n][n] += -2 / h[d] / h[d];
        
          for(int s=0; s<2; ++s)
            result[n][grid.getNeighborIndex(n,d,s*2-1)] += 1 / h[d] / h[d];
        }
      }
      else{
        // Get node world coordinates
        const Vector<number_type> x = grid.getCoordinates(n);

         // Determine whether this is a Dirichlet boundary condition
         // (otherwise it is assumed to be a Neumann condition).
        if(bf.isDirichlet(x)){
          result[n][n] = 1;

        } 
        else{
          for(int d=0; d<dim; ++d){
            for(int s=0; s<2; ++s){
              int side = s*2 -1;
              const int bneighbor = grid.getNeighborIndex(n,d,side);
              if(size_type(bneighbor) == grid.invalid_node){
                side *= -1;
                const int neighbor = grid.getNeighborIndex(n,d,side);
                result[n][n] += -1. / h[d] / h[d];
                result[n][neighbor] += 1. / h[d] / h[d];
              }
            }
          }
        } 

      } // Boundary

    } // n

  }

};

#endif
