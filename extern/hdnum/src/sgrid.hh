#ifndef HDNUM_SGRID_HH
#define HDNUM_SGRID_HH
#include <limits>
#include <assert.h>

namespace hdnum {
  /** @brief Structured Grid for Finite Differences

      \tparam N  A continuous type representing coordinate values.
      \tparam DF A boolean function which defines the domain.
      \tparam dimension The grid dimension.
  */
  template<class N, class DF, int dimension>
  class SGrid
  {
  public:

    //! Export size type.
    typedef std::size_t size_type;

    //! Export number type.
    typedef N number_type;

    //! Type of the function defining the domain.
    typedef DF DomainFunction;

    enum { dim = dimension };

    //! Side definitions for usage in getNeighborIndex(..)
    static const int positive = 1;
    static const int negative = -1;

    
  private:
  
    const Vector<number_type> extent;
    const Vector<size_type>   size;
    const DomainFunction & df;
    Vector<number_type> h;
    Vector<size_type> offsets;
    std::vector<size_type> node_map;
    std::vector<size_type> grid_map;
    std::vector<bool> inside_map;
    std::vector<bool> boundary_map;

    size_t n_nodes;
    
    inline Vector<size_type> index2grid(size_type index) const 
    {
      Vector<size_type> c(dim);
      for(int d=dim-1; d>=0; --d){
        c[d] = index / offsets[d];
        index -= c[d] * offsets[d];
      }
      return c;
    }

    inline Vector<number_type> grid2world(const Vector<size_type> & c) const
    {
      Vector<number_type> w(dim);
      for(int d=dim-1; d>=0; --d)
        w[d] = c[d] * h[d];
      return w;
    }

    inline Vector<number_type> index2world(size_type index) const
    {
      Vector<number_type> w(dim);
      Vector<size_type> c = index2grid(index);
      return grid2world(c);
    }


  public:

    //! The value which is returned to indicate an invalid node.
    const size_type invalid_node;
    
    /** \brief Constructor
     
        \param[in] extent_ The extent of the grid domain. The actual
        computational domain may be smaller and is defined by the domain
        function df_.

        \param[in] size_ The number of nodes in each grid dimension.

        \param[in] df_ The domain function. It has to provide a boolean
        function evaluate(Vector<number_type> x) which returns true if the
        node which is positioned at the coordinates of x is within the
        computational domain.

    */
    SGrid(const Vector<number_type> extent_, 
          const Vector<size_type> size_, 
          const DomainFunction & df_)
      : extent(extent_), size(size_), df(df_), 
        h(dim), offsets(dim),
        invalid_node(std::numeric_limits<size_type>::max())
    {
      // Determine total number of nodes, increment offsets, and cell
      // widths.
      n_nodes = 1;
      offsets.resize(dim);
      h.resize(dim);
      for(int d=0; d<dim; ++d){
        n_nodes *= size[d];
        offsets[d] = d==0 ? 1 : size[d-1] * offsets[d-1];
        h[d] = extent[d] / number_type(size[d]-1);
      }

      // Initialize maps.
      node_map.resize(0);
      inside_map.resize(n_nodes);
      grid_map.resize(n_nodes);
      boundary_map.resize(0);
      boundary_map.resize(n_nodes,false);

      for(size_type n=0; n<n_nodes; ++n){
        Vector<size_type> c = index2grid(n);
        Vector<number_type> x = grid2world(c);
      
        inside_map[n] = df.evaluate(x);
        if(inside_map[n]){
          node_map.push_back(n);
          grid_map[n] = node_map.size()-1;
        }
        else
          grid_map[n] = invalid_node;
      }

      // Find boundary nodes
      for(size_type n=0; n<node_map.size(); ++n){
        for(int d=0; d<dim; ++d){
          for(int s=0; s<2; ++s){
            const int side = s*2-1;
            const size_type neighbor = getNeighborIndex(n,d,side,1);
            if(neighbor == invalid_node)
              boundary_map[node_map[n]] = true;
          }
        }
      }

    }

    /** \brief Provides the index of the k-th neighbor of the node with index ln.
     
        \param[in] ln Index of the node whose neighbor is to be
        determined.  

        \param[in] n_dim The axes which connects the node and its neighbor
        (e.g. n_dim = 0 for a neighbor in the direction of the x-axes

        \param[in] n_side Determines whether the neighbor is in positive
        of negative direction of the given axes. Should be either
        SGrid::positive or SGrid::negative .

        \param[in] k For k=1 it will return the direct neighbor. Higher
        values will give distant nodes in the given direction. If the
        indicated node is not within the grid any more, then invalid_node
        will be returned. For k=0 it will simply return ln.

        \return size_type The index of the neighbor node.
    */
    size_type getNeighborIndex(const size_type ln, const size_type n_dim, const int n_side, const int k = 1) const
    {
      const size_type n = node_map[ln];
      const Vector<size_type> c = index2grid(n);
      size_type neighbors[2];
      neighbors[0] = c[n_dim];
      neighbors[1] = size[n_dim]-c[n_dim]-1;

      assert(n_side == 1 || n_side == -1);
      if(size_type(k) > neighbors[(n_side+1)/2])
        return invalid_node;

      const size_type neighbor = n + offsets[n_dim] * n_side * k;

      if(!inside_map[neighbor])
        return invalid_node;

      return grid_map[neighbor];
    }

    /** \brief Returns true if the node is on the boundary of the
        discrete compuational domain.
    */
    bool isBoundaryNode(const size_type ln) const
    {
      return boundary_map[node_map[ln]];
    }

    /** \brief Returns the number of nodes which are in the
        compuational domain.
    */
    size_type getNumberOfNodes() const
    {
      return node_map.size();
    }

    Vector<size_type> getGridSize() const
    {
      return size;
    }

    /** \brief Returns the cell width h of the structured grid.
    */
    Vector<number_type> getCellWidth() const
    {
      return h;
    }

    /** \brief Returns the world coordinates of the node with the
        given node index.
    */
    Vector<number_type> getCoordinates(const size_type ln) const
    {
      return index2world(node_map[ln]);
    }

    std::vector<Vector<number_type> > getNodeCoordinates() const
    {
      std::vector<Vector<number_type> > coords;
      for(size_type n=0; n<node_map.size(); ++n){
        coords.push_back(Vector<number_type>(dim));
        coords.back() = index2world(node_map[n]);
      }
      return coords;
    }
  
  };

};

#endif // HDNUM_SGRID_HH
