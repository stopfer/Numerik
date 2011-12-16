#include <iostream>
#include <math.h>

#include "hdnum.hh"
#include "initial_value_problem.h"

template<class Model>
class ImplicitTrapezSolver {
public:
    typedef typename Model::time_type time_type;
    typedef typename Model::number_type number_type;
    //! constructor stores reference to the model

    ImplicitTrapezSolver(const Model & model_)
    : model(model_), u(model.size()), f(model.size()) ,x(model.size()), s(model.size()),w(model.size()),p(model.size()), q(model.size()){
        model.initialize(t, u);
        dt = 0.1;
    }

    //! set time step for subsequent steps

    void set_dt(time_type dt_) {
        dt = dt_;
    }

    //! do one step

    void step() {

        hdnum::Vector<number_type> iterValue = u;
        number_type tol = 0.001;
        number_type error = 100000;
        //newton solve
        while (error > tol) {
            hdnum::DenseMatrix<number_type> A;
            hdnum::Vector<number_type> b;
            this->evaluateNewtonJacobiFunction(iterValue, A);
            this->evaluateNewtonFunction(iterValue, b);
            for(size_t bb=0;bb<b.size();++bb)
                b[bb]*=(number_type(-1.0));
            //do the Lr decompositon
            hdnum::row_equilibrate(A, s);
            hdnum::lr_fullpivot(A, p, q);
            hdnum::apply_equilibrate(s, b);
            hdnum::permute_forward(p, b);
            hdnum::solveL(A, b, b);
            x = 0;
            hdnum::solveR(A, x, b);
            hdnum::permute_backward(q, x);
            z = x;

            //update the value         
            newValue = iterValue + z;

            //compute error
            tmp = newValue - iterValue;
            error = tmp.two_norm() / iterValue.two_norm();
            iterValue = newValue;
        }
        u = iterValue;
        t += dt;
    }

    void evaluateNewtonFunction
    (
        const hdnum::Vector<number_type>& inputVector,
        hdnum::Vector<number_type>& resultVectors
    ) {
        // "sub" function which is passed to the
        // newton solver
        
        hdnum::Vector<number_type> t1;
        hdnum::Vector<number_type> t2;
        model.f(t + dt, inputVector, t1);
        model.f(t, inputVector, t2);
        for(size_t tt=0;tt<t1.size();++tt){
            t1*=(dt / 2.0);
            t2*=(dt / 2.0);
        }
        resultVectors = inputVector-t1-u-t2;
    }

    void evaluateNewtonJacobiFunction
    (
            const hdnum::Vector<number_type>& inputVector,
            hdnum::DenseMatrix<number_type>& resultMatrix
    ) {
        // jacobi matrix for the "sub" model
        // which is passed to the newton solver
        
        hdnum::DenseMatrix<number_type> temp;
        hdnum::DenseMatrix<number_type> identity(inputVector.size(), inputVector.size(), 0);
        // identity
        for (size_t i = 0; i < inputVector.size(); ++i)
            identity(i, i) = 1;
        // evaluate the jacobi matrix of the model
        model.f_x(t + dt, u, temp);
        for(size_t c=0;c<temp.colsize();++c)
        for(size_t r=0;r<temp.rowsize();++r)
            temp(r,c)*=(dt / 2.0) ;
        resultMatrix = identity - temp;
    }
    //! set current state
    void set_state(time_type t_, const hdnum::Vector<number_type>& u_) {
        t = t_;
        u = u_;
    }
    //! get current state
    const hdnum::Vector<number_type>& get_state() const {
        return u;
    }
    //! get current time
    time_type get_time() const {
        return t;
    }
    //! get dt used in last step (i.e. to compute current state)
    time_type get_dt() const {
        return dt;
    }
private:
    const Model & model;
    time_type t, dt;
    hdnum::Vector<number_type> u;
    hdnum::Vector<number_type> f;
    hdnum::Vector<number_type> tmp, s, w, x, z, newValue;
    hdnum::Array<size_t> p,q;
    number_type maxerroror_;
};

int main() {
    typedef double Number; // define a number type
    const Number t0 = 0.0; // initial time
    Number tStep = 0.005; // delta t
    const Number tMax = 100000; // end time
    hdnum::Vector<Number> u0(2); // initial state
    //set values of u0 and A
    u0[0] = 3;
    u0[1] = 4;
    typedef InitialValueProblem<Number> Model; // Model type
    Model model(u0, t0); // instantiate model
    typedef ImplicitTrapezSolver<Model> Solver; // solver
    Solver solver(model); // instantiate solver
    // initialize solver
    solver.set_dt(tStep);
    solver.set_state(t0, u0);
    size_t iter=0;
    while (solver.get_time() <= tMax) // the time loop
    {
        solver.step(); // advance model by one time step
        // cout for convergence information
        if(iter%100000==0)
            std::cout<<"u0 for t="<<solver.get_time()<<" ="<<solver.get_state()[0]<<" u1  ="<<solver.get_state()[1]<<"\n";
        ++iter;
    }
    return 0;
 }
