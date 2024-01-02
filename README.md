The solver presented here consists in a discretization of the equations of the model Uyttendaele and Shambaugh (1990) using the finite volume method (FVM) with Matlab. The first order upwind discretization of the convection term is employed, 
and the source terms are calculated using inputs from single-phase CFD calculation of the airflow.
The boundary conditions are the fiber flow rate, the initial fiber temperature and diameter.
The fiber properties are set by the user and the non-newtonian viscosity is calculated using the Cross-law model.
The results of the solver are the fiber and diameter temperature.