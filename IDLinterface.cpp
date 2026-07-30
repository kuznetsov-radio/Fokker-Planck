#include "main.h"

#ifndef LINUX
extern "C" __declspec(dllexport) int CREATE_SOLVER(int argc, void **argv)
#else
extern "C" int CREATE_SOLVER(int argc, void **argv)
#endif
{
 int *general_param=(int*)argv[0];
 double *grid_param=(double*)argv[1];
 double *inf_param=(double*)argv[2];
 double *loop_param=(double*)argv[3];
 double *turb_param=(double*)argv[4];

 create_solver(general_param, grid_param, inf_param, loop_param, turb_param);

 return 0;
}

#ifndef LINUX
extern "C" __declspec(dllexport) int DELETE_SOLVER(int argc, void **argv)
#else
extern "C" int DELETE_SOLVER(int argc, void **argv)
#endif
{
 delete_solver();

 return 0;
}

#ifndef LINUX
extern "C" __declspec(dllexport) int GET_GRIDS(int argc, void **argv)
#else
extern "C" int GET_GRIDS(int argc, void **argv)
#endif
{
 double *E_out=(double*)argv[0];
 double *mu_out=(double*)argv[1];
 double *z_out=(double*)argv[2];
 double *B_loop=(double*)argv[3];
 double *n_loop=(double*)argv[4];

 get_grids(E_out, mu_out, z_out, B_loop, n_loop);

 return 0;
}

#ifndef LINUX
extern "C" __declspec(dllexport) int SOLVE(int argc, void **argv)
#else
extern "C" int SOLVE(int argc, void **argv)
#endif
{
 double *time=(double*)argv[0];
 double *S_in=(double*)argv[1];
 double *J_out=(double*)argv[2];
 double *f_out=(double*)argv[3];
 double *n_fast_out=(double*)argv[4];

 solve(time, S_in, J_out, f_out, n_fast_out);

 return 0;
}