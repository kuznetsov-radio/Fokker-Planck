#pragma once

#ifndef __MAIN_H__
#define __MAIN_H__

#include "solver.h"

#ifndef LINUX
#include <windows.h>

#ifdef BUILD_DLL
    #define DLL_EXPORT __declspec(dllexport)
#else
    #define DLL_EXPORT __declspec(dllimport)
#endif

extern "C" __declspec(dllexport) void get_grids(double *E_out, double *mu_out,  double *z_out, double *B_loop,double* n_loop);
extern "C" __declspec(dllexport) void solve(double* time, double* S_in, double *J_out, double *f_out, double *n_fast_out);
extern "C" __declspec(dllexport) void create_solver(int const *general_param, double const *grid_param, double const *inf_param, double const *loop_param, double const * turb_param);
extern "C" __declspec(dllexport) void delete_solver();
#else
void get_grids(double *E_out, double *mu_out,  double *z_out, double *B_loop,double* n_loop);
void solve(double* time, double* S_in, double *J_out, double *f_out, double *n_fast_out);
void create_solver(int const *general_param, double const *grid_param, double const *inf_param, double const *loop_param, double const * turb_param);
void delete_solver();
#endif
#endif // __MAIN_H__
