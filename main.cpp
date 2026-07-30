#include "main.h"
#include "MyTypes.h"
#include "solver.h"
#ifndef LINUX
#include <iostream>
#endif

Solver* p_sol=NULL;

//sign of current
//add curl B

#ifndef LINUX
BOOL APIENTRY DllMain(HINSTANCE hinstDLL, DWORD fdwReason, LPVOID lpvReserved)
{
    switch (fdwReason){
        case DLL_PROCESS_ATTACH:{
            break;
        }
        case DLL_PROCESS_DETACH:{
            if (p_sol!=NULL){delete p_sol;}
            break;
        }
    }
    return TRUE; // succesful
}
#endif

#ifndef LINUX
extern "C" __declspec(dllexport) void create_solver(int const *general_param, double const *grid_param, double const *inf_param,double const *loop_param, double const * turb_param)
#else
void create_solver(int const *general_param, double const *grid_param, double const *inf_param,double const *loop_param, double const * turb_param)
#endif
{
    ///general_param = [NE, Nmu, Nz_output, Nz_loop_param, switch_reverse_current, switch_turb]
    ///grid_param=[Emin,zmax,tmax,dt0,dt_max,gamma, r_gridE, r_gridz]
    ///inf_param=[n_inf, log_inf, conductivity_inf]
    ///loop_param=[*z_loop,*B_loop,*n_loop,*coulomb_log,*conductivity]
    ///turb_param=[*E_turb, *lamb_turb]
    p_sol =new Solver(general_param, grid_param, inf_param, loop_param, turb_param);
}

#ifndef LINUX
extern "C" __declspec(dllexport) void delete_solver()
#else
void delete_solver()
#endif
{
    if (p_sol!=NULL)
    {
        delete p_sol;
        p_sol=NULL;
    }
#ifndef LINUX
    else MessageBoxA(0, "First you must to create solver!", "DLL Error in function \"get_grids\"", MB_OK | MB_ICONINFORMATION);
#endif

}

#ifndef LINUX
extern "C" __declspec(dllexport) void get_grids(double *E_out, double *mu_out,  double *z_out, double *dlnB_out,double* n_loop_out)
#else
void get_grids(double *E_out, double *mu_out,  double *z_out, double *dlnB_out,double* n_loop_out)
#endif
{
    if (p_sol!=NULL)
    {
        p_sol->get_E(E_out);
        p_sol->get_mu(mu_out);
        p_sol->get_z(z_out);
        p_sol->get_dlnB(dlnB_out);
        p_sol->get_n_background(n_loop_out);
    }
#ifndef LINUX
    else MessageBoxA(0, "First you must to create solver!", "DLL Error in function \"get_grids\"", MB_OK | MB_ICONINFORMATION);
#endif
}

#ifndef LINUX
extern "C" __declspec(dllexport) void solve(double* time, double * S_in, double *J_in, double *f_in, double *n_fast_in)
#else
void solve(double* time, double * S_in, double *J_in, double *f_in, double *n_fast_in)
#endif
{
    if (p_sol!=NULL)
    {
        Solver & sol=*p_sol;
        sol.define_initial_f(S_in, f_in,J_in,n_fast_in);
        double time0=sol.get_time();
        while(1) {
            if (sol.switch_reverse_current){
                for (int cycle=0; cycle<2; cycle+=1){
                    sol.f_step();
                    sol.J_step();
                }
            }else{
                sol.f_step();
            }
            if (sol.change_step_time_or_break()){break;}
        }
        sol.calc_n_fast();
        double dt=sol.get_time()-time0;
        *time+=dt;
        sol.clear_step();
    }
#ifndef LINUX
    else MessageBoxA(0, "First you must to create solver!", "DLL Error in function \"solve\"", MB_OK | MB_ICONINFORMATION);
#endif
}

