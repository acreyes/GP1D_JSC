
! slope limiters
#define MINMOD  1
#define MC      2
#define VANLEER 3

! Riemann solvers
#define HLL     1
#define HLLC    2
#define ROE     3

! primitive vars
#define DENS_VAR 1
#define VELX_VAR 2
#define VELY_VAR 3
#define VELZ_VAR 4
#define MAGY_VAR 5
#define MAGZ_VAR 6
#define PRES_VAR 7
#define EINT_VAR 8
#define GAMC_VAR 9
#define GAME_VAR 10
#define SHOK_VAR 11
#define NUMB_VAR 11
#define NSYS_VAR 7 /* total number of equations of the conservative system */

! conservative vars
#define MOMX_VAR 2
#define MOMY_VAR 3
#define MOMZ_VAR 4
#define ENER_VAR 7

! waves
#define SHOCKLEFT 1
#define ALFVNLEFT 2
#define SLOWWLEFT 3
#define CTENTROPY 4
#define SLOWWRGHT 5
#define ALFVNRGHT 6
#define SHOCKRGHT 7
#define NUMB_WAVE 7

! setup parameters
#define MAX_STRING_LENGTH 80

! BC
#define OUTFLOW  1
#define PERIODIC 2
#define REFLECT  3
#define USER     4

!Kernel Quadratures
#define EXACT 1
#define MIDPT 2
#define TRAPZ 3
#define SIMPS 4

#define PI 4.*ATAN(1.)
