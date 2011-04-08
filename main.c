
#include <math.h>
#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE_struct_ls.h"

#ifdef M_PI
  #define PI M_PI
#else
  #define PI 3.14159265358979
#endif

int main (int argc, char *argv[])
{

   int i, j, k;

   int myid, num_procs;

   int n, N, pi, pj;
   double h, h2;
   int ilower[2], iupper[2];

   int solver_id;
   int n_pre, n_post;
   int rap, relax, skip, sym;
   int time_index;

   int num_iterations;
   double final_res_norm;

   int print_solution;

  HYPRE_StructGrid     grid;
  HYPRE_StructStencil  stencil;
  HYPRE_StructMatrix   A;
  HYPRE_StructVector   b;
  HYPRE_StructVector   x;
  HYPRE_StructSolver   solver;
  HYPRE_StructSolver   precond;

/* Initialize MPI */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

/* 1. Set up a grid */
/* Create an empty 2D grid object */
  HYPRE_StructGridCreate(MPI_COMM_WORLD, 2, &grid);

/* Add a new box to the grid */
  HYPRE_StructGridSetExtents(grid, ilower, iupper);

/* This is a collective call finalizing the grid assembly.
The grid is now ``ready to be used'' */
  HYPRE_StructGridAssemble(grid);

/* 2. Define the discretization stencil */
/* Define the geometry of the stencil */
  int offsets[5][2] = {{0,0}, {-1,0}, {1,0}, {0,-1}, {0,1}};

/* Create an empty 2D, 5-pt stencil object */
  HYPRE_StructStencilCreate(2, 5, &stencil);

/* Assign stencil entries */
  for (i = 0; i < 5; i++){
    HYPRE_StructStencilSetElement(stencil, i, offsets[i]);
    }

/* 3. Set up Struct Vectors for b and x */
  double *values;

/* Create an empty vector object */
  HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &b);
  HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &x);

/* Indicate that the vector coefficients are ready to be set */
  HYPRE_StructVectorInitialize(b);
  HYPRE_StructVectorInitialize(x);

  values = calloc((n*n), sizeof(double));


/* Free memory */
  HYPRE_StructGridDestroy(grid);
  HYPRE_StructStencilDestroy(stencil);
  HYPRE_StructMatrixDestroy(A);
  HYPRE_StructVectorDestroy(b);
  HYPRE_StructVectorDestroy(x);

/* Finalize MPI */
  MPI_Finalize();

  return (0);
}
