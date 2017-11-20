/*
 * Benchmakr Cholesky factorization and inverse using different libraries.
 * - Lapack (serial)
 * - ScaLapack (parallel)
 * - Elemental (parallel)
 * - Elpa (parallel)
 */

#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/multithread_info.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>

#include <deal.II/lac/vector.h>

#include <deal.II/lac/scalapack.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_solver.h>

#ifdef WITH_ELEMENTAL
#include <El.hpp>
#endif

#include <fstream>
#include <iostream>


using namespace dealii;

template <typename T, typename T2>
double SPD(const T i, const T j, const T2 size, boost::random::mt19937 &gen, boost::random::uniform_01<> &dist)
{
  const double val = dist(gen);
  Assert (val >= 0. && val <= 1.,
          ExcInternalError());
  if (i==j)
    // since A(i,j) < 1 and
    // a symmetric diagonally dominant matrix is SPD
    return val + size;
  else
    return val;
}

template <typename NumberType>
void test(const unsigned int size, const unsigned int block_size)
{
  MPI_Comm mpi_communicator(MPI_COMM_WORLD);
  const unsigned int this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator));

  ConditionalOStream pcout (std::cout, (this_mpi_process ==0));
  TimerOutput  computing_timer(mpi_communicator,
                               pcout,
                               TimerOutput::never,
                               TimerOutput::wall_times);

  // test multiplication with random vectors
  boost::random::mt19937 gen;
  boost::random::uniform_01<> dist;

  {
    const unsigned int n_mpi_processes(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD));
    TimerOutput::Scope t(computing_timer, "Lapack");

    // do Lapack only for 1 core and single block size (which is irrelevant)
    // otherwise we simply won't have enough memory
    if (n_mpi_processes == 1 &&
        block_size == 16)
      {
        FullMatrix<NumberType> A(size), invA(size);
        for (unsigned int i = 0; i < size; ++i)
          for (unsigned int j = i; j < size; ++j)
            {
              const double val = SPD(i,j,size,gen,dist);
              A(i,j) = val;
              A(j,i) = val;
            }

        invA.invert(A);
      }
  }

  int sc_grid_rows = 0, sc_grid_cols = 0;
  {
    TimerOutput::Scope t(computing_timer, "ScaLapack");
    std::pair<int,int> sizes = std::make_pair(size,size), block_sizes = std::make_pair(block_size,block_size);
    std::shared_ptr<ProcessGrid> grid = std::make_shared<ProcessGrid>(mpi_communicator,sizes,block_sizes);
    ScaLAPACKMatrix<NumberType> scalapack_A (size, grid, block_size);

    sc_grid_rows = grid->get_process_grid_rows();
    sc_grid_cols = grid->get_process_grid_columns();

    const unsigned int n_local_rows = scalapack_A.local_m();
    const unsigned int n_local_cols = scalapack_A.local_n();

    for (unsigned int i = 0; i < n_local_rows; ++i)
      {
        const int glob_i = scalapack_A.global_row(i);
        for (unsigned int j = 0; j < n_local_cols; ++j)
          {
            const int glob_j = scalapack_A.global_column(j);
            const double val = SPD(glob_i,glob_j,size,gen,dist);
            if (glob_i == glob_j)
              scalapack_A.local_el(i,j) = val;
            // FIXME: hard-code lower triangular:
            else if (glob_i < glob_j)
              scalapack_A.local_el(i,j) = val;
            else
              scalapack_A.local_el(i,j) = 0.;
          }
      }

    scalapack_A.compute_cholesky_factorization();
    scalapack_A.invert();
  }

  int el_grid_width = 0, el_grid_height = 0;
  {
    TimerOutput::Scope t(computing_timer, "Elemental");
#ifdef WITH_ELEMENTAL
    // El::PrintVersion(std::ostream &os = std::cout)

    // Change the algorithmic blocksize
    El::SetBlocksize(block_size);

    // Create default column-major process grid:
    El::Grid grid(mpi_communicator);
    el_grid_width  = grid.Width();
    el_grid_height = grid.Height();

    // Build a size x size distributed matrix over
    // for various combinations MC/MR, see http://libelemental.org/documentation/dev/core/dist_matrix/Element/MC_MR.html
    El::DistMatrix<double,El::MC,El::MR> el_A(size,size,grid);
    El::UpperOrLower uplo = El::CharToUpperOrLower('L');

    El::HermitianUniformSpectrum(el_A, size, 1e-9, 10 );

    // custom algorithm for SPD inverse, namely
    // Bientinesi et al "Families of algorithms related to the inversion of a symmetric positive definite matrix", in particular, the variant 2 algorithm from Fig. 9.
    El::HPDInverse(uplo, el_A);

    // const double cond = El::Condition(el_A); // 2-norm
    // const double cond = El::OneCondition(el_A);

    /* can modify non-local entries:
    if( A.DistRank() == 0 )
    {
        A.Reserve( 3 );
        A.QueueUpdate( 50, 50, 1. );
        A.QueueUpdate( 51, 51, 2. );
        A.QueueUpdate( 52, 52, 3. );
    }
    else if( A.DistRank() == 1 )
    {
        A.Reserve( 2 );
        A.QueueUpdate( 0, 0, 17. );
        A.QueueUpdate( 1, 0, 18. );
    }
    A.ProcessQueues();
     */

    /*
    // graphs:
    const int n = 20;
    El::Graph graph(n); // graph with n vertices
    graph.Reserve( 2*n );
    for ( El::Int i=0; i<n; ++i )
      {
        if ( i != 0 )
          graph.QueueConnection( i, i-1 );
        if ( i != n-1 )
          graph.QueueConnection( i, i+1 );
      }
    graph.ProcessQueues();
    */

#endif
  }

  // print wall clock times and some data
  const std::map<std::string,double> summary = computing_timer.get_summary_data(TimerOutput::total_wall_time);
  pcout << size << " "
        << block_size << " "
        << sc_grid_rows << " "
        << sc_grid_cols << " "
        << summary.at("Lapack") << " "
        << summary.at("ScaLapack") << " "
        << summary.at("Elemental") << " "
        << el_grid_width << " "
        << el_grid_height
        << std::endl;


}

/*
extern "C" {
  void openblas_set_num_threads(int num_threads);
}
*/

int main (int argc,char **argv)
{

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      // openblas_set_num_threads(1);
#ifdef WITH_ELEMENTAL
      El::Environment env( argc, argv );
#endif

      const unsigned int n_mpi_processes(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD));
      const unsigned int this_mpi_process(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD));
      const int n_threads = dealii::MultithreadInfo::n_threads();
      if (this_mpi_process==0)
        {
          std::cout << "-------------------------------------------------------------------------" << std::endl
#ifdef DEBUG
                    << "--     . running in DEBUG mode" << std::endl
#else
                    << "--     . running in OPTIMIZED mode" << std::endl
#endif
                    << "--     . running with " << n_mpi_processes << " MPI process" << (n_mpi_processes == 1 ? "" : "es") << std::endl;

          if (n_threads>1)
            std::cout << "--     . using " << n_threads << " threads " << (n_mpi_processes == 1 ? "" : "each") << std::endl;
          std::cout << "-------------------------------------------------------------------------" << std::endl;
        }

      const std::vector<unsigned int> sizes = {{1000, 5000, 10000, 20000, 30000}};
      const std::vector<unsigned int> blocks = {{16,32,64}};

      for (const auto &s : sizes)
        for (const auto &b : blocks)
          test<double>(s,b);

    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    };
}
