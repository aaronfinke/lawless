#include <scitbx/math/eigensystem.h>
#include <scitbx/array_family/shared.h>
#include <tnt_array1d.h>
#include <tnt_array2d.h>
#include <jama_eig.h>

namespace phaser {

// solve symmetric eigenvalue problem (return eigenvalues and eigenvectors)
template <class T>
class Upper_symmetric_eigenvector_solve
{
  public:
  Upper_symmetric_eigenvector_solve(TNT::Fortran_Matrix<T>& A,TNT::Vector<T>& eigenvals)
  {
  //input matrix A replaced with eigenvectors
    TNT::Subscript N = A.num_rows();
    scitbx::af::versa<T, scitbx::af::c_grid<2> >  Aversa(scitbx::af::c_grid<2>(N,N));
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) Aversa(i,j) = A(i+1,j+1);
    scitbx::math::eigensystem::real_symmetric<T> eigens(Aversa.ref());
    A = TNT::Fortran_Matrix<T>(N,N,eigens.vectors().begin());
    eigenvals = TNT::Vector<T>(N,eigens.values().begin());
  }
};

template <class T>
class SymmetricPseudoinverse
{
  private:
  TNT::Fortran_Matrix<T> inverse_matrix;
  public:
  SymmetricPseudoinverse(TNT::Fortran_Matrix<T>& M, int& filtered, bool return_matrix=false, int min_to_filter=0)
  {
    TNT::Subscript N(M.num_rows());
    TNT::Vector<T> eigenvals(N);
//#define __PHASER_DEBUG_PSEUDOINVERSE__
#ifdef __PHASER_DEBUG_PSEUDOINVERSE__
    TNT::Fortran_Matrix<T> Mcopy = M;
#endif
    Upper_symmetric_eigenvector_solve<T>(M, eigenvals); //M now eigenvectors in place
  // Set condition number for range of eigenvalues, depending on variable type
    T evCondition(std::min(1.e9,0.01/std::numeric_limits<T>::epsilon()));
    T evMax(eigenvals(1)); //scitbx, largest first
    T evMin = evMax/evCondition;
    T ZERO(0);
    TNT::Fortran_Matrix<T> lambdaInv(N,N,ZERO);
    for (int i = 1; i <= N; i++)
      lambdaInv(i,i) = (eigenvals(i) > evMin) ? 1/eigenvals(i) : ZERO;
    if (min_to_filter)
    {
      PHASER_ASSERT (min_to_filter < N);
      for (int i=N-min_to_filter+1; i <= N; i++)
        lambdaInv(i,i) = ZERO;
    }
    inverse_matrix = M * lambdaInv * TNT::transpose(M);
    //count number filtered
    filtered = 0;
    for (int i = 1; i <= N; i++) 
      if (eigenvals(i) <= evMin || i > N-min_to_filter) filtered++;
    if (return_matrix)
    {
      TNT::Fortran_Matrix<T> lambda(N,N,ZERO);
      for (int i = 1; i <= N; i++)
        lambda(i,i) = (eigenvals(i) > evMin) ? eigenvals(i) : ZERO;
      M = M * lambda * TNT::transpose(M);
    }
#ifdef __PHASER_DEBUG_PSEUDOINVERSE__
    // Check that pseudoinverse obeys Moore-Penrose conditions
    // (Golub & van Loan, p. 257)
    TNT::Fortran_Matrix<T> testmat = Mcopy * inverse_matrix * Mcopy;
    T maxdev(0),maxelement(0);
    for (int i = 1; i <= N; i++)
      for (int j = 1; j <= N; j++)
      {
        maxelement = std::max(maxelement,std::abs(Mcopy(i,j)));
        maxdev = std::max(maxdev,std::abs(Mcopy(i,j)-testmat(i,j)));
      }
    std::cout << "Moore-Penrose condition 1, maximum relative error: " << maxdev/maxelement << "\n";
    testmat = inverse_matrix * Mcopy * inverse_matrix;
    maxdev = 0.; maxelement = 0;
    for (int i = 1; i <= N; i++)
      for (int j = 1; j <= N; j++)
      {
        maxelement = std::max(maxelement,std::abs(inverse_matrix(i,j)));
        maxdev = std::max(maxdev,std::abs(inverse_matrix(i,j)-testmat(i,j)));
      }
    std::cout << "Moore-Penrose condition 2, maximum relative error: " << maxdev/maxelement << "\n";
    testmat = Mcopy * inverse_matrix;
    maxdev = 0.; maxelement = 0;
    for (int i = 1; i <= N; i++)
      for (int j = 1; j <= N; j++)
      {
        maxelement = std::max(maxelement,std::abs(testmat(i,j)));
        maxdev = std::max(maxdev,std::abs(testmat(i,j)-testmat(j,i)));
      }
    std::cout << "Moore-Penrose condition 3, maximum relative error: " << maxdev/maxelement << "\n";
    testmat = inverse_matrix * Mcopy;
    maxdev = 0.; maxelement = 0;
    for (int i = 1; i <= N; i++)
      for (int j = 1; j <= N; j++)
      {
        maxelement = std::max(maxelement,std::abs(testmat(i,j)));
        maxdev = std::max(maxdev,std::abs(testmat(i,j)-testmat(j,i)));
      }
    std::cout << "Moore-Penrose condition 4, maximum relative error: " << maxdev/maxelement << "\n";
#endif
  }
  TNT::Fortran_Matrix<T> getInv()
  { return inverse_matrix; }
};

}
