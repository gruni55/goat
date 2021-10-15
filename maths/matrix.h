/*****************************************************************//**
 * \file   matrix.h
 * \brief  This file contains the definition of a template for 3x3 matrices, which can be used together with the Vector class, defined in this project. 
 * 
 * \author Thomas
 * \date   July 2021
 *********************************************************************/
#ifndef MATRIX_H
#define MATRIX_H
#include "vector.h"
#include <iostream>
/**
 * @brief This class represents a threedimensional (numeric) Matrix as Template
 * 
 * The Matrix class represents a template class for 3x3 numeric matrices. Operators for standard matrix operation are provided. 
 * This class is intended to work with the Vector template class defined in this SDK.
 */
template<class T> class Matrix 
{
 public :
 Matrix () ///< Standard constructor, initializes the Matrix with zero.
 {  
 for (int i=0; i<3; i++)
  for (int j=0; j<3; j++)
   M[i][j]=0.0; 
 }

 void binWrite (std::ofstream &os) ///< writes the Matrix into a binary file, represented by os
 {
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
  os.write ((char *) &M[i][j], sizeof (M[i][j])); 
 }

 void binRead (std::ifstream &is) ///< reads the Matrix from a binary file, represented by is
 {
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
  is.read ((char *) &M[i][j], sizeof (M[i][j])); 
 }
 

 /**
  * @brief Main constructor for the template class Matrix.
  * This constructor uses three vectors as parameters, which represents the three columns of the matrix
  * \param a first column
  * \param b second column
  * \param c third column
  */
 Matrix (Vector<T> a, Vector<T> b, Vector<T> c)  
 {
  for (int i=0; i<3; i++)
     M[i][0]=a[i];
  for (int i=0; i<3; i++)
     M[i][1]=b[i]; 
  for (int i=0; i<3; i++)
     M[i][2]=c[i]; 
 }
 
 Matrix& operator = (const Matrix& A) ///< assignment operator 
 {
 if (this==&A) return *this;
 for (int i=0; i<3; i++)
  for (int j=0; j<3; j++) 
    M[i][j]=A.M[i][j];
 return *this;
 }

 bool operator == (const Matrix& A) ///< comparison operator (two matrices are equal, if all components are exactly the same)
 {
  bool Erg=true;
  for (int i=0; i<3; i++) 
   for (int j=0; j<3; j++)
      Erg=Erg && (M[i][j]==A.M[i][j]);
   return Erg;
 } 
 
 bool operator != (const Matrix& M) ///< negative comparison operator, returns false, if at least one of the components are not equal
 {
  bool Erg=false;
  for (int i=0; i<3; i++) 
   for (int j=0; j<3; j++)
      Erg=Erg || (M[i][j]!=M.M[i][j]);
   return Erg;
 } 
  
 /**
  * @brief This operator returns the value of the j-th element in the i-th column
  * @param i index of the column
  * @param j index of the row
  */
 T& operator () (int i, int j) { return M[i][j]; } 

 /** @name special operators
  * @brief These operators are the common definitions of the corresponding matrix operations
  */
 ///@{
 Matrix operator + (const Matrix& A) ///< adding two matrices (same type)
 {
 Matrix<T> Erg;
  for (int i=0; i<3; i++)
   for (int j=0; j<3; j++)
    Erg.M[i][j]=M[i][j]+A.M[i][j];
  return Erg;
 }

 Matrix& operator +=(const Matrix& A) 
 {
  for (int i=0; i<3; i++)
   for (int j=0; j<3; j++)
    M[i][j]+=A.M[i][j];
  return *this;
 }

Matrix operator - ()  ///< minus sign operator
 { 
  Matrix Erg;
  for (int i=0; i<3; i++)
   for (int j=0; j<3; j++)
    Erg.M[i][j]=-M[i][j];
  return Erg;  
 }

 Matrix& operator -=(const Matrix& A) ///< subtraction operator
 {
  for (int i=0; i<3; i++)
   for (int j=0; j<3; j++)
    M[i][j]-=A.M[i][j];
  return *this;
 }
 
 Matrix& operator *= (const Matrix& A) ///< multiplication with matrix A
 {
  T m[3][3];
  for (int i=0; i<3; i++)
   for (int j=0; j<3; j++)
   {
     m[i][j]=0;
     for (int l=0; l<3; l++)
     m[i][j]+=M[i][l]*A.M[l][j];
    }
   for (int i=0; i<3; i++)
     for (int j=0; j<3; j++)
      M[i][j]=m[i][j];
   return *this;
 }

 Matrix& operator /=(const T& x) ///< division by a scalar
 {
  for (int i=0; i<3; i++) 
   for (int j=0; j<3; j++)
    M[i][j]/=x; 
   return *this; 
 }

 Matrix& operator *=(const T& x) ///< multiplication with a scalar
 {
  for (int i=0; i<3; i++) 
   for (int j=0; j<3; j++)
    M[i][j]*=x; 
   return *this; 
 }
 
 
 friend Matrix operator - (const Matrix& A, const Matrix& B) ///< subtraction of two matrices (same type)
 {
 Matrix Erg;
  for (int i=0; i<3; i++)
   for (int j=0; j<3; j++)
    Erg.M[i][j]=A.M[i][j]-B.M[i][j];
  return Erg;
 }
 
 friend Matrix operator * (const Matrix& A, const Matrix& B) ///< multiplication of two matrices (same type)
 {
  Matrix Erg;
  for (int i=0; i<3; i++)
   for (int j=0; j<3; j++)
     for (int l=0; l<3; l++)
      Erg.M[i][j]+=A.M[i][l]*B.M[l][j];
  return Erg; 
 }
 
 friend Matrix operator * (const Matrix& A, const T& x) ///< multiplication of matrix with a scalar
 {
  Matrix Erg;
  for (int i=0; i<3; i++)
   for (int j=0; j<3; j++)
      Erg.M[i][j]=A.M[i][j]*x;
  return Erg; 
 }
 
 friend Matrix operator * (const T& x, const Matrix& A) ///< multiplication of scalar with a matrix
 {
  Matrix Erg;
  for (int i=0; i<3; i++)
   for (int j=0; j<3; j++)
      Erg.M[i][j]=x*A.M[i][j];
  return Erg; 
 }
  
 friend Vector<T> operator * (const Matrix& A, const Vector<T>& r) ///< multiplication of a matrix with a vector
 {
  Vector<T> Erg;
  for (int i=0; i<3; i++) 
   for (int j=0; j<3; j++)
    Erg[i]+=A.M[i][j]*r[j];
  return Erg;  
 }
 
 friend Vector<T> operator * (const Vector<T>& r, const Matrix& A) ///< multiplication of a vector with a matrix
 {
  Vector<T> Erg;
  for (int i=0; i<3; i++) 
   for (int j=0; j<3; j++)
    Erg[i]+=A.M[i][j]*r[j];
  return Erg;  
 }
  
friend Matrix operator / (const Matrix& A, const T& x) ///< division of a matrix with a scalar
 {
  Matrix Erg;
  for (int i=0; i<3; i++)
   for (int j=0; j<3; j++)
      Erg.M[i][j]=A.M[i][j] / x;
  return Erg; 
 }
 ///@} 

/**
 * @brief Output of the matrix into an ostream
 * 
 * \param os stream 
 * \param A  matrix
 * \return changed stream
 */
 friend std::ostream& operator << (std::ostream& os, const Matrix& A) 
 {
  for (int i=0; i<3; i++)
   os << A.M[i][0] << "  " << A.M[i][1] << "  " << A.M[i][2] << std::endl;
  return os; 
 }  

 /**
  * @brief Input of the matrix from an istream
  *
  * \param is stream
  * \param A  matrix
  * \return changed stream
  */

 friend std::istream& operator >> (std::istream& is, Matrix& A)
 {
  for (int i=0; i<3; i++)
   for (int j=0; j<3; j++)
    is >> A.M[i][j];
  return is; 
 }
 

 /**
  * @brief Calculates the transpose of the matrix A.
  * 
  * \param A Matrix to be transposed (represented by Aij)
  * \return transposed matrix (=Aji)
  */
 friend Matrix transpose (const Matrix& A)
 { 
  Matrix Erg;
  for (int i=0; i<3; i++)
   for (int j=0; j<3; j++)
    Erg.M[i][j]=A.M[j][i];   
  return Erg;
 }  
  
 /**
  * @brief Calculates the trace of matrix A, i.e. the sum over all diagonal elements.
  * 
  * \param A Matrix
  * \return trace
  */
 friend T trace(const Matrix& A)
 {
  T Erg=0;
  for (int i=0; i<3; i++)
   Erg+=A.M[i][i];
  return Erg;
 }
 
 /**
  * @brief Calculates the determinant of the matrix A.
  * 
  * \param A Matrix
  * \return Determinant of A
  */
 friend T det(const Matrix& A)
 {
  T Erg; 
  Erg=A.M[0][0]*(A.M[1][1]*A.M[2][2]-A.M[2][1]*A.M[1][2])
     -A.M[0][1]*(A.M[1][0]*A.M[2][2]-A.M[2][0]*A.M[1][2])
     +A.M[0][2]*(A.M[1][0]*A.M[2][1]-A.M[2][0]*A.M[1][1]);
  return Erg;
 }

 /**
  * @brief This function returns a matrix where all elements of the i-th column and the j-th row set to zero
  * 
  * \param M Matrix
  * \param i column
  * \param j row
  * \return Exchanged matrix
  */
 friend Matrix cutMatrix (const Matrix& M,int i, int j)
 {
  Matrix Erg;
  for (int k=0; k<3; k++)
   for (int l=0; l<3; l++)
    {
     if ((k!=i)&&(l!=j)) Erg.M[k][l]=M.M[k][l];
     if ((l==j)&&(k==i)) Erg.M[k][l]=1; 
     else
     {
       if (l==j) Erg.M[k][l]=0;
       if (k==i) Erg.M[k][l]=0;
     }
    }
  return Erg; 
 }
 /**
  * @brief Calculates the inverse of a matrix, if possible.
  * 
  * \param A Matrix to invert 
  * \param invertable true, if A is invertable otherwise false
  * \return Inverse of matrix A, if A is invertable, otherwise the function returns a zero matrix.
  */
 friend Matrix invert (const Matrix<T>& A,bool &invertable)
 {
  Matrix Erg;
  T C;
  C=det(A);
  if (C==0.0) invertable=false;
  else
  {
   Erg.M[0][0]=A.M[1][1]*A.M[2][2] - A.M[2][1]*A.M[1][2];
   Erg.M[1][0]=A.M[2][0]*A.M[1][2] - A.M[1][0]*A.M[2][2];
   Erg.M[2][0]=A.M[1][0]*A.M[2][1] - A.M[2][0]*A.M[1][1];

   Erg.M[0][1]=A.M[2][1]*A.M[0][2] - A.M[0][1]*A.M[2][2];
   Erg.M[1][1]=A.M[0][0]*A.M[2][2] - A.M[2][0]*A.M[0][2];
   Erg.M[2][1]=A.M[2][0]*A.M[0][1] - A.M[0][0]*A.M[2][1];

   Erg.M[0][2]=A.M[0][1]*A.M[1][2] - A.M[1][1]*A.M[0][2];
   Erg.M[1][2]=A.M[1][0]*A.M[0][2] - A.M[0][0]*A.M[1][2];
   Erg.M[2][2]=A.M[0][0]*A.M[1][1] - A.M[0][1]*A.M[1][0];
   invertable=true;
  }
 return Erg/C; 
 }
 /**
  * @brief Calculates the inverse of a matrix, if possible.
  *
  * \param A Matrix to invert
  * \return Inverse of matrix A, if A is invertable, otherwise the function returns a zero matrix.
  */

 friend Matrix invert (const Matrix<T>& A)
 {
  bool dummy;
  return invert (A,dummy); 
 } 

 void clear() ///< clears the Matrix (all values are zero)
 {
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      M[i][j]=0;
 } 

 T M[3][3]; 
};
/** @name special operators for mixed types
 * @brief Operators for calculations with mixed types.
 */
 ///@{
// Matrix-Matrix Operatoren mit gemischten Typen (std::complex<double> /double)
Matrix<std::complex<double> > operator * (const Matrix<double>&, const Matrix<std::complex<double> >&);
Matrix<std::complex<double> > operator * (const Matrix<std::complex<double> >&, const Matrix<double>&);
Matrix<std::complex<double> > operator + (const Matrix<double>&, const Matrix<std::complex<double> >&);
Matrix<std::complex<double> > operator + (const Matrix<std::complex<double> >&, const Matrix<double>&);
Matrix<std::complex<double> > operator - (const Matrix<double>&, const Matrix<std::complex<double> >&);
Matrix<std::complex<double> > operator - (const Matrix<std::complex<double> >&, const Matrix<double>&);

// Matrix-Vektor Operatoren mit gemischten Typen (std::complex<double> /double)
Vector<std::complex<double> > operator * (const Matrix<double>&, const Vector<std::complex<double> >&);
Vector<std::complex<double> > operator * (const Matrix<std::complex<double> >&, const Vector<double>&);

// Matrix-Skalar Operatoren mit gemischten Typen (std::complex<double> /double)
Matrix<std::complex<double> > operator * (const Matrix<double>&, const std::complex<double> &);
Matrix<std::complex<double> > operator * (const Matrix<std::complex<double> >&, const double&);
Matrix<std::complex<double> > operator * (const double&, const Matrix<std::complex<double> >&);
Matrix<std::complex<double> > operator * (const std::complex<double> &, const Matrix<double>&);
Matrix<std::complex<double> > operator / (const Matrix<double>&, const std::complex<double> &);
Matrix<std::complex<double> > operator / (const Matrix<std::complex<double> >&, const double&);
///@}

// Drehmatrizen um die 3 Raumachsen
/** @name Rotation matrices
 *  @brief Rotation matrices around the principal axis of the cartesian coordinate system (angles are given in radiants)
 */
 ///@{
Matrix<double> Dx(double phi); ///< Rotation around x axis
Matrix<double> Dy(double phi); ///< Rotation around y axis
Matrix<double> Dz(double phi); ///< Rotation around z axis
Matrix<double> rotMatrix(const Vector<double> a, double gamma); ///<   Rotation matrix for rotation around the axis a by the angle gamma
/**
 * @brief Rotation matrix calculated from rotation around x-axis, y-axis and z-axis 
 * @param alpha Rotation angle around x-axis
 * @param  beta Rotation angle around y-axis
 * @param gamma Rotation angle around z-axis
 * @return Rotation Matrix D=Dx(alpha)*Dy(beta)*Dz(gamma)
 */
Matrix<double> rotMatrix(double alpha, double beta, double gamma); 

/**
 * @brief Calculates rotation matrix for a rotation around an axis passing through a given point and pointing in a defined direction. 
 * @param n Point on the axis 
 * @param k Direction of the axis
 * @param gamma Rotation angle
 * @return Matrix, which belongs to the given rotation.
 */
Matrix <double> rotMatrixA(Vector<double> n, Vector<double> k, double gamma); 
Matrix<double> rotMatrix(Vector<double> P, double dtheta, double dphi);
///@}
Matrix<double> unity();  ///< unity matrix (double precision)
Matrix<std::complex<double> > cunity (); ///< unity matrix (complex)
Matrix<double> null(); ///< Null-matrix (all components are zero) 
/**
 * @brief Calculates the transformation matrices from the laboratory coordinate system into a local system and vice versa.
 * \param e0 direction of the first local axis
 * \param e1 direction of the second local axis
 * \param e2 direction of the third local axis
 * \return H transformation matrix from laboratory system into local coordinate system
 * \return R transformation matrix from local system into laboratory coordinate system
 */
void trafo (const Vector<double> & e0, 
            const Vector<double> & e1, 
            const Vector<double> & e2, 
            Matrix<double> & H,
            Matrix<double> & R);

const Matrix<double> UNITY=Matrix<double>(Vector<double>(1,0,0),Vector<double>(0,1,0),Vector<double>(0,0,1)); ///< Unity matrix (double-valued)
const Matrix<std::complex<double> > CUNITY=Matrix<std::complex<double> >(Vector<std::complex<double> >(1,0,0),Vector<std::complex<double> >(0,1,0),Vector<std::complex<double> >(0,0,1));///< (complex-valued)

#endif // MATRIX_H
