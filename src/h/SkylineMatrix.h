/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include <string>
#include <climits>

#ifdef _DEBUG_
#include "Outputter.h"
#endif

//! CSkylineMatrix class is used to store the FEM stiffness matrix in skyline storage
template <class T_>
class CSkylineMatrix
{
    
private:
//! Store the stiffness matrkix in skyline storage
    T_* data_;
    
//! Dimension of the stiffness matrix
    unsigned int NEQ_;

//! Maximum half bandwith
    unsigned int MK_;

//! Size of the storage used to store the stiffness matrkix in skyline
    unsigned int NWK_;

//! Column hights
    unsigned int* ColumnHeights_;
    
//! Diagonal address of all columns in data_
    unsigned int* DiagonalAddress_;
    
public:

//! constructors
    inline CSkylineMatrix();
    inline CSkylineMatrix(unsigned int N);
    
//! destructor
    inline ~CSkylineMatrix();

//! operator (i,j) where i and j numbering from 1
//! For the sake of efficiency, the index bounds are not checked
    inline T_& operator()(unsigned int i, unsigned int j);

#ifdef _DEBUG_
//! operator (i) where i numbering from 1
//! For the sake of efficiency, the index bounds are not checked
    inline T_& operator()(unsigned int i);
#endif
    
//! Allocate storage for the skyline matrix
    inline void Allocate();
    
//! Calculate the column height, used with the skyline storage scheme
    void CalculateColumnHeight(unsigned int* LocationMatrix, size_t ND);

//! Calculate the maximum half bandwidth ( = max(ColumnHeights) + 1 )
    void CalculateMaximumHalfBandwidth();

//! Calculate address of diagonal elements in banded matrix
//! Caution: Address is numbered from 1 !
    void CalculateDiagnoalAddress();

//! Assemble the element stiffness matrix to the global stiffness matrix
    void Assembly(double* Matrix, unsigned int* LocationMatrix, size_t ND);

//! Return pointer to the ColumnHeights_
    inline unsigned int* GetColumnHeights();

//! Return the maximum half bandwidth
    inline unsigned int GetMaximumHalfBandwidth() const;

//! Return pointer to the DiagonalAddress_
    inline unsigned int* GetDiagonalAddress();

//! Return the dimension of the stiffness matrix
    inline unsigned int dim() const;
    
//! Return the size of the storage used to store the stiffness matrkix in skyline
    inline unsigned int size() const;

}; /* class definition */

//! constructor functions
template <class T_>
inline CSkylineMatrix<T_>::CSkylineMatrix()
{
    NEQ_ = 0;
    MK_  = 0;
    NWK_ = 0;
    
    data_ = nullptr;
    ColumnHeights_ = nullptr;
    DiagonalAddress_ = nullptr;
}

template <class T_>
inline CSkylineMatrix<T_>::CSkylineMatrix(unsigned int N)
{
    NEQ_ = N;
    MK_  = 0;
    NWK_ = 0;

    data_ = nullptr;
    
    ColumnHeights_ = new unsigned int [NEQ_];
    for (unsigned int i = 0; i < NEQ_; i++)
        ColumnHeights_[i] = 0;

    DiagonalAddress_ = new unsigned int [NEQ_ + 1];
    for (unsigned int i = 0; i < NEQ_ + 1; i++)
        DiagonalAddress_[i] = 0;
}

//! destructor function
template <class T_>
inline CSkylineMatrix<T_>::~CSkylineMatrix<T_>()
{
    if (ColumnHeights_)
        delete[] ColumnHeights_;
    
    if (DiagonalAddress_)
        delete[] DiagonalAddress_;
    
    if (data_)
        delete[] data_;
}

//! operator function (i,j) where i and j numbering from 1
template <class T_>
inline T_& CSkylineMatrix<T_>::operator()(unsigned int i, unsigned int j)
{
    if (j >= i)
        return data_[DiagonalAddress_[j - 1] + (j - i) - 1];
    else
        return data_[DiagonalAddress_[i - 1] + (i - j) - 1];
}

#ifdef _DEBUG_
//! operator function (i) where i numbering from 1
template <class T_>
inline T_& CSkylineMatrix<T_>::operator()(unsigned int i)
{
    return data_[i];
}
#endif

//! Allocate storage for the matrix
template <class T_>
inline void CSkylineMatrix<T_>::Allocate()
{
    NWK_ = DiagonalAddress_[NEQ_] - DiagonalAddress_[0];

    data_ = new T_[NWK_];
    for (unsigned int i = 0; i < NWK_; i++)
        data_[i] = T_(0);
}

//! Return pointer to the ColumnHeights_
template <class T_>
inline unsigned int* CSkylineMatrix<T_>::GetColumnHeights()
{
    return ColumnHeights_;
}

//! Return the maximum half bandwidth
template <class T_>
inline unsigned int CSkylineMatrix<T_>::GetMaximumHalfBandwidth() const
{
    return(MK_);
}

//! Return pointer to the DiagonalAddress_
template <class T_>
inline unsigned int* CSkylineMatrix<T_>::GetDiagonalAddress()
{
    return DiagonalAddress_;
}

//! Return the dimension of the stiffness matrix
template <class T_>
inline unsigned int CSkylineMatrix<T_>::dim() const
{
    return(NEQ_);
}

//! Return the size of the storage used to store the stiffness matrkix in skyline
template <class T_>
inline unsigned int CSkylineMatrix<T_>::size() const
{
   return(NWK_);
}

//  Calculate the column height, used with the skyline storage scheme
template <class T_>
void CSkylineMatrix<T_>::CalculateColumnHeight(unsigned int* LocationMatrix, size_t ND)
{    
//  Look for the row number of the first non-zero element
    unsigned int nfirstrow = INT_MAX;
    for (unsigned int i = 0; i < ND; i++)
        if (LocationMatrix[i] && LocationMatrix[i] < nfirstrow)
            nfirstrow = LocationMatrix[i];
    
//  Calculate the column height contributed by this element
    for (unsigned int i = 0; i < ND; i++)
    {
        unsigned int column = LocationMatrix[i];
        if (!column)
            continue;
        
        unsigned int Height = column - nfirstrow;
        if (ColumnHeights_[column-1] < Height) ColumnHeights_[column-1] = Height;
    }
}

// Maximum half bandwidth ( = max(ColumnHeights) + 1 )
template <class T_>
void CSkylineMatrix<T_>::CalculateMaximumHalfBandwidth()
{
    MK_ = ColumnHeights_[0];
    
    for (unsigned int i=1; i<NEQ_; i++)
        if (MK_ < ColumnHeights_[i])
            MK_ = ColumnHeights_[i];
    
    MK_ = MK_ + 1;
}

//    Assemble the banded global stiffness matrix (skyline storage scheme)
template <class T_>
void CSkylineMatrix<T_>::Assembly(double* Matrix, unsigned int* LocationMatrix, size_t ND)
{
//  Assemble global stiffness matrix
    for (unsigned int j = 0; j < ND; j++)
    {
        unsigned int Lj = LocationMatrix[j];    // Global equation number corresponding to jth DOF of the element
        if (!Lj) continue;
        
//      Address of diagonal element of column j in the one dimensional element stiffness matrix
        unsigned int DiagjElement = (j+1)*j/2;
        
        for (unsigned int i = 0; i <= j; i++)
        {
            unsigned int Li = LocationMatrix[i];    // Global equation number corresponding to ith DOF of the element
            
            if (!Li) continue;
            
            (*this)(Li,Lj) += Matrix[DiagjElement + j - i];
        }
    }
    
    return;
}

//    Calculate address of diagonal elements in banded matrix
//    Caution: Address is numbered from 1 !
template <class T_>
void CSkylineMatrix<T_>::CalculateDiagnoalAddress()
{
    //    Calculate the address of diagonal elements
    //    M(0) = 1;  M(i+1) = M(i) + H(i) + 1 (i = 0:NEQ)
    DiagonalAddress_[0] = 1;
    for (unsigned int col = 1; col <= NEQ_; col++)
        DiagonalAddress_[col] = DiagonalAddress_[col - 1] + ColumnHeights_[col-1] + 1;
    
#ifdef _DEBUG_
    COutputter* Output = COutputter::GetInstance();
    Output->PrintDiagonalAddress();
#endif
    
}

