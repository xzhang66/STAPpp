/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.02, October 27, 2017                                        */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include <string>

//! CSkylineMatrix class is used to store the FEM stiffness matrix in skyline storage
template <class T_>
class CSkylineMatrix
{
//! Store the stiffness matrkix in skyline storage
    T_* data_;
    
//! Dimension of the stiffness matrix
    unsigned int NEQ_;
    
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

//! operator (i,j) where i and j number from 1
//! For the sake of efficiency, the index bounds are not checked
    inline T_& operator()(unsigned int i, unsigned int j);
    
//! operator (i) where i numbers from 1
//! For the sake of efficiency, the index bounds are not checked
    inline T_ operator()(unsigned int i);

//! Allocate storage for the skyline matrix
    inline void Allocate();

//! Return pointer to the ColumnHeights_
    inline unsigned int* GetColumnHeights();
    
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
    
    data_ = nullptr;
    ColumnHeights_ = nullptr;
    DiagonalAddress_ = nullptr;
}

template <class T_>
inline CSkylineMatrix<T_>::CSkylineMatrix(unsigned int N)
{
    NEQ_ = N;

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

//! operator function (i,j) where i and j number from 1
template <class T_>
inline T_& CSkylineMatrix<T_>::operator()(unsigned int i, unsigned int j)
{
    if (j >= i)
        return data_[DiagonalAddress_[j - 1] + (j - i) - 1];
    else
        return data_[DiagonalAddress_[i - 1] + (i - j) - 1];
}

//! operator function (i) where i numbers from 1
template <class T_>
inline T_ CSkylineMatrix<T_>::operator()(unsigned int i)
{
    return data_[i];
}

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
