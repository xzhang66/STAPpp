/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

template <class T_>
class CSkylineMatrix
{
    T_* data_;  // Store the skyline matrix
    
    unsigned int NEQ_;  // Dimension of the matrix
    unsigned int NWK_;  // Size of data_
    
    unsigned int* ColumnHeights_;   // Column hights
    unsigned int* DiagonalAddress_; // Diagonal address of all columns in data_

public:

    // constructor
    inline CSkylineMatrix(int N);
    
    // destructor
    inline ~CSkylineMatrix();

    // operator []
    inline T_& operator()(int i, int j);

    // Allocate storage for the skyline matrix
    inline void Allocate();   // Allocate storage for the matrix
    
    // Set the diagonal address of ith column
    inline void DiagonalAddress(unsigned int i, unsigned int address);
    
    // Set the hight of ith column
    inline void ColumnHeight(unsigned int i, unsigned int height);

    // Set the element (i,j) to value
    inline void set(int i, int j, T_ value);

    // Return the address of element (i,j) in data_ (numbering from 0)
    inline unsigned int address(int i, int j) const;

    // Return the dimension of the matrix (NEQ_)
    inline unsigned int dim() const;
    
    // Return the size of the matrix (NWK_)
    inline unsigned int size() const;

}; /* class definition */

// constructor functions
template <class T_>
inline CSkylineMatrix<T_>::CSkylineMatrix(int N)
{
    NEQ_ = N;

    data_ = NULL;
    ColumnHeights_ = new unsigned int [NEQ_];
    DiagonalAddress_ = new unsigned int [NEQ_ + 1];
}

// destructor function
template <class T_>
inline CSkylineMatrix<T_>::~CSkylineMatrix<T_>()
{
    delete[] ColumnHeights_;
    delete[] DiagonalAddress_;
    
    if (data_)
        delete[] data_;
}

// operator function
template <class T_>
inline T_& CSkylineMatrix<T_>::operator()(int i, int j)
{
    return data_[address(i,j)];
}

// Allocate storage for the matrix
template <class T_>
inline void CSkylineMatrix<T_>::Allocate()
{
    NWK_ = DiagonalAddress_[NEQ_] - DiagonalAddress_[0];
    data_ = new T_[NWK_];
}

// Set the diagonal address of ith column
template <class T_>
inline void CSkylineMatrix<T_>::DiagonalAddress(unsigned int i, unsigned int address)
{
    DiagonalAddress_[i-1] = address;
}

// Set the hight of ith column
template <class T_>
inline void CSkylineMatrix<T_>::ColumnHeight(unsigned int i, unsigned int height)
{
    ColumnHeights_[i-1] = height;
}

// Return the address of element (i,j) in data_ (numbering from 0)
template <class T_>
inline unsigned int CSkylineMatrix<T_>::address(int i, int j) const
{
    if (j >= i)
        return DiagonalAddress_[j - 1] + (j - i) - 1;
    else
        return DiagonalAddress_[i - 1] + (i - j) - 1;
}

// Set the element (i,j) to value
template <class T_>
inline void CSkylineMatrix<T_>::set(int i, int j, T_ value)
{
    data_[address(i,j)] = value;
}

// Return the dimension of the matrix (NEQ_)
template <class T_>
inline unsigned int CSkylineMatrix<T_>::dim() const
{
    return(NEQ_);
}

// Return the size of the matrix (NWK_)
template <class T_>
inline unsigned int CSkylineMatrix<T_>::size() const
{
   return(NWK_);
}
