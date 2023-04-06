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

#include "Node.h"
#include "Material.h"

using namespace std;

template <class type> void clear( type* a, unsigned int N );	// Clear an array

//!	Element base class
/*!	All type of element classes should be derived from this base class */
class CElement
{
protected:

//!	Number of nodes per element
	unsigned int NEN_;

//!	Nodes of the element
	CNode** nodes_;

//!	Material of the element
	CMaterial* ElementMaterial_;	//!< Pointer to an element of MaterialSetList[][]
    
//! Location Matrix of the element
    unsigned int* LocationMatrix_;

//! Dimension of the location matrix
    unsigned int ND_;

public:

//!	Constructor
	CElement() : NEN_(0), nodes_(nullptr), ElementMaterial_(nullptr) {}

//! Virtual deconstructor
    virtual ~CElement() {
        if (nodes_)
            delete [] nodes_;
        
        if (ElementMaterial_)
            delete [] ElementMaterial_;
        
        if (LocationMatrix_)
            delete [] LocationMatrix_;
    }

//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList) = 0;

//!	Write element data to stream
	virtual void Write(COutputter& output) = 0;

//! Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
    virtual void GenerateLocationMatrix()
    {
        unsigned int i = 0;
        for (unsigned int N = 0; N < NEN_; N++)
            for (unsigned int D = 0; D < CNode::NDF; D++)
                LocationMatrix_[i++] = nodes_[N]->bcode[D];
    }

//! Return the size of the element stiffness matrix (stored as an array column by column)
    virtual unsigned int SizeOfStiffnessMatrix()
    {
        unsigned int size = 0;
        for (int i=1; i<= ND_; i++)
            size += i;
        
        return size;
    }

//!	Calculate element stiffness matrix (Upper triangular matrix, stored as an array column by colum)
	virtual void ElementStiffness(double* stiffness) = 0; 

//!	Calculate element stress 
	virtual void ElementStress(double* stress, double* Displacement) = 0;

//! Return number of nodes per element
    inline unsigned int GetNEN() { return NEN_; }
    
//!	Return nodes of the element
	inline CNode** GetNodes() { return nodes_; }

//!	Return material of the element
	inline CMaterial* GetElementMaterial() { return ElementMaterial_; }
    
    //! Return the Location Matrix of the element
    inline unsigned int* GetLocationMatrix() { return LocationMatrix_; }
    
    //! Return the dimension of the location matrix
    inline unsigned int GetND() { return ND_; }
};
