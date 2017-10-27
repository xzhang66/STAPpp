/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.0, October 14, 2017                                         */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include <fstream>

#include "Element.h"
#include "Material.h"
#include "Node.h"

using namespace std;

//! Element group class
class CElementGroup
{
private:
    
//! List of all nodes in the domain, obtained from CDomain object
    static CNode* NodeList_;
    
//! Element type of this group
    unsigned int ElementType_;

//! Number of elements in this group
    unsigned int NUME_;

//! Element List in this group
    CElement* ElementList_;

//! Number of material/section property sets in this group
    unsigned int NUMMAT_;

//! Material list in this group
    CMaterial* MaterialList_;
    
public:
    
//! Constructor
    CElementGroup();
    
//! Deconstructor
    ~CElementGroup();

//! Read element group data from stream Input
    bool Read(ifstream& Input);

//! Read bar element data from the input data file
    bool ReadBarElementData(ifstream& Input);
    
//! Return element type of this group
    unsigned int GetElementType() { return ElementType_; }
    
//! Return the number of elements in the group
    unsigned int GetNUME() { return NUME_; }

//! Return element List in this group
    CElement* GetElementList() { return ElementList_; }
    
//! Return the number of material/section property setss in this element group
    unsigned int GetNUMMAT() { return NUMMAT_; }
    
//! Return material list in this group
    CMaterial* GetMaterialList() { return MaterialList_; }

};
