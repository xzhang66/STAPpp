/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.0, October 14, 2017                                         */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "ElementGroup.h"
#include "Bar.h"
#include "Domain.h"

CNode* CElementGroup::NodeList_ = nullptr;

//! Constructor
CElementGroup::CElementGroup()
{
    if (!NodeList_)
    {
        CDomain* FEMData = CDomain::Instance();
        NodeList_ = FEMData->GetNodeList();
    }
    
    ElementType_ = 0;
    
    NUME_ = 0;
    ElementList_ = nullptr;
    
    NUMMAT_ = 0;
    MaterialList_ = nullptr;
}

//! Deconstructor
CElementGroup::~CElementGroup()
{
    if (ElementList_)
        delete [] ElementList_;
    
    if (MaterialList_)
        delete [] MaterialList_;
}

//! Read element group data from stream Input
bool CElementGroup::Read(ifstream& Input)
{
    Input >> ElementType_ >> NUME_ >> NUMMAT_;
    
    switch (ElementType_)
    {
        case 1:    // Bar element
            if (!ReadBarElementData(Input))
                return false;
            break;
            
        default:    // Invalid element type
            cout << "*** Error *** Elment type " << ElementType_ <<  " has not been implemented.\n\n";
            return false;
    }

    return true;
}

//  Read bar element data from the input data file
bool CElementGroup::ReadBarElementData(ifstream& Input)
{
//  Read material/section property lines
    MaterialList_ = new CBarMaterial[NUMMAT_];    // Materials for group EleGrp
    CBarMaterial* mlist = dynamic_cast<CBarMaterial*>(MaterialList_);
    
//  Loop over for all material property sets in this element group
    for (unsigned int mset = 0; mset < NUMMAT_; mset++)
        if (!mlist[mset].Read(Input, mset))
            return false;
    
//  Read element data lines
    ElementList_ = new CBar[NUME_];    // Elements of gorup EleGrp
    CBar* elist = dynamic_cast<CBar*>(ElementList_);
    
//  Loop over for all elements in this element group
    for (unsigned int Ele = 0; Ele < NUME_; Ele++)
        if (!elist[Ele].Read(Input, Ele, MaterialList_, NodeList_))
            return false;
    
    return true;
}
