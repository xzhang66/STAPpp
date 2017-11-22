/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.02, October 27, 2017                                        */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include <algorithm>

#include "Element.h"

//! Virtual deconstructor
CElement::~CElement()
{
    if (!nodes)
        delete [] nodes;
    
    if (!ElementMaterial)
        delete [] ElementMaterial;

    if (!LocationMatrix)
        delete [] LocationMatrix;
}
