// getSubGroups.cpp
//

#include "getsubgroups.hh"

#include <assert.h>
#define ASSERT assert

namespace scala
{
  //--------------------------------------------------------------
  std::vector<CCtbxSym::PointGroup> GetSubGroups(const hkl_symmetry& symm)
  // Make list of subgroups
  {
    int Nel = symm.Nelement();
    char LatticeType = symm.lattice_type();

    std::vector<CCtbxSym::PointGroup> subgroups;
    // Always put P-1 (or C1 etc) into list
    subgroups.push_back(CCtbxSym::PointGroup(LatticeType));

    // Loop pairs of symmetry elements, constructing pointgroup
    // from each pair
    for (int i=0;i<Nel;i++) {
      ASSERT(symm.NopInElement(i) >= 0);
      // 1st Symop from symmetry element i
      std::vector<double> Ri = symm.SymopInElement(0,i);
      for (int j=i;j<Nel;j++) {
	ASSERT(symm.NopInElement(j) >= 0);
	// 1st Symop from symmetry element j
	std::vector<double> Rj = symm.SymopInElement(0,j);
	// Generate point-group from symops i & j
	CCtbxSym::PointGroup PG(i, Ri, j, Rj, LatticeType);
	// Store if we don't have it already
	if (std::find(subgroups.begin(),subgroups.end(),
		      CCtbxSym::PointGroup(PG)) == subgroups.end()) {
	  subgroups.push_back(CCtbxSym::PointGroup(PG));
	  //^
	  //	  std::cout << "GetSubGroups subgroup added "
	  //		    <<  subgroups.back().RefPGname() << " reindex to standard "
	  //		    << subgroups.back().RefSGreindex().as_hkl() << "\n";
	  //^-
	} else {
	  //^
	  //	  std::cout << "GetSubGroups subgroup NOT added (same) "
	  //		    <<  subgroups.back().RefPGname() << " reindex to standard "
	  //		    << CCtbxSym::PointGroup(PG).RefSGreindex().as_hkl() << "\n";
	  //^-
	}
      }
    }
    //^
    //    std::cout << "Number of subgroups = " << subgroups.size() << "\n\n";
    //^-

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
    // We have a list of "subgroups" constructed from
    // symmetry elements - loop round subgroups adding any additional
    // elements which belong
    for (size_t k=0; k<subgroups.size();k++) {
      // Test all elements to see if they are present
      for (int j=0;j<Nel;j++) {
	if ( ! subgroups[k].HasElement(j)) {
	  // Element not yet there, test if it should be added
	  std::vector<double> Rj = symm.SymopInElement(0,j);
	  bool added = subgroups[k].AddElement(j, Rj);
	  added = added;
	}
      }
    }
    return subgroups;
  }
  //--------------------------------------------------------------
  void RemoveImplicitSubgroups(std::vector<CCtbxSym::PointGroup>& subgroups,
			       const hkl_symmetry& mergeSymm, const ReindexOp& reindex)
  // symops from mergeSymm can be reindexed into lattice frame by reindex
  {
    std::vector<bool> keepsubgroup(subgroups.size(), true); // flags to keep subgroup
    bool anyremoved = false;

    for (size_t igl=0;igl<subgroups.size();++igl) { // Loop lattice subgroups
      //^      std::cout << "\nGroup " << igl << " LGname " << subgroups[igl].LGname() << "\n";      
      // loop symmetry elements in mergeSymm
      for (int iel=0;iel<mergeSymm.Nelement();++iel) {
	// 1st op in element
	clipper::Symop invsymop = mergeSymm.ClipperSymopInElement(0,iel);
	invsymop = reindex.Symop(invsymop); // reindex to lattice cell
	std::vector<double> Rmatrix =  MVutil::SetVMat33(invsymop.rot()); // ... as matrix
	// Is this element in the subgroup?
	if (!subgroups[igl].HasElement(Rmatrix)) {
	  // fail, all elements must be present
	  keepsubgroup[igl] = false;
	  anyremoved = true;
	  break;
	}
      }
      //^
      //      if (keepsubgroup[igl]) {std::cout << "keep\n";
      //      } else {std::cout << "remove\n"; } //^-
    }
    if (anyremoved) {
      int ng = subgroups.size();
      int k = -1;
      for (int igl=0;igl<ng;++igl) { // Loop lattice subgroups
	if (keepsubgroup[igl]) {
          if (++k != igl) { // don't self-copy 
            subgroups[k] = subgroups[igl];  
	  }
	  //^
	  //	  std::cout << "Keeping group " << subgroups[igl].LGname() << "\n"; //^
	  //	} else {
	  //^
	  //	  std::cout << "Removing group " << subgroups[igl].LGname() << "\n"; //^
	}
      }
      subgroups.resize(++k);
      //^      std::cout << "\nRemaining groups:\n";
      //      for (int igl=0;igl<subgroups.size();++igl) { // Loop lattice subgroups
      //      	std::cout << "Lattice subgroup " << subgroups[igl].LGname() << "\n"; //^
      //      }
      //^-
    }
    if (subgroups.size() <= 0) {
      subgroups.assign(1,CCtbxSym::PointGroup(mergeSymm.symbol_xHM()));
    }      
  }
} // scala
