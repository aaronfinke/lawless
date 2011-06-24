// cctbx_utils.cpp
// Utility routines using cctbx calls

#include "cctbx_utils.hh"
#include "util.hh"
#include "scala_util.hh"
#include "string_util.hh"

namespace MVutil
{
  // * *  * *  * *  * *  * *  * *  * *  * *  * *  * *  * *  * * 
  // cctbx
  //--------------------------------------------------------------
 scitbx::mat3<double> SetSMat33(const std::vector<double>& Rmatrix)
    // Make rt_max object from matrix as vector, translation 000
  {
    return scitbx::mat3<double>(Rmatrix[0],Rmatrix[1],Rmatrix[2],
			     Rmatrix[3],Rmatrix[4],Rmatrix[5],
			     Rmatrix[6],Rmatrix[7],Rmatrix[8]);
  }
  //--------------------------------------------------------------
  scitbx::vec3<double> SetSvec3(const clipper::Vec3<double>& v3)
  // Make scitbx::vec3 from clipper::Vec3 <double>
  {
    return scitbx::vec3<double>(v3[0],v3[1],v3[2]);
  }
  //--------------------------------------------------------------
  sgtbx::rt_mx SetRtMx(scitbx::mat3<double>& cdmat,
		       const int r_den, const int t_den)
    // Make rt_max object from matrix, translation 000
  {
    return sgtbx::rt_mx(cdmat, scitbx::vec3<double>(0.0,0.0,0.0),
			r_den, t_den);

  }
  //--------------------------------------------------------------
  sgtbx::rt_mx SetRtMx(const std::vector<double>& Rmatrix,
		       const int r_den, const int t_den)
    // Make rt_max object from matrix as vector, translation 000
  {
    return sgtbx::rt_mx(SetSMat33(Rmatrix), scitbx::vec3<double>(0.0,0.0,0.0),
			r_den, t_den);
  }
  //--------------------------------------------------------------
  sgtbx::rt_mx SetRtMx(const clipper::RTop<double>& RT,
		       const int r_den, const int t_den)
  // Make rt_max object from RTop, matrix & translation
  {
    return sgtbx::rt_mx(SetSMat33(SetVMat33(RT.rot())),
			SetSvec3(RT.trn()),
			r_den, t_den);

  }
  // * *  * *  * *  * *  * *  * *  * *  * *  * *  * *  * *  * * 
  //--------------------------------------------------------------
  std::string FormatReindex_as_hkl(const scitbx::mat3<double>& op)
  {
    sgtbx::rt_mx R(op.transpose(), scitbx::vec3<double>(0.0,0.0,0.0),
		   sgtbx::cb_r_den,sgtbx::cb_t_den);
    std::string chb = R.as_xyz(false, false, "hkl");
    chb = StringUtil::Strip(chb, '*');
    return "["+chb+"]";
  }
  //--------------------------------------------------------------
}
