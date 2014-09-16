// selectsdcorrreflections.cpp
//
// Mark reflections with accept/reject flags to be used
// for SD correction optimisation
// The idea is to reduce the number of observations to speed up optimisation
//
// 1) no singletons
// 2) Roughly evenly distributed of |E^2| in lower intensity bins
// 3) all higher intensities E^2 > E2min
// 4) no ice rings
//
// Note that analysis is done on intensity bins, not E^2 bins, so the
// distribution will not be even on intensity bins: this doesn't matter
//

#include "selectsdcorrreflections.hh"
#include "normalise.hh"
#include "selectedobservations.hh"
#include "scala_util.hh"
#include "intensitybin.hh"

namespace scala {

  // ------------------------------------------------------------
  int SelectSDcorrReflections(hkl_unmerge_list& hkl_list,
                              const all_controls& controls,
                              const int& Nbintarget)
  // On entry:
  //   hkl_list    reflection list, scales applied if needed
  //   Nbintarget  target minimum number of reflections / intensity bin
  //
  // On exit:
  //   hkl_list    reflection list, reflection accept flags updated
  //
  // returns number of reflections accepted
  //
  {
    hkl_list.ResetReflAccept();  // set to accept everything
    // Overall Normalisation
    double MinIsigRatio = -1.0;  // no resolution cutoff
    bool Overall = true;
    Rings NoRings;
    ResoRange ResRangeN = hkl_list.ResLimRange();
    Normalise NormRes = SetNormalise(hkl_list, MinIsigRatio, Overall,
                                     ResRangeN, NoRings, 0);

    // Intensity bins etc
    int NintBin = controls.analysis.NiBins();
    //   Number of bins, number of "reference" bin,
    //   intensity at "reference" bin, maximum intensity
    float Iav = NormRes.Imean();
    float Jmax = NormRes.Imax();
    IntensityBin Irange(NintBin, NintBin/2, Iav, Jmax);

    std::vector<int> countsa(NintBin, 0);  // N accepted
    std::vector<int> countsr(NintBin, 0);  // N rejected

    // Omit ice rings
    Rings icerings;
    icerings.DefaultIceRings();

    hkl_list.rewind();
    reflection this_refl;
    int nvrefl = 0;  // count number of valid reflections with at least 2 observations
    while (hkl_list.next_reflection(this_refl) >= 0)  {
      if ((icerings.InRing(this_refl.invresolsq()) < 0) &&
        (this_refl.NvalidObservations() > 1)) {  // no singletons or in icering
        nvrefl++;
      }
    }

    // From Wilson statistics p(E^2 > Q) = exp(-Q)
    double E2min = -log(float(Nbintarget)/float(nvrefl));
    // E2min shouldn't be too large
    //    E2min = Min(E2min, 1.5);
    E2min = Min(E2min, 0.5);
    int nreflarge = Nint(float(nvrefl)* exp(-E2min)); // rough number > E2min
    nreflarge = nreflarge;
    // bin at E^2 = 1 (Emidbin) is NintBin/2
    double Emidbin = 1.0;
    // rough number / bin
    int nrefbin = (float(nvrefl)*(1.0 - exp(-Emidbin))/(0.5*float(NintBin)));
    double p0 = 1.0 - exp(-E2min);  // p = 1 at E2min
    double frac = Min(1.0, float(Nbintarget)/float(nrefbin)); // average fraction to accept
    //^
    //    std::cout << "nvrefl, Nbintarget, nreflarge, nrefbin "
    //        << nvrefl<<" "<< Nbintarget<<" " << nreflarge << " " << nrefbin << " E2min " << E2min <<"\n";
    //    std::cout << "Iav, Jmax " << Iav <<" "<<Jmax <<"\n";
    //    std::cout << "Frac "<< frac  << "\n";
    //^-


    hkl_list.rewind();
    int refAccepted = 0;
    int mint;
    int raccept = 0; // accept
    int rindex;
    int inring = 0;
    int rejf = 0;    // rejected on frac

    // * * * * Loop all reflections
    while ((rindex = hkl_list.next_reflection(this_refl)) >= 0)  {
      if (icerings.InRing(this_refl.invresolsq()) >= 0) {
        // reject
        raccept = +1;
        inring++;
        mint = 0;
      } else if (this_refl.NvalidObservations() > 1) {  // no singletons
        // Select all (I+ & I-) accepted observations for all datasets
        SelectedObservations allobs(this_refl, -1, ALL);
        IsigI AvIsig = allobs.Average();  // average I, 1/variance weight
        raccept = 0; // accept
        double E2 = NormRes.applyAvg(AvIsig.I(), this_refl.invresolsq());
        if (E2 < E2min) {
          // Accept a fraction frac of these scaled by p(E^2)
          double p = exp(-E2); // p(true(E2) > E2)
          double acc = 1.0;
          if (p > 0.0) {acc = (1-frac)*(1-p)/p0 + frac;}
          if (p >= 1.0 || FRandom(1.0) > acc) { // reject all negative <I>
            // reject
            raccept = +1;
            rejf++;
            //^
            //      std::cout << frac <<" "<<AvIsig.I()<<" "<<E2<<" "<< p<<" "<<acc<<"  **\n"; //^-
          }
        }
        mint = Irange.bin(AvIsig.I());
      } else { // singletons
        mint = 0;
        raccept = +1;
      }
      if (raccept != 0) { // not accepted
        this_refl.SetStatus(raccept);
        hkl_list.replace_reflection(this_refl); // store updated reflection
        countsr.at(mint)++;
      } else {
        // Count reflections accepted
        countsa.at(mint)++;
        refAccepted++;
      }
    }  // end loop reflections

    //^
    //    std::cout << "\n *** Number rejected in ice rings " << inring <<"\n";
    //
    //    std::cout << "\n *** Number rejected on frac " << rejf <<"\n";
    //    //^
    //    for (int i=0;i<countsa.size();++i) {
    //      std::cout <<"Count  Ibin acc rej " << i
    //          << " " << countsa[i]<< " " << countsr[i] <<"\n";
    //    }
    //^-
    return refAccepted;
  }
}
