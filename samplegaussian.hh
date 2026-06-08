class SampleGaussian
{
  // Class to generate random numbers (double) with Gaussian distribution
  // Box-Muller algorithm from Numerical Recipes
public:
  SampleGaussian() : gotone(false) {}

  // Get number with zero mean & unit variance
  double Get();
  double Get(const double& Mean, const double& SD);

private:
  double previous; // algorithm generates numbers in pairs, 
  bool gotone;    //  so buffer one of them
};
