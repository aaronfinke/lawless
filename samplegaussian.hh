class SampleGaussian
{
  // Class to generate random numbers (float) with Gaussian distribution
  // Box-Muller algorithm from Numerical Recipes
public:
  SampleGaussian() : gotone(false) {}

  // Get number with zero mean & unit variance
  float Get();
  float Get(const float& Mean, const float& SD);

private:
  float previous; // algorithm generates numbers in pairs, 
  bool gotone;    //  so buffer one of them
};
