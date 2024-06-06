#ifndef SPEC_PK
#define SPEC_PK


class SpectrumPeak {
 public:
  SpectrumPeak(double mass, double intens) {
    mass_ = mass;
    intens_ = intens;
  }

  double mass_;
  double intens_;


};


#endif
