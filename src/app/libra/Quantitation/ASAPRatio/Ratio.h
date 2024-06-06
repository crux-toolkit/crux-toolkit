#ifndef RATIO_H
#define RATIO_H


class Ratio {

 public:

 Ratio(double mean, double std, int numpeps) : mean_(mean), standard_dev_(std), number_peptides_(numpeps) { }

  double mean_;
  double standard_dev_;
  int number_peptides_;

};


#endif
