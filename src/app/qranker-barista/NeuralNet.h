#ifndef NEURAL_NET_H
#define NEURAL_NET_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <assert.h>
#include <cstring>
#include <math.h>
#include <vector>
using namespace std;


class State
{
 public:
 State() 
   : len_x(0), 
   len_dx(0),
   x((double*)0),  
   dx((double*)0)  {}
  inline void resize(int n) {clear(); x = new double[n]; len_x=n; dx = new double[n]; len_dx = n;}
  inline void resize_dx(int n) {clear(); dx = new double[n]; len_dx = n;}
  inline void resize_x(int n) {clear(); x = new double[n]; len_x = n;}
  inline void clear(){if(len_x) delete[] x; if(len_dx) delete[] dx; len_x = 0; len_dx=0;}
  ~State(){clear();}
 
  int len_x;
  int len_dx;
  double *x;
  double *dx;
 
};

class Sigmoid
{
public:
  Sigmoid():num_neurons(0){}
  ~Sigmoid() {}
  void clear(){num_neurons = 0;}
  void resize(int m);
  Sigmoid& operator=(Sigmoid &S);
  void clone(Sigmoid &S);

  inline int get_num_neurons(){return num_neurons;}
 
  void write_to_file(ofstream &outfile);
  void read_from_file(ifstream &infile);
 
  void fprop(State& down, State &up);
  void bprop(State &down, State &up);
   
 protected:
  int num_neurons;
};



class Linear
{
 public:
 Linear(): 
   num_neurons(0),
   num_features(0),
   has_bias(1),
   num_refs(0),
   w(0),
   bias(0),
   dw(0),
   dbias(0) {}
  ~Linear() {clear();}
  void init(int m, int n, int has_b);
  void resize(int m, int n, int has_b);
  void clear();
  Linear& operator=(Linear &L);
  void copy(Linear &L);
  void clone(Linear &L);
  void make_random();
  
  inline void set_num_features(int nf) {num_features = nf;}
  inline int get_num_features() const {return num_features;}
  inline void set_num_neurons(int nr) {num_neurons = nr;}
  inline int get_num_neurons() const{return num_neurons;}
   
  inline double* get_weights() {return w;}
  inline double* get_dweights() {return dw;}
  inline double* get_bias() {return bias;}
  inline double* get_dbias() {return dbias;}

  void write_to_file(ofstream &outfile);
  void read_from_file(ifstream &infile);
 
  void fprop(State &down, State &up);
  void bprop(State &down, State &up);
  void clear_gradients();
  void update(double mu, double weight_decay=0.0);
  void update1(double mu, double weight_decay = 0.0);
  
 protected:
  int num_neurons;
  int num_features;
  int has_bias;

  int *num_refs;
  double *w;
  double *bias;
  double *dw;
  double *dbias;
  
};

/************************NeuralNet******************************************/
class NeuralNet {
 public:
 NeuralNet():is_linear(0){}
  ~NeuralNet(){}
  void clear();
  void resize_states();
  void initialize(int nfeatures, int num_hu, int is_lin, int has_bias);
  NeuralNet& operator=(NeuralNet &N);
  void clone(NeuralNet &N);
  void copy(NeuralNet &N);
  void make_random();

  double* fprop(double *down);
  void clear_gradients();
  double* bprop(double *up);
  void update(double mu, double weight_decay=0.0);
  void update1(double mu, double weight_decay=0.0);

  double* get_weights(int mod);
  double* get_dweights(int mod);
  double* get_bias(int mod);
  double* get_dbias(int mod);
  inline int get_num_features(){return lin1.get_num_features();}
  inline int get_num_hu(){return lin1.get_num_neurons();}
  inline int get_num_layers(){if (is_linear) return 1; else return 2;}

 protected:
  int is_linear;
  
  State start;
  Linear lin1;
  State s1;
  Sigmoid sigm1;
  State s2;
  Linear lin2;
  State finish;
};


#endif /*NEURAL_NET_H*/
