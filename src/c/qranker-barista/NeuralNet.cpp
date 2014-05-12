#include "NeuralNet.h"
#include "utils.h"

/*****************Sigmoid**********************/

void Sigmoid :: resize(int m)
{
  num_neurons = m;
}
 

Sigmoid& Sigmoid :: operator=(Sigmoid &S)
{
  num_neurons = S.num_neurons;
  return *this;
}

void Sigmoid :: clone(Sigmoid &S)
{
  num_neurons = S.num_neurons;
}

void Sigmoid :: fprop(State& down, State &up)
{
  for(int k = 0; k < num_neurons; k++)
      up.x[k] = 1.0/(1.0+exp(-down.x[k]));
}
 
void Sigmoid :: bprop(State &down, State &up)
{
  for(int k = 0; k < num_neurons; k++)
    down.dx[k] = up.dx[k]*(1.0-up.x[k])*up.x[k];
}


/******************** Linear*************************/
void Linear :: make_random()
{
  if(num_neurons > 0 && num_features > 0)
    {
      for(int k = 0; k < num_neurons; k++)
	{
	  for(int j = 0; j < num_features; j++)
	    w[k*num_features+j] = ((double)myrandom()/UNIFORM_INT_DISTRIBUTION_MAX - 0.5)/(num_features*num_neurons);
	  if(has_bias)
	    bias[k] = ((double)myrandom()/UNIFORM_INT_DISTRIBUTION_MAX - 0.5)/(num_features*num_neurons);;
	}
    }
}

void Linear :: resize(int m, int n, int has_b)
{
  if(num_neurons != m || num_features != n)
    {
      clear();
      num_neurons = m;
      num_features = n;
      has_bias = has_b;
      w = new double[num_neurons*num_features];
      bias = new double[num_neurons];
      dw = new double[num_neurons*num_features];
      dbias = new double[num_neurons];
      num_refs = new int[1];
      num_refs[0] = 1;
  }
  for(int k = 0; k < num_neurons; k++)
    {
      for(int j = 0; j < num_features; j++)
	w[k*num_features+j] = ((double)myrandom()/UNIFORM_INT_DISTRIBUTION_MAX - 0.5)/(num_features*num_neurons);
      if(has_bias)
	bias[k] = ((double)myrandom()/UNIFORM_INT_DISTRIBUTION_MAX - 0.5)/(num_features*num_neurons);;
    }
  memset(dw,0,sizeof(double)*num_neurons*num_features);
  if(has_bias)
    memset(dbias, 0, sizeof(double)*num_neurons);
  
}


void Linear :: clear()
{
  if(num_refs)
    {
      num_refs[0]--;
      if(num_refs[0] < 1)
	{
	  delete[] w; w = 0;
	  delete[] bias; bias = 0;
	  delete[] dw; dw = 0;
	  delete[] dbias; dbias = 0;
	  delete[] num_refs; num_refs = 0;
	  num_neurons = 0;
	  num_features = 0;
	}
    }
}

Linear& Linear :: operator=(Linear &L)
{
  resize(L.num_neurons, L.num_features, L.has_bias);
  for(int k = 0; k < num_neurons; k++)
    {
      for(int j = 0; j < num_features; j++)
	w[k*num_features+j] = L.w[k*num_features+j];
      if(has_bias)
	bias[k] = L.bias[k];
    }
  return *this;
}

void Linear :: copy(Linear &L)
{
  assert(L.num_neurons == num_neurons);
  assert(L.num_features == num_features);
  for(int k = 0; k < num_neurons; k++)
    {
      for(int j = 0; j < num_features; j++)
	w[k*num_features+j] = L.w[k*num_features+j];
      if(has_bias)
	bias[k] = L.bias[k];
    }
}


void Linear :: clone(Linear &L)
{
  clear();
  num_neurons = L.num_neurons;
  num_features = L.num_features;
  has_bias = L.has_bias;
  w = L.w;
  dw = L.dw;
  bias = L.bias;
  dbias = L.dbias;
  num_refs = L.num_refs;
  num_refs[0]++;
}


void Linear :: clear_gradients()
{
  memset(dw,0,sizeof(double)*num_neurons*num_features);
  if(has_bias)
    memset(dbias,0,sizeof(double)*num_neurons);
}

void Linear :: fprop(State &down, State &up)
{
  for(int k = 0; k < num_neurons; k++)
    {
      double d = 0.0;
      for(int j = 0; j < num_features; j++)
	d += w[k*num_features+j]*down.x[j];
      //if there is a bias
      if(has_bias)
	d += bias[k];
      up.x[k] = d;
    }
}


void Linear :: bprop(State &down, State &up)
{
  memset(down.dx,0,sizeof(double)*num_features);
  for(int j = 0; j < num_features; j++)
    for(int k = 0; k < num_neurons; k++)
      down.dx[j] += up.dx[k]*w[k*num_features+j];
 
  for(int k = 0; k < num_neurons; k++)
    {
      for(int j = 0; j < num_features; j++)
	dw[k*num_features+j] += up.dx[k]*down.x[j];
      //if there is a bias
      if(has_bias)
	dbias[k] += up.dx[k];
    }
}

void Linear :: update(double mu, double weight_decay)
{
  for(int k = 0; k < num_neurons; k++)
    {
      for(int j = 0; j < num_features; j++)
	{
	  if(weight_decay > 0.0)
	    dw[k*num_features+j] += weight_decay*w[k*num_features+j];
	  w[k*num_features+j] -= mu*dw[k*num_features+j];
	}
      //if there is a bias
      if(has_bias)
	{
	  if(weight_decay > 0.0)
            dbias[k] += weight_decay*bias[k];
	  bias[k] -= mu*dbias[k];
	}
    }
}

void Linear :: update1(double mu, double weight_decay)
{
  for(int k = 0; k < num_neurons; k++)
    {
      for(int j = 0; j < num_features; j++)
	{
	  if(weight_decay > 0.0)
	    dw[k*num_features+j] += weight_decay*w[k*num_features+j];
	  w[k*num_features+j] -= (mu/(num_features*num_neurons))*dw[k*num_features+j];
	}
      //if there is a bias
      if(has_bias)
	{
	   if(weight_decay > 0.0)
            dbias[k] += weight_decay*bias[k];
	   bias[k] -= (mu/(num_features*num_neurons))*dbias[k];
	}
    }
}

/**************** Neural Net ********************/
void NeuralNet :: make_random()
{
  lin1.make_random();
  if(!is_linear)
    lin2.make_random();
}

void NeuralNet :: resize_states()
{
  int num_features = lin1.get_num_features();
  int num_neurons1 = lin1.get_num_neurons();
    
  //first state
  start.resize_dx(num_features);
  if(!is_linear)
    {
      assert(num_neurons1 == sigm1.get_num_neurons());
      s1.resize(num_neurons1);
      assert(num_neurons1 == lin2.get_num_features());
      s2.resize(num_neurons1);
    }
    //last state
    finish.resize_x(1);
}

void NeuralNet :: initialize(int nfeatures, int num_hu, int is_lin, int has_bias)
{
  is_linear = is_lin;
  //cout << "initializing net " << nfeatures << " " << num_hu << endl;
  lin1.resize(num_hu, nfeatures, has_bias);
  if(!is_linear)
    {
      sigm1.resize(num_hu);
      lin2.resize(1,num_hu, has_bias);
    }
  resize_states();
}

NeuralNet& NeuralNet :: operator=(NeuralNet &N)
{
  is_linear = N.is_linear;
  lin1 = N.lin1;
  if(!is_linear)
    {
      sigm1 = N.sigm1;
      lin2 = N.lin2;
    }
  resize_states();
  return *this;
}

void NeuralNet :: copy(NeuralNet &N)
{
  assert(is_linear == N.is_linear);
  lin1.copy(N.lin1);
  if(!is_linear)
    {
      lin2.copy(N.lin2);
    }
}


void NeuralNet :: clone(NeuralNet &N)
{
  is_linear = N.is_linear;
  lin1.clone(N.lin1);
  if(!is_linear)
    {
      sigm1.clone(N.sigm1);
      lin2.clone(N.lin2);
    }
  resize_states();
}


double* NeuralNet :: fprop(double *x)
{
  start.x = x;
  if(is_linear)
    lin1.fprop(start,finish);
  else
    {
      lin1.fprop(start,s1);
      sigm1.fprop(s1,s2);
      lin2.fprop(s2,finish);
    }
  return finish.x;

}


void NeuralNet :: clear_gradients()
{
  lin1.clear_gradients();
  if(!is_linear)
    lin2.clear_gradients();
}


double* NeuralNet :: bprop(double *dx)
{
  finish.dx = dx;
  if(is_linear)
    lin1.bprop(start,finish);
  else
    {
      lin2.bprop(s2,finish);
      sigm1.bprop(s1,s2);
      lin1.bprop(start,s1);
    }
    
  return start.dx;
}

void NeuralNet :: update(double mu, double weight_decay)
{
  lin1.update(mu, weight_decay);
  if(!is_linear)
    lin2.update(mu, weight_decay);
}

void NeuralNet :: update1(double mu, double weight_decay)
{
  lin1.update1(mu, weight_decay);
  if(!is_linear)
    lin2.update1(mu, weight_decay);
}


double* NeuralNet :: get_weights(int mod)
{
  if(mod == 1)
    return lin1.get_weights();
  else if(mod == 2)
    {
      if(is_linear)
	{
	  cout << "neural net was initialized with only 1 layer, layer " << mod << " does not exits" << endl;
	  return 0;
	}
      return lin2.get_weights();
    }
  else
    {
      if(is_linear)
	cout << "neural net was initialized with only 1 layer, ";
      else
	cout << "neural net was initialized with 2 layers, ";
      cout << "layer " << mod << "does not exist" << endl;
      return 0;
    }

}

double* NeuralNet :: get_dweights(int mod)
{
  if(mod == 1)
    return lin1.get_dweights();
  else if(mod == 2)
    {
      if(is_linear)
	{
	  cout << "neural net was initialized with only 1 layer, layer " << mod << " does not exits" << endl;
	  return 0;
	}
      return lin2.get_dweights();
    }
  else
    {
      if(is_linear)
	cout << "neural net was initialized with only 1 layer, ";
      else
	cout << "neural net was initialized with 2 layers, ";
      cout << "layer " << mod << "does not exist" << endl;
      return 0;
    }
}

double* NeuralNet :: get_bias(int mod)
{
  if(mod == 1)
    return lin1.get_bias();
  else if(mod == 2)
    {
      if(is_linear)
	{
	  cout << "neural net was initialized with only 1 layer, layer " << mod << " does not exits" << endl;
	  return 0;
	}
    return lin2.get_bias();
    }
  else
    {
      if(is_linear)
	cout << "neural net was initialized with only 1 layer, ";
      else
	cout << "neural net was initialized with 2 layers, ";
      cout << "layer " << mod << "does not exist" << endl;
      return 0;
    }
}

double* NeuralNet :: get_dbias(int mod)
{
  if(mod == 1)
    return lin1.get_dbias();
  if(mod == 2)
    {
      if(is_linear)
	{
	  cout << "neural net was initialized with only 1 layer, layer " << mod << " does not exits" << endl;
	  return 0;
	}
      return lin2.get_dbias();
    }
  else
    {
      if(is_linear)
	cout << "neural net was initialized with only 1 layer, ";
      else
	cout << "neural net was initialized with 2 layers, ";
      cout << "layer " << mod << "does not exist" << endl;
      return 0;
    }
}

/********************************************************************/

double check_gradients(NeuralNet &n, double *x)
{
  double sum = 0.0;
  
  double h = 10e-6;
    
  double *c = n.fprop(x);
  double cbefore = c[0];

  double *gc = new double[1];
  gc[0] = 1;
  n.clear_gradients();
  n.bprop(gc);

  double *w = n.get_weights(1);
  double *dw = n.get_dweights(1);
  for (int j = 0; j < n.get_num_hu(); j++)
    for (int i = 0; i < n.get_num_features(); i++)
      {
	int ii = j*n.get_num_features()+i;
	double ww = w[ii];
	w[ii] += h;
	c = n.fprop(x);
	double d = c[0]-cbefore;
	sum += fabs(dw[ii]-d/h);
	w[ii]=ww;
    }
    
  double *bias = n.get_bias(1);
  double *dbias = n.get_dbias(1);
  for (int i = 0; i < n.get_num_hu();i++)
    {
      double bb = bias[i];
      bias[i] +=h;
      c = n.fprop(x);
      double d = c[0]-cbefore;
      sum += fabs(dbias[i]-d/h);
      bias[i] =bb;
    }
    
  if(n.get_num_layers() > 1)
    {
      double *w = n.get_weights(2);
      double *dw = n.get_dweights(2);
      for (int i = 0; i < n.get_num_hu(); i++)
	{
	  double ww = w[i];
	  w[i] += h;
	  c = n.fprop(x);
	  double d = c[0]-cbefore;
	  sum += fabs(dw[i]-d/h);
	  w[i]=ww;
	}

      double *bias = n.get_bias(2);
      double *dbias = n.get_dbias(2);
      for (int i = 0; i < 1;i++)
	{
	  double bb = bias[i];
	  bias[i] +=h;
	  c = n.fprop(x);
	  double d = c[0]-cbefore;
	  sum += fabs(dbias[i]-d/h);
	  bias[i] =bb;
	}
    }
  
  n.update(0.005);
  return sum;
}

double check_gradients_hinge(NeuralNet &n, double *x, int y)
{
  double sum = 0.0;
  
  double h = 10e-6;
    
  double *c = n.fprop(x);
  double cbefore = 1.0-c[0]*y;

  double *gc = new double[1];
  gc[0] = -1*y;
  n.clear_gradients();
  n.bprop(gc);

  double *w = n.get_weights(1);
  double *dw = n.get_dweights(1);
  for (int j = 0; j < n.get_num_hu(); j++)
    for (int i = 0; i < n.get_num_features(); i++)
      {
	int ii = j*n.get_num_features()+i;
	double ww = w[ii];
	w[ii] += h;
	c = n.fprop(x);
	double d = (1.0-c[0]*y)-cbefore;
	sum += fabs(dw[ii]-d/h);
	w[ii]=ww;
    }

  
  double *bias = n.get_bias(1);
  double *dbias = n.get_dbias(1);
  for (int i = 0; i < n.get_num_hu();i++)
    {
      double bb = bias[i];
      bias[i] +=h;
      c = n.fprop(x);
      double d = (1.0-c[0]*y)-cbefore;
      sum += fabs(dbias[i]-d/h);
      bias[i] =bb;
    }
    
  if(n.get_num_layers() > 1)
    {
      double *w = n.get_weights(2);
      double *dw = n.get_dweights(2);
      for (int i = 0; i < n.get_num_hu(); i++)
	{
	  double ww = w[i];
	  w[i] += h;
	  c = n.fprop(x);
	  double d = (1.0-c[0]*y)-cbefore;
	  sum += fabs(dw[i]-d/h);
	  w[i]=ww;
	}

      double *bias = n.get_bias(2);
      double *dbias = n.get_dbias(2);
      for (int i = 0; i < 1;i++)
	{
	  double bb = bias[i];
	  bias[i] +=h;
	  c = n.fprop(x);
	  double d = (1.0-c[0]*y)-cbefore;
	  sum += fabs(dbias[i]-d/h);
	  bias[i] =bb;
	}
    }
  
  //n.update(0.005);
  return sum;
}


double check_gradients_sigmoid(NeuralNet &n, double *x, int y)
{
  double sum = 0.0;
  
  double h = 10e-6;
    
  double *c = n.fprop(x);
  double a = exp(y*c[0]);
  double cbefore = (1/(1+a));

  double *gc = new double[1];
  gc[0] = -a/((1+a)*(1+a))*y;
  n.clear_gradients();
  n.bprop(gc);

  double *w = n.get_weights(1);
  double *dw = n.get_dweights(1);
  for (int j = 0; j < n.get_num_hu(); j++)
    for (int i = 0; i < n.get_num_features(); i++)
      {
	int ii = j*n.get_num_features()+i;
	double ww = w[ii];
	w[ii] += h;
	c = n.fprop(x);
	double aa = exp(y*c[0]);
	double cafter = (1/(1+aa));
	double d = cafter-cbefore;
	sum += fabs(dw[ii]-d/h);
	w[ii]=ww;
    }

  
  double *bias = n.get_bias(1);
  double *dbias = n.get_dbias(1);
  for (int i = 0; i < n.get_num_hu();i++)
    {
      double bb = bias[i];
      bias[i] +=h;
      c = n.fprop(x);
      double aa = exp(y*c[0]);
      double cafter = (1/(1+aa));
      double d = cafter-cbefore;
      sum += fabs(dbias[i]-d/h);
      bias[i] =bb;
    }
    
  if(n.get_num_layers() > 1)
    {
      double *w = n.get_weights(2);
      double *dw = n.get_dweights(2);
      for (int i = 0; i < n.get_num_hu(); i++)
	{
	  double ww = w[i];
	  w[i] += h;
	  c = n.fprop(x);
	  double aa = exp(y*c[0]);
	  double cafter = (1/(1+aa));
	  double d = cafter-cbefore;
	  sum += fabs(dw[i]-d/h);
	  w[i]=ww;
	}

      double *bias = n.get_bias(2);
      double *dbias = n.get_dbias(2);
      for (int i = 0; i < 1;i++)
	{
	  double bb = bias[i];
	  bias[i] +=h;
	  c = n.fprop(x);
	  double aa = exp(y*c[0]);
	  double cafter = (1/(1+aa));
	  double d = cafter-cbefore;
	  sum += fabs(dbias[i]-d/h);
	  bias[i] =bb;
	}
    }
  
  //n.update();
  return sum;
}


/*
void check_gradients_rank(NeuralNet &n, NeuralNet *nets, double *x1, int y1, double *x2, int y2)
{

  double diff = 0;
  double h = 0.0001;
    
  double *c1 = nets[0].fprop(x1);
  double *c2 = nets[1].fprop(x2);
 
  double d = c1[0]-c2[0]; 
  diff -= d;

  double *gc = new double[1];
  
  n.clear_gradients();
  gc[0] = 1.0;
  nets[0].bprop(gc);
  gc[0] = -1.0;
  nets[1].bprop(gc);

  double *w = n.get_weights(1);
  double *dw = n.get_dweights(1);
  for (int j = 0; j < n.get_num_hu(); j++)
    {
      for (int i = 0; i < n.get_num_features(); i++)
	{
	  int ii = j*n.get_num_features()+i;
	  double ww = w[ii];
	  w[ii] += h;
	  c1 = nets[0].fprop(x1);
	  c2 = nets[1].fprop(x2);
	  d = c1[0]-c2[0];
	  
	  diff += d;
	  cout << dw[i] << " " << diff/h << " " << dw[i]-diff/h << endl;
	  w[ii]=w;
	  diff -= d;
	}
    }
    
}
*/
/******************************************************/
/*
int main()
{
  int num_features = 17;
  int num_hu = 5;
  int has_bias = 1;
  int is_lin = 0;
  
  NeuralNet n;
  n.initialize(num_features, num_hu, is_lin, has_bias);

  double sm = 0;
  int nexamples = 10000;
  
  for( int k = 0; k < nexamples; k++)
    {
      double *x = new double[num_features];
      for (int i = 0; i < num_features; i++)
	{
	  x[i] = (double)rand()/(double)RAND_MAX;
	}
      int label = 1;
      sm += check_gradients_hinge(n, x,label);
    }
  cout << sm/(double)nexamples << endl;

  sm = 0;
  nexamples = 100;
  for( int k = 0; k < nexamples; k++)
    {
      double *x = new double[num_features];
      for (int i = 0; i < num_features; i++)
	{
	  x[i] = (double)rand()/(double)RAND_MAX;
	}
      int label = 1;
      sm += check_gradients_sigmoid(n, x,label);
    }
  cout << sm/(double)nexamples << endl;

*/

  /*
  NeuralNet *nets = new NeuralNet[2];
  
  nets[0].clone(n);
  nets[1].clone(n);
  
  double *x1 = new double[num_features];
  for (int i = 0; i < num_features; i++)
    {
      x1[i] = (double)rand()/(double)RAND_MAX;
    }
  int label1 = 1;
  
  double *x2 = new double[num_features];
  for (int i = 0; i < num_features; i++)
    {
      x2[i] = (double)rand()/(double)RAND_MAX;
    }
  int label2 = -1;

  check_gradients_rank(n ,nets, x1, label1, x2, label2);
  
  x1 = new double[num_features];
  for (int i = 0; i < num_features; i++)
    {
      x1[i] = (double)rand()/(double)RAND_MAX;
    }
  label1 = 1;
  
  
  for (int i = 0; i < num_features; i++)
    {
      x2[i] = (double)rand()/(double)RAND_MAX;
    }
  label2 = -1;

  check_gradients_rank(n, nets, x1, label1, x2, label2);
  */




  //double *out = n.fprop(x);
  //cout << out[0] << endl;
  //double *gc = new double(1);
  //gc[0] = -1*(double)label;
  //n.bprop(gc);

  //}

