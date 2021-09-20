//#include </usr/include/python2.7/Python.h>
//#include </usr/include/python3.5/Python.h>
//#include </usr/include/python3.7m/Python.h>
//#include </usr/local/include/Eigen/Dense>
#define _USE_MATH_DEFINES
#include <Python.h>
#include "Eigen/Dense"

#include <cmath>



using Eigen::MatrixXd;
using namespace Eigen;
using namespace std;

//Declare paraArrayFunc
double paraX(double originIn[2], double defIn[2], double r, double T, double alph, double th)
{
  double Pdx, Pdy, a, Vd;
  Vd = 1.0;
  a = alph;
  Pdx = defIn[0] - originIn[0];
  Pdy = defIn[1] - originIn[1];
  double X = -(cos(th)*(a*r + Pdx*cos(th) + Pdy*sin(th) - sqrt(pow(T,2)*pow(Vd,2) - pow(Pdx,2)*pow(sin(th),2) - pow(Pdy,2)*pow(cos(th),2) + pow(Pdx,2)*pow(a,2) + pow(Pdy,2)*pow(a,2) + pow(r,2) - 2*T*1*r + 2*Pdx*a*r*cos(th) + 2*Pdy*a*r*sin(th) + 2*Pdx*Pdy*cos(th)*sin(th) - 2*Pdx*T*a*cos(th) - 2*Pdy*T*a*sin(th)) - T*Vd*a))/(pow(a,2) - 1);
  return X + originIn[0];

}

double paraY(double originIn[2], double defIn[2], double r, double T, double alph, double th)
{
  double Pdx, Pdy, a, Vd;
  Vd = 1.0;
  a = alph;
  Pdx = defIn[0] - originIn[0];
  Pdy = defIn[1] - originIn[1];
  double Y = -(sin(th)*(a*r + Pdx*cos(th) + Pdy*sin(th) - sqrt(pow(T,2)*pow(Vd,2) - pow(Pdx,2)*pow(sin(th),2) - pow(Pdy,2)*pow(cos(th),2) + pow(Pdx,2)*pow(a,2) + pow(Pdy,2)*pow(a,2) + pow(r,2) - 2*T*1*r + 2*Pdx*a*r*cos(th) + 2*Pdy*a*r*sin(th) + 2*Pdx*Pdy*cos(th)*sin(th) - 2*Pdx*T*a*cos(th) - 2*Pdy*T*a*sin(th)) - T*Vd*a))/(pow(a,2) - 1);
  return Y+originIn[1];


}


double * para_path3(double originIn[2], double defIn[2], double r, double T, double alph)
{
  const int size = 200;

  Vector2d origin, def;
  origin << originIn[0], originIn[1];
  def << defIn[0], defIn[1];
  double pi = 3.141592653589793238462643383279;

  ArrayXXd  out(size,2);

  VectorXd  angles;
  VectorXd  paraOut(2);
  static double doubArrayOut[size];
  angles.setLinSpaced(size/2,-pi,pi);
  //out.setZero((int)size, 2);
  out.col(0) = angles.array().cos();
  out.col(1) = angles.array().sin();


  for(int i = 0; i < size / 2; i++){

    //doubArrayOut[i] = (double)cos(angles(i));
    doubArrayOut[i] = paraX(originIn, defIn, r, T, alph, angles(i));

  }

  int q = 0;
  for(int i = size / 2; i < size; i++){

    doubArrayOut[i] = paraY(originIn, defIn, r, T, alph, angles(q));
    q++;

  }

  return doubArrayOut;
};

////////////////  NEW FUNCTIONS

//Define Convex Function
double convex_function(double originIn[2], double defIn[2], double r, double T, double alph, double th, double G_Circ[3])
{

  double Xp = paraX(originIn, defIn, r, T, alph, th);
  double Yp = paraY(originIn, defIn, r, T, alph, th);

  double Xg = G_Circ[0];
  double Yg = G_Circ[1];
  double Rad = G_Circ[2];

  return pow((Xp - Xg),2) + pow((Yp - Yg),2) - pow(Rad,2);
}

//checks if invader passes through target region
int CisInRegion(double Invader[2], double P_star[2], double G_Circ[3]){

  double Xp, Yp;
  double v[2];
  double Xg = G_Circ[0];
  double Yg = G_Circ[1];
  double Rad = G_Circ[2];
  double G;

  v[0] = P_star[0] - Invader[0];
  v[1] = P_star[1] - Invader[1];

  for(double t = 0; t <= 1; t+=0.05){
    Xp = Invader[0] + t*v[0];
    Yp = Invader[1] + t*v[1];
    G = pow((Xp - Xg),2) + pow((Yp - Yg),2) - pow(Rad,2);
    if(G < 0){

      return 1;

    }


  }

  return 0;
}

//////  Minimum Index of Array

int min_Idx(double originIn[2], double defIn[2], double r, double T, double alph, double G_Circ[3], int size)
{

  VectorXd  angles;
  angles.setLinSpaced(size,-M_PI,M_PI);
  int minIdx = 0;
  double minValue = convex_function(originIn, defIn, r, T, alph, angles(0), G_Circ);
  double current;

  for(int i = 1; i < size; i++){

    current = convex_function(originIn, defIn, r, T, alph, angles(i), G_Circ);
    if (current < minValue){
      minIdx = i;
      minValue = current;
    }
  }

  return minIdx;
}

/// Numerical Gradient
double numGrad_convex(double originIn[2], double defIn[2], double r, double T, double alph, double G_Circ[3], double init)
{
  double h = 0.0001;
  double init_ang = init;
  double f1 = convex_function(originIn, defIn, r, T, alph, init_ang+(h/2), G_Circ);
  double f2 = convex_function(originIn, defIn, r, T, alph, init_ang-(h/2), G_Circ);

  return (f1 - f2) / (h);
}

/// minimization function IMPLIMENT THIS

double CminiConvex(double originIn[2], double defIn[2], double r, double T, double alph, double G_Circ[3], int size)
{

  VectorXd  angles;
  angles.setLinSpaced(size,-M_PI,M_PI);
  int courseMin_id = min_Idx(originIn, defIn, r, T, alph, G_Circ, size);
  double init = angles(courseMin_id);
  double minAng = init;
  double prevAng = 5;
  double lr = 0.001;
  double grad;

  for(int i = 0; i < 1500; i++){

    grad = numGrad_convex(originIn, defIn, r, T, alph, G_Circ, minAng);
    minAng = minAng - lr*grad;

    if (abs(prevAng - minAng) <= 0.00001){
      //std::cout << "Hellow" << std::endl;
      //std::cout << "Hellow: ";
      //std::cout << i << std::endl;
      return minAng;

    }

    prevAng = minAng;
  }



  return minAng;


}



/////////////////////////////  PYTHON STUFF


PyObject *Convert_Big_Array(double array[], int length)
  { PyObject *pylist, *item;
    int i;
    pylist = PyList_New(length);
    for (i=0; i<length; i++) {
      //item = PyInt_FromLong(array[i]);
      item = PyFloat_FromDouble(array[i]);
      PyList_SetItem(pylist, i, item);
    }
    return pylist;
  }

static PyObject* miniConvex(PyObject* self, PyObject* args)
  {
    int size = 100;
    double r,T, alph;
    double originIn[2];
    double defIn[2];
    double G_Circ[3];

    if(!PyArg_ParseTuple(args,"(dd)(dd)ddd(ddd)",&originIn[0], &originIn[1], &defIn[0], &defIn[1] , &r, &T, &alph, &G_Circ[0], &G_Circ[1], &G_Circ[2]))
      return NULL;

    double ansReturn = CminiConvex(originIn, defIn, r, T, alph, G_Circ, size);
    //double ansReturn = 1;
    //return Py_BuildValue("d", ansReturn );
    return PyFloat_FromDouble(ansReturn);

  }

static PyObject* isInRegion(PyObject* self, PyObject* args)
  {
    double Invader[2];
    double P_star[2];
    double G_Circ[3];

    if(!PyArg_ParseTuple(args,"(dd)(dd)(ddd)",&Invader[0], &Invader[1], &P_star[0], &P_star[1], &G_Circ[0], &G_Circ[1], &G_Circ[2]))
      return NULL;

   
    //double ansReturn = 1;
    //return Py_BuildValue("d", ansReturn );
    return Py_BuildValue("i", CisInRegion(Invader, P_star, G_Circ));

  }



// Para Function Python Object
static PyObject* paraPathFunc(PyObject* self, PyObject* args)
{

  double r,T, alph;
  double originIn[2];
  double defIn[2];
  double arr3[4];
  arr3[0] = 0;
  arr3[1] = 1;
  arr3[2] = 0;
  arr3[3] = 0;



  //create [item]
  std::string s = "[";
  char ch = 'd';
  for(int i = 0; i < 2*2; i++)
  {
    s.push_back(ch);
  }
  char chend = ']';
  s.push_back(chend);
  const char * c = s.c_str();

  if(!PyArg_ParseTuple(args,"(dd)(dd)ddd",&originIn[0], &originIn[1], &defIn[0], &defIn[1] , &r, &T, &alph))
    return NULL;
  double *dubout;
  dubout = para_path3(originIn, defIn, r, T, alph);


  return Convert_Big_Array(dubout, 200);
  //return Py_BuildValue("(dd)", (arr3[0],arr3[1]));

}

double CdirectionalCostFunction(double P0[2], double P1[2], double P2[2])
{
  double Vector0[2], Vector1[2];

  Vector0[0] = P1[0] - P0[0];
  Vector0[1] = P1[1] - P0[1];

  Vector1[0] = P2[0] - P1[0];
  Vector1[1] = P2[1] - P1[1];

  double dotProd = Vector0[0] * Vector1[0] + Vector0[1] * Vector1[1];
  double A = sqrt(pow(Vector0[0],2) +pow(Vector0[1],2) );
  double B = sqrt(pow(Vector1[0],2) +pow(Vector1[1],2) );
  double ang = acos(dotProd / (A*B));

  return sqrt(pow((1-cos(ang)),2) + pow(sin(ang),2));


}

static PyObject* dirCostFunc(PyObject* self, PyObject* args)
  {

    double P_0[2];
    double P_1[2];
    double P_2[2];

    if(!PyArg_ParseTuple(args,"(dd)(dd)(dd)",&P_0[0], &P_0[1], &P_1[0], &P_1[1], &P_2[0], &P_2[1]))
      return NULL;

    double ansReturn = CdirectionalCostFunction(P_0, P_1, P_2);
    //double ansReturn = 1;
    //return Py_BuildValue("d", ansReturn );
    return PyFloat_FromDouble(ansReturn);

  }



//declare Test Fib function and relevant PY stuff////////////////
int Cfib(int n)
{
  if (n < 2)
    return n;

  else

    return Cfib(n-1) + Cfib(n-2);
}



static PyObject* fib(PyObject* self, PyObject* args)
{
  int n;

  if(!PyArg_ParseTuple(args,"i", &n))
    return NULL;

  return Py_BuildValue("i", Cfib(n));

}


static PyObject* version(PyObject* self)
{
  return Py_BuildValue("s", "Version 1.0");
}

///////////////////////////////////////////////////////

static PyMethodDef myMethods[] = {

  {"fib", fib, METH_VARARGS, "Calc Fib Seq"},
  {"paraPathFunc",paraPathFunc , METH_VARARGS, "Calc para path"},
  {"miniConvex", miniConvex , METH_VARARGS, "Calc minimum"},
  {"dirCostFunc", dirCostFunc , METH_VARARGS, "Calc directionalCostFunction"},
  {"isInRegion", isInRegion , METH_VARARGS, "Calc directionalCostFunction"},
  {"version", (PyCFunction)version, METH_NOARGS, "Returns Versions"},
  {NULL, NULL, 0, NULL}


};


static struct PyModuleDef zeppNumLib2 = {

  PyModuleDef_HEAD_INIT,
  "zeppNumLib2",
  "Fibonacci Module",
  -1,
  myMethods

};

PyMODINIT_FUNC PyInit_zeppNumLib2(void)
{

  return PyModule_Create(&zeppNumLib2);

}
