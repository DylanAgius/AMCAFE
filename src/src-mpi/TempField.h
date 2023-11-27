// temperature class (obtained from Moose simulation)

#ifndef TEMPFIELD_H
#define TEMPFIELD_H

#include "Grid.h"
#include <string>
#include <vector>
#include "Partition.h"
#include "BasePlate.h"

class TempField
{
 public:
  // define constructor
  TempField(Grid &g, Partition &, BasePlate &);
  unsigned int simpleHash(unsigned int input);
  std::pair<double, int> tuplify(int, int);
  double randomGaussian();
  void InitializeSchwalbach();
  void InitializeAnalytic();
  void SchwalbachTempCurr();
  void AnalyticTempCurr(double tcurr,std::vector<double> &TempOut,std::vector<int> &icellid,int Ntot);
  void AnalyticTempCurrAct(double tcurr,std::vector<double> &TempOut,std::vector<int> &icellid,int Ntot);    
  void SchwalbachTempCurr(double tcurr, std::vector<double> &TempOut);
  std::vector<std::vector<double>> Temp;
  std::vector<double> DDtTemp,dXM,TempCurr,TempCurrAct,lamXYZ,bmSTD,bmLx,bmX0,bmDX,bmPeriod,offset,shiftL;
  std::vector<int> nXM,nTTemp,ispvec;
  
  int NtM,indexM,patternID,nSource,tInd;
  double dtM,bmV,bmP,bmEta,rcut,tcut,T0,Ci,DelT,alpha,T0targ,zlaserOff,tBeg0;
  std::vector<double> a1,a2,Qp,tBeg; // 
 private:
  Grid *_xyz;
  Partition *_part;
  BasePlate *_bp;
  std::string *filnambaseptr;
}; // end class TempField

#endif
