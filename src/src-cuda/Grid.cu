// member function definitions for Grid.C
#include "SetPrecision.cuh"
#include "Grid.cuh"
#include "math.h"
#include "iostream"
#include "string"
#include <algorithm>
#include <vector>
#include "fstream"
#include "sstream"
#include "numeric"

__global__ void UpdateLaserGlobal(Grid *dg, real *laser_coor, real *laser_coor2){
  //********************************************************
  // NOTE: ONLY RUN WITH 1 BLOCK AND 1 THREAD (THI IS SERIAL 
  // CODE THAT UPDATES GLOBAL VARIABLES)
  //********************************************************
  int tid=threadIdx.x + blockDim.x*blockIdx.x;
  if (tid==0){dg->UpdateLaser(laser_coor,laser_coor2);}
} // end UpdateLaserGlobal

__global__ void UpdateTime2Global(Grid *dg, const real dt)
{
  dg->time+=dt;
  dg->tInd+=1;
}



//constructor
Grid::Grid(std::string &filInput)
{
  // initialize default settings which will change if specified input file
  time =0.0;
  tInd =0;
  nX[0]=128; nX[1]=128; nX[2]=64;
  dX[0]=1e-6; dX[1]=1e-6;dX[2]=1e-6;
  nDim = 3;
  outint = 1;
  patternID = 1;
  tL = 1620;
  tS = 1531.5;
  ntype = 0; // "Moore" = 0 and  "VonNeumann" =1
  nnodePerCell = pow(2,nDim);
  beamSTD[0]=7.5e-5;  beamSTD[1]=7.5e-5;  beamSTD[2]=7.5e-5;
  layerT = 25e-6;
  bpH = 0.0;
  mu = .01; // note that this is rate/ (\mu m)^3 for base plate tessellation
  T0 = 300;
  ictrl = 3;
  gsize[0]=0; gsize[1]=0;
  meltparam[0]=75e-6; meltparam[1]=162.75e-6;meltparam[2]=75e-6; meltparam[4]=75e-6;
  bhatch = 1.53*meltparam[2];
  bmDelT = 4.0/3.0*meltparam[0]/bmV;
  bmV = 500e-3;
  rNmax = .002;
  offset[0]=0.;offset[1]=0.;offset[2]=0.; // positive value means starting outside of domain
  outNL = 0;
  gth0=0.;
  gth=0.;
  Avel=0.;
  nvel=0.;
  //read data from input file
  readInputFile(filInput);
  //if (gsize[0]==0){gsize={nX[0]*dX[0]*2,nX[1]*dX[1]*2};}
  lX[0] = nX[0]*dX[0];lX[1]=nX[1]*dX[1];lX[2]=nX[2]*dX[2];
  bpH< std::numeric_limits<real>::epsilon() ? bpH=layerT : bpH=bpH;
  nZlayer = int(round(layerT/dX[2]));
  isp=0;
  inewscanflg=1;
  inewlayerflg=1;
  Nzhg = std::min(nX[2],int(ceil(bpH*.9999/dX[2])));
  indlayer=1;
  ilaserLoc= Nzhg + indlayer*nZlayer;
  Ntd=int(ceil(gsize[1]/bhatch))+1;
  Nsd=int(ceil(gsize[0]/(bmV*bmDelT)))+1;
  NpT=Nsd*Ntd;
  lcoor = (real*)malloc(2*NpT*sizeof(real));
  lcoor2 = (real*)malloc(2*NpT*sizeof(real));
  real xlmin=(lX[0]-gsize[0])/2.,ylmin=(lX[1]-gsize[1])/2.;
  int k;
  real gmid[2];
  gmid[0]=lX[0]/2.;
  gmid[1]=lX[1]/2.;
  gth+=gth0;
  for (int j2=0;j2<Ntd;++j2){
    for (int j1=0;j1<Nsd;++j1){
      k=Nsd*j2+j1;
      lcoor[2*k]=j1*bmDelT*bmV+xlmin;
      lcoor[2*k+1]=j2*bhatch+ylmin;
      lcoor2[2*k] = cos(gth)*(lcoor[2*k]-gmid[0])-
	sin(gth)*(lcoor[2*k+1]-gmid[1])+gmid[0];
      lcoor2[2*k+1] = sin(gth)*(lcoor[2*k]-gmid[0])+
	cos(gth)*(lcoor[2*k+1]-gmid[1])+gmid[1];
    } // j1
  } // j2   
  //grid box
  gbox[0]=-bhatch/2.;
  gbox[1]=lX[0]+bhatch/2.;
  gbox[2]=-bhatch/2.;
  gbox[3]=lX[1]+bhatch/2.;               
  nlayerTot=int(ceil( (real)(nX[2]-Nzhg)/(real)nZlayer));
} // end constructor

__device__ void Grid::UpdateLaser(real *laser_coor, real *laser_coor2){
  //********************************************************
  // NOTE: ONLY RUN WITH 1 BLOCK AND 1 THREAD (THI IS SERIAL 
  // CODE THAT UPDATES GLOBAL VARIABLES)
  //********************************************************
  int itmp,iflg=0,irep=0,k;
  real x,y,gmid[2];


  itmp = (isp+1) - ((isp+1)/Nsd)*Nsd;
  if (itmp==0){inewscanflg=1;}

  if (isp==(NpT-1)){
      inewscanflg=1;
      inewlayerflg=1;
      isp=0;
      indlayer+=1;
  }
  while(irep==0 || isp==0){
    irep+=1;
    itmp=isp;
    if (inewscanflg==0 && itmp<(NpT-1)){
      isp+=1;
    } else {
      while (itmp<(NpT-1) && iflg==0){
	itmp+=1;
	x=laser_coor2[2*itmp];
	y=laser_coor2[2*itmp+1];
	if (x>gbox[0] && x<gbox[1] && y>gbox[2] && y<gbox[3]){iflg=1;}
      } // while (itmp< ...
      if (itmp<(NpT-1)){
	isp=itmp;
	inewscanflg=0;
      } else {
	inewlayerflg=1;
	inewscanflg=1;
	indlayer+=1;
	ilaserLoc= Nzhg + indlayer*nZlayer;	
	isp=0;	
	bcheck=indlayer>nlayerTot;
	// update grid 
	gmid[0]=lX[0]/2.;
	gmid[1]=lX[1]/2.;	
	gth+=gth0;		
	for (int j2=0;j2<Ntd;++j2){
	  for (int j1=0;j1<Nsd;++j1){
	    k=Nsd*j2+j1;
	    laser_coor2[2*k] = cos(gth)*(laser_coor[2*k]-gmid[0])-
	      sin(gth)*(laser_coor[2*k+1]-gmid[1])+gmid[0];
	    laser_coor2[2*k+1] = sin(gth)*(laser_coor[2*k]-gmid[0])+
	      cos(gth)*(laser_coor[2*k+1]-gmid[1])+gmid[1];
	  } // j1
	} // j2	

      } // if/else (itmp<NpT-1...
    } // if (inewscanflg==0...
  } // while(irep==0...

  if (irep==1){inewlayerflg=0;} 
} // end UpdateLaser 

__device__ void Grid::GetNeighbors(int &jvox, int *ineigh)
{
  // ineigh is static array length 27 where 
  // ineigh[26] is # of neighbors for voxel jvox
  int j3,j2,j1,jst, itmp[3]={-1,0,1},nx0=nX[0],nx1=nX[1],nx2=nX[2],cc;
  j3 = jvox /( nx0*nx1) ;
  j2 = (jvox - nx0*nx1*j3)/nx0 ;
  j1 = jvox - nx0*nx1*j3 - nx0*j2;
  cc=0;
  for (int i3 =0;i3<3;++i3){
    if ( (j3+itmp[i3]<0) || (j3+itmp[i3]>=nx2)){continue;}
    for (int i2 =0;i2<3;++i2){
      if ( (j2+itmp[i2]<0) || (j2+itmp[i2]>=nx1)){continue;}
      for (int i1 =0;i1<3;++i1){
	if ( (j1+itmp[i1]<0) || (j1+itmp[i1]>=nx0)){continue;}
	jst = nx0*nx1*(j3+itmp[i3])+nx0*(j2+itmp[i2])+j1+itmp[i1];
	if (jst !=jvox){
	  ineigh[cc] = jst;
	  cc+=1;
	}
      }
    }
  }
  ineigh[26]=cc;
}
void Grid::readInputFile(std::string &filInput)
{
  std::ifstream filIn;
  std::string inputData,keyword;
  int k=0,n=0;
  filIn.open(filInput.c_str());
  while (!filIn.eof()){
    char c;
    filIn.get(c);
    if (c == '#') {
      // skip comments indicated by "#" in input file
      while (filIn && c != '\n') {
	filIn.get(c);
      }
    }
    if (k >= n) {
      n = 2*n;
      const int m = 10000;
      if (n < m) n = m;
      inputData.resize(n);
    }
    inputData[k++] = c;
  }
  filIn.close();
  inputData.resize(k);  
  std::istringstream simInput(inputData);
  simInput >> keyword;
  while(simInput){
    std::transform(keyword.begin(),keyword.end(),keyword.begin(),
                   [](unsigned char s){return std::tolower(s);});
    if (keyword=="dx") {
      simInput >> keyword;
      dX[0]=std::stod(keyword);
      dX[1]=dX[0];
      dX[2]=dX[0];
    }
    if (keyword=="bmv") {
      simInput >> keyword;
      bmV=std::stod(keyword);
    }
    if (keyword=="bmp") {
      simInput >> keyword;
      bmP=std::stod(keyword);
    }
    if (keyword=="bmdelt") {
      simInput >> keyword;
      bmDelT=std::stod(keyword);
    }
    if (keyword=="nx") {
      simInput >> keyword;
      nX[0]=std::stoi(keyword);
      simInput >> keyword;
      nX[1]=std::stoi(keyword);
      simInput >> keyword;
      nX[2]=std::stoi(keyword);
    }
    if (keyword=="r") {
      simInput >> keyword;
      rNmax=std::stod(keyword);
    }
    if (keyword=="lt") {
      simInput >> keyword;
      layerT=std::stod(keyword);
    }
    if (keyword=="tl") {
      simInput >> keyword;
      tL=std::stod(keyword);
    }
    if (keyword=="ts") {
      simInput >> keyword;
      tS=std::stod(keyword);
    }
    if (keyword=="bhatch") {
      simInput >> keyword;
      bhatch=std::stod(keyword);
    }
    if (keyword=="mu") {
      simInput >> keyword;
      mu=std::stod(keyword);
    }
    if (keyword=="bph") {
      simInput >> keyword;
      bpH=std::stod(keyword);
    }
    if (keyword=="lrate") {
      simInput >> keyword;
      lrate=std::stod(keyword);
    }
    if (keyword=="outnl") {
      simInput >> keyword;
      outNL=std::stoi(keyword);
    }
    if (keyword=="ntsd") {
      simInput >> keyword;
      nTsd=std::stoi(keyword);
    }
    if (keyword=="offset") {
      simInput >> keyword;
      offset[0]=std::stod(keyword);
      simInput >> keyword;
      offset[1]=std::stod(keyword);
      simInput >> keyword;
      offset[2]=std::stod(keyword);
    }
    if (keyword=="meltparam") {
      simInput >> keyword;
      meltparam[0]=std::stod(keyword);
      simInput >> keyword;
      meltparam[1]=std::stod(keyword);
      simInput >> keyword;
      meltparam[2]=std::stod(keyword);
      simInput >> keyword;
      meltparam[3]=std::stod(keyword);
    }
    if (keyword=="patternid") {
      simInput >> keyword;
      patternID=std::stoi(keyword);
    }
    if (keyword=="outint") {
      simInput >> keyword;
      outint=std::stoi(keyword);
    }
    if (keyword=="gridsize") {
      simInput >> keyword;
      gsize[0]=std::stod(keyword);
      simInput >> keyword;
      gsize[1]=std::stod(keyword);
    }
    if (keyword=="gridtheta") {
      simInput >> keyword;
      gth0=std::stod(keyword);
    }
    if (keyword=="dendritegrowth2") {
      simInput >> keyword;
      Avel=std::stod(keyword);
      simInput >> keyword;
      nvel=std::stod(keyword);
    }

    simInput >> keyword;
  } // while(simInput)
} // readInputFile
void Grid::UpdateTime2(const real &dtIn)
{
  time +=dtIn;
  tInd +=1;
} // end UpdateTime2

