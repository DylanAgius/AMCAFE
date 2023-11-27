// member functions for TempField

#include "Grid.h"
#include "TempField.h"
#include "fstream"
// #include "iostream"
#include "math.h"
#include "Partition.h"
#include "BasePlate.h"
#include "numeric"
#include <algorithm>
#include "mpi.h"
#include <ctime>
#include <random>
#include <functional>



// constructor
TempField::TempField(Grid &g, Partition & part, BasePlate &bp)
{
  _xyz = &g;
  _part = &part;
  _bp = &bp;
  TempCurr.resize(_part->ncellLoc+_part->nGhost,0.0);
  DDtTemp.assign(_part->ncellLoc+_part->nGhost,0.0);
  tInd = 0;
  bmV = _xyz->bmV;
  patternID = _xyz->patternID;
  T0 = _xyz->T0;
  zlaserOff=1.0; // 1.0 (this specifies where laser z value is - see SchwalbachTempCurr)
  DelT = _xyz->bmDelT;
  bmDX = {DelT*bmV,_xyz->bhatch,_xyz->layerT}; // (SD, TD, BD) per layer
  offset=_xyz->offset;
 // ispvec.assign(_xyz->NpT,0);
  
  std::random_device rd;
  std::mt19937 gen(rd());
  

    std::mt19937 rng(rd());

  

  


	
  //It is here that the addtional patterns must be added.
  if (patternID==1 || patternID==3){
   ispvec.assign(_xyz->NpT,0);
    std::iota(ispvec.begin(),ispvec.end(),0);
    
    for (int value : ispvec) {
       std::cout << value << " ";
   }
  } // if (patternID==1...
  if (patternID==2 || patternID==4){
   ispvec.assign(_xyz->NpT,0);
    int k;
    for (int j=0;j<_xyz->Ntd;++j){
      k=_xyz->Nsd*j;
      if (fmod(j,2)==0){
	std::iota(ispvec.begin()+k,ispvec.begin()+k+_xyz->Nsd,k);
      } // if (fmod(j,2)==0...
      else {
	for (int j1=0;j1<_xyz->Nsd;++j1){
	  ispvec[k+j1]=k+_xyz->Nsd-1-j1;
	} // for (int j1=0...
      } // else (fmod(j,2)==0...
    } // for (int j=0...
  } // if (patternID==2...
  
  // pattern 5 is now a random selection
  if (patternID==5){
  // A simple hash function (not secure for cryptographic use)



    std::vector<unsigned int> numbers;
    for (int i = 0; i < _xyz->NpT; ++i) {
        ispvec.push_back(i);
    }

    // Use the hash function to reorder the numbers
    for (int i = 0; i < _xyz->NpT; ++i) {
        unsigned int hashValue = simpleHash(i + 1);
        unsigned int newIndex = hashValue % (i + 1);

       // Swap numbers[i] and numbers[newIndex]
        unsigned int temp = ispvec[i];
        ispvec[i] = ispvec[newIndex];
        ispvec[newIndex] = temp;
    }
    }
    
  if (patternID==8){
  //std::srand(static_cast<unsigned int>(std::time(nullptr)));
 
    //ispvec=_xyz->ispvec;
    
   std::vector<int> values(_xyz->NpT);
   std::iota(values.begin(), values.end(), 0);


    std::vector<std::pair<double, int>> pairs;
    pairs.reserve(values.size());

   for (int i = 0; i < values.size(); ++i) {
       pairs.push_back(tuplify(values[i], i));
    }

    std::sort(pairs.begin(), pairs.end());

  
    ispvec.reserve(pairs.size());

   for (const auto& p : pairs) {
        ispvec.push_back(p.second);
  }
    
  //   for (int value : ispvec) {
  //      std::cout << value << " ";
  //  }

   
  }
    
    
  if (patternID==7){
     int k = 0;
     std::vector<int> test;
     int Nsd1=_xyz->Nsd/2;
     int Nsd2=_xyz->Nsd/2;
    
   for (int j = 0; j < _xyz->Ntd*2; ++j) {
           k = j * Nsd2; // Incrementally increase k by Nsd2
        
        if (j % 2 == 0) {
            std::vector<int> val;
            for (int idx = k; idx < k + Nsd1; ++idx) {
                val.push_back(idx);
            }
            
            std::srand(static_cast<unsigned int>(std::time(nullptr)));
            std::random_shuffle(val.begin(), val.end());
            
            ispvec.insert(ispvec.end(), val.begin(), val.end());
        } else {
            test.clear();
            for (int j1 = 0; j1 < Nsd1; ++j1) {
                test.push_back(k + Nsd1 - 1 - j1);
            }
            
            std::srand(static_cast<unsigned int>(std::time(nullptr)));
            std::random_shuffle(test.begin(), test.end());
            
            ispvec.insert(ispvec.end(), test.begin(), test.end());
        }
        }
        }
    
   // Pattern 6 is now the Dehoff filling method
 if (patternID==6){
    // Initialise values single values
    int lines = 5; // this will have be user defined
    int skip = 10; //this will have be user defined
    int intv=skip; 
    int init=0;
    int extraval2=0;
    int start=0;
    int start2=1;
    int ispvsum=0;
    
    FILE *fptr;
	fptr=fopen("check.txt","a");
      
			fprintf(fptr,"%d, %d",_xyz->Nsd,_xyz->Ntd);
			fprintf(fptr,"\n");
			
			fclose(fptr);
   
    
    // Initialise additional vectors for Dehoff fill
    std::vector<std::vector<int> > ispvecc(_xyz->Ntd,std::vector<int>(_xyz->Nsd));
    std::vector<std::vector<int> > ispvecn(_xyz->Ntd,std::vector<int>(_xyz->Nsd));
    std::vector<int> test;


    // Calculate the total sum expected in the ispvec array
    // This is used to control the while loop
    for (int v=0; v<_xyz->Nsd; v++){
        ispvsum += v;
     }
    // Outside loop which iterates across the expected number of laser passes
    for (int j=0; j<_xyz->Ntd; j++){
         // Need to update loop values since these are used within the while and for loops below
	 int check=0;
 	 int check2=0;
	 int extraval2=0;
	 int sum=0;
	 int size=_xyz->Nsd;
	 int totalv=0;
	 int totalv2=0;/////added to see the time
	 
	 

        
	 while(sum <= ispvsum-_xyz->Nsd){
	      for (int i=j; i<_xyz->Ntd; i+=lines){
	          if (i%2==0){
	             std::vector<int> test;
		     // define indeces of the first pass
		     for (int k = init; k < _xyz->Nsd; k += skip) {
                        test.push_back(k);
			totalv +=1;  // this is used to define the size of the array that is formed
			} // for int k
			if (start==0){
			    if (totalv< _xyz->Nsd/skip){
				check2=1;
				//This adds an extra few points which is based on the location selection being larger the final pass value and the indeces 
				// selection restarts
				test.push_back(extraval2);
			    } //if totalv
			 } // if start
			 // if start is not zero (this is used to determine if the pass starts from zero or at the 'skip' location
			 else {
			     if ((totalv<(_xyz->Nsd/skip)) || (totalv<totalv2)){
				check2=1;
				test.push_back(extraval2);
			      }//if totalv
			  }//else 
			  //adds to the vector which defines the passes along the beam length
			  ispvecc[i].insert(ispvecc[i].end(),test.begin(),test.end());
	                  totalv2=0; // updates the total value to be used again in the next pass
			  //update array to remove leading 0s
			  for (int del=0; del<_xyz->Nsd;del++){
			      ispvecn[i][del]=ispvecc[i][_xyz->Nsd+del];
			   }
		    }//if i%2
		    //These are the odd value lines which must also alternate from starting from zero to the "skip" value
	            else {
			std::vector<int> test;
			for (int x=intv; x<_xyz->Nsd; x += skip){
			    test.push_back(x);
			    totalv2 += 1;
			 }//for  int x=intv
			 if (start2==0){
                            if (totalv2<(_xyz->Nsd/skip)){
			       check2=1;
			       test.push_back(extraval2);
			     } //if totalv2
		          } // if start2=0
			  else {
                             if ((totalv2< (_xyz->Nsd/skip)) || (totalv2 <totalv)){
				check2=1;
				// Add the extra value if the pass is larger than the are area.
				test.push_back(extraval2);
			      } //if totalv2 ||
			   } //else
			   //adds to the vector which defines the passes along the beam length
			   ispvecc[i].insert(ispvecc[i].end(),test.begin(),test.end());
			   totalv=0;
			   //update array to remove leading 0s
			   for (int del=0; del<_xyz->Nsd;del++){
			       ispvecn[i][del]=ispvecc[i][_xyz->Nsd+del];
			    } // for del =0
		     }//else

	         }//for i=j
                 //calculate the sum 
		 sum=accumulate(ispvecn[j].begin(),ispvecn[j].end(),0);
                 // update iteration values for the starting points of each pass
		 init += 1;
		 intv += 1;
		 //check to see the next vector to ensure length is correct
		 if (check2==0){
		    extraval2=0;
		    check2=1;
		 }//if check2
		 else {
                    extraval2++;
		  }//else 
		

	      } //while loop
	      
	    
	      //change order because now the starting point is an odd
	      if ((j+1)%2==0) {
		intv=skip;
		init=0;
		start=0;
		start2=1;
	       }//if ((j+1)
	      else {
		intv=0;
		init=skip;
		start=1;
		start2=0;
	      }//else	
         }//for j=0


	//Now need to flatten the array
	int value=0;

	fptr=fopen("check.txt","a");
	for (int i=0; i<_xyz->Ntd; i++){
		for (int j=0; j<_xyz->Nsd; j++){
			ispvec[value]=ispvecn[i][j]+(_xyz->Nsd)*i;
			
			fprintf(fptr,"%d",ispvec[value]);
			fprintf(fptr,"\n");
			value++;
		  
			
		} // for int j=0
	} // for int i=0 
	
    fclose(fptr);
  ////////////////

 } // if patternID==6
   
 ////// Location of additional patterns
} // end TempField

unsigned int TempField::simpleHash(unsigned int input) {
    return (input * 2654435761u) % 4294967296u;
}


  std::pair<double, int> TempField::tuplify(int x, int y) {
    //std::normal_distribution<double> normalDist(0, 1);
      double orderliness = 0.5;
    double first = orderliness * y + randomGaussian();
    return std::make_pair(first, x);
}

// Function to generate normally distributed random numbers using Box-Muller transform
double TempField::randomGaussian() {
    double u1 = static_cast<double>(rand()) / RAND_MAX;
    double u2 = static_cast<double>(rand()) / RAND_MAX;
    return sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
}




void TempField::InitializeAnalytic()
{
  /*
    This this 2 double ellipsoids (one encompassing another) to represent the temperature
    field. This approach is used in Rogers_Madison_Tikare, Comp Mat Sci, 2017.
    SCAN INFORMATION
    patternID=0: scan in +X direction only
    patternID=1: scan in alternating +/- X direction (same direction in
		 layer above and below
    patternID=2: scan in alternating +/- X direction (different direction 
		 in layer above and below
    patternID=3: scan in +X direction for layer i and +Y direction for 
		 layer i+1 ...
    patternID=4: scan in alternating +/- X direction for layer i and 
		 alternating +/- Y direction for layer i+1
   */
  a1.resize(6);
  Qp.resize(2);
  a1[0] = _xyz->meltparam[0];
  a1[1] = _xyz->meltparam[2];
  a1[2] = _xyz->meltparam[3];
  a1[3] = _xyz->meltparam[1];
  a1[4] = a1[1];
  a1[5] = a1[2];
  

} // end InitializeAnalytic 

void TempField::AnalyticTempCurr(double tcurr,std::vector<double> & TempOut, std::vector<int> &icellid, int Ntot)
{
  //computes temp field based on Schwalbach et al
  int j1,j2,j3,iplay;
  double x0,y0,x,y,z,dsq,dsq2,bx,by,xi,xp,yp,dirp,zp;
  std::vector<double> rij1(3),xs1(3),xs2(3);
  iplay=_xyz->nX[0]*_xyz->nX[1]*_xyz->ilaserLoc;
  TempOut.assign(Ntot,T0);
  xi = _xyz->tL*(1+std::numeric_limits<double>::epsilon() );
  int js1, js2;
  std::vector<double> a1m(6);
  for (int j=0;j<6;++j){a1m[j]=a1[j];}
  // x,y,z spatial location of source (in grid reference frame)
    x=_xyz->lcoor[2*ispvec[_xyz->isp]]-offset[0];
  y=_xyz->lcoor[2*ispvec[_xyz->isp]+1]-offset[1];
  z=_xyz->ilaserLoc*_xyz->dX[2]-offset[2];
  for (int j=0;j<Ntot;++j){
    j3 = floor(icellid[j]/(_xyz->nX[0]*_xyz->nX[1]));
    j2 = floor( (icellid[j]- _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
    j1 = icellid[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
    x0 = (double(j1)+.5)*(_xyz->dX[0]);
    y0 = (double(j2)+.5)*(_xyz->dX[1]);
    zp = (double(j3)+.5)*(_xyz->dX[2]);
    if (zp>z){continue;}
    xp=cos(_xyz->gth)*(x0-_xyz->LX[0]/2.)+
      sin(_xyz->gth)*(y0-_xyz->LX[1]/2.)+_xyz->LX[0]/2.;
    yp=-sin(_xyz->gth)*(x0-_xyz->LX[0]/2.) + 
      cos(_xyz->gth)*(y0-_xyz->LX[1]/2.)+_xyz->LX[1]/2.;
    if (fmod(_xyz->isp,_xyz->Nsd)==0){
      dirp=(_xyz->lcoor[2*ispvec[_xyz->isp+1]]-_xyz->lcoor[2*ispvec[_xyz->isp]]);
    } else {
      dirp=(_xyz->lcoor[2*ispvec[_xyz->isp]]-_xyz->lcoor[2*ispvec[_xyz->isp-1]]);
    }
    rij1[0] = xp-x;
    rij1[1] = yp-y;
    rij1[2] = zp-z;
    if (dirp*rij1[0]>0){ 
      //xp,yp is in front of laser
      dsq = pow(rij1[1]/a1m[1],2.0)+pow(rij1[2]/a1m[2],2.0);
      if (dsq<1.0 && (fabs(rij1[0])<bmDX[0]) ){
	TempOut[j]=xi;
      } else {
	TempOut[j] = _xyz->tS;
      }
    } else {
      dsq = pow(rij1[0]/a1m[3],2.0)+pow(rij1[1]/a1m[4],2.0)+pow(rij1[2]/a1m[5],2.0);
      if (dsq<1.0){
	TempOut[j]=xi;
      } else {
	TempOut[j] = _xyz->tS;
      }
    } // if (dirp*rij1[0]>0...
  } // for (int j=0...
  // update to new scan if nothing melted (i.e. laser well outside domain)
  double tmelt=_xyz->tL;
  int n1=_part->ncellLoc,icheck,ichecktmp;
  x=_xyz->lcoor2[2*ispvec[_xyz->isp]]-offset[0];
  y=_xyz->lcoor2[2*ispvec[_xyz->isp]+1]-offset[1];
  if (x<_xyz->gbox[0] || x>_xyz->gbox[1] || y<_xyz->gbox[2] || y>_xyz->gbox[3]){
    ichecktmp=std::any_of(TempOut.begin(),TempOut.begin()+n1,[&tmelt]
			  (double tchk){return tchk >= tmelt;});
    MPI_Allreduce(&ichecktmp,&icheck,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if (icheck==0){ 
      _xyz->inewscanflg=1;
    } // if (icheck==0...
  } // if (x<box[0...
} // end AnalyticTempCurr()          

void TempField::AnalyticTempCurrAct(double tcurr,std::vector<double> & TempOut, std::vector<int> &icellid, int Ntot)
{
  //computes temp field based on Schwalbach et al
  int j1,j2,j3,iplay;
  double x0,y0,x,y,z,dsq,dsq2,bx,by,xi,xp,yp,dirp,zp,rij,tc;
  std::vector<double> rij1(3),xs1(3),xs2(3);
  iplay=_xyz->nX[0]*_xyz->nX[1]*_xyz->ilaserLoc;
  
  
  //added definitions
  std::vector<double> lam(3);
  bmSTD = _xyz->beamSTD;
  bmV = _xyz->bmV;
  T0targ = _xyz->T0targ;
  bmEta = _xyz->beamEta;
  T0 = _xyz->T0;
  bmDX = {DelT*bmV,2.7*bmSTD[1],_xyz->layerT};
  double x1 = pow( pow(bmSTD[0],2.0)*pow(bmSTD[1],2.0)*pow(bmSTD[2],2.0),.5);
  //bmP = T0targ*_xyz->cP*_xyz->rho*pow(2.0,.5)*pow(M_PI,1.5)*x1/DelT;
  bmP=_xyz->bmP;
  rcut = pow( -2* pow(*std::min_element(bmSTD.begin(),bmSTD.end()),2.0)*log(.001),.5);
  Ci = bmEta*bmP*DelT/(_xyz->rho*_xyz->cP*pow(2.0,.5)*pow(M_PI,1.5));
  
  TempOut.assign(Ntot,T0);
  xi = _xyz->tL*(1+std::numeric_limits<double>::epsilon() );
  int js1, js2;
  std::vector<double> a1m(6);
  for (int j=0;j<6;++j){a1m[j]=a1[j];}
  // x,y,z spatial location of source (in grid reference frame)
    x=_xyz->lcoor[2*ispvec[_xyz->isp]]-offset[0];
  y=_xyz->lcoor[2*ispvec[_xyz->isp]+1]-offset[1];
  z=_xyz->ilaserLoc*_xyz->dX[2]-offset[2];
  for (int j=0;j<Ntot;++j){
    j3 = floor(icellid[j]/(_xyz->nX[0]*_xyz->nX[1]));
    j2 = floor( (icellid[j]- _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
    j1 = icellid[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
    x0 = (double(j1)+.5)*(_xyz->dX[0]);
    y0 = (double(j2)+.5)*(_xyz->dX[1]);
    zp = (double(j3)+.5)*(_xyz->dX[2]);
    if (zp>z){continue;}
    xp=cos(_xyz->gth)*(x0-_xyz->LX[0]/2.)+
      sin(_xyz->gth)*(y0-_xyz->LX[1]/2.)+_xyz->LX[0]/2.;
    yp=-sin(_xyz->gth)*(x0-_xyz->LX[0]/2.) + 
      cos(_xyz->gth)*(y0-_xyz->LX[1]/2.)+_xyz->LX[1]/2.;
    if (fmod(_xyz->isp,_xyz->Nsd)==0){
      dirp=(_xyz->lcoor[2*ispvec[_xyz->isp+1]]-_xyz->lcoor[2*ispvec[_xyz->isp]]);
    } else {
      dirp=(_xyz->lcoor[2*ispvec[_xyz->isp]]-_xyz->lcoor[2*ispvec[_xyz->isp-1]]);
    }
    rij1[0] = xp-x;
    rij1[1] = yp-y;
    rij1[2] = zp-z;
    if (dirp*rij1[0]>0){ 
      //xp,yp is in front of laser
      dsq = pow(rij1[1]/a1m[1],2.0)+pow(rij1[2]/a1m[2],2.0);
      //if (dsq<1.0 && (fabs(rij1[0])<bmDX[0]) ){
	//TempOut[j]=xi;
      //} else {
      
        //Added the following
          tc = _xyz->time - DelT;
         
          lam = {pow(bmSTD[0],2.0)+2*alpha*(_xyz->time - tc), 
	       pow(bmSTD[1],2.0)+2*alpha*(_xyz->time - tc),
	       pow(bmSTD[2],2.0)+2*alpha*(_xyz->time - tc)};
	       
	  rij = pow(xp-x,2.0)/2/lam[0] +pow(yp-y,2.0)/2/lam[1] +
	  pow(zp-z,2.0)/2/lam[2];
          TempOut[j] += Ci/pow(lam[0]*lam[1]*lam[2],.5)*exp(-rij);
	//TempOut[j] = _xyz->tS;
      //}
    } else {
      dsq = pow(rij1[0]/a1m[3],2.0)+pow(rij1[1]/a1m[4],2.0)+pow(rij1[2]/a1m[5],2.0);
      //if (dsq<1.0){
	//TempOut[j]=xi;
      //} else {
          //Added the following
          tc = _xyz->time - DelT;
         
          lam = {pow(bmSTD[0],2.0)+2*alpha*(_xyz->time - tc), 
	       pow(bmSTD[1],2.0)+2*alpha*(_xyz->time - tc),
	       pow(bmSTD[2],2.0)+2*alpha*(_xyz->time - tc)};
	       
	  rij = pow(xp-x,2.0)/2/lam[0] +pow(yp-y,2.0)/2/lam[1] +
	  pow(zp-z,2.0)/2/lam[2];
          TempOut[j] += Ci/pow(lam[0]*lam[1]*lam[2],.5)*exp(-rij);
	//TempOut[j] = _xyz->tS;
      //}
    } // if (dirp*rij1[0]>0...
  } // for (int j=0...
 
	
  // update to new scan if nothing melted (i.e. laser well outside domain)
  double tmelt=_xyz->tL;
  int n1=_part->ncellLoc,icheck,ichecktmp;
  x=_xyz->lcoor2[2*ispvec[_xyz->isp]]-offset[0];
  y=_xyz->lcoor2[2*ispvec[_xyz->isp]+1]-offset[1];
  if (x<_xyz->gbox[0] || x>_xyz->gbox[1] || y<_xyz->gbox[2] || y>_xyz->gbox[3]){
    ichecktmp=std::any_of(TempOut.begin(),TempOut.begin()+n1,[&tmelt]
			  (double tchk){return tchk >= tmelt;});
    MPI_Allreduce(&ichecktmp,&icheck,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if (icheck==0){ 
      _xyz->inewscanflg=1;
    } // if (icheck==0...
  } // if (x<box[0...
} // end AnalyticTempCurrAct()            
void TempField::InitializeSchwalbach()
{
  /*
    This follows the model: Schwalbach, Edwin J., et al. "A discrete source model of powder 
    bed fusion additive manufacturing thermal history." Additive Manufacturing 25 (2019): 485-498.
    SCAN INFORMATION
    patternID=0: scan in +X direction only
    patternID=1: scan in alternating +/- X direction (same direction in
		 layer above and below
    patternID=2: scan in alternating +/- X direction (different direction 
		 in layer above and below
    patternID=3: scan in +X direction for layer i and +Y direction for 
		 layer i+1 ...
    patternID=4: scan in alternating +/- X direction for layer i and 
		 alternating +/- Y direction for layer i+1
   */
  bmSTD = _xyz->beamSTD;
  bmV = _xyz->bmV;
  T0targ = _xyz->T0targ;
  bmEta = _xyz->beamEta;
  patternID = _xyz->patternID;
  T0 = _xyz->T0;
  zlaserOff=1.0; // 1.0 (this specifies where laser z value is - see SchwalbachTempCurr)
  double tb,minTemp;
  if (patternID==0 || patternID==1 || patternID==2 || patternID==3 || patternID==4){
    DelT = 4.0/3.0*bmSTD[0]/bmV;
    //bmDX = {DelT*bmV,4.0/3.0*bmSTD[1],_xyz->layerT};
    bmDX = {DelT*bmV,2.7*bmSTD[1],_xyz->layerT};
    double x1 = pow( pow(bmSTD[0],2.0)*pow(bmSTD[1],2.0)*pow(bmSTD[2],2.0),.5);
    bmP = T0targ*_xyz->cP*_xyz->rho*pow(2.0,.5)*pow(M_PI,1.5)*x1/DelT;
    rcut = pow( -2* pow(*std::min_element(bmSTD.begin(),bmSTD.end()),2.0)*log(.001),.5);
    Ci = bmEta*bmP*DelT/(_xyz->rho*_xyz->cP*pow(2.0,.5)*pow(M_PI,1.5));
    alpha = _xyz->kappa/_xyz->cP/_xyz->rho;
    minTemp = 0.250; // cut off change in temperature for tcut
    tcut = pow(Ci/(2*alpha*minTemp),2.0/3.0);		       
    nSource = ceil(tcut/DelT);
    offset=_xyz->offset;
    bmLx={_xyz->LX[0]+1*bmDX[0],_xyz->LX[1]+1*bmDX[1],_xyz->LX[2]};
    if (patternID==1){
      offset={0.0,-_xyz->nX[1]*_xyz->dX[1]/2.0,0.0}; // positive value means starting outside domain
      shiftL={3*bmDX[0],0.0,0.0};
      bmLx={_xyz->LX[0]+shiftL[0],_xyz->LX[1]+1*bmDX[1],_xyz->LX[2]};
    } // if (patternID==1...
    nTTemp = {int(floor(bmLx[0]/bmDX[0] ))+1,int(floor(bmLx[1]/bmDX[1]))+1,
            int(floor(bmLx[2]/bmDX[2]))};
    bmLx={(nTTemp[0]-1)*bmDX[0],(nTTemp[1]-1)*bmDX[1],(nTTemp[2]-1)*bmDX[2]};
  }

} // end InitializeSchwalbach

void TempField::SchwalbachTempCurr()
{
  //computes temp field based on Schwalbach et al
  int Ntot = _part->nGhost+_part->ncellLoc, j1,j2,j3,iplay;
  double x0,y0,z0,rij,x,y,z,tc,xst;
  std::vector<double> lam(3);
  iplay=_xyz->nX[0]*_xyz->nX[1]*_xyz->ilaserLoc;
  TempCurr.assign(Ntot,T0);
  if (patternID==0){
    for (int j=0;j<Ntot;++j){
      //if (_part->icellidLoc[j] >= iplay){continue;}
      j3 = floor(_part->icellidLoc[j]/(_xyz->nX[0]*_xyz->nX[1]));
      j2 = floor( (_part->icellidLoc[j]- _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
      j1 = _part->icellidLoc[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
      x0 = (double(j1)+.5)*(_xyz->dX[0]);
      y0 = (double(j2)+.5)*(_xyz->dX[1]);
      z0 = (double(j3)+.5)*(_xyz->dX[2]);
      for (int jt=0;jt<nSource;++jt){
	// tc,x,y,z space-time location of source
	tc = _xyz->time - jt*DelT;
	x = fmod(round(tc/DelT),(nTTemp[0]))*bmDX[0] - offset[0];
	y = floor(fmod(round(tc/DelT),(nTTemp[0]*nTTemp[1]))/nTTemp[0])*bmDX[1]-offset[1];
	z = (floor((tc/DelT+DelT*1e-6)/(nTTemp[0]*nTTemp[1]))+zlaserOff)*bmDX[2] + _bp->height - 
	  ceil(offset[2]/_xyz->dX[2])*_xyz->dX[2];
	lam = {pow(bmSTD[0],2.0)+2*alpha*(_xyz->time - tc), 
	       pow(bmSTD[1],2.0)+2*alpha*(_xyz->time - tc),
	       pow(bmSTD[2],2.0)+2*alpha*(_xyz->time - tc)};
	rij = pow(x-x0,2.0)/2/lam[0] +pow(y-y0,2.0)/2/lam[1] +
	  pow(z-z0,2.0)/2/lam[2];
	TempCurr[j] += Ci/pow(lam[0]*lam[1]*lam[2],.5)*exp(-rij);
      } // for (int j1...
    } // for (int j...
  } // if (patternID==0 ...
  if (patternID==1){
    int js1, js2;
    y = floor(fmod(round(_xyz->time/DelT),(nTTemp[0]*nTTemp[1]))/nTTemp[0])*bmDX[1]-offset[1];
    z = (floor((_xyz->time/DelT+DelT*1e-6)/(nTTemp[0]*nTTemp[1]))+zlaserOff)*bmDX[2] + _bp->height - 
      ceil(offset[2]/_xyz->dX[2])*_xyz->dX[2];
    for (int j=0;j<Ntot;++j){
      j3 = floor(_part->icellidLoc[j]/(_xyz->nX[0]*_xyz->nX[1]));
      j2 = floor( (_part->icellidLoc[j]- _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
      j1 = _part->icellidLoc[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
      x0 = (double(j1)+.5)*(_xyz->dX[0]);
      y0 = (double(j2)+.5)*(_xyz->dX[1]);
      z0 = (double(j3)+.5)*(_xyz->dX[2]);
      for (int jt=0;jt<nSource;++jt){
	// tc,x,y,z space-time location of source
	tc = jt*DelT;
	js1 = fmod( round(_xyz->time/DelT),nTTemp[0]*nTTemp[1]);
	js2 = fmod(floor(js1/nTTemp[0])+1,2);
	x = (.5*pow(-1,js2)+.5)*bmLx[0] - 
	  pow(-1,js2)*fmod(js1,nTTemp[0])*bmDX[0] + pow(-1,js2)*tc*bmV - 
	  (pow(-1,js2)+1)/2.0*shiftL[0] + pow(-1,js2)*offset[0] ;
	//if (_part->myid==0){if (j==0 & jt==0){std::cout<< tInd<<","<<x<<","<<y<<","<<z<<std::endl;}}
	lam = {pow(bmSTD[0],2.0)+2*alpha*(tc), 
	       pow(bmSTD[1],2.0)+2*alpha*(tc),
	       pow(bmSTD[2],2.0)+2*alpha*(tc)};
	rij = pow(x-x0,2.0)/2/lam[0] +pow(y-y0,2.0)/2/lam[1] +
	  pow(z-z0,2.0)/2/lam[2];
	TempCurr[j] += Ci/pow(lam[0]*lam[1]*lam[2],.5)*exp(-rij);
      } // for (int j1...
    } // for (int j...
  } // if (patternID==1 ...
  if (patternID==2){
    int js1, js2,js3;
    for (int j=0;j<Ntot;++j){
      j3 = floor(_part->icellidLoc[j]/(_xyz->nX[0]*_xyz->nX[1]));
      j2 = floor( (_part->icellidLoc[j]- _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
      j1 = _part->icellidLoc[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
      x0 = (double(j1)+.5)*(_xyz->dX[0]);
      y0 = (double(j2)+.5)*(_xyz->dX[1]);
      z0 = (double(j3)+.5)*(_xyz->dX[2]);
      for (int jt=0;jt<nSource;++jt){
	// tc,x,y,z space-time location of source
	tc = _xyz->time - jt*DelT;
	js1 = fmod( round(tc/DelT),nTTemp[0]*nTTemp[1]);
	js3 = fmod( floor( round(tc/DelT)/(nTTemp[0]*nTTemp[1]) ),2);
	js2 = fmod(floor(js1/nTTemp[0])+1+js3,2);
	x = (.5*pow(-1,js2)+.5)*bmLx[0] - 
	  pow(-1,js2)*fmod(js1,nTTemp[0])*bmDX[0] - offset[0];
	y = floor(fmod(round(tc/DelT),(nTTemp[0]*nTTemp[1]))/nTTemp[0])*bmDX[1]-offset[1];
	z = (floor((tc/DelT+DelT*1e-6)/(nTTemp[0]*nTTemp[1]))+zlaserOff)*bmDX[2] + _bp->height - 
	  ceil(offset[2]/_xyz->dX[2])*_xyz->dX[2];
	//if (_part->myid==0){if (j==0 & jt==0){std::cout<< tInd<<","<<x<<","<<y<<","<<z<<std::endl;}}
	lam = {pow(bmSTD[0],2.0)+2*alpha*(_xyz->time - tc), 
	       pow(bmSTD[1],2.0)+2*alpha*(_xyz->time - tc),
	       pow(bmSTD[2],2.0)+2*alpha*(_xyz->time - tc)};
	rij = pow(x-x0,2.0)/2/lam[0] +pow(y-y0,2.0)/2/lam[1] +
	  pow(z-z0,2.0)/2/lam[2];
	TempCurr[j] += Ci/pow(lam[0]*lam[1]*lam[2],.5)*exp(-rij);
      } // for (int j1...
    } // for (int j...
  } // if (patternID==2 ...
  if (patternID==3){
    int js;
    for (int j=0;j<Ntot;++j){
      //if (_part->icellidLoc[j] >= iplay){continue;}
      j3 = floor(_part->icellidLoc[j]/(_xyz->nX[0]*_xyz->nX[1]));
      j2 = floor( (_part->icellidLoc[j]- _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
      j1 = _part->icellidLoc[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
      x0 = (double(j1)+.5)*(_xyz->dX[0]);
      y0 = (double(j2)+.5)*(_xyz->dX[1]);
      z0 = (double(j3)+.5)*(_xyz->dX[2]);
      for (int jt=0;jt<nSource;++jt){
	// tc,x,y,z space-time location of source
	tc = _xyz->time - jt*DelT;
	js = fmod( floor( round(tc/DelT)/(nTTemp[0]*nTTemp[1])),2);
	if (js==0){
	  x = fmod(round(tc/DelT),(nTTemp[0]))*bmDX[0] - offset[0];
	  y = floor(fmod(round(tc/DelT),(nTTemp[0]*nTTemp[1]))/nTTemp[0])*bmDX[1]-offset[1];
	} else {
	  y = fmod(round(tc/DelT),(nTTemp[0]))*bmDX[0] - offset[0];
	  x = floor(fmod(round(tc/DelT),(nTTemp[0]*nTTemp[1]))/nTTemp[0])*bmDX[1]-offset[1];
	}
	z = (floor((tc/DelT+DelT*1e-6)/(nTTemp[0]*nTTemp[1]))+zlaserOff)*bmDX[2] + _bp->height - 
	  ceil(offset[2]/_xyz->dX[2])*_xyz->dX[2];
	lam = {pow(bmSTD[0],2.0)+2*alpha*(_xyz->time - tc), 
	       pow(bmSTD[1],2.0)+2*alpha*(_xyz->time - tc),
	       pow(bmSTD[2],2.0)+2*alpha*(_xyz->time - tc)};
	rij = pow(x-x0,2.0)/2/lam[0] +pow(y-y0,2.0)/2/lam[1] +
	  pow(z-z0,2.0)/2/lam[2];
	TempCurr[j] += Ci/pow(lam[0]*lam[1]*lam[2],.5)*exp(-rij);
      } // for (int j1...
    } // for (int j...
  } // if (patternID==3 ...
  if (patternID==4){
    int js1, js2,js3;
    for (int j=0;j<Ntot;++j){
      j3 = floor(_part->icellidLoc[j]/(_xyz->nX[0]*_xyz->nX[1]));
      j2 = floor( (_part->icellidLoc[j]- _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
      j1 = _part->icellidLoc[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
      x0 = (double(j1)+.5)*(_xyz->dX[0]);
      y0 = (double(j2)+.5)*(_xyz->dX[1]);
      z0 = (double(j3)+.5)*(_xyz->dX[2]);
      for (int jt=0;jt<nSource;++jt){
	// tc,x,y,z space-time location of source
	tc = _xyz->time - jt*DelT;
	js1 = fmod( round(tc/DelT),nTTemp[0]*nTTemp[1]);
	js2 = fmod(floor(js1/nTTemp[0])+1,2);
	js3 = fmod( floor( round(tc/DelT)/(nTTemp[0]*nTTemp[1])),2);
	if (js3==0){
	  x = (.5*pow(-1,js2)+.5)*bmLx[0] - 
	    pow(-1,js2)*fmod(js1,nTTemp[0])*bmDX[0] - offset[0];
	  y = floor(fmod(round(tc/DelT),(nTTemp[0]*nTTemp[1]))/nTTemp[0])*bmDX[1]-offset[1];
	} else {
	  y = (.5*pow(-1,js2)+.5)*bmLx[0] - 
	    pow(-1,js2)*fmod(js1,nTTemp[0])*bmDX[0] - offset[0];
	  x = floor(fmod(round(tc/DelT),(nTTemp[0]*nTTemp[1]))/nTTemp[0])*bmDX[1]-offset[1];
	}
	z = (floor((tc/DelT+DelT*1e-6)/(nTTemp[0]*nTTemp[1]))+zlaserOff)*bmDX[2] + _bp->height - 
	  ceil(offset[2]/_xyz->dX[2])*_xyz->dX[2];
	//if (_part->myid==0){if (j==0 & jt==0){std::cout<< tInd<<","<<x<<","<<y<<","<<z<<std::endl;}}
	lam = {pow(bmSTD[0],2.0)+2*alpha*(_xyz->time - tc), 
	       pow(bmSTD[1],2.0)+2*alpha*(_xyz->time - tc),
	       pow(bmSTD[2],2.0)+2*alpha*(_xyz->time - tc)};
	rij = pow(x-x0,2.0)/2/lam[0] +pow(y-y0,2.0)/2/lam[1] +
	  pow(z-z0,2.0)/2/lam[2];
	TempCurr[j] += Ci/pow(lam[0]*lam[1]*lam[2],.5)*exp(-rij);
      } // for (int j1...
    } // for (int j...
  } // if (patternID==4 ...
} // end SchwalbachTempCurr()

void TempField::SchwalbachTempCurr(double tcurr,std::vector<double> & TempOut )
{
  //computes temp field based on Schwalbach et al
  int Ntot = _part->nGhost+_part->ncellLoc, j1,j2,j3,iplay;
  double x0,y0,z0,rij,x,y,z,tc;
  std::vector<double> lam(3);
  if (patternID==0){
    for (int j=0;j<Ntot;++j){
      j3 = floor(_part->icellidLoc[j]/(_xyz->nX[0]*_xyz->nX[1]));
      j2 = floor( (_part->icellidLoc[j]- _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
      j1 = _part->icellidLoc[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
      x0 = (double(j1)+.5)*(_xyz->dX[0]);
      y0 = (double(j2)+.5)*(_xyz->dX[1]);
      z0 = (double(j3)+.5)*(_xyz->dX[2]);
      for (int jt=0;jt<nSource;++jt){
	// tc,x,y,z space-time location of source
	tc = tcurr  - jt*DelT;
	x = fmod((tc/DelT),(nTTemp[0]))*bmDX[0] - offset[0];
	y = floor(fmod((tc/DelT),(nTTemp[0]*nTTemp[1]))/nTTemp[0])*bmDX[1]-offset[1];
	z = (floor((tc/DelT)/(nTTemp[0]*nTTemp[1]))+zlaserOff)*bmDX[2] + _bp->height - 
	  ceil(offset[2]/_xyz->dX[2])*_xyz->dX[2];
	lam = {pow(bmSTD[0],2.0)+2*alpha*(tcurr - tc), 
	       pow(bmSTD[1],2.0)+2*alpha*(tcurr - tc),
	       pow(bmSTD[2],2.0)+2*alpha*(tcurr - tc)};
	rij = pow(x-x0,2.0)/2/lam[0] +pow(y-y0,2.0)/2/lam[1] +
	  pow(z-z0,2.0)/2/lam[2];
	TempOut[j] += Ci/pow(lam[0]*lam[1]*lam[2],.5)*exp(-rij);
      } // for (int j1...
    } // for (int j...
  } // if (patternID==0 ...
  if (patternID==1){
    int js1, js2;
    y = floor(fmod(round(_xyz->time/DelT),(nTTemp[0]*nTTemp[1]))/nTTemp[0])*bmDX[1]-offset[1];
    z = (floor((_xyz->time/DelT+DelT*1e-6)/(nTTemp[0]*nTTemp[1]))+zlaserOff)*bmDX[2] + _bp->height - 
      ceil(offset[2]/_xyz->dX[2])*_xyz->dX[2];
    for (int j=0;j<Ntot;++j){
      j3 = floor(_part->icellidLoc[j]/(_xyz->nX[0]*_xyz->nX[1]));
      j2 = floor( (_part->icellidLoc[j]- _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
      j1 = _part->icellidLoc[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
      x0 = (double(j1)+.5)*(_xyz->dX[0]);
      y0 = (double(j2)+.5)*(_xyz->dX[1]);
      z0 = (double(j3)+.5)*(_xyz->dX[2]);
      for (int jt=0;jt<nSource;++jt){
	// tc,x,y,z space-time location of source
	tc = jt*DelT;
	js1 = fmod( round(_xyz->time/DelT),nTTemp[0]*nTTemp[1]);
	js2 = fmod(floor(js1/nTTemp[0])+1,2);
	x = (.5*pow(-1,js2)+.5)*bmLx[0] - 
	  pow(-1,js2)*fmod(js1,nTTemp[0])*bmDX[0] + pow(-1,js2)*tc*bmV - 
	  (pow(-1,js2)+1)/2.0*shiftL[0] + pow(-1,js2)*offset[0] ;
	lam = {pow(bmSTD[0],2.0)+2*alpha*(tc), 
	       pow(bmSTD[1],2.0)+2*alpha*(tc),
	       pow(bmSTD[2],2.0)+2*alpha*(tc)};
	rij = pow(x-x0,2.0)/2/lam[0] +pow(y-y0,2.0)/2/lam[1] +
	  pow(z-z0,2.0)/2/lam[2];
	TempOut[j] += Ci/pow(lam[0]*lam[1]*lam[2],.5)*exp(-rij);
      } // for (int j1...
    } // for (int j...
  } // if (patternID==1 ...
  if (patternID==2){
    int js1, js2,js3;
    for (int j=0;j<Ntot;++j){
      j3 = floor(_part->icellidLoc[j]/(_xyz->nX[0]*_xyz->nX[1]));
      j2 = floor( (_part->icellidLoc[j]- _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
      j1 = _part->icellidLoc[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
      x0 = (double(j1)+.5)*(_xyz->dX[0]);
      y0 = (double(j2)+.5)*(_xyz->dX[1]);
      z0 = (double(j3)+.5)*(_xyz->dX[2]);
      for (int jt=0;jt<nSource;++jt){
	// tc,x,y,z space-time location of source
	tc = _xyz->time - jt*DelT;
	js1 = fmod( round(tc/DelT),nTTemp[0]*nTTemp[1]);
	js3 = fmod( floor( round(tc/DelT)/(nTTemp[0]*nTTemp[1]) ),2);
	js2 = fmod(floor(js1/nTTemp[0])+1+js3,2);
	x = (.5*pow(-1,js2)+.5)*bmLx[0] - 
	  pow(-1,js2)*fmod(js1,nTTemp[0])*bmDX[0] - offset[0];
	y = floor(fmod(round(tc/DelT),(nTTemp[0]*nTTemp[1]))/nTTemp[0])*bmDX[1]-offset[1];
	z = (floor((tc/DelT+DelT*1e-6)/(nTTemp[0]*nTTemp[1]))+zlaserOff)*bmDX[2] + _bp->height - 
	  ceil(offset[2]/_xyz->dX[2])*_xyz->dX[2];
	//if (_part->myid==0){if (j==0 & jt==0){std::cout<< tInd<<","<<x<<","<<y<<","<<z<<std::endl;}}
	lam = {pow(bmSTD[0],2.0)+2*alpha*(_xyz->time - tc), 
	       pow(bmSTD[1],2.0)+2*alpha*(_xyz->time - tc),
	       pow(bmSTD[2],2.0)+2*alpha*(_xyz->time - tc)};
	rij = pow(x-x0,2.0)/2/lam[0] +pow(y-y0,2.0)/2/lam[1] +
	  pow(z-z0,2.0)/2/lam[2];
	TempOut[j] += Ci/pow(lam[0]*lam[1]*lam[2],.5)*exp(-rij);
      } // for (int j1...
    } // for (int j...
  } // if (patternID==2 ...
  if (patternID==3){
    int js;
    for (int j=0;j<Ntot;++j){
      j3 = floor(_part->icellidLoc[j]/(_xyz->nX[0]*_xyz->nX[1]));
      j2 = floor( (_part->icellidLoc[j]- _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
      j1 = _part->icellidLoc[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
      x0 = (double(j1)+.5)*(_xyz->dX[0]);
      y0 = (double(j2)+.5)*(_xyz->dX[1]);
      z0 = (double(j3)+.5)*(_xyz->dX[2]);
      for (int jt=0;jt<nSource;++jt){
	// tc,x,y,z space-time location of source
	tc = _xyz->time - jt*DelT;
	js = fmod( floor( round(tc/DelT)/(nTTemp[0]*nTTemp[1])),2);
	if (js==0){
	  x = fmod(round(tc/DelT),(nTTemp[0]))*bmDX[0] - offset[0];
	  y = floor(fmod(round(tc/DelT),(nTTemp[0]*nTTemp[1]))/nTTemp[0])*bmDX[1]-offset[1];
	} else {
	  y = fmod(round(tc/DelT),(nTTemp[0]))*bmDX[0] - offset[0];
	  x = floor(fmod(round(tc/DelT),(nTTemp[0]*nTTemp[1]))/nTTemp[0])*bmDX[1]-offset[1];
	}
	z = (floor((tc/DelT+DelT*1e-6)/(nTTemp[0]*nTTemp[1]))+zlaserOff)*bmDX[2] + _bp->height - 
	  ceil(offset[2]/_xyz->dX[2])*_xyz->dX[2];
	lam = {pow(bmSTD[0],2.0)+2*alpha*(_xyz->time - tc), 
	       pow(bmSTD[1],2.0)+2*alpha*(_xyz->time - tc),
	       pow(bmSTD[2],2.0)+2*alpha*(_xyz->time - tc)};
	rij = pow(x-x0,2.0)/2/lam[0] +pow(y-y0,2.0)/2/lam[1] +
	  pow(z-z0,2.0)/2/lam[2];
	TempOut[j] += Ci/pow(lam[0]*lam[1]*lam[2],.5)*exp(-rij);
      } // for (int j1...
    } // for (int j...
  } // if (patternID==3 ...
  if (patternID==4){
    int js1, js2,js3;
    for (int j=0;j<Ntot;++j){
      j3 = floor(_part->icellidLoc[j]/(_xyz->nX[0]*_xyz->nX[1]));
      j2 = floor( (_part->icellidLoc[j]- _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
      j1 = _part->icellidLoc[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
      x0 = (double(j1)+.5)*(_xyz->dX[0]);
      y0 = (double(j2)+.5)*(_xyz->dX[1]);
      z0 = (double(j3)+.5)*(_xyz->dX[2]);
      for (int jt=0;jt<nSource;++jt){
	// tc,x,y,z space-time location of source
	tc = _xyz->time - jt*DelT;
	js1 = fmod( round(tc/DelT),nTTemp[0]*nTTemp[1]);
	js2 = fmod(floor(js1/nTTemp[0])+1,2);
	js3 = fmod( floor( round(tc/DelT)/(nTTemp[0]*nTTemp[1])),2);
	if (js3==0){
	  x = (.5*pow(-1,js2)+.5)*bmLx[0] - 
	    pow(-1,js2)*fmod(js1,nTTemp[0])*bmDX[0] - offset[0];
	  y = floor(fmod(round(tc/DelT),(nTTemp[0]*nTTemp[1]))/nTTemp[0])*bmDX[1]-offset[1];
	} else {
	  y = (.5*pow(-1,js2)+.5)*bmLx[0] - 
	    pow(-1,js2)*fmod(js1,nTTemp[0])*bmDX[0] - offset[0];
	  x = floor(fmod(round(tc/DelT),(nTTemp[0]*nTTemp[1]))/nTTemp[0])*bmDX[1]-offset[1];
	}
	z = (floor((tc/DelT+DelT*1e-6)/(nTTemp[0]*nTTemp[1]))+zlaserOff)*bmDX[2] + _bp->height - 
	  ceil(offset[2]/_xyz->dX[2])*_xyz->dX[2];
	lam = {pow(bmSTD[0],2.0)+2*alpha*(_xyz->time - tc), 
	       pow(bmSTD[1],2.0)+2*alpha*(_xyz->time - tc),
	       pow(bmSTD[2],2.0)+2*alpha*(_xyz->time - tc)};
	rij = pow(x-x0,2.0)/2/lam[0] +pow(y-y0,2.0)/2/lam[1] +
	  pow(z-z0,2.0)/2/lam[2];
	TempOut[j] += Ci/pow(lam[0]*lam[1]*lam[2],.5)*exp(-rij);
      } // for (int j1...
    } // for (int j...
  } // if (patternID==4 ...
} // end SchwalbachTempCurr()

