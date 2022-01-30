//  by Sompob Saralamba sompob@tropmedres.ac
// 
functions{
  
  real PDF_normal(real mu, real sig, real x){
    //return exp(normal_lpdf(x| mu, sig));
    return (1/(exp( ((x-mu)*(x-mu))/(2*sig*sig))*sqrt(2*pi())*sig));    
  }
  
 // based on White et al. 1992 
 // The effects of multiplication and synchronicity on the vascular distribution of parasites in falciparum malaria 
  vector InitDist(real N0, real mu, real sigma, real PMR){
    real sumP;
    vector[48] p;
    vector[48] ls;
    
    p=rep_vector(0,48);
    //based on White et al. model
    for(age in 1:48)
      p[age] = (PDF_normal(mu,sigma,age-48)/PMR)+PDF_normal(mu,sigma,age)+(PDF_normal(mu,sigma,age+48)*PMR);
    
    sumP = sum(p);
    for(age in 1:48)
      ls[age] = pow(10,N0)*p[age]/sumP;
    
    return ls;
  }
  
  // for the parasite aging 
  vector RotateRight(vector ls, real PMR){
    //size of the input vector is fixed at 48
    vector[48] tmpvec;
    tmpvec[1] = ls[48]*PMR;
    for(age in 2:48)
      tmpvec[age]=ls[age-1];
      
    return tmpvec;
  }
  
  // Hill function of drug concentration
  real Eff(real conc, real gamma, real ec50, real emax){
    return emax*pow(conc,gamma)/(pow(conc,gamma) + pow(ec50,gamma));
  }
  
  // sequestration function ; 
  // based on Saralamba et al.'s model 
  // probability of the parasites at age to be Ring stage
  real PRingfunc(int age, int a1,int a2){
    return ((age < a1)?1.0:exp(log(0.5)*(age-a1)/(a2-a1)));
  }

  // the function for counting the circulate parasites (rings)
  // based on Saralamba et al.'s model 
  real CountRing(vector ls,int a1, int a2){
    vector[48] pls;
    for(i in 1:48)
      pls[i]=PRingfunc(i,a1,a2);
    
    return sum(pls .* ls);
  }  
  
  // calculate the concentration @ time t
  real Concfn(real xm, real ym, real ke, real t){
    return ((t<=xm)?((ym/xm)*t):ym*exp(-ke*(t-xm)));
  }

  // Recovery function
  // assumed to be a sigmoid function with a lag time
  // rate - the damaged rate per hour
  // lag - the lag-time (in hours)
  // return recovery rate per hour
  real Recovfunc(real t, real rm, real rate, real lag){
    return rm/(1.0 +exp(-rate*(t-lag)));
  }

  // calculate the concentration and its effect from t=1 to MAX_ROW hours
  matrix DoseResponse(real xm, real ym, real ke, real gamma, real ec50, real emax, int everyH, int Ndrug){
    int MAX_ROW = 28*24; //4 weeks
    
    matrix[MAX_ROW,2] ec;
    int j;    //for counting the time interval between each dose
    int nd;   // number of doses
    int drugperiod;
    real conc;
    
    
    drugperiod = everyH*Ndrug;  
    ec=rep_matrix(0,MAX_ROW,2);   //maximum time that the model can run is MAX_ROW hrs
    
    conc=0.0;
    j = 0;  
		nd = 1;  
		for(i in 1:drugperiod){
		  if(j < everyH){
		    conc = Concfn(xm,ym,ke,j);
		    //print(conc);
		    ec[i,1]=conc;
			  ec[i,2]=Eff(ec[i,1],gamma,ec50,emax);
			}
			if((j==everyH)&&(nd<Ndrug)){
			  nd=nd+1;
			  conc=conc+Concfn(xm,ym,ke,j);
			  //print(conc);
			  ec[i,1]=conc;
			  ec[i,2]=Eff(ec[i,1],gamma,ec50,emax);
			  j=0;
			}
			j=j+1;
		}
		return ec;    
  }

  // calculate the concentration and the recovery rate for each parasite stage
  matrix DoseRecovery(real xm, real ym, real ke, real gamma, real ec50, real emax,vector recovermax, real lag ,int everyH, int Ndrug){
    int MAX_ROW = 28*24; //4 weeks
    
    matrix[MAX_ROW,5] ec;
    //matrix[240,2] ec;
    //int	i;    //time from 1 to drugperiod
    int j;    //for counting the time interval between each dose
    int nd;   // number of doses
    int drugperiod;
    real conc;
    
    
    drugperiod = everyH*Ndrug;  
    ec=rep_matrix(0,MAX_ROW,5);   //maximum time that the model can run is MAX_ROW hrs
    
    conc=0.0;
    j = 0;  
		nd = 1;  
		for(i in 1:drugperiod){
		  if(j < everyH){
		    conc = Concfn(xm,ym,ke,j);
		    //print(conc);
		    ec[i,1]=conc;
			  ec[i,2]=Eff(ec[i,1],gamma,ec50,emax);
        // Recovfunc(real t, real rm, real rate, real lag)
        ec[i,3]=Recovfunc(i, recovermax[1], Eff(ec[i,1],gamma,ec50,emax),lag);  
        ec[i,4]=Recovfunc(i, recovermax[2], Eff(ec[i,1],gamma,ec50,emax),lag);
        ec[i,5]=Recovfunc(i, recovermax[3], Eff(ec[i,1],gamma,ec50,emax),lag);
			}
			if((j==everyH)&&(nd<Ndrug)){
			  nd=nd+1;
			  conc=conc+Concfn(xm,ym,ke,j);
			  //print(conc);
			  ec[i,1]=conc;
			  ec[i,2]=Eff(ec[i,1],gamma,ec50,emax);
        ec[i,3]=Recovfunc(i, recovermax[1], Eff(ec[i,1],gamma,ec50,emax),lag);
        ec[i,4]=Recovfunc(i, recovermax[2], Eff(ec[i,1],gamma,ec50,emax),lag);
        ec[i,5]=Recovfunc(i, recovermax[3], Eff(ec[i,1],gamma,ec50,emax),lag);
			  j=0;
			}
			j=j+1;
		}
		return ec;    
  }



  //observed time for Pailin data
  //based on the data used in 
  int[] ObservedTimeP(int np){
    int tmp[31];  //obseved time in hours
    int time[np];
    
    	tmp[1]=0;
			tmp[2]=2;
			tmp[3]=4;
			tmp[4]=6;
			tmp[5]=8;
			tmp[6]=12;
			tmp[7]=18;
			tmp[8]=24;
			tmp[9]=30;
			tmp[10]=36;
			tmp[11]=42;
			tmp[12]=48;
			tmp[13]=54;
			tmp[14]=60;
			tmp[15]=66;
			tmp[16]=72;
			tmp[17]=78;
			tmp[18]=84;
			tmp[19]=90;
			tmp[20]=96;
			tmp[21]=102;
			tmp[22]=108;
			tmp[23]=114;
			tmp[24]=120;
			tmp[25]=126;
			tmp[26]=132;
			tmp[27]=138;
			tmp[28]=144;
			tmp[29]=150;
			tmp[30]=156;
			tmp[31]=162;
			
			for(i in 1:np){
			  time[i]=tmp[i];
			}
      return time;
  }
  
  //changing between the states
  matrix ChangeState(vector state1, vector state2, vector rate){
    vector[48] outfromstate1;
    vector[48] remaininstate1;
    vector[48] remaininstate2;
    
    outfromstate1 = state1 .* rate;
    //out must not greater than what has left
    for(i in 1:num_elements(state1)){
      if(outfromstate1[i] > state1[i]){
        outfromstate1[i] = state1[i];
      }
    }
    remaininstate1 = state1 - outfromstate1;
    remaininstate2 = state2 + outfromstate1;
    return append_col(remaininstate1,remaininstate2);
  }  
  
  
  //return the rate vector for all sensitive zones
  vector vecRates(vector SensitiveZones, vector rates){
    
    int RingsZone[2] = {6,26};
    int TrophZone[2] = {27,38};
    int SchiZone[2] = {39,44};
    vector[48] temp;
    
    temp = SensitiveZones;
    temp[RingsZone[1]:RingsZone[2]] = SensitiveZones[RingsZone[1]:RingsZone[2]] * rates[1];
    temp[TrophZone[1]:TrophZone[2]] = SensitiveZones[TrophZone[1]:TrophZone[2]] * rates[2];
    temp[SchiZone[1]:SchiZone[2]] = SensitiveZones[SchiZone[1]:SchiZone[2]] * rates[3];
    
    return temp;
  }
  

  // maing funciton of the model
  // the parasites from the susceptible state will be marked as injuted based on the Eff function
  // and then will be moved to the Injured state.
  // At the injurd state, the parasites can be retunred either to the susceptible state or the dead state.
  // 
  vector DM(real N0,real mu,real sigma,real PMR, real xm,real ym,real ke,int everyH,int Ndrug,real gamma, real ec50, real emax, vector recovermax, vector deadrates, real lag ,int runmax){
      
       int MAX_ROW = 28*24;    //4 weeks
      matrix[MAX_ROW,2] ecls; //for concentration and damage rates
      matrix[MAX_ROW,5] recls; //for concentration and recovery rates
    
      //matrix[240,2] ecls; //for concentration and response
      

      vector[48] suscp = rep_vector(0.0,48);
      vector[48] injured = rep_vector(0.0,48);
      vector[48] dead = rep_vector(0.0,48); 
      vector[48] SensitiveZones = rep_vector(0.0,48);
      vector[48] RecoverRateVector;
      vector[48] DeadRateVector;
      vector[48] InjuryRateVector;
      vector[3] recoverrates_tmp;

      matrix[48,2] tmpMatrix; // keeping the results from ChangeState funciton
      vector[(runmax+1)] outputvector;  //include t=0

      
      // recovery rates
      recls = DoseRecovery(xm,ym,ke,gamma,ec50,emax,recovermax,lag,everyH,Ndrug);

      SensitiveZones[6:44] = rep_vector(1.0,39);
      
      // rate vectors
      // RecoverRateVector = vecRates(SensitiveZones, recoverates);
      DeadRateVector = vecRates(SensitiveZones, deadrates);
      
      //percentage of damaged parasites
			ecls = DoseResponse(xm,ym,ke,gamma,ec50,emax,everyH,Ndrug);
			
			//initial age distribution of the parasites
      suscp = InitDist(N0,mu,sigma,PMR);
      
      outputvector[1]=log10(CountRing(suscp,11,14)); //count at t=0 ; before taking drug
      
      for(t in 1:runmax){
        recoverrates_tmp = to_vector(recls[t,3:5]);   
        RecoverRateVector = vecRates(SensitiveZones, recoverrates_tmp);
        
        //from susceptible -> injured (Injured)
        InjuryRateVector = SensitiveZones*(ecls[t,2]*0.01);
        tmpMatrix = ChangeState(suscp,injured, InjuryRateVector);
        suscp = col(tmpMatrix,1);
        injured = col(tmpMatrix,2);
        
        //from injured -> susceptible (Recover)
        tmpMatrix = ChangeState(injured, suscp, RecoverRateVector);
        injured = col(tmpMatrix,1);
        suscp = col(tmpMatrix,2);
        
        //from injured -> dead (Dead)
        tmpMatrix = ChangeState(injured, dead, DeadRateVector);
        injured = col(tmpMatrix,1);
        dead = col(tmpMatrix,2);

        //aging
        suscp = RotateRight(suscp, PMR);
        injured = RotateRight(injured, PMR);
        suscp[1] = suscp[1] + injured[1];
        injured[1] = 0.0;

        // //observed  parasites  
        outputvector[t+1] = log10(CountRing(suscp + injured,11,14)); //count at time t
 
      }

      return outputvector;
    }

    //find the clearance time from the output vector
    real ClearanceTime(vector countvector, real deteclim){
      real ct;
      int veclen;
      int indx;
      int run;
      run = 1;
      veclen = num_elements(countvector);
      indx=1;
      ct = 240;
      
      while(indx <= veclen && run==1){
        if(countvector[indx] <= deteclim){
          ct = indx;
          run = 0;
        }
        indx=indx+1;
      }
      return ct;
    }
  
  
} //end function block

data{
  int Nid;  // number of patients
  int NPPARA[Nid]; //number of data points
  real y_obs[Nid,22];  
  int Nparms;
  
  real xm[Nid];  //time at max conc
  real ym[Nid];  //maximum conc
  real ke[Nid];   //elimination rate of the drug
  real DetecLim[Nid];  //detection limit
  vector[Nparms] lb;  //lower bound
  vector[Nparms] ub;  //upper bound
}
parameters{

  real<lower=lb[1],upper=ub[1]> N0_mean;
  real<lower=lb[2],upper=ub[2]> mu_mean;
  real<lower=lb[3],upper=ub[3]> sig_mean;
  real<lower=lb[4],upper=ub[4]> pmr_mean;
  real<lower=lb[5],upper=ub[5]> gamma_mean;
  real<lower=lb[6],upper=ub[6]> ec50_mean;
  real<lower=lb[7],upper=ub[7]> deadrateR_mean;
  real<lower=lb[8],upper=ub[8]> deadrateT_mean;
  real<lower=lb[9],upper=ub[9]> deadrateS_mean;  
  real<lower=lb[10],upper=ub[10]> recoverrateR_mean;
  real<lower=lb[11],upper=ub[11]> recoverrateT_mean;
  real<lower=lb[12],upper=ub[12]> recoverrateS_mean;
  real<lower=lb[13],upper=ub[13]> recoverlag_mean;
  
  cholesky_factor_corr[Nparms] L;  
  vector<lower=0>[Nparms]  omega;
  matrix[Nparms,Nid] etaStd;
  real<lower=0> SD[Nid];

}
transformed parameters {
  // real y[Nid,22];
  matrix[Nparms,Nid] phiInd;
  matrix<lower=0>[Nparms,Nid] thetaInd;
  vector[Nparms] phiPop;
  vector[Nparms] thetaPop;

  matrix[Nparms,Nparms] Cor;
  matrix[Nparms,Nparms] Cov;

  thetaPop[1] = N0_mean;
  thetaPop[2] = mu_mean;
  thetaPop[3] = sig_mean;
  thetaPop[4] = pmr_mean;
  thetaPop[5] = gamma_mean;
  thetaPop[6] = ec50_mean;
  thetaPop[7] = deadrateR_mean;
  thetaPop[8] = deadrateT_mean;
  thetaPop[9] = deadrateS_mean;	
  thetaPop[10] = recoverrateR_mean;
  thetaPop[11] = recoverrateT_mean;
  thetaPop[12] = recoverrateS_mean;
  thetaPop[13] = recoverlag_mean;

  phiPop = log((thetaPop - lb) ./ (ub - thetaPop));
  phiInd = rep_matrix(phiPop, Nid) + diag_pre_multiply(omega, L*etaStd);
  thetaInd = rep_matrix(lb,Nid) .* (1 - inv_logit(phiInd)) + rep_matrix(ub, Nid) .* inv_logit(phiInd);
  
  Cor = L * L';
  Cov = quad_form_diag(Cor,omega);

}

model
{
  matrix[Nid,22] mod_pred;
  matrix[Nid,241] mod_tmp;
  matrix[Nid,3] RecoverRate;
  matrix[Nid,3] DeadRate;
  
  int timewant[Nid, 22];
  
  L ~ lkj_corr_cholesky(2);
  omega ~ normal(0,1);
  to_row_vector(etaStd) ~ normal(0,1);
  SD ~ cauchy(0,1);

  for(id in 1:Nid){

    timewant[id,1:NPPARA[id]] = ObservedTimeP(NPPARA[id]);   //using Painlin observed times
            
    // Dead Rates
    DeadRate[id,1] = thetaInd[7,id];
    DeadRate[id,2] = thetaInd[8,id];
    DeadRate[id,3] = thetaInd[9,id];	
    // Recovery Rates
    RecoverRate[id,1] = thetaInd[10,id];
    RecoverRate[id,2] = thetaInd[11,id];
    RecoverRate[id,3] = thetaInd[12,id];

	//the proposed model output
    mod_tmp[id,1:241] = to_row_vector(DM(thetaInd[1,id],thetaInd[2,id],thetaInd[3,id],thetaInd[4,id],xm[id],ym[id],ke[id],24,7,thetaInd[5,id],thetaInd[6,id],99.99,
                          to_vector(RecoverRate[id,1:3]),to_vector(DeadRate[id,1:3]),thetaInd[13,id],240));
        
    for(i in 1:NPPARA[id]){
      mod_pred[id, i] = mod_tmp[id,timewant[id,i]+1];
    }
    // print(mod_pred[id,1:NPPARA[id]]);
    y_obs[id,1:NPPARA[id]] ~ normal(mod_pred[id,1:NPPARA[id]], SD[id]);

  }
  
}
generated quantities{
  matrix[Nid,22] y_pred;
  matrix[Nid,22] ySplit_pred;
  matrix[Nid,22] log_lik;
  matrix[Nid,22] mod_pred;
  matrix[Nid,22] modSplit_pred;
  matrix[Nid,240+1] mod_tmp;  //mod_24_7
  matrix[Nid,240+1] mod_tmpSplit; //mod_12_14

  matrix[Nid,240+1] mod_24_1;
  matrix[Nid,240+1] mod_24_3;
  matrix[Nid,240+1] mod_24_5;
  
  matrix[Nid,240+1] mod_12_2;
  matrix[Nid,240+1] mod_12_6;
  matrix[Nid,240+1] mod_12_10;
  
  

  matrix[Nid,3] RecoverRate;
  matrix[Nid,3] DeadRate;
  real AS7CT[Nid];
  real SplitCT[Nid];
  int timewant[Nid,22];
  
  for(id in 1:Nid){

    timewant[id, 1:NPPARA[id]] = ObservedTimeP(NPPARA[id]);
    // Dead Rates
    DeadRate[id,1] = thetaInd[7,id];
    DeadRate[id,2] = thetaInd[8,id];
    DeadRate[id,3] = thetaInd[9,id];	
    // Recovery Rates
    RecoverRate[id,1] = thetaInd[10,id];
    RecoverRate[id,2] = thetaInd[11,id];
    RecoverRate[id,3] = thetaInd[12,id];
  
    //once a day
    //every 24 Ndrug = 7
    mod_tmp[id,1:241] = to_row_vector(DM(thetaInd[1,id],thetaInd[2,id],thetaInd[3,id],thetaInd[4,id],xm[id],ym[id],ke[id],24,7,thetaInd[5,id],thetaInd[6,id],99.99,
                          to_vector(RecoverRate[id,1:3]),to_vector(DeadRate[id,1:3]),thetaInd[13,id],240));
    //twice a day
    //every 12 Ndrug = 14
    mod_tmpSplit[id,1:241] = to_row_vector(DM(thetaInd[1,id],thetaInd[2,id],thetaInd[3,id],thetaInd[4,id],xm[id],ym[id],ke[id],12,14,thetaInd[5,id],thetaInd[6,id],99.99,
                          to_vector(RecoverRate[id,1:3]),to_vector(DeadRate[id,1:3]),thetaInd[13,id],240));

    AS7CT[id]=ClearanceTime(to_vector(mod_tmp[id,1:241]),DetecLim[id]);
    SplitCT[id]=ClearanceTime(to_vector(mod_tmpSplit[id,1:241]),DetecLim[id]);

    for(i in 1:NPPARA[id]){
      mod_pred[id, i] = mod_tmp[id,timewant[id,i]+1];
      modSplit_pred[id, i] = mod_tmpSplit[id,timewant[id,i]+1];
      y_pred[id,i] = normal_rng(mod_pred[id,i], SD[id]);
      ySplit_pred[id,i] = normal_rng(modSplit_pred[id,i], SD[id]);

      log_lik[id,i] = normal_lpdf(y_pred[id,i]| mod_pred[id,i],SD[id]);
    }

  }
  
}
