#include "TOV_library.h"

//TOV equations
double dPdr(double r, double m, double P, double rho){
    if(P<0) {printf("P<0 in dP/dr %f %f\n",P,rho);exit(1);}
    return -Ggrav*(P+rho)*(m+4*M_PI*r*r*r*P)/(r*(r-2*Ggrav*m/c2));
}
double dmdr(double r, double rho){
    return 4*M_PI*r*r*rho;
}
double dphidr(double dPdr, double P, double rho){
    return -dPdr/(P+rho);
}
//Tidal deformability equations
double dHdr(double beta){
    return beta;
}
double dbetadr(double r, double m, double P, double H, double beta, double dens, double ddens){
    double Gc2=Ggrav/c2;
    double factor=1.0-2*m*Gc2/r; factor=2.0/factor; 
    double term1=-2*M_PI*Gc2*(5*dens+9*P+(dens+P)*ddens)+3.0/r/r+factor*Gc2*(m/r/r+4*M_PI*r*P)*Gc2*(m/r/r+4*M_PI*r*P);
    term1*=factor*H;
    double term2=-1+Gc2*m/r+Gc2*2*M_PI*r*r*(dens-P);
    term2*=factor*beta/r;
    return term1+term2;
}
//Derivatives related to EoS
double ddensdPPP(double P, struct PP eos){
  int i,j;
  double rho=densityPP(P, eos);
  for(i=eos.Nrmd-1; i>=0;i--){
      if(rho>=eos.rmd[i]){
          j=eos.pieces[i+1];
          return (1.0+epsilonpoly(rho,eos.K[j],eos.G[j],eos.cteps[j],eos.rmd[i])+eos.K[j]*pow(rho, eos.G[j]-1))/(eos.K[j]*eos.G[j]*pow(rho, eos.G[j]-1)); 
      }
  }
  j=eos.pieces[0];
  return (1.0+epsilonpoly(rho,eos.K[j],eos.G[j],eos.cteps[j],0.0)+eos.K[j]*pow(rho, eos.G[j]-1))/(eos.K[j]*eos.G[j]*pow(rho, eos.G[j]-1)); 

}
double ddensdPTASPP(double P, struct TASPP eos){
  double rho;
  rho=densityTASPP(P, eos);
  return 1.0/csTASPP(rho, eos);
}
double ddensdPTab(double P, struct Tab eos){
    double y1,y2, PL, PR;
    int i=locateQuantity(eos.Ptab, P, eos.lenEos);

    if(i==eos.lenEos-1){
        y1=(eos.rhotab[i]*(1+eos.epstab[i])-eos.rhotab[i-1]*(1+eos.epstab[i-1]))/(eos.Ptab[i]-eos.Ptab[i-1]);
        y2=(eos.rhotab[i-1]*(1+eos.epstab[i-1])-eos.rhotab[i-2]*(1+eos.epstab[i-2]))/(eos.Ptab[i-1]-eos.Ptab[i-2]);
        PL=eos.Ptab[i-2];
        PR=eos.Ptab[i];
    }
    else if(i==0){
        y1=(eos.rhotab[i+2]*(1+eos.epstab[i+2])-eos.rhotab[i+1]*(1+eos.epstab[i+1]))/(eos.Ptab[i+2]-eos.Ptab[i+1]);
        y2=(eos.rhotab[i+1]*(1+eos.epstab[i+1])-eos.rhotab[i]*(1+eos.epstab[i]))/(eos.Ptab[i+1]-eos.Ptab[i]);
        PL=eos.Ptab[i];
        PR=eos.Ptab[i+2];
    }
    else{
        y1=(eos.rhotab[i+1]*(1+eos.epstab[i+1])-eos.rhotab[i]*(1+eos.epstab[i]))/(eos.Ptab[i+1]-eos.Ptab[i]);
        y2=(eos.rhotab[i]*(1+eos.epstab[i])-eos.rhotab[i-1]*(1+eos.epstab[i-1]))/(eos.Ptab[i]-eos.Ptab[i-1]);
        PL=eos.Ptab[i-1];
        PR=eos.Ptab[i+1];
    }

    y1=log10(y1); y2=log10(y2); PL=log10(PL); PR=log10(PR);
    return pow(10,y2+(y1-y2)/(PR-PL)*(log10(P)-PL));

}

double solveStarMRP_PP(double *Mmax, double rho_c, int saveRadialProfile, int checkAccuracyMetric, double *Maxerror, struct PP eos){
    double *k1, *k2, *k3, *k4, *sol; FILE *file; FILE *temp;
    if((k1=(double*) malloc(3*sizeof(double)))==NULL){printf("error reserving memory\n"); exit(1);}
    if((k2=(double*) malloc(3*sizeof(double)))==NULL){printf("error reserving memory\n"); exit(1);}
    if((k3=(double*) malloc(3*sizeof(double)))==NULL){printf("error reserving memory\n"); exit(1);}
    if((k4=(double*) malloc(3*sizeof(double)))==NULL){printf("error reserving memory\n"); exit(1);}
    if((sol=(double*) malloc(3*sizeof(double)))==NULL){printf("error reserving memory\n"); exit(1);} 
    //Pressure, mass, phi

    if(saveRadialProfile){
        char buffer[100];
        char filename[104]="PPprofile_";
        sprintf(buffer, "%g", rho_c);
        strcat(filename, buffer);
        strcat(filename, ".txt");
        file=fopen(filename, "w");
        fprintf(file, "#1-r(km)\t2-m(Ms)\t3-rho(g/cm3)\t4-phi(not rescaled)\n");
    }
    if(checkAccuracyMetric){
        temp=fopen("temp_file.txt", "w");
    }

    double dr=1.0, rho, r, mprev, eps;
    //Initial condition
    r=dr; rho=rho_c; eps=epsilonPP(rho_c,eos);
    sol[0]=pressurePP(rho_c, eos); 
    sol[1]=4*M_PI*r*r*r*rho*(1+eps);  
    sol[2]=0;

    //First iteration out of the loop to deal with the singularity
    k1[0]=0;
    k1[1]=0;
    k1[2]=0;

    k2[0]=dPdr(r+dr*0.5, sol[1], sol[0], rho*(1.0+eps))/c2;
    k2[1]=dmdr(r+dr*0.5, rho*(1.0+eps));
    k2[2]=dphidr(k2[0], sol[0], rho*(1.0+eps));

    rho=densityPP(sol[0]+0.5*dr*k2[0],eos); eps=epsilonPP(rho,eos);
    k3[0]=dPdr(r+dr*0.5, sol[1]+dr*0.5*k2[1], sol[0]+dr*0.5*k2[0], rho*(1.0+eps))/c2;
    k3[1]=dmdr(r+dr*0.5, rho*(1.0+eps));
    k3[2]=dphidr(k3[0], sol[0]+dr*0.5*k2[0], rho*(1.0+eps));

    rho=densityPP(sol[0]+dr*k3[0],eos); eps=epsilonPP(rho,eos);
    k4[0]=dPdr(r+dr, sol[1]+dr*k3[1], sol[0]+dr*k3[0], rho*(1.0+eps))/c2;
    k4[1]=dmdr(r+dr, rho*(1.0+eps));
    k4[2]=dphidr(k4[0], sol[0]+dr*k3[0], rho*(1.0+eps));

    sol[0]+=dr*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6;
    sol[1]+=dr*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6;
    sol[2]+=dr*(k1[2]+2*k2[2]+2*k3[2]+k4[2])/6;

    r+=dr;
    do{
        rho=densityPP(sol[0],eos); eps=epsilonPP(rho,eos);
        k1[0]=dPdr(r,sol[1],sol[0],rho*(1.0+eps))/c2;
        k1[1]=dmdr(r,rho*(1.0+eps));
        k1[2]=dphidr(k1[0],sol[0],rho*(1.0+eps));

        rho=densityPP(sol[0]+0.5*dr*k1[0],eos); eps=epsilonPP(rho,eos);
        k2[0]=dPdr(r+dr*0.5, sol[1]+0.5*dr*k1[1], sol[0]+0.5*dr*k1[0], rho*(1.0+eps))/c2;
        k2[1]=dmdr(r+dr*0.5, rho*(1.0+eps));
        k2[2]=dphidr(k2[0], sol[0]+0.5*dr*k1[0], rho*(1.0+eps));

        rho=densityPP(sol[0]+0.5*dr*k2[0],eos); eps=epsilonPP(rho,eos);
        k3[0]=dPdr(r+dr*0.5, sol[1]+dr*0.5*k2[1], sol[0]+dr*0.5*k2[0], rho*(1.0+eps))/c2;
        k3[1]=dmdr(r+dr*0.5, rho*(1.0+eps));
        k3[2]=dphidr(k3[0], sol[0]+dr*0.5*k2[0], rho*(1.0+eps));

        rho=densityPP(sol[0]+dr*k3[0],eos); eps=epsilonPP(rho,eos);
        k4[0]=dPdr(r+dr, sol[1]+dr*k3[1], sol[0]+dr*k3[0], rho*(1.0+eps))/c2;
        k4[1]=dmdr(r+dr, rho*(1.0+eps));
        k4[2]=dphidr(k4[0], sol[0]+dr*k3[0], rho*(1.0+eps));

        sol[0]+=dr*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6; mprev=sol[1];
        sol[1]+=dr*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6;
        sol[2]+=dr*(k1[2]+2*k2[2]+2*k3[2]+k4[2])/6;

        if(saveRadialProfile){
            fprintf(file, "%e %.15e %.15e %.15e\n",r*1e-5, sol[1]*0.001/1.98847e30, densityPP(sol[0],eos), sol[2]);
        }
        if(checkAccuracyMetric){
            fprintf(temp, "%.15e\t%.15e\n",densityPP(sol[0],eos), sol[2]);
        }
       
        if(r>2 && (sol[1]-mprev)/mprev<1e-12) break; //stopping condition by saturation of mass
        r+=dr;
        
    }while (sol[0]>5e4); //Safety stopping condition: very low pressure
    (*Mmax)=sol[1]*0.001/1.98847e30;
    if(saveRadialProfile) fclose(file);
    if(checkAccuracyMetric){ 
        fclose(temp);
        
        double R, Rsch, phiR, phiSch, renorm;
        R=r;
        Rsch=2*Ggrav*sol[1]/c2;
        phiR=sol[2];
        phiSch=0.5*log(1-Rsch/R);
        renorm=-phiR+phiSch;
   
        const unsigned MAX_LENGTH = 256;
        char buffer[MAX_LENGTH], line2[MAX_LENGTH];
        int element=0;
        char *token;
        double phi, P, quantity;
        double maxCte=0, minCte=1e30;
        FILE *fp = fopen("temp_file.txt", "r");
        while (fgets(buffer, MAX_LENGTH, fp)){
            memcpy(line2, buffer, MAX_LENGTH);
            token = strtok(line2, "\t");
            element=0;
            while( token != NULL ) {
                if(element==0) rho=strtod(token,NULL);
                if(element==1) phi=strtod(token,NULL);
                token = strtok(NULL, "\t");
                element++;
            }
            phi+=renorm;
            eps=epsilonPP(rho, eos);
            P=pressurePP(rho, eos);
            quantity=(rho*(1+eps)+P)/rho*exp(phi);

            if (quantity>maxCte)
                maxCte=quantity;
            else if (quantity<minCte)
                minCte=quantity;
        }
        fclose(fp);
        remove("temp_file.txt");
        if(saveRadialProfile){
            printf("An accurate integration keeps quantity mu_B*exp(phi) constant\n");
            printf("The maximum absolute difference in the constant quantity was %.15e\n", maxCte-minCte);
            printf("The maximum relative difference in the constant quantity was %.15e\n", (maxCte-minCte)/minCte); 
        }
        (*Maxerror)=(maxCte-minCte)/minCte;
    }

    free(k1); free(k2); free(k3); free(k4); free(sol);
    return r*1e-5;
}

double solveStarMRPL_PP(double *Mmax, double *LambdaNS, double rho_c, int saveRadialProfile, int checkAccuracyMetric, double *Maxerror, struct PP eos){
    double *k1, *k2, *k3, *k4, *sol; FILE *file; FILE *temp;
    if((k1=(double*) malloc(5*sizeof(double)))==NULL){printf("error reserving memory\n"); exit(1);}
    if((k2=(double*) malloc(5*sizeof(double)))==NULL){printf("error reserving memory\n"); exit(1);}
    if((k3=(double*) malloc(5*sizeof(double)))==NULL){printf("error reserving memory\n"); exit(1);}
    if((k4=(double*) malloc(5*sizeof(double)))==NULL){printf("error reserving memory\n"); exit(1);}
    if((sol=(double*) malloc(5*sizeof(double)))==NULL){printf("error reserving memory\n"); exit(1);} 
    //Pressure, mass, phi, H, beta

    if(saveRadialProfile){
        char buffer[100];
        char filename[104]="PPprofile_";
        sprintf(buffer, "%g", rho_c);
        strcat(filename, buffer);
        strcat(filename, ".txt");
        file=fopen(filename, "w");
        fprintf(file, "#1-r(km)\t2-m(Ms)\t3-rho(g/cm3)\t4-phi(not rescaled)\t5-Lambda\n");
    }
    if(checkAccuracyMetric){
        temp=fopen("temp_file.txt", "w");
    }

    double R, M, H, beta, y, C, C5, lovek2, Lambda;
    double dr=1.0, rho, r, mprev, eps;
    //Initial condition
    r=dr; rho=rho_c; eps=epsilonPP(rho_c,eos);
    sol[0]=pressurePP(rho_c, eos); 
    sol[1]=4*M_PI*r*r*r*rho*(1+eps);  
    sol[2]=0;
    sol[3]=r*r; 
    sol[4]=2*r;

    //First iteration out of the loop to deal with the singularity
    k1[0]=0;
    k1[1]=0;
    k1[2]=0;
    k1[3]=dHdr(sol[4]); 
    k1[4]=dbetadr(r,sol[1],sol[0],sol[3],sol[4],rho*(1.0+eps),ddensdPPP(sol[0],eos));
    

    k2[0]=dPdr(r+dr*0.5, sol[1], sol[0], rho*(1.0+eps))/c2;
    k2[1]=dmdr(r+dr*0.5, rho*(1.0+eps));
    k2[2]=dphidr(k2[0], sol[0], rho*(1.0+eps));
    k2[3]=dHdr(sol[4]);
    k2[4]=dbetadr(r+dr*0.5, sol[1],sol[0],sol[3],sol[4],rho*(1+eps),ddensdPPP(sol[0],eos));

    rho=densityPP(sol[0]+0.5*dr*k2[0],eos); eps=epsilonPP(rho,eos);
    k3[0]=dPdr(r+dr*0.5, sol[1]+dr*0.5*k2[1], sol[0]+dr*0.5*k2[0], rho*(1.0+eps))/c2;
    k3[1]=dmdr(r+dr*0.5, rho*(1.0+eps));
    k3[2]=dphidr(k3[0], sol[0]+dr*0.5*k2[0], rho*(1.0+eps));
    k3[3]=dHdr(sol[4]+dr*0.5*k2[4]);
    k3[4]=dbetadr(r+dr*0.5, sol[1]+dr*0.5*k2[1],sol[0]+dr*0.5*k2[0],sol[3]+dr*0.5*k2[3],sol[4]+dr*0.5*k2[4],rho*(1+eps),ddensdPPP(sol[0]+dr*0.5*k2[0], eos));

    rho=densityPP(sol[0]+dr*k3[0],eos); eps=epsilonPP(rho,eos);
    k4[0]=dPdr(r+dr, sol[1]+dr*k3[1], sol[0]+dr*k3[0], rho*(1.0+eps))/c2;
    k4[1]=dmdr(r+dr, rho*(1.0+eps));
    k4[2]=dphidr(k4[0], sol[0]+dr*k3[0], rho*(1.0+eps));
    k4[3]=dHdr(sol[4]+dr*k3[4]);
    k4[4]=dbetadr(r+dr, sol[1]+dr*k3[1], sol[0]+dr*k3[0],sol[3]+dr*k3[3], sol[4]+dr*k3[4], rho*(1+eps),ddensdPPP(sol[0]+dr*k3[0],eos));

    sol[0]+=dr*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6;
    sol[1]+=dr*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6;
    sol[2]+=dr*(k1[2]+2*k2[2]+2*k3[2]+k4[2])/6;
    sol[3]+=dr*(k1[3]+2*k2[3]+2*k3[3]+k4[3])/6;
    sol[4]+=dr*(k1[4]+2*k2[4]+2*k3[4]+k4[4])/6;

    r+=dr;
    do{
        rho=densityPP(sol[0],eos); eps=epsilonPP(rho,eos);
        k1[0]=dPdr(r,sol[1],sol[0],rho*(1.0+eps))/c2;
        k1[1]=dmdr(r,rho*(1.0+eps));
        k1[2]=dphidr(k1[0],sol[0],rho*(1.0+eps));
        k1[3]=dHdr(sol[4]); 
        k1[4]=dbetadr(r,sol[1],sol[0],sol[3],sol[4],rho*(1.0+eps),ddensdPPP(sol[0],eos));

        rho=densityPP(sol[0]+0.5*dr*k1[0],eos); eps=epsilonPP(rho,eos);
        k2[0]=dPdr(r+dr*0.5, sol[1]+0.5*dr*k1[1], sol[0]+0.5*dr*k1[0], rho*(1.0+eps))/c2;
        k2[1]=dmdr(r+dr*0.5, rho*(1.0+eps));
        k2[2]=dphidr(k2[0], sol[0]+0.5*dr*k1[0], rho*(1.0+eps));
        k2[3]=dHdr(sol[4]);
        k2[4]=dbetadr(r+dr*0.5, sol[1],sol[0],sol[3],sol[4],rho*(1+eps),ddensdPPP(sol[0],eos));

        rho=densityPP(sol[0]+0.5*dr*k2[0],eos); eps=epsilonPP(rho,eos);
        k3[0]=dPdr(r+dr*0.5, sol[1]+dr*0.5*k2[1], sol[0]+dr*0.5*k2[0], rho*(1.0+eps))/c2;
        k3[1]=dmdr(r+dr*0.5, rho*(1.0+eps));
        k3[2]=dphidr(k3[0], sol[0]+dr*0.5*k2[0], rho*(1.0+eps));
        k3[3]=dHdr(sol[4]+dr*0.5*k2[4]);
        k3[4]=dbetadr(r+dr*0.5, sol[1]+dr*0.5*k2[1],sol[0]+dr*0.5*k2[0],sol[3]+dr*0.5*k2[3],sol[4]+dr*0.5*k2[4],rho*(1+eps),ddensdPPP(sol[0]+dr*0.5*k2[0], eos));

        rho=densityPP(sol[0]+dr*k3[0],eos); eps=epsilonPP(rho,eos);
        k4[0]=dPdr(r+dr, sol[1]+dr*k3[1], sol[0]+dr*k3[0], rho*(1.0+eps))/c2;
        k4[1]=dmdr(r+dr, rho*(1.0+eps));
        k4[2]=dphidr(k4[0], sol[0]+dr*k3[0], rho*(1.0+eps));
        k4[3]=dHdr(sol[4]+dr*k3[4]);
        k4[4]=dbetadr(r+dr, sol[1]+dr*k3[1], sol[0]+dr*k3[0],sol[3]+dr*k3[3], sol[4]+dr*k3[4], rho*(1+eps),ddensdPPP(sol[0]+dr*k3[0],eos));

        sol[0]+=dr*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6; mprev=sol[1];
        sol[1]+=dr*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6;
        sol[2]+=dr*(k1[2]+2*k2[2]+2*k3[2]+k4[2])/6;
        sol[3]+=dr*(k1[3]+2*k2[3]+2*k3[3]+k4[3])/6;
        sol[4]+=dr*(k1[4]+2*k2[4]+2*k3[4]+k4[4])/6;

        if(saveRadialProfile){     
            R=r; M=sol[1]; H=sol[3]; beta=sol[4];  rho=densityPP(sol[0],eos); eps=epsilonPP(rho,eos);
            y=R*beta/H-4*M_PI*R*R*R*rho*(1+eps)/M; //added correction for nonzero density at surface
            C=M/R*Ggrav/c2; C5=C*C*C*C*C;
            lovek2=8.0/5*C5*(1-2*C)*(1-2*C)*(2+2*C*(y-1)-y)/((2*C*(6-3*y+3*C*(5*y-8)))+4*C*C*C*(13-11*y+C*(3*y-2)+2*C*C*(1+y))+3*(1-2*C)*(1-2*C)*(2-y+2*C*(y-1))*log(1-2*C));   
            Lambda=2.0/3*lovek2/C5;
            fprintf(file, "%e %.15e %.15e %.15e %.15e\n",r*1e-5, sol[1]*0.001/1.98847e30, densityPP(sol[0],eos), sol[2], Lambda);
        }
        if(checkAccuracyMetric){
            fprintf(temp, "%.15e\t%.15e\n",densityPP(sol[0],eos), sol[2]);
        }
       
        if(r>2 && (sol[1]-mprev)/mprev<1e-12) break; //stopping condition by saturation of mass
        r+=dr;
        
    }while (sol[0]>5e4); //Safety stopping condition: very low pressure
    (*Mmax)=sol[1]*0.001/1.98847e30;
    R=r; M=sol[1]; H=sol[3]; beta=sol[4];  rho=densityPP(sol[0],eos); eps=epsilonPP(rho,eos);
    y=R*beta/H-4*M_PI*R*R*R*rho*(1+eps)/M; //added correction for nonzero density at surface
    C=M/R*Ggrav/c2; C5=C*C*C*C*C;
    lovek2=8.0/5*C5*(1-2*C)*(1-2*C)*(2+2*C*(y-1)-y)/((2*C*(6-3*y+3*C*(5*y-8)))+4*C*C*C*(13-11*y+C*(3*y-2)+2*C*C*(1+y))+3*(1-2*C)*(1-2*C)*(2-y+2*C*(y-1))*log(1-2*C));   
    Lambda=2.0/3*lovek2/C5;
    (*LambdaNS)=Lambda;

    if(saveRadialProfile) fclose(file);

    if(checkAccuracyMetric){ 
        fclose(temp);
        
        double R, Rsch, phiR, phiSch, renorm;
        R=r;
        Rsch=2*Ggrav*sol[1]/c2;
        phiR=sol[2];
        phiSch=0.5*log(1-Rsch/R);
        renorm=-phiR+phiSch;
   
        const unsigned MAX_LENGTH = 256;
        char buffer[MAX_LENGTH], line2[MAX_LENGTH];
        int element=0;
        char *token;
        double phi, P, quantity;
        double maxCte=0, minCte=1e30;
        FILE *fp = fopen("temp_file.txt", "r");
        while (fgets(buffer, MAX_LENGTH, fp)){
            memcpy(line2, buffer, MAX_LENGTH);
            token = strtok(line2, "\t");
            element=0;
            while( token != NULL ) {
                if(element==0) rho=strtod(token,NULL);
                if(element==1) phi=strtod(token,NULL);
                token = strtok(NULL, "\t");
                element++;
            }
            phi+=renorm;
            eps=epsilonPP(rho, eos);
            P=pressurePP(rho, eos);
            quantity=(rho*(1+eps)+P)/rho*exp(phi);

            if (quantity>maxCte)
                maxCte=quantity;
            else if (quantity<minCte)
                minCte=quantity;
        }
        fclose(fp);
        remove("temp_file.txt");
        if(saveRadialProfile){
            printf("An accurate integration keeps quantity mu_B*exp(phi) constant\n");
            printf("The maximum absolute difference in the constant quantity was %.15e\n", maxCte-minCte);
            printf("The maximum relative difference in the constant quantity was %.15e\n", (maxCte-minCte)/minCte); 
        }
        (*Maxerror)=(maxCte-minCte)/minCte;
    }

    free(k1); free(k2); free(k3); free(k4); free(sol);
    return r*1e-5;
}

double solveStarMRP_TASPP(double *Mmax, double rho_c, int saveRadialProfile, int checkAccuracyMetric, double *Maxerror, struct TASPP eos){
    double *k1, *k2, *k3, *k4, *sol; FILE *file; FILE *temp;
    if((k1=(double*) malloc(3*sizeof(double)))==NULL){printf("error reserving memory\n"); exit(1);}
    if((k2=(double*) malloc(3*sizeof(double)))==NULL){printf("error reserving memory\n"); exit(1);}
    if((k3=(double*) malloc(3*sizeof(double)))==NULL){printf("error reserving memory\n"); exit(1);}
    if((k4=(double*) malloc(3*sizeof(double)))==NULL){printf("error reserving memory\n"); exit(1);}
    if((sol=(double*) malloc(3*sizeof(double)))==NULL){printf("error reserving memory\n"); exit(1);} 
    //Pressure, mass, phi

    if(saveRadialProfile){
        char buffer[100];
        char filename[104]="TASPPprofile_";
        sprintf(buffer, "%g", rho_c);
        strcat(filename, buffer);
        strcat(filename, ".txt");
        file=fopen(filename, "w");
        fprintf(file, "#1-r(km)\t2-m(Ms)\t3-rho(g/cm3)\t4-phi(not rescaled)\n");
    }
    if(checkAccuracyMetric){
        temp=fopen("temp_file.txt", "w");
    }

    double dr=1.0, rho, r, mprev, eps;
    //Initial condition
    r=dr; rho=rho_c; eps=epsilonTASPP(rho_c,eos);
    sol[0]=pressureTASPP(rho_c, eos); 
    sol[1]=4*M_PI*r*r*r*rho*(1+eps);  
    sol[2]=0;

    //First iteration out of the loop to deal with the singularity
    k1[0]=0;
    k1[1]=0;
    k1[2]=0;

    k2[0]=dPdr(r+dr*0.5, sol[1], sol[0], rho*(1.0+eps))/c2;
    k2[1]=dmdr(r+dr*0.5, rho*(1.0+eps));
    k2[2]=dphidr(k2[0], sol[0], rho*(1.0+eps));

    rho=densityTASPP(sol[0]+0.5*dr*k2[0],eos); eps=epsilonTASPP(rho,eos);
    k3[0]=dPdr(r+dr*0.5, sol[1]+dr*0.5*k2[1], sol[0]+dr*0.5*k2[0], rho*(1.0+eps))/c2;
    k3[1]=dmdr(r+dr*0.5, rho*(1.0+eps));
    k3[2]=dphidr(k3[0], sol[0]+dr*0.5*k2[0], rho*(1.0+eps));

    rho=densityTASPP(sol[0]+dr*k3[0],eos); eps=epsilonTASPP(rho,eos);
    k4[0]=dPdr(r+dr, sol[1]+dr*k3[1], sol[0]+dr*k3[0], rho*(1.0+eps))/c2;
    k4[1]=dmdr(r+dr, rho*(1.0+eps));
    k4[2]=dphidr(k4[0], sol[0]+dr*k3[0], rho*(1.0+eps));

    sol[0]+=dr*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6;
    sol[1]+=dr*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6;
    sol[2]+=dr*(k1[2]+2*k2[2]+2*k3[2]+k4[2])/6;

    r+=dr;
    do{
        rho=densityTASPP(sol[0],eos); eps=epsilonTASPP(rho,eos);
        k1[0]=dPdr(r,sol[1],sol[0],rho*(1.0+eps))/c2;
        k1[1]=dmdr(r,rho*(1.0+eps));
        k1[2]=dphidr(k1[0],sol[0],rho*(1.0+eps));

        rho=densityTASPP(sol[0]+0.5*dr*k1[0],eos); eps=epsilonTASPP(rho,eos);
        k2[0]=dPdr(r+dr*0.5, sol[1]+0.5*dr*k1[1], sol[0]+0.5*dr*k1[0], rho*(1.0+eps))/c2;
        k2[1]=dmdr(r+dr*0.5, rho*(1.0+eps));
        k2[2]=dphidr(k2[0], sol[0]+0.5*dr*k1[0], rho*(1.0+eps));

        rho=densityTASPP(sol[0]+0.5*dr*k2[0],eos); eps=epsilonTASPP(rho,eos);
        k3[0]=dPdr(r+dr*0.5, sol[1]+dr*0.5*k2[1], sol[0]+dr*0.5*k2[0], rho*(1.0+eps))/c2;
        k3[1]=dmdr(r+dr*0.5, rho*(1.0+eps));
        k3[2]=dphidr(k3[0], sol[0]+dr*0.5*k2[0], rho*(1.0+eps));

        rho=densityTASPP(sol[0]+dr*k3[0],eos); eps=epsilonTASPP(rho,eos);
        k4[0]=dPdr(r+dr, sol[1]+dr*k3[1], sol[0]+dr*k3[0], rho*(1.0+eps))/c2;
        k4[1]=dmdr(r+dr, rho*(1.0+eps));
        k4[2]=dphidr(k4[0], sol[0]+dr*k3[0], rho*(1.0+eps));

        sol[0]+=dr*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6; mprev=sol[1];
        sol[1]+=dr*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6;
        sol[2]+=dr*(k1[2]+2*k2[2]+2*k3[2]+k4[2])/6;

        if(saveRadialProfile){
            fprintf(file, "%e %.15e %.15e %.15e\n",r*1e-5, sol[1]*0.001/1.98847e30, densityTASPP(sol[0],eos), sol[2]);
        }
        if(checkAccuracyMetric){
            fprintf(temp, "%.15e\t%.15e\n",densityTASPP(sol[0],eos), sol[2]);
        }
       
        if(r>2 && (sol[1]-mprev)/mprev<1e-12) break; //stopping condition by saturation of mass
        r+=dr;
        
    }while (sol[0]>5e4); //Safety stopping condition: very low pressure
    (*Mmax)=sol[1]*0.001/1.98847e30;

    if(saveRadialProfile) fclose(file);

    if(checkAccuracyMetric){ 
        fclose(temp);
        
        double R, Rsch, phiR, phiSch, renorm;
        R=r;
        Rsch=2*Ggrav*sol[1]/c2;
        phiR=sol[2];
        phiSch=0.5*log(1-Rsch/R);
        renorm=-phiR+phiSch;
   
        const unsigned MAX_LENGTH = 256;
        char buffer[MAX_LENGTH], line2[MAX_LENGTH];
        int element=0;
        char *token;
        double phi, P, quantity;
        double maxCte=0, minCte=1e30;
        FILE *fp = fopen("temp_file.txt", "r");
        while (fgets(buffer, MAX_LENGTH, fp)){
            memcpy(line2, buffer, MAX_LENGTH);
            token = strtok(line2, "\t");
            element=0;
            while( token != NULL ) {
                if(element==0) rho=strtod(token,NULL);
                if(element==1) phi=strtod(token,NULL);
                token = strtok(NULL, "\t");
                element++;
            }
            phi+=renorm;
            eps=epsilonTASPP(rho, eos);
            P=pressureTASPP(rho, eos);
            quantity=(rho*(1+eps)+P)/rho*exp(phi);
            if (quantity>maxCte)
                maxCte=quantity;
            else if (quantity<minCte)
                minCte=quantity;
        }
        fclose(fp);
        remove("temp_file.txt");
        if(saveRadialProfile){
            printf("An accurate integration keeps quantity mu_B*exp(phi) constant\n");
            printf("The maximum absolute difference in the constant quantity was %.15e\n", maxCte-minCte);
            printf("The maximum relative difference in the constant quantity was %.15e\n", (maxCte-minCte)/minCte); 
        }
        (*Maxerror)=(maxCte-minCte)/minCte;
    }

    free(k1); free(k2); free(k3); free(k4); free(sol);
    return r*1e-5;
}

double solveStarMRPL_TASPP(double *Mmax, double *LambdaNS, double rho_c, int saveRadialProfile, int checkAccuracyMetric, double *Maxerror, struct TASPP eos){
    double *k1, *k2, *k3, *k4, *sol; FILE *file; FILE *temp;
    if((k1=(double*) malloc(5*sizeof(double)))==NULL){printf("error reserving memory\n"); exit(1);}
    if((k2=(double*) malloc(5*sizeof(double)))==NULL){printf("error reserving memory\n"); exit(1);}
    if((k3=(double*) malloc(5*sizeof(double)))==NULL){printf("error reserving memory\n"); exit(1);}
    if((k4=(double*) malloc(5*sizeof(double)))==NULL){printf("error reserving memory\n"); exit(1);}
    if((sol=(double*) malloc(5*sizeof(double)))==NULL){printf("error reserving memory\n"); exit(1);} 
    //Pressure, mass, phi, H, beta

    if(saveRadialProfile){
        char buffer[100];
        char filename[104]="TASPPprofile_";
        sprintf(buffer, "%g", rho_c);
        strcat(filename, buffer);
        strcat(filename, ".txt");
        file=fopen(filename, "w");
        fprintf(file, "#1-r(km)\t2-m(Ms)\t3-rho(g/cm3)\t4-phi(not rescaled)\t5-Lambda\n");
    }
    if(checkAccuracyMetric){
        temp=fopen("temp_file.txt", "w");
    }

    double R, M, H, beta, y, C, C5, lovek2, Lambda;
    double dr=1.0, rho, r, mprev, eps;
    //Initial condition
    r=dr; rho=rho_c; eps=epsilonTASPP(rho_c,eos);
    sol[0]=pressureTASPP(rho_c, eos); 
    sol[1]=4*M_PI*r*r*r*rho*(1+eps);  
    sol[2]=0;
    sol[3]=r*r; 
    sol[4]=2*r;

    //First iteration out of the loop to deal with the singularity
    k1[0]=0;
    k1[1]=0;
    k1[2]=0;
    k1[3]=dHdr(sol[4]); 
    k1[4]=dbetadr(r,sol[1],sol[0],sol[3],sol[4],rho*(1.0+eps),ddensdPTASPP(sol[0],eos));
    

    k2[0]=dPdr(r+dr*0.5, sol[1], sol[0], rho*(1.0+eps))/c2;
    k2[1]=dmdr(r+dr*0.5, rho*(1.0+eps));
    k2[2]=dphidr(k2[0], sol[0], rho*(1.0+eps));
    k2[3]=dHdr(sol[4]);
    k2[4]=dbetadr(r+dr*0.5, sol[1],sol[0],sol[3],sol[4],rho*(1+eps),ddensdPTASPP(sol[0],eos));

    rho=densityTASPP(sol[0]+0.5*dr*k2[0],eos); eps=epsilonTASPP(rho,eos);
    k3[0]=dPdr(r+dr*0.5, sol[1]+dr*0.5*k2[1], sol[0]+dr*0.5*k2[0], rho*(1.0+eps))/c2;
    k3[1]=dmdr(r+dr*0.5, rho*(1.0+eps));
    k3[2]=dphidr(k3[0], sol[0]+dr*0.5*k2[0], rho*(1.0+eps));
    k3[3]=dHdr(sol[4]+dr*0.5*k2[4]);
    k3[4]=dbetadr(r+dr*0.5, sol[1]+dr*0.5*k2[1],sol[0]+dr*0.5*k2[0],sol[3]+dr*0.5*k2[3],sol[4]+dr*0.5*k2[4],rho*(1+eps),ddensdPTASPP(sol[0]+dr*0.5*k2[0], eos));

    rho=densityTASPP(sol[0]+dr*k3[0],eos); eps=epsilonTASPP(rho,eos);
    k4[0]=dPdr(r+dr, sol[1]+dr*k3[1], sol[0]+dr*k3[0], rho*(1.0+eps))/c2;
    k4[1]=dmdr(r+dr, rho*(1.0+eps));
    k4[2]=dphidr(k4[0], sol[0]+dr*k3[0], rho*(1.0+eps));
    k4[3]=dHdr(sol[4]+dr*k3[4]);
    k4[4]=dbetadr(r+dr, sol[1]+dr*k3[1], sol[0]+dr*k3[0],sol[3]+dr*k3[3], sol[4]+dr*k3[4], rho*(1+eps),ddensdPTASPP(sol[0]+dr*k3[0],eos));

    sol[0]+=dr*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6;
    sol[1]+=dr*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6;
    sol[2]+=dr*(k1[2]+2*k2[2]+2*k3[2]+k4[2])/6;
    sol[3]+=dr*(k1[3]+2*k2[3]+2*k3[3]+k4[3])/6;
    sol[4]+=dr*(k1[4]+2*k2[4]+2*k3[4]+k4[4])/6;

    r+=dr;
    do{
        rho=densityTASPP(sol[0],eos); eps=epsilonTASPP(rho,eos);
        k1[0]=dPdr(r,sol[1],sol[0],rho*(1.0+eps))/c2;
        k1[1]=dmdr(r,rho*(1.0+eps));
        k1[2]=dphidr(k1[0],sol[0],rho*(1.0+eps));
        k1[3]=dHdr(sol[4]); 
        k1[4]=dbetadr(r,sol[1],sol[0],sol[3],sol[4],rho*(1.0+eps),ddensdPTASPP(sol[0],eos));

        rho=densityTASPP(sol[0]+0.5*dr*k1[0],eos); eps=epsilonTASPP(rho,eos);
        k2[0]=dPdr(r+dr*0.5, sol[1]+0.5*dr*k1[1], sol[0]+0.5*dr*k1[0], rho*(1.0+eps))/c2;
        k2[1]=dmdr(r+dr*0.5, rho*(1.0+eps));
        k2[2]=dphidr(k2[0], sol[0]+0.5*dr*k1[0], rho*(1.0+eps));
        k2[3]=dHdr(sol[4]);
        k2[4]=dbetadr(r+dr*0.5, sol[1],sol[0],sol[3],sol[4],rho*(1+eps),ddensdPTASPP(sol[0],eos));

        rho=densityTASPP(sol[0]+0.5*dr*k2[0],eos); eps=epsilonTASPP(rho,eos);
        k3[0]=dPdr(r+dr*0.5, sol[1]+dr*0.5*k2[1], sol[0]+dr*0.5*k2[0], rho*(1.0+eps))/c2;
        k3[1]=dmdr(r+dr*0.5, rho*(1.0+eps));
        k3[2]=dphidr(k3[0], sol[0]+dr*0.5*k2[0], rho*(1.0+eps));
        k3[3]=dHdr(sol[4]+dr*0.5*k2[4]);
        k3[4]=dbetadr(r+dr*0.5, sol[1]+dr*0.5*k2[1],sol[0]+dr*0.5*k2[0],sol[3]+dr*0.5*k2[3],sol[4]+dr*0.5*k2[4],rho*(1+eps),ddensdPTASPP(sol[0]+dr*0.5*k2[0], eos));

        rho=densityTASPP(sol[0]+dr*k3[0],eos); eps=epsilonTASPP(rho,eos);
        k4[0]=dPdr(r+dr, sol[1]+dr*k3[1], sol[0]+dr*k3[0], rho*(1.0+eps))/c2;
        k4[1]=dmdr(r+dr, rho*(1.0+eps));
        k4[2]=dphidr(k4[0], sol[0]+dr*k3[0], rho*(1.0+eps));
        k4[3]=dHdr(sol[4]+dr*k3[4]);
        k4[4]=dbetadr(r+dr, sol[1]+dr*k3[1], sol[0]+dr*k3[0],sol[3]+dr*k3[3], sol[4]+dr*k3[4], rho*(1+eps),ddensdPTASPP(sol[0]+dr*k3[0],eos));

        sol[0]+=dr*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6; mprev=sol[1];
        sol[1]+=dr*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6;
        sol[2]+=dr*(k1[2]+2*k2[2]+2*k3[2]+k4[2])/6;
        sol[3]+=dr*(k1[3]+2*k2[3]+2*k3[3]+k4[3])/6;
        sol[4]+=dr*(k1[4]+2*k2[4]+2*k3[4]+k4[4])/6;

        if(saveRadialProfile){     
            R=r; M=sol[1]; H=sol[3]; beta=sol[4];  rho=densityTASPP(sol[0],eos); eps=epsilonTASPP(rho,eos);
            y=R*beta/H-4*M_PI*R*R*R*rho*(1+eps)/M; //added correction for nonzero density at surface
            C=M/R*Ggrav/c2; C5=C*C*C*C*C;
            lovek2=8.0/5*C5*(1-2*C)*(1-2*C)*(2+2*C*(y-1)-y)/((2*C*(6-3*y+3*C*(5*y-8)))+4*C*C*C*(13-11*y+C*(3*y-2)+2*C*C*(1+y))+3*(1-2*C)*(1-2*C)*(2-y+2*C*(y-1))*log(1-2*C));   
            Lambda=2.0/3*lovek2/C5;
            fprintf(file, "%e %.15e %.15e %.15e %.15e\n",r*1e-5, sol[1]*0.001/1.98847e30, densityTASPP(sol[0],eos), sol[2], Lambda);
        }
        if(checkAccuracyMetric){
            fprintf(temp, "%.15e\t%.15e\n",densityTASPP(sol[0],eos), sol[2]);
        }
       
        if(r>2 && (sol[1]-mprev)/mprev<1e-12) break; //stopping condition by saturation of mass
        r+=dr;
        
    }while (sol[0]>5e4); //Safety stopping condition: very low pressure
    (*Mmax)=sol[1]*0.001/1.98847e30;
    R=r; M=sol[1]; H=sol[3]; beta=sol[4];  rho=densityTASPP(sol[0],eos); eps=epsilonTASPP(rho,eos);
    y=R*beta/H-4*M_PI*R*R*R*rho*(1+eps)/M; //added correction for nonzero density at surface
    C=M/R*Ggrav/c2; C5=C*C*C*C*C;
    lovek2=8.0/5*C5*(1-2*C)*(1-2*C)*(2+2*C*(y-1)-y)/((2*C*(6-3*y+3*C*(5*y-8)))+4*C*C*C*(13-11*y+C*(3*y-2)+2*C*C*(1+y))+3*(1-2*C)*(1-2*C)*(2-y+2*C*(y-1))*log(1-2*C));   
    Lambda=2.0/3*lovek2/C5;
    (*LambdaNS)=Lambda;

    if(saveRadialProfile) fclose(file);

    if(checkAccuracyMetric){ 
        fclose(temp);
        
        double R, Rsch, phiR, phiSch, renorm;
        R=r;
        Rsch=2*Ggrav*sol[1]/c2;
        phiR=sol[2];
        phiSch=0.5*log(1-Rsch/R);
        renorm=-phiR+phiSch;
   
        const unsigned MAX_LENGTH = 256;
        char buffer[MAX_LENGTH], line2[MAX_LENGTH];
        int element=0;
        char *token;
        double phi, P, quantity;
        double maxCte=0, minCte=1e30;
        FILE *fp = fopen("temp_file.txt", "r");
        while (fgets(buffer, MAX_LENGTH, fp)){
            memcpy(line2, buffer, MAX_LENGTH);
            token = strtok(line2, "\t");
            element=0;
            while( token != NULL ) {
                if(element==0) rho=strtod(token,NULL);
                if(element==1) phi=strtod(token,NULL);
                token = strtok(NULL, "\t");
                element++;
            }
            phi+=renorm;
            eps=epsilonTASPP(rho, eos);
            P=pressureTASPP(rho, eos);
            quantity=(rho*(1+eps)+P)/rho*exp(phi);
            if (quantity>maxCte)
                maxCte=quantity;
            else if (quantity<minCte)
                minCte=quantity;
        }
        fclose(fp);
        remove("temp_file.txt");
        if(saveRadialProfile){
            printf("An accurate integration keeps quantity mu_B*exp(phi) constant\n");
            printf("The maximum absolute difference in the constant quantity was %.15e\n", maxCte-minCte);
            printf("The maximum relative difference in the constant quantity was %.15e\n", (maxCte-minCte)/minCte); 
        }
        (*Maxerror)=(maxCte-minCte)/minCte;
    }

    free(k1); free(k2); free(k3); free(k4); free(sol);
    return r*1e-5;
}

double solveStarMRP_TAB(double *Mmax, double rho_c, int saveRadialProfile, int checkAccuracyMetric, double *Maxerror, struct Tab eos){
    double *k1, *k2, *k3, *k4, *sol; FILE *file; FILE *temp;
    if((k1=(double*) malloc(3*sizeof(double)))==NULL){printf("error reserving memory\n"); exit(1);}
    if((k2=(double*) malloc(3*sizeof(double)))==NULL){printf("error reserving memory\n"); exit(1);}
    if((k3=(double*) malloc(3*sizeof(double)))==NULL){printf("error reserving memory\n"); exit(1);}
    if((k4=(double*) malloc(3*sizeof(double)))==NULL){printf("error reserving memory\n"); exit(1);}
    if((sol=(double*) malloc(3*sizeof(double)))==NULL){printf("error reserving memory\n"); exit(1);} 
    //Pressure, mass, phi

    if(saveRadialProfile){
        char buffer[100];
        char filename[104]="TABprofile_";
        sprintf(buffer, "%g", rho_c);
        strcat(filename, buffer);
        strcat(filename, ".txt");
        file=fopen(filename, "w");
        fprintf(file, "#1-r(km)\t2-m(Ms)\t3-rho(g/cm3)\t4-phi(not rescaled)\n");
    }
    if(checkAccuracyMetric){
        temp=fopen("temp_file.txt", "w");
    }

    double dr=1.0, rho, r, mprev, eps;
    //Initial condition
    r=dr; rho=rho_c; eps=epsilonTab(rho_c,eos);
    sol[0]=pressureTab(rho_c, eos); 
    sol[1]=4*M_PI*r*r*r*rho*(1+eps);  
    sol[2]=0;

    //First iteration out of the loop to deal with the singularity
    k1[0]=0;
    k1[1]=0;
    k1[2]=0;

    k2[0]=dPdr(r+dr*0.5, sol[1], sol[0], rho*(1.0+eps))/c2;
    k2[1]=dmdr(r+dr*0.5, rho*(1.0+eps));
    k2[2]=dphidr(k2[0], sol[0], rho*(1.0+eps));

    rho=densityTab(sol[0]+0.5*dr*k2[0],eos); eps=epsilonTab(rho,eos);
    k3[0]=dPdr(r+dr*0.5, sol[1]+dr*0.5*k2[1], sol[0]+dr*0.5*k2[0], rho*(1.0+eps))/c2;
    k3[1]=dmdr(r+dr*0.5, rho*(1.0+eps));
    k3[2]=dphidr(k3[0], sol[0]+dr*0.5*k2[0], rho*(1.0+eps));

    rho=densityTab(sol[0]+dr*k3[0],eos); eps=epsilonTab(rho,eos);
    k4[0]=dPdr(r+dr, sol[1]+dr*k3[1], sol[0]+dr*k3[0], rho*(1.0+eps))/c2;
    k4[1]=dmdr(r+dr, rho*(1.0+eps));
    k4[2]=dphidr(k4[0], sol[0]+dr*k3[0], rho*(1.0+eps));

    sol[0]+=dr*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6;
    sol[1]+=dr*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6;
    sol[2]+=dr*(k1[2]+2*k2[2]+2*k3[2]+k4[2])/6;

    r+=dr;
    do{
        rho=densityTab(sol[0],eos); eps=epsilonTab(rho,eos);
        k1[0]=dPdr(r,sol[1],sol[0],rho*(1.0+eps))/c2;
        k1[1]=dmdr(r,rho*(1.0+eps));
        k1[2]=dphidr(k1[0],sol[0],rho*(1.0+eps));

        rho=densityTab(sol[0]+0.5*dr*k1[0],eos); eps=epsilonTab(rho,eos);
        k2[0]=dPdr(r+dr*0.5, sol[1]+0.5*dr*k1[1], sol[0]+0.5*dr*k1[0], rho*(1.0+eps))/c2;
        k2[1]=dmdr(r+dr*0.5, rho*(1.0+eps));
        k2[2]=dphidr(k2[0], sol[0]+0.5*dr*k1[0], rho*(1.0+eps));

        rho=densityTab(sol[0]+0.5*dr*k2[0],eos); eps=epsilonTab(rho,eos);
        k3[0]=dPdr(r+dr*0.5, sol[1]+dr*0.5*k2[1], sol[0]+dr*0.5*k2[0], rho*(1.0+eps))/c2;
        k3[1]=dmdr(r+dr*0.5, rho*(1.0+eps));
        k3[2]=dphidr(k3[0], sol[0]+dr*0.5*k2[0], rho*(1.0+eps));

        rho=densityTab(sol[0]+dr*k3[0],eos); eps=epsilonTab(rho,eos);
        k4[0]=dPdr(r+dr, sol[1]+dr*k3[1], sol[0]+dr*k3[0], rho*(1.0+eps))/c2;
        k4[1]=dmdr(r+dr, rho*(1.0+eps));
        k4[2]=dphidr(k4[0], sol[0]+dr*k3[0], rho*(1.0+eps));

        sol[0]+=dr*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6; mprev=sol[1];
        sol[1]+=dr*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6;
        sol[2]+=dr*(k1[2]+2*k2[2]+2*k3[2]+k4[2])/6;

        if(saveRadialProfile){
            fprintf(file, "%e %.15e %.15e %.15e\n",r*1e-5, sol[1]*0.001/1.98847e30, densityTab(sol[0],eos), sol[2]);
        }
        if(checkAccuracyMetric){
            fprintf(temp, "%.15e\t%.15e\n",densityTab(sol[0],eos), sol[2]);
        }
       
        if(r>2 && (sol[1]-mprev)/mprev<1e-12) break; //stopping condition by saturation of mass
        r+=dr;
        
    }while (sol[0]>5e4); //Safety stopping condition: very low pressure
    (*Mmax)=sol[1]*0.001/1.98847e30;

    if(saveRadialProfile) fclose(file);

    if(checkAccuracyMetric){ 
        fclose(temp);
        
        double R, Rsch, phiR, phiSch, renorm;
        R=r;
        Rsch=2*Ggrav*sol[1]/c2;
        phiR=sol[2];
        phiSch=0.5*log(1-Rsch/R);
        renorm=-phiR+phiSch;
   
        const unsigned MAX_LENGTH = 256;
        char buffer[MAX_LENGTH], line2[MAX_LENGTH];
        int element=0;
        char *token;
        double phi, P, quantity;
        double maxCte=0, minCte=1e30;
        FILE *fp = fopen("temp_file.txt", "r");
        while (fgets(buffer, MAX_LENGTH, fp)){
            memcpy(line2, buffer, MAX_LENGTH);
            token = strtok(line2, "\t");
            element=0;
            while( token != NULL ) {
                if(element==0) rho=strtod(token,NULL);
                if(element==1) phi=strtod(token,NULL);
                token = strtok(NULL, "\t");
                element++;
            }
            phi+=renorm;
            eps=epsilonTab(rho, eos);
            P=pressureTab(rho, eos);
            quantity=(rho*(1+eps)+P)/rho*exp(phi);
            if (quantity>maxCte)
                maxCte=quantity;
            else if (quantity<minCte)
                minCte=quantity;
        }
        fclose(fp);
        remove("temp_file.txt");
        if(saveRadialProfile){
            printf("An accurate integration keeps quantity mu_B*exp(phi) constant\n");
            printf("The maximum absolute difference in the constant quantity was %.15e\n", maxCte-minCte);
            printf("The maximum relative difference in the constant quantity was %.15e\n", (maxCte-minCte)/minCte); 
        }
        (*Maxerror)=(maxCte-minCte)/minCte;
    }

    free(k1); free(k2); free(k3); free(k4); free(sol);
    return r*1e-5;
}

double solveStarMRPL_TAB(double *Mmax, double *LambdaNS, double rho_c, int saveRadialProfile, int checkAccuracyMetric, double *Maxerror, struct Tab eos){
    double *k1, *k2, *k3, *k4, *sol; FILE *file; FILE *temp;
    if((k1=(double*) malloc(5*sizeof(double)))==NULL){printf("error reserving memory\n"); exit(1);}
    if((k2=(double*) malloc(5*sizeof(double)))==NULL){printf("error reserving memory\n"); exit(1);}
    if((k3=(double*) malloc(5*sizeof(double)))==NULL){printf("error reserving memory\n"); exit(1);}
    if((k4=(double*) malloc(5*sizeof(double)))==NULL){printf("error reserving memory\n"); exit(1);}
    if((sol=(double*) malloc(5*sizeof(double)))==NULL){printf("error reserving memory\n"); exit(1);} 
    //Pressure, mass, phi, H, beta

    if(saveRadialProfile){
        char buffer[100];
        char filename[104]="TABprofile_";
        sprintf(buffer, "%g", rho_c);
        strcat(filename, buffer);
        strcat(filename, ".txt");
        file=fopen(filename, "w");
        fprintf(file, "#1-r(km)\t2-m(Ms)\t3-rho(g/cm3)\t4-phi(not rescaled)\t5-Lambda\n");
    }
    if(checkAccuracyMetric){
        temp=fopen("temp_file.txt", "w");
    }

    double R, M, H, beta, y, C, C5, lovek2, Lambda;
    double dr=1.0, rho, r, mprev, eps;
    //Initial condition
    r=dr; rho=rho_c; eps=epsilonTab(rho_c,eos);
    sol[0]=pressureTab(rho_c, eos); 
    sol[1]=4*M_PI*r*r*r*rho*(1+eps);  
    sol[2]=0;
    sol[3]=r*r; 
    sol[4]=2*r;

    //First iteration out of the loop to deal with the singularity
    k1[0]=0;
    k1[1]=0;
    k1[2]=0;
    k1[3]=dHdr(sol[4]); 
    k1[4]=dbetadr(r,sol[1],sol[0],sol[3],sol[4],rho*(1.0+eps),ddensdPTab(sol[0],eos));
    

    k2[0]=dPdr(r+dr*0.5, sol[1], sol[0], rho*(1.0+eps))/c2;
    k2[1]=dmdr(r+dr*0.5, rho*(1.0+eps));
    k2[2]=dphidr(k2[0], sol[0], rho*(1.0+eps));
    k2[3]=dHdr(sol[4]);
    k2[4]=dbetadr(r+dr*0.5, sol[1],sol[0],sol[3],sol[4],rho*(1+eps),ddensdPTab(sol[0],eos));

    rho=densityTab(sol[0]+0.5*dr*k2[0],eos); eps=epsilonTab(rho,eos);
    k3[0]=dPdr(r+dr*0.5, sol[1]+dr*0.5*k2[1], sol[0]+dr*0.5*k2[0], rho*(1.0+eps))/c2;
    k3[1]=dmdr(r+dr*0.5, rho*(1.0+eps));
    k3[2]=dphidr(k3[0], sol[0]+dr*0.5*k2[0], rho*(1.0+eps));
    k3[3]=dHdr(sol[4]+dr*0.5*k2[4]);
    k3[4]=dbetadr(r+dr*0.5, sol[1]+dr*0.5*k2[1],sol[0]+dr*0.5*k2[0],sol[3]+dr*0.5*k2[3],sol[4]+dr*0.5*k2[4],rho*(1+eps),ddensdPTab(sol[0]+dr*0.5*k2[0], eos));

    rho=densityTab(sol[0]+dr*k3[0],eos); eps=epsilonTab(rho,eos);
    k4[0]=dPdr(r+dr, sol[1]+dr*k3[1], sol[0]+dr*k3[0], rho*(1.0+eps))/c2;
    k4[1]=dmdr(r+dr, rho*(1.0+eps));
    k4[2]=dphidr(k4[0], sol[0]+dr*k3[0], rho*(1.0+eps));
    k4[3]=dHdr(sol[4]+dr*k3[4]);
    k4[4]=dbetadr(r+dr, sol[1]+dr*k3[1], sol[0]+dr*k3[0],sol[3]+dr*k3[3], sol[4]+dr*k3[4], rho*(1+eps),ddensdPTab(sol[0]+dr*k3[0],eos));

    sol[0]+=dr*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6;
    sol[1]+=dr*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6;
    sol[2]+=dr*(k1[2]+2*k2[2]+2*k3[2]+k4[2])/6;
    sol[3]+=dr*(k1[3]+2*k2[3]+2*k3[3]+k4[3])/6;
    sol[4]+=dr*(k1[4]+2*k2[4]+2*k3[4]+k4[4])/6;

    r+=dr;
    do{
        rho=densityTab(sol[0],eos); eps=epsilonTab(rho,eos);
        k1[0]=dPdr(r,sol[1],sol[0],rho*(1.0+eps))/c2;
        k1[1]=dmdr(r,rho*(1.0+eps));
        k1[2]=dphidr(k1[0],sol[0],rho*(1.0+eps));
        k1[3]=dHdr(sol[4]); 
        k1[4]=dbetadr(r,sol[1],sol[0],sol[3],sol[4],rho*(1.0+eps),ddensdPTab(sol[0],eos));

        rho=densityTab(sol[0]+0.5*dr*k1[0],eos); eps=epsilonTab(rho,eos);
        k2[0]=dPdr(r+dr*0.5, sol[1]+0.5*dr*k1[1], sol[0]+0.5*dr*k1[0], rho*(1.0+eps))/c2;
        k2[1]=dmdr(r+dr*0.5, rho*(1.0+eps));
        k2[2]=dphidr(k2[0], sol[0]+0.5*dr*k1[0], rho*(1.0+eps));
        k2[3]=dHdr(sol[4]);
        k2[4]=dbetadr(r+dr*0.5, sol[1],sol[0],sol[3],sol[4],rho*(1+eps),ddensdPTab(sol[0],eos));

        rho=densityTab(sol[0]+0.5*dr*k2[0],eos); eps=epsilonTab(rho,eos);
        k3[0]=dPdr(r+dr*0.5, sol[1]+dr*0.5*k2[1], sol[0]+dr*0.5*k2[0], rho*(1.0+eps))/c2;
        k3[1]=dmdr(r+dr*0.5, rho*(1.0+eps));
        k3[2]=dphidr(k3[0], sol[0]+dr*0.5*k2[0], rho*(1.0+eps));
        k3[3]=dHdr(sol[4]+dr*0.5*k2[4]);
        k3[4]=dbetadr(r+dr*0.5, sol[1]+dr*0.5*k2[1],sol[0]+dr*0.5*k2[0],sol[3]+dr*0.5*k2[3],sol[4]+dr*0.5*k2[4],rho*(1+eps),ddensdPTab(sol[0]+dr*0.5*k2[0], eos));

        rho=densityTab(sol[0]+dr*k3[0],eos); eps=epsilonTab(rho,eos);
        k4[0]=dPdr(r+dr, sol[1]+dr*k3[1], sol[0]+dr*k3[0], rho*(1.0+eps))/c2;
        k4[1]=dmdr(r+dr, rho*(1.0+eps));
        k4[2]=dphidr(k4[0], sol[0]+dr*k3[0], rho*(1.0+eps));
        k4[3]=dHdr(sol[4]+dr*k3[4]);
        k4[4]=dbetadr(r+dr, sol[1]+dr*k3[1], sol[0]+dr*k3[0],sol[3]+dr*k3[3], sol[4]+dr*k3[4], rho*(1+eps),ddensdPTab(sol[0]+dr*k3[0],eos));

        sol[0]+=dr*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6; mprev=sol[1];
        sol[1]+=dr*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6;
        sol[2]+=dr*(k1[2]+2*k2[2]+2*k3[2]+k4[2])/6;
        sol[3]+=dr*(k1[3]+2*k2[3]+2*k3[3]+k4[3])/6;
        sol[4]+=dr*(k1[4]+2*k2[4]+2*k3[4]+k4[4])/6;

        if(saveRadialProfile){     
            R=r; M=sol[1]; H=sol[3]; beta=sol[4];  rho=densityTab(sol[0],eos); eps=epsilonTab(rho,eos);
            y=R*beta/H-4*M_PI*R*R*R*rho*(1+eps)/M; //added correction for nonzero density at surface
            C=M/R*Ggrav/c2; C5=C*C*C*C*C;
            lovek2=8.0/5*C5*(1-2*C)*(1-2*C)*(2+2*C*(y-1)-y)/((2*C*(6-3*y+3*C*(5*y-8)))+4*C*C*C*(13-11*y+C*(3*y-2)+2*C*C*(1+y))+3*(1-2*C)*(1-2*C)*(2-y+2*C*(y-1))*log(1-2*C));   
            Lambda=2.0/3*lovek2/C5;
            fprintf(file, "%e %.15e %.15e %.15e %.15e\n",r*1e-5, sol[1]*0.001/1.98847e30, densityTab(sol[0],eos), sol[2], Lambda);
        }
        if(checkAccuracyMetric){
            fprintf(temp, "%.15e\t%.15e\n",densityTab(sol[0],eos), sol[2]);
        }
       
        if(r>2 && (sol[1]-mprev)/mprev<1e-12) break; //stopping condition by saturation of mass
        r+=dr;
        
    }while (sol[0]>5e4); //Safety stopping condition: very low pressure
    (*Mmax)=sol[1]*0.001/1.98847e30;
    R=r; M=sol[1]; H=sol[3]; beta=sol[4];  rho=densityTab(sol[0],eos); eps=epsilonTab(rho,eos);
    y=R*beta/H-4*M_PI*R*R*R*rho*(1+eps)/M; //added correction for nonzero density at surface
    C=M/R*Ggrav/c2; C5=C*C*C*C*C;
    lovek2=8.0/5*C5*(1-2*C)*(1-2*C)*(2+2*C*(y-1)-y)/((2*C*(6-3*y+3*C*(5*y-8)))+4*C*C*C*(13-11*y+C*(3*y-2)+2*C*C*(1+y))+3*(1-2*C)*(1-2*C)*(2-y+2*C*(y-1))*log(1-2*C));   
    Lambda=2.0/3*lovek2/C5;
    (*LambdaNS)=Lambda;

    if(saveRadialProfile) fclose(file);

    if(checkAccuracyMetric){ 
        fclose(temp);
        
        double R, Rsch, phiR, phiSch, renorm;
        R=r;
        Rsch=2*Ggrav*sol[1]/c2;
        phiR=sol[2];
        phiSch=0.5*log(1-Rsch/R);
        renorm=-phiR+phiSch;
   
        const unsigned MAX_LENGTH = 256;
        char buffer[MAX_LENGTH], line2[MAX_LENGTH];
        int element=0;
        char *token;
        double phi, P, quantity;
        double maxCte=0, minCte=1e30;
        FILE *fp = fopen("temp_file.txt", "r");
        while (fgets(buffer, MAX_LENGTH, fp)){
            memcpy(line2, buffer, MAX_LENGTH);
            token = strtok(line2, "\t");
            element=0;
            while( token != NULL ) {
                if(element==0) rho=strtod(token,NULL);
                if(element==1) phi=strtod(token,NULL);
                token = strtok(NULL, "\t");
                element++;
            }
            phi+=renorm;
            eps=epsilonTab(rho, eos);
            P=pressureTab(rho, eos);
            quantity=(rho*(1+eps)+P)/rho*exp(phi);
            if (quantity>maxCte)
                maxCte=quantity;
            else if (quantity<minCte)
                minCte=quantity;
        }
        fclose(fp);
        remove("temp_file.txt");
        if(saveRadialProfile){
            printf("An accurate integration keeps quantity mu_B*exp(phi) constant\n");
            printf("The maximum absolute difference in the constant quantity was %.15e\n", maxCte-minCte);
            printf("The maximum relative difference in the constant quantity was %.15e\n", (maxCte-minCte)/minCte); 
        }
        (*Maxerror)=(maxCte-minCte)/minCte;
    }

    free(k1); free(k2); free(k3); free(k4); free(sol);
    return r*1e-5;
}

void starSequence_PP(double rho_ini, double rho_end, double rho_step, int checkAccuracy, struct PP eos){
    FILE *f;
    f=fopen("PPstarSeq.txt", "w");
    if(checkAccuracy)
        fprintf(f, "#1-rho_c(g/cm3)\t2-M(Ms)\t3-R(km)\t4-Lambda()\t5-MaxDev\n");
    else
        fprintf(f, "#1-rho_c(g/cm3)\t2-M(Ms)\t3-R(km)\t4-Lambda()\n");
    double rho_c, Lambda, M, R, error;
    int s=SIGN(rho_end-rho_ini);
    for(rho_c=rho_ini; s*rho_c<=s*rho_end; rho_c+=rho_step){ //supports rho_end>rho_ini and viceversa
        R=solveStarMRPL_PP(&M, &Lambda, rho_c, 0, checkAccuracy, &error, eos);
        if(checkAccuracy)
            fprintf(f, "%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n", rho_c, M, R, Lambda, error);
        else
            fprintf(f, "%.15e\t%.15e\t%.15e\t%.15e\n", rho_c, M, R, Lambda);
            
    }
    fclose(f);
}
void starSequence_TASPP(double rho_ini, double rho_end, double rho_step, int checkAccuracy, struct TASPP eos){
    FILE *f;
    f=fopen("TASPPstarSeq.txt", "w");
    if(checkAccuracy)
        fprintf(f, "#1-rho_c(g/cm3)\t2-M(Ms)\t3-R(km)\t4-Lambda()\t5-MaxDev\n");
    else
        fprintf(f, "#1-rho_c(g/cm3)\t2-M(Ms)\t3-R(km)\t4-Lambda()\n");
    double rho_c, Lambda, M, R, error;
    int s=SIGN(rho_end-rho_ini);
    for(rho_c=rho_ini; s*rho_c<=s*rho_end; rho_c+=rho_step){ //supports rho_end>rho_ini and viceversa
        R=solveStarMRPL_TASPP(&M, &Lambda, rho_c, 0, checkAccuracy, &error, eos);
        if(checkAccuracy)
            fprintf(f, "%.15e\t%.15e\t%.15e\t%.15e\n", rho_c, M, R, Lambda);
        else
            fprintf(f, "%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n", rho_c, M, R, Lambda, error);
            
    }
    fclose(f);
}
void starSequence_TAB(double rho_ini, double rho_end, double rho_step, int checkAccuracy, struct Tab eos){
    FILE *f;
    f=fopen("TABstarSeq.txt", "w");
    if(checkAccuracy)
        fprintf(f, "#1-rho_c(g/cm3)\t2-M(Ms)\t3-R(km)\t4-Lambda()\t5-MaxDev\n");
    else
        fprintf(f, "#1-rho_c(g/cm3)\t2-M(Ms)\t3-R(km)\t4-Lambda()\n");
    double rho_c, Lambda, M, R, error;
    int s=SIGN(rho_end-rho_ini);
    for(rho_c=rho_ini; s*rho_c<=s*rho_end; rho_c+=rho_step){ //supports rho_end>rho_ini and viceversa
        R=solveStarMRPL_TAB(&M, &Lambda, rho_c, 0, checkAccuracy, &error, eos);
        if(checkAccuracy)
            fprintf(f, "%.15e\t%.15e\t%.15e\t%.15e\n", rho_c, M, R, Lambda);
        else
            fprintf(f, "%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n", rho_c, M, R, Lambda, error);
            
    }
    fclose(f);
}