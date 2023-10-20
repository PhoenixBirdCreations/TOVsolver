#include "read_par_tov.h"

void read_par_tov(char *parfile, int *saveRadialProfile, int *checkAccuracy, int *sequence, double **seqdetails){
    FILE *fp = fopen(parfile, "r");
    const unsigned MAX_LENGTH = 256;
    char buffer[MAX_LENGTH], line2[MAX_LENGTH];

    fgets(buffer, MAX_LENGTH, fp); //skip header
    
    //star or sequence
    fgets(buffer, MAX_LENGTH, fp);
    memcpy(line2, buffer, MAX_LENGTH);
    if (strstr(line2, "STAR") != NULL) {
        (*sequence)=0;
        (*seqdetails)=(double*)malloc(1*sizeof(double));
    }
    else{
        (*sequence)=1;
        (*seqdetails)=(double*)malloc(3*sizeof(double));
    } 

    //accuracy
    fgets(buffer, MAX_LENGTH, fp); //skip line
    fgets(buffer, MAX_LENGTH, fp); //skip header
    fgets(buffer, MAX_LENGTH, fp);
    memcpy(line2, buffer, MAX_LENGTH);
    (*checkAccuracy)=atoi(line2);

    //radial profile
    fgets(buffer, MAX_LENGTH, fp); //skip line
    fgets(buffer, MAX_LENGTH, fp); //skip header
    fgets(buffer, MAX_LENGTH, fp);
    memcpy(line2, buffer, MAX_LENGTH);
    (*saveRadialProfile)=atoi(line2);

    //central density
    fgets(buffer, MAX_LENGTH, fp); //skip line
    fgets(buffer, MAX_LENGTH, fp); //skip header
    fgets(buffer, MAX_LENGTH, fp);
    memcpy(line2, buffer, MAX_LENGTH);
    (*seqdetails)[0]=strtod(line2,NULL);

    //sequence details
    if (*sequence){
        fgets(buffer, MAX_LENGTH, fp); //skip line
        fgets(buffer, MAX_LENGTH, fp); //skip header
        fgets(buffer, MAX_LENGTH, fp);
        memcpy(line2, buffer, MAX_LENGTH);
        (*seqdetails)[0]=strtod(line2,NULL);
        fgets(buffer, MAX_LENGTH, fp);
        memcpy(line2, buffer, MAX_LENGTH);
        (*seqdetails)[1]=strtod(line2,NULL);
        fgets(buffer, MAX_LENGTH, fp);
        memcpy(line2, buffer, MAX_LENGTH);
        (*seqdetails)[2]=strtod(line2,NULL);
        printf("Calculating a sequence of %d stars\n", (int)(fabs(((*seqdetails)[0]-(*seqdetails)[1])/(*seqdetails)[2])));
    }
}


