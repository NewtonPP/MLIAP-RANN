#define MLIAP_RANN
#ifdef MLIAP_RANN

#include "mliap_descriptor_rann.h"
#include "mliap_descriptor.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "memory.h"
#include "mliap_data.h"
#include "neigh_list.h"
#include "pair_mliap.h"
#include "tokenizer.h"
#include "rann_fingerprint.h"
#include "pair_rann.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

static constexpr int MAXLINE = 1024;

MLIAPDescriptorRANN::MLIAPDescriptorRANN(LAMMPS *_lmp, char * paramfilename) : Pointers(_lmp),  MLIAPDescriptor(_lmp){
    read_paramfile(paramfilename);
}

double compute_cutoff (double rij, double rc, double dr){
int out; 
if (rij < (rc-dr)){
    out = 1;
}
else if (rij > rc ){
    out = 0;
}
else {
    int numerator = rc - rij;
    int denominator = dr;
    out = 1 -  numerator/denominator;
    out = pow(out,4);
    out = 1 - out;
    out = pow (out, 2);
}
return out;
}


void MLIAPDescriptorRANN::compute_descriptors (class MLIAPData *data){
    int n = data->nlistatoms;
    int* numneighs  = data->numneighs; // Number of neighbors for an atom
    int* iatoms = data->iatoms; //Index of atom
    int* jatoms = data->jatoms; //Index of each neighbor
    int* jelems = data->jelems; //Elements of each neighbor;
    int re;
    double Fn;
    double descriptor;
    double euler = 2.71828;
    double alpha;   
    double rc, dr;
    double sij;

    for (int ii = 0; ii<n ; ii++){
        int atomi = iatoms[ii]; //Global atom index for the atom is now stored in atomi
        int neighbori = numneighs[ii]; //Neighbor of atomi
        int *jelems = data->lmp_firstneigh[ii]; // This contains the neighbors for atomi

    for (int jj =0; jj<neighbori;jj++){
        int jindex = jelems[jj];
        double rij = data->rij[ii][jj];

         if (ii != jj){
             descriptor += pow((rij/re),n) * pow(euler,alpha * (rij/re)) * compute_cutoff(rij,rc, dr) * sij; 
         }
    }
    }
}

void MLIAPDescriptorRANN::read_paramfile(char *paramfilename){
// Initializing the flags  to zero. This is the checkmark for the keywords in the file

    int atomtypes_flag = 0;
    int mass_flag = 0;
    int fingerprintsperelement_flag = 0;
    int fingerprints_flag = 0;
    int fingerprintconstants_flag = 0;

for (int i =0; i<nelements; i++) delete[] elements[i];
delete [] elements;
memory->destroy(radelem);
memory->destroy(cutsq);
memory->destroy(wjelem);


//Open rann param file on proc 0;

FILE *fpparam; 

if (comm->me == 0 ){
    fpparam = utils::open_potential(paramfilename,lmp,nullptr);
    if(fpparam == nullptr){
        error->one(FLERR, "Cannot open SO3 parameter file {}: {}", paramfilename, utils::getsyserror);
    }
}

char line [MAXLINE] = {'\0'};
char *ptr;
int eof = 0;
int n , nwords;


while(true){
    if(comm->me == 0){
    ptr = utils::fgets_trunc(line,MAXLINE, fpparam);
    if(ptr == nullptr){
        eof = 1;
        fclose(fpparam);
    }
    n = strlen(line)+1;
    }

MPI_Bcast(&eof,1,MPI_INT,0, world);
if (eof) break;
MPI_Bcast(&n,1,MPI_INT,0,world);
MPI_Bcast(line,1,MPI_CHAR,0,world);

if(ptr == strchr(line, '#')) *ptr = '\0';
 nwords = utils::count_words(line);
 if(nwords == 0) continue;

 Tokenizer words(line, "' \t\n\r\f");
 std::string keyword = words.next();

if (comm->me != 0){
    utils::logmesg(lmp, "RANN keyword {}{} \n", keyword);
}

if (keyword == "atomtypes"){
    for (int iatom = 0; iatom < nelements; iatom++){ //NEED TO CHECK THIS
        elements[iatom] = utils::strdup(words.next());
        atomtypes_flag = 1;
    }
}
if(keyword == "mass"){
    for (int iatom = 0; iatom < nelements ; iatom++){
        wjelem[iatom] = utils::numeric(FLERR,words.next(), false, lmp); // NEED TO CHECK IF wj_elem IS THE VECTOR OF MASS
        mass_flag = 1 ;
    }
}
//HOLD FINGERPRINTS PER ELEMENT
if (keyword == "fingerprintsperelement"){
    for (int i = 0; i<nelements; i++){
        fingerprintsperelement[i] = utils::inumeric(FLERR, words.next(), false, lmp);
       }
    words.next();
    for(int i=0; i<nelements; i++){
        for (int j=0; j<fingerprintsperelement[i]; j++){
            fingerprints[i][j] = utils::strdup(words.next()); 
           keyword = words.next();
           if(keyword != "fingerprints"){
            break;
           }
        }
    } 
    }

if(keyword == "fingerprintconstants"){
        
}

}


}



#endif