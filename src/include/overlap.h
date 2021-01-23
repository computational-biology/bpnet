#ifndef OVERLAPLIB_OVERLAP_C
#define OVERLAPLIB_OVERLAP_C

#include <iostream>
//#include <iomanip>
#include <fstream>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <cstdlib>
#include <cstdio>
#include <sstream>

//#include <set>
#include <iterator>
#include <ctime>

//#include <gsl/gsl_spmatrix.h>
//#include <igraph/igraph.h>






//#include "cblas.h"


#define PI 3.141592654

#define MAX_ATOM_IN_A_BASE 1000
#define ATOM_TYPE 13
#define NUM_ATOM_IN_RESIDUE_CLASS

using namespace std;
using std::showpoint;




//enum ResidueName {ALA = 1, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS, PHE, PRO, SER, THR, TRP, TYR, VAL };

enum OvlpNitrogenousBase {Adenine=1,Guanine,Cytosine,Uracil};  // Now it is only Canonical.

enum OvlpAtomVariant {C1_p=1,C2,C4,C5,C6,C8,N1,N2,N3,N4,N6,N7,N9,O2,O4,O6};  // C1_p is either C1' or C1*

enum OvlpAtomType {Carbon = 6, Nitrogen =7, Oxygen = 8,Unknown_Atom=190};

enum OvlpFlag {Init=0,FullSkip,OneSkip,FullOverlap,SingleInclude,WholeSurfaceSkip};



namespace OvlpGen {
      const double epsilon = 0.000001;
      const char* MONTHS[] = {
              "January", "February", "March", "April", "May", "June",
              "July", "August", "September", "October", "November", "December"
      };
      //const double PI = 3.141592654;
      void check_sysconfig(){
      }
      string get_base_file_name(std::string in){
            return in.substr(0, in.find_last_of("."));
      }
      int diff(int x,int y){
            if(x > y){
                  return x - y;
            }else{
                  return y - x;
            }
      }
      std::string today(){
            std::ostringstream oss;
            time_t     rawtime;
            struct tm* timeinfo;

            time( &rawtime );
            timeinfo = localtime( &rawtime );
            //string(13)
            oss<<timeinfo->tm_mday <<" "<< MONTHS[ timeinfo->tm_mon ] << " "<< (timeinfo->tm_year + 1900);
            return oss.str();
      }


      //enum OvlpAtomVariant {C1_p=1,C2,C4,C5,C6,C8,N1,N2,N3,N4,N6,N7,N9,O2,O4,O6};
}


namespace OvlpParameters {
      const std::string file_version = "00.00.02";
      std::string pdb_accn="";
      const double rad_tolerance   = 0.20;

      const double carbon_radius   = 1.908;
      const double nitrogen_radius = 1.824;
      const double oxygen_radius   = 1.6612;

      const double carbon_tolerant_radius   = carbon_radius   + OvlpParameters::rad_tolerance;
      const double nitrogen_tolerant_radius = nitrogen_radius + OvlpParameters::rad_tolerance;
      const double oxygen_tolerant_radius   = oxygen_radius   + OvlpParameters::rad_tolerance;
      const int num_atom_types = 3;

      const OvlpAtomType valid_atoms[] = {Carbon,Nitrogen,Oxygen};

      const double epsilon = 0.0001;
      const int dots_per_angstrom_sqr  = 30;
      const int cluster_cardinality    = 300;
      double ovlp_cutoff = 0.0001;
};

bool Ovlp_is_equal_radius(double x,double y){
      return fabs(x - y) < OvlpParameters::epsilon ? true : false;
}


char get_char_from_nucleobase(OvlpNitrogenousBase nb){
      if(nb == Adenine)  return 'A';
      if(nb == Guanine)  return 'G';
      if(nb == Cytosine) return 'C';
      if(nb == Uracil)   return 'U';
      else{
            cerr<<"Error..... Invalid Nitrogeneous base    "<<nb<<endl;
            exit(1);
      }

}

OvlpAtomType get_base_atom_from_radius(double radius){
      if(Ovlp_is_equal_radius(radius, OvlpParameters::carbon_radius))  return Carbon;
      if(Ovlp_is_equal_radius(radius, OvlpParameters::nitrogen_radius))  return Nitrogen;
      if(Ovlp_is_equal_radius(radius, OvlpParameters::oxygen_radius)) return Oxygen;
      else{
            cerr<<"Error..... Invalid radius    "<<radius<<endl;
            exit(1);
      }

}




double ovlp_get_atom_radius(char atom_type){
      if(atom_type == 'C') return OvlpParameters::carbon_radius;
      if(atom_type == 'N') return OvlpParameters::nitrogen_radius;
      if(atom_type == 'O') return OvlpParameters::oxygen_radius;
      else{
            cerr<<"Error..... No such atom type found    "<<atom_type<<endl;
            exit(1);
      }
}


OvlpAtomType get_atom_base_type(char atom_type[]){
      if(strlen(atom_type)>=3) return Unknown_Atom;
      switch(atom_type[0]){
            case 'C' : switch(atom_type[1]){
                        case '2' :
                        case '4' :
                        case '5' :
                        case '6' :
                        case '8' : return Carbon;
                              //case '1' : if(atom_type[2] == '\'' || atom_type[2] == '*') return C1_p;
                        default  : return Unknown_Atom;
                  }
            case 'N' : switch(atom_type[1]){
                        case '1' :
                        case '2' :
                        case '3' :
                        case '4' :
                        case '6' :
                        case '7' :
                        case '9' : return Nitrogen;
                        default  : return Unknown_Atom;
                  }
            case 'O' : switch(atom_type[1]){
                        case '2' :
                        case '4' :
                        case '6' : return Oxygen;
                        default  : return Unknown_Atom;
                  }
            default : return Unknown_Atom;
      }
}

OvlpAtomVariant get_atom_type(char atom_type[]){
      switch(atom_type[0]){
            case 'C' : switch(atom_type[1]){
                        case '2' : return C2;
                        case '4' : return C4;
                        case '5' : return C5;
                        case '6' : return C6;
                        case '8' : return C8;
                        case '1' : if(atom_type[2] == '\'' || atom_type[2] == '*') return C1_p;
                        default  : cerr<<"Error.... No atom type found with name "<<atom_type<<endl;
                              exit(1);
                  }
            case 'N' : switch(atom_type[1]){
                        case '1' : return N1;
                        case '2' : return N2;
                        case '3' : return N3;
                        case '4' : return N4;
                        case '6' : return N6;
                        case '7' : return N7;
                        case '9' : return N9;
                        default  : cerr<<"Error.... No atom type found with name "<<atom_type<<endl;
                              exit(1);
                  }
            case 'O' : switch(atom_type[1]){
                        case '2' : return O2;
                        case '4' : return O4;
                        case '6' : return O6;
                        default  : cerr<<"Error.... No atom type found with name "<<atom_type<<endl;
                              exit(1);
                  }
            default : cerr<<"Error.... No atom type found with name "<<atom_type<<endl;
                  exit(1);
      }
}




class OverlapResidueClass;
class OverlapAtom;

template <typename T>
class OverlapStack {
private:


      T* elem;
      int top;
      int size;
public:
      OverlapStack(int size){
            this->elem = new T[size];
            this->size = size;
            this->top  = -1;
      }
      ~OverlapStack(){
            //cout<<"Stack deleted"<<endl;
            delete[] elem;
      }

      bool isFull(){
            return top == size-1? true:false;
      }
      bool isEmpty(){
            return this->top == -1? true : false;
      }
      void push(T elem){
            assert(this->top < size-1);
            this->elem[++this->top] = elem;
      }
      T pop(){
            assert(this->top > -1);
            return this->elem[this->top--];
      }
      int get_no_elements(){
            return top + 1;
      }
};

template <typename T>
class OverlapLinkedList {
public:
      T data;
      OverlapLinkedList * next;
public:
      OverlapLinkedList(T data){
            this->data = data;
            this->next = NULL;
      }
      OverlapLinkedList * insert_at_begining(T data){
            OverlapLinkedList<T>* tmp = new OverlapLinkedList<T>(data);
            tmp->next = this;
            return tmp;
      }
      void display(){
            OverlapLinkedList * tmp=this;;
            while(tmp != NULL){
                  cout<<tmp->data;
                  tmp = tmp->next;
            }
      }
      friend class FileData;
      friend class ResidueClass;
      friend class NUtil;
};




class OvlpOutFileBaseInfo {
      int mCorBaseSerial;
      int mPdbBaseSerial;
      std::string mBaseName;
      std::string mChainName;
      std::string mInsCode;

      std::string mBasePairName;
      std::string mBasePairType;
      double      mPairValue;

public:
      ~OvlpOutFileBaseInfo(){
            ;//cout<<"OvlpOutFileBaseInfo deleted"<<endl;
      }
      int get_cor_serial(){
            return mCorBaseSerial;
      }
      int get_pdb_Serial(){
            return mPdbBaseSerial;
      }
      std::string get_base_name(){
            return mBaseName;
      }
      std::string get_chain_name(){
            return mChainName;
      }
      std::string get_pair_name(){
            return mBasePairName;
      }
      std::string get_pair_type(){
            return mBasePairType;
      }
      double get_pair_value(){
            return mPairValue;
      }
      OvlpOutFileBaseInfo(int tmp_cor_srial,int tmp_pdb_serial,std::string tmp_base_name,std::string tmp_ins_code, std::string tmp_chain_name,std::string tmp_base_pair_name,std::string tmp_base_pair_type,double tmp_value){
            mCorBaseSerial  = tmp_cor_srial;
            mPdbBaseSerial  = tmp_pdb_serial;
            mBaseName       = tmp_base_name;
            mInsCode        = tmp_ins_code;
            mChainName      = tmp_chain_name;
            mBasePairName   = tmp_base_pair_name;
            mBasePairType   = tmp_base_pair_type;
            mPairValue      = tmp_value;
      }
      void fprint(FILE* fp){
            fprintf(fp, "  %6d  %6d  %3s  %3s  %4s  %2s  %6.2lf",
                    mCorBaseSerial,
                    mPdbBaseSerial,
                    mBaseName.c_str(),
                    mChainName.c_str(),
                    mBasePairName.c_str(),
                    mBasePairType.c_str(),
                    mPairValue);
      }

};

class OvlpOutFileRow {
      int mCorBaseSerial;
      int mPdbBaseSerial;
      std::string mBaseName;
      std::string mChainName;
      std::string mInsCode;
      OvlpOutFileBaseInfo ** mPairs;
      int mNumPairs;
public:
      std::string get_pair_name(int from_cor_serial, int to_cor_srtial){
            assert(from_cor_serial == mCorBaseSerial);
            std::string pair_name;
            for(int i=0;i<mNumPairs;i++){
                  if(mPairs[i]->get_cor_serial() == to_cor_srtial){
                        return mPairs[i]->get_pair_name();
                  }
            }
            return "NOT_FOUND";
      }

      std::string get_pair_type(int from_cor_serial, int to_cor_srtial){
            assert(from_cor_serial == mCorBaseSerial);
            std::string pair_type;
            for(int i=0;i<mNumPairs;i++){
                  if(mPairs[i]->get_cor_serial() == to_cor_srtial){
                        return mPairs[i]->get_pair_type();
                  }
            }
            return "NOT_FOUND";
      }
      int get_cor_serial(){
            return mCorBaseSerial;
      }
      int get_pdb_serial(){
            return mPdbBaseSerial;
      }
      std::string get_base_name(){
            return mBaseName;
      }
      std::string get_chain_name(){
            return mChainName;
      }
      int get_number_of_pairs(){
            return mNumPairs;
      }
      OvlpOutFileBaseInfo * get_base_pair_info(int index){
            return mPairs[index];
      }

      ~OvlpOutFileRow(){
            //cout <<"From OvlpOutFileRow OvlpOutFileBaseInfo deleted"<<endl;
            for(int i=0; i<mNumPairs; ++i){
                  delete mPairs[i];
            }
            //cout <<"From OvlpOutFileRow mPairs deleted"<<endl;
            delete [] mPairs;
            //cout<<"OvlpOutFileRow deletd"<<endl;
      }

      OvlpOutFileRow(std::string line){
            stringstream ss;
            ss<<line;
            ss>>mCorBaseSerial;
            ss>>mPdbBaseSerial;
            ss>>mBaseName;
            ss>>mInsCode;
            ss>>mChainName;
            mNumPairs = (line.length() - 25)/36;  //34 is the numor characters in the outfile for every base pair info.
            if(mNumPairs == 0){
                  mPairs = NULL;
            }else{
                  mPairs = new OvlpOutFileBaseInfo *[mNumPairs];
                  for(int i=0;i<mNumPairs;i++){
                        int start = 26+ i*36;
                        std::string base_line = line.substr(start, start+36);

                        int tmp_cor_srial;
                        int tmp_pdb_serial;
                        std::string tmp_base_name;
                        std::string tmp_chain_name;
                        std::string tmp_ins_code;
                        std::string tmp_base_pair_name;
                        std::string tmp_base_pair_type;
                        double tmp_value;
                        ss<<base_line;
                        ss>>tmp_cor_srial;
                        ss>>tmp_pdb_serial;
                        ss>>tmp_base_name;
                        ss>>tmp_ins_code;
                        ss>>tmp_chain_name;
                        ss>>tmp_base_pair_name;
                        ss>>tmp_base_pair_type;
                        ss>>tmp_value;
                        mPairs[i]= new OvlpOutFileBaseInfo(tmp_cor_srial,tmp_pdb_serial,tmp_base_name,tmp_ins_code, tmp_chain_name,tmp_base_pair_name,tmp_base_pair_type,tmp_value);
                  }
            }
      }


      void fprint(FILE* fp){
            fprintf(fp,"%6d  %6d  %3s  %3s",mCorBaseSerial,mPdbBaseSerial,mBaseName.c_str(),mChainName.c_str());
            for(int i=0;i<mNumPairs;i++){
                  mPairs[i]->fprint(fp);
            }
            fprintf(fp,"\n");
      }
};

class OvlpOutFileRowArray {
      OvlpOutFileRow ** mOutFileRowArray;
      int mNumRows;
public:
      ~OvlpOutFileRowArray(){
            for(int i=0; i<mNumRows; ++i){
                  delete mOutFileRowArray[i];
            }
            //cout<<"mOutFileRowArray deleted"<<endl;
            delete[] mOutFileRowArray;
            //cout<<"OvlpOutFileRowArray deleted"<<endl;
      }
      OvlpOutFileRowArray(OverlapStack<OvlpOutFileRow *> *stack){
            mNumRows            = stack->get_no_elements();
            mOutFileRowArray    = new OvlpOutFileRow *[mNumRows];
            for(int i=mNumRows-1;i>=0;i--){
                  mOutFileRowArray[i] = stack->pop();
                  //mOutFileRowArray[i]->fprint(stdout);
            }
            assert(stack->isEmpty() == true);
      }
      int get_num_rows(){
            return mNumRows;
      }
      OvlpOutFileRow * get_outrow(int array_index){
            return mOutFileRowArray[array_index];
      }
};


OvlpOutFileRowArray * populate_array_from_out_file(std::string out_filename){
      std::string line;
      ifstream infile_out;
      stringstream ss;
      infile_out.open(out_filename.c_str(), ios_base::in);
      //except(infile_out.is_open() == false, "Unable to open file "+out_filename);
      int num_residues = 0;
      while(infile_out.eof() != true){
            getline(infile_out,line);
            //cout<<line<<endl;
            if(line.substr(0,7) == "#HEADER"){
                  if(line.substr(0, 37) == "#HEADER   Cleaned number of residues:"){
                        ss << line.substr(37, 7);
                        ss>>num_residues;
                        //cout<< "Cleaned number of Residues in .out file = "<<num_residues<<endl;
                  }
            }
            if(line[0] != '#'){
                  break;
            }
      }
      OverlapStack<OvlpOutFileRow *> stack = OverlapStack<OvlpOutFileRow *>(num_residues);

      for(int i=0;i<num_residues;i++){
            //cout<<line.length()- 21<<" : "<<line<<endl;
            OvlpOutFileRow * out_row = new OvlpOutFileRow(line);
            stack.push(out_row);
            //out_row->fprint(stdout);
            //out_entry[i] = new OvlpOutFileRow(line);
            //out_entry[i]->print();

            getline(infile_out,line);
      }
      //cout<<"OverlapStack size "<<stack.get_no_elements()<<endl;
      OvlpOutFileRowArray * out_array = new OvlpOutFileRowArray(&stack);
      return out_array;
}


class OvlpFileData {
public:
      int base1;
      int base1_pdb;
      std::string base1_name;
      std::string base1_chain;
      int base2;
      int base2_pdb;
      std::string base2_name;
      std::string base2_chain;
      std::string pair_name;
      std::string pair_type;
      double overlap;
      ~OvlpFileData(){
            cout<<"OvlpFileData deleted"<<endl;
      }
      OvlpFileData(int base1,int base1_pdb,std::string base1_name,std::string base1_chain,
               int base2,int base2_pdb,std::string base2_name,std::string base2_chain,
               std::string pair_name,
               std::string pair_type,
               double overlap
      ){
            this->base1         = base1;
            this->base1_pdb     = base1_pdb;
            this->base1_name    = base1_name;
            this->base1_chain   = base1_chain;

            this->base2         = base2;
            this->base2_pdb     = base2_pdb;
            this->base2_name    = base2_name;
            this->base2_chain   = base2_chain;
            this->pair_name     = pair_name;
            this->pair_type     = pair_type;

            this->overlap       = overlap;



      }
      void fprint(FILE* fp){
            fprintf(fp,"OVLP    %6d:%-6d  %6d:%-6d   %3s:%-3s   %3s-%-3s  %4s  %2s  : %8.2lf",
                    base1,
                    base2,
                    base1_pdb,
                    base2_pdb,
                    base1_name.c_str(),
                    base2_name.c_str(),
                    base1_chain.c_str(),
                    base2_chain.c_str(),
                    pair_name.c_str(),
                    pair_type.c_str(),
                    overlap);
            /*    fprintf(fp,"      (");
                for(int i=0;i<mNumBase1_atoms;i++){
                    fprintf(fp,"%s ",this->base1_atom_names[i].c_str());
                }
                fprintf(fp,":");
                for(int i=0;i<mNumBase2_atoms;i++){
                    fprintf(fp,"%s ",this->base2_atom_names[i].c_str());
                }
                fprintf(fp,")");*/
            /*    map<string,OverlapLinkedList<string>*> :: iterator itr;
                fprintf(fp, "    [ ");
                for(itr = overlap_atom_map->begin(); itr != overlap_atom_map->end(); ++itr){
                    fprintf(fp," ( %s:",itr->first.c_str());
                    OverlapLinkedList<string>* sll = itr->second;
                    while(sll != NULL){
                        fprintf(fp,"%s ",sll->data.c_str());
                        sll=sll->next;
                    }
                    fprintf(fp,")");
                }
                fprintf(fp,"]");*/
            fprintf(fp,"\n");
      }


};






class OvlpPoint3D {
protected:
      double x;
      double y;
      double z;
public:
      OvlpPoint3D(double x, double y, double z){
            this->x = x;
            this->y = y;
            this->z = z;
      }

      double dist(OvlpPoint3D * p){
            return sqrt ( (x-p->x) * (x-p->x) + (y-p->y) * (y-p->y) + (z-p->z) * (z-p->z) );

      }
      double dist_sqr(OvlpPoint3D * p){
            return (x-p->x) * (x-p->x) + (y-p->y) * (y-p->y) + (z-p->z) * (z-p->z) ;

      }
      friend ostream& operator << (ostream &out, const OvlpPoint3D * p){
            out <<"("<<p->x<<","<<p->y<<","<<p->z<<")";
            //cout<<"pointer version\n";
            return out;
      }
      OvlpPoint3D * operator + (OvlpPoint3D * p){
            return new OvlpPoint3D(x+p->x,y+p->y,z+p->z);
      }
      friend ostream& operator << (ostream &out, const OvlpPoint3D & p){
            out <<"("<<p.x<<","<<p.y<<","<<p.z<<")";
            cout<<"reference version\n";
            return out;
      }
      double getX(){
            return x;
      }
      double getY(){
            return y;
      }
      double getZ(){
            return z;
      }
};

class OvlpSphere :public OvlpPoint3D {
protected:
      double r;
public:
      OvlpSphere(double x, double y, double z, double r): OvlpPoint3D(x,y,z) {
            this->r = r;
      }
      double getR(){
            return r;
      }
      bool isFullyNonTouching1(OvlpSphere * sph){
            if(this->dist_sqr(sph) > (this->r + sph->r) * (this->r + sph->r) ) return true;
            return false;
      }
      bool isFullyInscribedIn1(OvlpSphere * sph){
            if(this->dist_sqr(sph) < (sph->r - this->r)*(sph->r - this->r)) return true;
            return false;
      }
      OvlpSphere * scale_sphere(OvlpPoint3D * pt){
            return new OvlpSphere(x+pt->getX(),y+pt->getY(),z+pt->getZ(),r);
      }
      bool isIntersecting(OvlpSphere * sph){
            double redii_sum = this->r + sph->r;
            return this->dist_sqr(sph) <= redii_sum * redii_sum ? true : false ;
      }
      double surface_gap(OvlpSphere * sph){
            double center_dist = this->dist(sph);
            double redii_sum   = this->r + sph->r;
            return center_dist - redii_sum;
      }
      double dist_between_centers(OvlpSphere * sph){
            return this->dist(sph);
      }
      double surface_area(){
            return 4*PI*r*r;
      }
      double volume(){
            return 4.0 * (PI * r * r * r ) / 3.0 ;
      }
      double spherical_cap_curved_surface_area(double h){    //h is cap height. Ref: https://en.wikipedia.org/wiki/Spherical_cap
            return 2 * PI * r * h;
      }
      double spherical_cap_volume(double h){
            return PI * h * h * (3.0*r - h) / 3.0;
      }
      double sphere_sphere_intersection_plane_dist(OvlpSphere * sph){  // http://mathworld.wolfram.com/Sphere-SphereIntersection.html
            assert(this->isIntersecting(sph) == true);
            double d = this->dist_between_centers(sph);
            double x = (d * d - sph->r * sph->r + this->r * this->r )/ (2.0 * d);
            return x;
      }
      double sphere_sphere_intersection_cap_height(OvlpSphere * sph){
            double x = this->sphere_sphere_intersection_plane_dist(sph);
            return this->r - x ;
      }
      double sphere_sphere_union_surface_area(OvlpSphere * sph){
            //double area;
            double this_h        = this->sphere_sphere_intersection_cap_height(sph);
            double sph_h         = sph->sphere_sphere_intersection_cap_height(this);
            double this_cap_area = this->spherical_cap_curved_surface_area(this_h);
            double sph_cap_area  = sph->spherical_cap_curved_surface_area(sph_h);
            double this_area     = this->surface_area();
            double sph_area      = sph->surface_area();
            return (this_area+sph_area)-(this_cap_area+sph_cap_area) ;
      }
};








class OvlpSurfaceCluster {
public:
      OvlpSphere * smallestSphere;
      int num_points;
      //double max_spread;
      OverlapLinkedList<OvlpPoint3D *>* point_cluster;

      OvlpSurfaceCluster(OvlpPoint3D * point[],int num_points,int position, double max_spread){
            this->num_points = num_points;
            //this->max_spread = max_spread;
            smallestSphere = new OvlpSphere(point[position]->getX(),point[position]->getY(),point[position]->getZ(),max_spread);// = point[position];
            //cout<<"...  num ...points"<<num_points<<endl;
            point_cluster = new OverlapLinkedList<OvlpPoint3D *>(point[num_points-1]);
            for(int i=num_points-2;i>=0;i--){
                  point_cluster = point_cluster->insert_at_begining(point[i]);
            }
      }
      friend ostream& operator << (ostream &out, const OvlpSurfaceCluster * s){
            //out <<"[mid ="<<s->mid_point<<", no pts="<<s->num_points<<", max spread ="<<s->max_spread<<"]\n";
            //s->point_cluster->display();
            cout<<"OvlpSphere cout not properly overwriten...";
            return out;
      }
      int get_num_points(){
            return num_points;
      }
};




class OvlpUtil {
public:
      static string get_base_file_name(std::string in){
            return in.substr(0, in.find_last_of("."));
      }

      static OvlpAtomType get_atom_with_tolerant_radius(double radius){
            if(fabs(radius - OvlpParameters::carbon_tolerant_radius)   < OvlpParameters::epsilon)   return Carbon;
            if(fabs(radius - OvlpParameters::nitrogen_tolerant_radius) < OvlpParameters::epsilon)   return Nitrogen;
            if(fabs(radius - OvlpParameters::oxygen_tolerant_radius)   < OvlpParameters::epsilon)   return Oxygen;
            else{
                  cerr<<"Error..... Invalid radius    "<<radius<<endl;
                  exit(1);
            }
      }
      static bool is_valid_radius_with_tolerance(double radius){
            if(fabs(radius - OvlpParameters::carbon_tolerant_radius)   < OvlpParameters::epsilon) return true;
            if(fabs(radius - OvlpParameters::nitrogen_tolerant_radius) < OvlpParameters::epsilon) return true;
            if(fabs(radius - OvlpParameters::oxygen_tolerant_radius)   < OvlpParameters::epsilon) return true;
            return false;
      }
      static double get_radius_with_tolerance(OvlpAtomType atom){
            if(atom == Carbon)   return OvlpParameters::carbon_tolerant_radius;
            if(atom == Nitrogen) return OvlpParameters::nitrogen_tolerant_radius;
            if(atom == Oxygen)   return OvlpParameters::oxygen_tolerant_radius;
            else{
                  cerr<<"Invalid atom type supplied in get_radius_with_tolerance() function."<<endl;
                  exit(15);
            }
      }
      static bool is_valid_atom(OvlpAtomType atom){
            if(atom == Carbon || atom == Nitrogen || atom == Oxygen) return true;
            else return false;
      }
      static void get_mid_point(OvlpPoint3D * point_array[], int size, int* position, double* radius){
            double matrix[size][size];
            for(int i=0;i<size;i++){
                  for(int j=i;j<size;j++){
                        if(i == j)
                              matrix[i][j] = 0.0;
                        else{
                              matrix[i][j] = point_array[i]->dist(point_array[j]);
                              matrix[j][i] = matrix[i][j];
                        }
                  }
            }

            double min_radius = 50000.0;

            int pos = -1;

            for(int i=0;i<size;i++){
                  double max_radius = 0.0;
                  for(int j=0;j<size;j++){
                        if(matrix[i][j]>=max_radius) max_radius = matrix[i][j];
                  }
                  if(max_radius<min_radius){
                        min_radius = max_radius;
                        pos = i;
                  }
            }
            *position = pos;
            *radius   = min_radius+ OvlpParameters::epsilon;
      }
      static char to_char_from_Nitrogenoue_base(OvlpNitrogenousBase nbase){
            if(nbase == Uracil)         return 'U';
            else if(nbase == Adenine)   return 'A';
            else if(nbase == Guanine)   return 'G';
            else if(nbase == Cytosine)  return 'C';
            else{
                  cerr<<"Error... Invalid Nitrogenous base number number"<<nbase;
                  exit(11);
            }
      }
};






class OvlpAtomSurfPoints {
private:
      double radius;
      int num_surf_pts;
      //OvlpPoint3D** surf_points;
      OvlpAtomType atom;
      OverlapLinkedList<OvlpSurfaceCluster *>* point_clusters;
      int num_nodes;
public:
      void display(){
            cout<<"OverlapAtom rad : "<<radius<<", no of surf pts : "<<num_surf_pts<<endl;
            point_clusters->display();
      }
      int get_total_surface_points(){
            return num_surf_pts;
      }
      ~OvlpAtomSurfPoints(){

            ; //cout<< "OvlpAtomSurfPoints destructor invoked"<<endl;
      }
      OvlpAtomSurfPoints(int radius, OvlpAtomType atom, OverlapStack<OvlpPoint3D *>* pointStack){
            this->radius        = radius;
            this->atom          = atom;
            this->num_surf_pts  = pointStack->get_no_elements();
            //cout<<"working... fine  "<<num_surf_pts<<endl;
            this->num_nodes = num_surf_pts / OvlpParameters::cluster_cardinality;
            OverlapStack<OvlpSurfaceCluster *> surf_stack = OverlapStack<OvlpSurfaceCluster *>(5000);

            int last_node_points = num_surf_pts % OvlpParameters::cluster_cardinality;

            if(last_node_points > 0){
                  //cout<<"> "<<last_node_points<<endl;
                  num_nodes ++;
                  OvlpPoint3D ** last_points = new OvlpPoint3D *[last_node_points];
                  for(int i=0;i<last_node_points;i++){
                        last_points[i] = pointStack->pop();
                        //cout<<"Working...1"<<endl;
                  }

                  int pos = -1;
                  double rad = 0.0;
                  OvlpUtil::get_mid_point(last_points,last_node_points,&pos,&rad);
                  if(pos==-1){
                        cout<<"position cannot be -1"<<endl;
                        exit(12);
                  }
                  OvlpSurfaceCluster * surf_clust = new OvlpSurfaceCluster(last_points,last_node_points,pos,rad);
                  surf_stack.push(surf_clust);
            }
            //cout<<"working... fine  "<<num_nodes<<endl;
            //cout<<"Working...2   num_nodes"<<num_nodes<<endl;
            int loop_var= last_node_points>0?num_nodes-1:num_nodes;
            for(int i=0;i<loop_var;i++){
                  //OvlpSurfaceCluster* = new OvlpSurfaceCluster();
                  OvlpPoint3D ** points = new OvlpPoint3D *[OvlpParameters::cluster_cardinality];
                  for(int j=0;j< OvlpParameters::cluster_cardinality;j++){
                        points[j] = pointStack->pop();
                  }

                  //assert(pointStack->isEmpty());
                  int pos = -1;
                  double rad = 0.0;
                  OvlpUtil::get_mid_point(points, OvlpParameters::cluster_cardinality,&pos,&rad);
                  if(pos==-1){
                        cout<<"Error... position cannot be -1"<<endl;
                        exit(13);
                  }
                  //cout<<"Working...3  points: "<<points[0]<<", "<<points[pos]<<endl;
                  OvlpSurfaceCluster * surf_clust = new OvlpSurfaceCluster(points, OvlpParameters::cluster_cardinality,pos,rad);
                  surf_stack.push(surf_clust);
            }
            point_clusters = new OverlapLinkedList<OvlpSurfaceCluster *>(surf_stack.pop());
            while(surf_stack.isEmpty() != true){
                  point_clusters = point_clusters->insert_at_begining(surf_stack.pop());
            }


            //cout<<"just before point cluster"<<endl;
            //point_clusters->display();

/*        this->surf_points   = new OvlpPoint3D*[num_surf_pts];
        for(int i=num_surf_pts-1; i>=0; i--){
            this->surf_points[i] = pointStack->pop();
        }*/


      }
      OvlpAtomType get_atom_type(){
            return this->atom;
      }
      OverlapLinkedList<OvlpSurfaceCluster *>* get_list(){
            return this->point_clusters;
      }
      void test(){
            cout<<"Testing of non- nullity...  atom"<<atom<<endl;
      }

      friend class ResidueClass;

};



class OverlapAtom :public OvlpSphere {
private:
      int serial;
      OvlpAtomVariant type;
      std::string atom_name;
      int occupancy;
      double temp_factor;
      double radius_with_tolerance;
      OvlpAtomSurfPoints * surf_points;
public:
      std::string get_atom_name(){
            return atom_name;
      }
      void print(){
            cout<<"ATOM "<<serial<<endl;
      }
      OverlapAtom(int serial, OvlpAtomVariant type, std::string atom_name, double x, double y, double z, double r, double occupancy, double temp_factor, OvlpAtomSurfPoints * surf_pts): OvlpSphere(x,y,z,r){
            this->serial        = serial;
            this->type          = type;
            this->atom_name = atom_name;
            this->occupancy     = occupancy;
            this->temp_factor   = temp_factor;
            this->surf_points   = surf_pts;
            this->radius_with_tolerance = r + OvlpParameters::rad_tolerance;
      }
      bool isFullyNonTouchingWith(OvlpSphere * sph){
            //if(this->dist_sqr(sph) > (this->radius_with_tolerance + sph->getR()) * (this->radius_with_tolerance + sph->getR()) ) return true;
            if(this->dist(sph)>this->radius_with_tolerance+sph->getR()) return true;
            return false;
      }
      bool FullyContains(OvlpSphere * sph){
            //if(this->dist_sqr(sph) < (this->radius_with_tolerance - sph->getR()) * (this->radius_with_tolerance - sph->getR() )) return true;
            if(this->dist(sph) < this->radius_with_tolerance - sph->getR()) return true;
            return false;
      }
      bool isFullyNonTouchingWith(OverlapAtom * a){
            if(this->dist_sqr(a) > (this->radius_with_tolerance + a->radius_with_tolerance) * (this->radius_with_tolerance + a->radius_with_tolerance) ) return true;
            return false;
      }
      bool FullyContains(OverlapAtom * a){
            if(this->dist_sqr(a) < (this->radius_with_tolerance- a->radius_with_tolerance)*(this->radius_with_tolerance-a->radius_with_tolerance )) return true;
            return false;
      }
      double get_radius_with_tolerance(){
            return radius_with_tolerance;
      }
      int get_total_surface_points(){
            return surf_points->get_total_surface_points();
      }
      OvlpAtomSurfPoints * get_sueface_points(){
            //cout<<"!@#"<<endl;
            //surf_points->display();
            //cout<<"!@#"<<endl;
            return surf_points;
      }
      OvlpPoint3D scale_point(OvlpPoint3D * p){
            return OvlpPoint3D(x+p->getX(),y+p->getY(),z+p->getZ());
      }

      int get_serial_no(){
            return serial;
      }
      bool overlapped(OverlapAtom * a){
            if(this->dist(a)<=this->radius_with_tolerance+a->radius_with_tolerance)return true;
            else return false;
      }
};

class OvlpStringArray {
private:
      std::string* array;
      int num_elem;
public:
      ~OvlpStringArray(){
            delete[] array;
      }
      OvlpStringArray(int elements){
            num_elem    = elements;
            array       = new std::string[num_elem];
      }
      void set_string_at(std::string str,int index){
            assert(index<num_elem);
            array[index] = str;
      }
      int search(std::string str){
            for(int i=0;i<num_elem;i++){
                  if(array[i] == str){
                        return i;
                  }
            }
            return -1;
      }
};


class OvlpNUtil {
public:
      static OvlpStringArray * read_words_from_file(std::string file_name){
            ifstream infile_stream;
            //infile_stream.open(const std::string &__s, optional ios_base::openmode __mode = ios_base::in)
            infile_stream.open(file_name.c_str(), std::ios_base::in);
            if(infile_stream.is_open() == false){
                  cerr<<"Error... unable to open file "<<file_name<<endl;
                  exit(1);
            }
            std::string line;
            OverlapStack<std::string> stack = OverlapStack<std::string>(100);
            while(infile_stream.eof() != true){
                  std::getline(infile_stream, line);
                  if(infile_stream.eof() == true){
                        break;
                  }
                  stack.push(line);
                  //cout<<"LINE  "<<line<<endl;
            }
            int num_elem = stack.get_no_elements();
            OvlpStringArray * str_array = new OvlpStringArray(num_elem);
            for(int i=num_elem-1;i>=0;i--){
                  str_array->set_string_at(stack.pop(),i);
            }
            infile_stream.close();
            return str_array;
      }
      static bool valid_base_atom(std::string atom, std::string base_atoms[],int size){
            for(int i=0;i<size;i++){
                  if(atom == base_atoms[i]) return true;
            }
            return false;
      }
      /*static bool valid_atom_for_base(std::string line){
          //std::string adenine_atoms[] = {"N9","C8","N7","C5","C4","C6","N6","C2","N3","N1"};
          stringstream ss;
          ss<<line;
          std::string tag;
          int id;
          std::string atom;
          std::string base;
          ss>>tag;
          ss>>id;
          ss>>atom;
          ss>>base;
          //cout<<"Working "<<tag<<", "<<id<<", "<<atom<<", "<<base<<", Len ="<<adenine_atoms->size()<<endl;

      }*/

      static void compute_surface(OverlapAtom * atom_array[], int array_size,int* total_points, int* total_overlap_points, int* total_surface_points){
            int surf_pts    = 0;
            int overlap_pts = 0;
            //OvlpPoint3D* t_point=NULL;
            //cout<<"This is working???   1"<<endl;

            //cout<<"start of single base pair\n";
            //cout<<"residue no = "<<residue_no<<" : ";
#pragma omp parallel for reduction(+:surf_pts, overlap_pts)
            for(int i=0;i<array_size;i++){
                  surf_pts = surf_pts + atom_array[i]->get_total_surface_points();
                  int overlapped_atom_index[200];
                  int overlapped_atom_count = 0;

                  for(int j=0;j<array_size;j++){
                        if(i != j){
                              if(atom_array[i]->overlapped(atom_array[j])==true){
                                    overlapped_atom_index[overlapped_atom_count] = j;
                                    overlapped_atom_count ++;
                                    //cout<<"No overlap  "<<atom_array[i]->get_serial_no()<<" And "<<atom_array[j]->get_serial_no()<<endl;
                              }
                        }
                  }
//cout<<"This of 12345\n";
                  if(overlapped_atom_count == 0){
                        //    cout<<"This of 12345\n";
                        overlap_pts = overlap_pts + 0;
                  }else{
                        OverlapAtom * jth_atom;
                        OverlapLinkedList<OvlpSurfaceCluster *>* ith_pts = atom_array[i]->get_sueface_points()->get_list();
                        //cout<<"Working 1 "<<endl;
                        while(ith_pts != NULL){
                              //cout<<"This of 2\n";
                              //ith_pts = ith_pts->next;
                              OvlpSphere mid_sphere = OvlpSphere( ith_pts->data->smallestSphere->getX()+atom_array[i]->getX(),
                                                          ith_pts->data->smallestSphere->getY()+atom_array[i]->getY(),
                                                          ith_pts->data->smallestSphere->getZ()+atom_array[i]->getZ(),
                                                          ith_pts->data->smallestSphere->getR()+atom_array[i]->getR()
                              );

                              //OvlpSphere* mid_sphere =  ith_pts->data->smallestSphere->scale_sphere(AtomArray[i]);// AtomArray[i]->scale_point(ith_pts->data->mid_point);
                              //cout<<"This of Suspect\n";
                              int single_point_index[200];
                              int single_point_count = 0;
                              //cout<<"Working 2\n";
                              OvlpFlag flag=Init;
                              for(int j=0;j<overlapped_atom_count;j++){ // i != j is unnecessary here. it is ommitted in the previous loop.
                                    jth_atom = atom_array[overlapped_atom_index[j]];

                                    //double dist         = mid_pt_of_ith->dist(jth_atom);
                                    //double radius_tol   = jth_atom->get_radius_with_tolerance();
                                    //double spread       = ith_pts->data->max_spread;
                                    //cout<<"This of 3\n";
                                    if(jth_atom->isFullyNonTouchingWith(&mid_sphere)){
                                          flag = FullSkip;
                                    }else if(jth_atom->FullyContains(&mid_sphere)){
                                          flag = FullOverlap;
                                          break;
                                    }else{
                                          single_point_index[single_point_count] = overlapped_atom_index[j];
                                          single_point_count ++;
                                    }
                              }
                              if(flag == FullOverlap){
                                    overlap_pts = overlap_pts + ith_pts->data->get_num_points();
                                    //cout<<"Overlap  "<<ith_pts->data->get_num_points()<<endl;

                              }else if(flag == FullSkip && single_point_count == 0){
                                    overlap_pts = overlap_pts + 0;
                                    //cout<<"Skip  "<<endl;
                              }else{
                                    OverlapLinkedList<OvlpPoint3D *>* tmp = ith_pts->data->point_cluster;

                                    OverlapAtom * current_atom;

                                    while(tmp!=NULL){
                                          OvlpPoint3D surface_point = atom_array[i]->scale_point(tmp->data);
                                          OvlpFlag final_flag=Init;
                                          //cout<<"This of 2222222222\n";

                                          for(int j=0;j<single_point_count;j++){
                                                //cout<<"Single  "<<current_atom->get_serial_no()<<endl;
                                                current_atom = atom_array[single_point_index[j]];
                                                if( current_atom->dist(&surface_point) <= current_atom->get_radius_with_tolerance()){
                                                      final_flag = SingleInclude;
                                                      break;
                                                }
                                          }
                                          if(final_flag == SingleInclude){
                                                overlap_pts = overlap_pts + 1;
                                          }
                                          tmp = tmp->next;
                                    }
                              }
                              ith_pts = ith_pts->next;
                        }
                  }

                  //cout<<"Working total overlap points for a single atlom = "<<overlap_pts<<endl;

            }
            *total_points           = surf_pts;
            *total_overlap_points   = overlap_pts;
            *total_surface_points   = surf_pts - overlap_pts;
      }
};








class OvlpAllSurfacePoints {
private:
      int num_atoms;
      OvlpAtomSurfPoints ** all_surf_points;
public:
      ~OvlpAllSurfacePoints(){
            //cout<<"OvlpAtomSurfPoints array deleted from OvlpAllSurfacePoints "<<endl;
            delete [] all_surf_points;
            //cout<<"OvlpAllSurfacePoints deleted"<<endl;
      }
      OvlpAllSurfacePoints(OverlapStack<OvlpAtomSurfPoints *>* stack){
            this->num_atoms = stack->get_no_elements();
            all_surf_points = new OvlpAtomSurfPoints *[num_atoms];
            for(int i=num_atoms-1;i>=0;i--){
                  all_surf_points[i] = stack->pop();
            }
      }
      OvlpAtomSurfPoints * get_surface_points(OvlpAtomType atom){
            for(int i=0;i<num_atoms;i++){
                  //cout<<"A = "<<this->all_surf_points[i]->get_atom_type()<<", B="<<atom<<endl;

                  if(this->all_surf_points[i]->get_atom_type() == atom){
                        return this->all_surf_points[i];
                  }
            }
            cerr<<"Wrong atom type supplied";
            exit(15);
            return NULL;
      }
      void display(){
            cout<<"Types of atom : "<<num_atoms<<endl;
            all_surf_points[0]->display();

      }


};

class OvlpSurfaceDataFile {
private:
      std::string file_name;
      std::string file_base;
      std::string file_ext;
      std::ifstream   surf_file;
      int num_basic_atom;


public:
      ~OvlpSurfaceDataFile(){
            ;//cout<<"OvlpSurfaceDataFile Deleted"<<endl;
      }
      OvlpSurfaceDataFile(std::string file_name, OvlpAtomType atom_type[],int num_basic_atoms){
            file_name = file_name+"/surface.xyz";
            this->file_name      = file_name;
            this->file_base      = file_name.substr(0, file_name.find_last_of("."));
            this->file_ext       = file_name.substr(file_name.find_last_of(".") + 1, file_name.length());
            this->num_basic_atom = num_basic_atom;
            if(strncasecmp(file_ext.c_str(), "xyz",6)) {
                  cerr<<"Error    Unable to generate surface.   ."<<file_ext<<" file has been supplied insted of .xyz file"<<endl;
                  exit(1);
            }
            //cout << file_base<<endl<<file_ext<<endl;
      }
      OvlpAllSurfacePoints * generate_surface_points(int atom_variants){
            OverlapStack<OvlpAtomSurfPoints *> atom_surf_stack = OverlapStack<OvlpAtomSurfPoints *>(1000);

            //OvlpAtomSurfPoints** all_surf_points = new OvlpAtomSurfPoints*[atom_variants];
            OverlapStack<OvlpPoint3D *> PointStack= OverlapStack<OvlpPoint3D *>(20000);
            std::string     line;
            surf_file.open(file_name.c_str());
            if(surf_file.is_open() == false){
                  cerr<<"Error... unable to open file "<<file_name<<endl;
                  exit(13);
            }
            double r, x, y, z;
            cout<<showpoint;
            int flag = 0;
            surf_file>>r;
            double curr_radius = r;
            OvlpAtomType atom;
            surf_file.seekg(0, ios_base::beg);
            while(surf_file.eof() != true){
                  surf_file>>r;
                  surf_file>>x;
                  surf_file>>y;
                  surf_file>>z;

                  if(OvlpUtil::is_valid_radius_with_tolerance(r) == true ){
                        if(flag == 0){
                              curr_radius = r;
                              flag = 1;
                              //cout<<"Abc 12345...........................................   "<<r<<endl;
                        }
                        //cout<<r<<" , "<<x<<", "<<y<<","<<z<<endl;
                        //cout<<"ppp "<<curr_radius<<", "<<r<<endl;

                        atom = OvlpUtil::get_atom_with_tolerant_radius(r);
                        assert(OvlpUtil::is_valid_atom(atom) == true);
                        if(Ovlp_is_equal_radius(curr_radius, r)){
//cout<<"ppp "<<curr_radius<<", "<<r<<endl;
                              PointStack.push(new OvlpPoint3D(x,y,z));
                        }else{
                              //cout<<"at popping"<<endl;
                              //int n = PointStack.get_no_elements();
                              //cout<<"n ="<<n<<endl;
                              OvlpAtomSurfPoints * surface_points = new OvlpAtomSurfPoints(curr_radius, OvlpUtil::get_atom_with_tolerant_radius(curr_radius),&PointStack);
                              atom_surf_stack.push(surface_points);

                              assert(PointStack.isEmpty() == true);
                              PointStack.push(new OvlpPoint3D(x,y,z));
                              flag = 0;
                        }
                  }
            }
            if(PointStack.isEmpty() != true){
                  //cout<<"Finally...     "<<endl;
                  //int n = PointStack.get_no_elements();
                  //cout<<"n ="<<n<<",r="<<curr_radius<<"atom ="<<OvlpUtil::get_atom_with_tolerant_radius(curr_radius)<<endl;
                  OvlpAtomSurfPoints * surface_points = new OvlpAtomSurfPoints(curr_radius, OvlpUtil::get_atom_with_tolerant_radius(curr_radius),&PointStack);
                  atom_surf_stack.push(surface_points);
            }
            OvlpAllSurfacePoints * all_surface_points = new OvlpAllSurfacePoints(&atom_surf_stack);
            return all_surface_points;
      }
};




















class OverlapResidueClass {
private:
      int residue_no;
      int num_atom;
      //OvlpNitrogenousBase nucleobase;
      std::string nucleobase;
      OverlapAtom ** AtomArray;
      int total_overlap_points;
      int total_points;
      int actual_surface_points;

      OvlpPoint3D* center_of_mass;

      //bool AtomArrayPopulated = false; // this is because AtomArray is not populated from constructor. So, every newly added function should check whether this is true before moving further.
public:
      ~OverlapResidueClass(){
            for(int i=0; i<num_atom; ++i){
                  delete AtomArray[i];
            }
            delete center_of_mass;
            delete[] AtomArray;
            //cout<<"OverlapResidueClass "<<endl;
      }
      bool get_atom_overlap_mapping_atoms(OverlapResidueClass * base2, std::string base1_atom[], int* num_base1_atom,std::string base2_atom[], int* num_base2_atom){
            int n1 = 0;
            int n2 = 0;
            for(int i=0;i<this->num_atom;i++){
                  std::string key = this->AtomArray[i]->get_atom_name();

                  for(int j=0;j<base2->num_atom;j++){
                        if(this->AtomArray[i]->isFullyNonTouchingWith(base2->AtomArray[j]) == false){
                              std::string elem = base2->AtomArray[j]->get_atom_name();

                              bool flag = false;
                              for(int i=0;i<n1;i++){
                                    if(key == base1_atom[i]){
                                          flag = true;
                                    }
                              }
                              if(flag == false){
                                    base1_atom[n1] = key;
                                    n1++;
                              }
                              flag = false;
                              for(int i=0;i<n2;i++){
                                    if(elem == base2_atom[i]){
                                          flag = true;
                                    }
                              }
                              if(flag == false){
                                    base2_atom[n2] = elem;
                                    n2++;
                              }
                        }
                  }

            }
            *num_base1_atom = n1;
            *num_base2_atom = n2;
            if(n2>0){
                  return true;
            }else{
                  return false;
            }

      }


      void print(){
            cout<<"Residue  "<<residue_no<<"      FIRST ATOM "<<AtomArray[0]<<endl;
      }
      OverlapAtom * get_atom(int index){
            return AtomArray[index];
      }
      int get_residue_number(){
            return this->residue_no;
      }
      OverlapResidueClass(int residue_number, std::string nucleo_base, OverlapStack<OverlapAtom *> *AtomStack){ // This constructor is non-conventional. This is created because atoms are first placed in a stack when created.
            this->residue_no    = residue_number;                                              // But this is the only way one can create a residue class.
            this->nucleobase    = nucleo_base;
            this->num_atom      = AtomStack->get_no_elements();
            this->AtomArray     = new OverlapAtom *[num_atom];
            double x_min        =  9999999.0;
            double x_max        = -9999999.0;
            double y_min        =  9999999.0;
            double y_max        = -9999999.0;
            double z_min        =  9999999.0;
            double z_max        = -9999999.0;




            for(int i=num_atom-1;i>=0;i--){
                  AtomArray[i] = AtomStack->pop();
                  if(AtomArray[i]->getX()<x_min) x_min = AtomArray[i]->getX();
                  if(AtomArray[i]->getX()>x_max) x_max = AtomArray[i]->getX();

                  if(AtomArray[i]->getY()<y_min) y_min = AtomArray[i]->getY();
                  if(AtomArray[i]->getY()>y_max) y_max = AtomArray[i]->getY();

                  if(AtomArray[i]->getZ()<z_min) z_min = AtomArray[i]->getZ();
                  if(AtomArray[i]->getZ()>z_max) z_max = AtomArray[i]->getZ();
            }
            this->center_of_mass = new OvlpPoint3D( (x_min+x_max)/2.0,(y_min+y_max)/2.0,(z_min+z_max)/2.0);
      }

/*    OverlapResidueClass(int residue_number, OvlpNitrogenousBase nucleo_base, OverlapStack<OverlapAtom*> *AtomStack){ // This constructor is non-conventional. This is created because atoms are first placed in a stack when created.
        this->residue_no    = residue_number;                                              // But this is the only way one can create a residue class.
        this->nucleobase    = nucleo_base;
        this->num_atom      = AtomStack->get_no_elements();
        this->AtomArray     = new OverlapAtom*[num_atom];
        double x_min        =  9999999.0;
        double x_max        = -9999999.0;
        double y_min        =  9999999.0;
        double y_max        = -9999999.0;
        double z_min        =  9999999.0;
        double z_max        = -9999999.0;




        for(int i=num_atom-1;i>=0;i--){
            AtomArray[i] = AtomStack->pop();
            if(AtomArray[i]->getX()<x_min) x_min = AtomArray[i]->getX();
            if(AtomArray[i]->getX()>x_max) x_max = AtomArray[i]->getX();

            if(AtomArray[i]->getY()<y_min) y_min = AtomArray[i]->getY();
            if(AtomArray[i]->getY()>y_max) y_max = AtomArray[i]->getY();

            if(AtomArray[i]->getZ()<z_min) z_min = AtomArray[i]->getZ();
            if(AtomArray[i]->getZ()>z_max) z_max = AtomArray[i]->getZ();
        }
        this->center_of_mass = new OvlpPoint3D( (x_min+x_max)/2.0,(y_min+y_max)/2.0,(z_min+z_max)/2.0);
    }*/
      double compute_overlap(OverlapResidueClass * res){
            int tot_pts;
            int surf_pts;
            int ovlp_pts;
            if(this->compute_joint_surface_area(res, &tot_pts, &ovlp_pts, &surf_pts) == true) {
                  //double this_surf = this->total_points - this->total_overlap_points;
                  //double res_surf  = res->total_points  - res->total_overlap_points;
                  //double tmp_surf  = tmp_res->total_points  - tmp_res->total_overlap_points;
                  double overlap = ((double)(this->actual_surface_points + res->actual_surface_points - surf_pts))/(2.0*(double) OvlpParameters::dots_per_angstrom_sqr);
                  //double overlap = (this->total_overlap_points + res->total_overlap_points - tmp_res->total_overlap_points)/(2.0*(double)OvlpParameters::dots_per_angstrom_sqr);
                  //if(this->residue_no<5 && this->residue_no<30)
                  /*if(overlap > 5.0){
                      cout<<"Details this, other, total "<<this->actual_surface_points<<", "<<res->actual_surface_points<<", "<<surf_pts<<endl;
                      cout<<"overlap between "<<this->residue_no<<" and "<<res->residue_no<<" is  "<<OvlpUtil::to_char_from_Nitrogenoue_base(this->nucleobase)<<":"<<OvlpUtil::to_char_from_Nitrogenoue_base(res->nucleobase)<<" = "<<overlap<<endl<<endl;
                  }*/
                  //cout<<"OVLP "<<overlap<<endl;
                  return overlap;
            }
            return -1.0;
      }
      std::string get_nucleobase(){
            return this->nucleobase;
            //return get_char_from_nucleobase(this->nucleobase);
      }

      bool is_enough_close_to(OverlapResidueClass * res,double dist){
            if(this->center_of_mass->dist_sqr(res->center_of_mass) <= dist * dist) return true;
            return false;
      }






      bool compute_joint_surface_area(OverlapResidueClass * res,int* total_pts, int* total_overlap_pts, int* act_surf_pts){
            if(this->center_of_mass->dist_sqr(res->center_of_mass)>100.0){
                  return false;
            }
            OverlapAtom * atom_array[50];
            int array_size = this->num_atom+res->num_atom;
            //OvlpPoint3D* t_point=NULL;
            //cout<<"This is working???   1"<<endl;

            //cout<<"start of single base pair\n";
            //cout<<"residue no = "<<residue_no<<" : ";

            int count = 0;
            for(int i=0;i<this->num_atom;i++){
                  atom_array[count] = this->AtomArray[i];
                  count ++;
            }
            for(int i=0;i<res->num_atom;i++){
                  atom_array[count] = res->AtomArray[i];
                  count ++;
            }
            OvlpNUtil::compute_surface(atom_array, count,total_pts,total_overlap_pts, act_surf_pts);
            return true;
      }

      void compute_surface_area(){
            int total_points;
            int total_surface_points;
            int total_overlap_points;
            OvlpNUtil::compute_surface(this->AtomArray,this->num_atom, &total_points, &total_overlap_points, &total_surface_points);
            this->total_points          = total_points;
            this->actual_surface_points = total_surface_points;
            this->total_overlap_points  = total_overlap_points;
      }













      void compute_surface_area_old(){
            int surf_pts    = 0;
            int overlap_pts = 0;
            //OvlpPoint3D* t_point=NULL;
            //cout<<"This is working???   1"<<endl;

            //cout<<"start of single base pair\n";
            //cout<<"residue no = "<<residue_no<<" : ";
            for(int i=0;i<num_atom;i++){
                  surf_pts = surf_pts + AtomArray[i]->get_total_surface_points();
                  int overlapped_atom_index[50];
                  int overlapped_atom_count = 0;

                  for(int j=0;j<num_atom;j++){
                        if(i != j){
                              if(AtomArray[i]->overlapped(AtomArray[j])==true){
                                    overlapped_atom_index[overlapped_atom_count] = j;
                                    overlapped_atom_count ++;
                                    //cout<<"No overlap  "<<AtomArray[i]->get_serial_no()<<endl;
                              }
                        }
                  }
//cout<<"This of 12345\n";
                  if(overlapped_atom_count == 0){
                        //    cout<<"This of 12345\n";
                        overlap_pts = overlap_pts + 0;
                  }else{
                        OverlapAtom * jth_atom;
                        OverlapLinkedList<OvlpSurfaceCluster *>* ith_pts = AtomArray[i]->get_sueface_points()->get_list();

                        while(ith_pts != NULL){
                              //cout<<"This of 2\n";
                              //ith_pts = ith_pts->next;
                              OvlpSphere mid_sphere = OvlpSphere(ith_pts->data->smallestSphere->getX()+AtomArray[i]->getX(),
                                                         ith_pts->data->smallestSphere->getY()+AtomArray[i]->getY(),
                                                         ith_pts->data->smallestSphere->getZ()+AtomArray[i]->getZ(),
                                                         ith_pts->data->smallestSphere->getR()+AtomArray[i]->getR());

                              //OvlpSphere* mid_sphere =  ith_pts->data->smallestSphere->scale_sphere(AtomArray[i]);// AtomArray[i]->scale_point(ith_pts->data->mid_point);
                              //cout<<"This of Suspect\n";
                              int single_point_index[50];
                              int single_point_count = 0;
                              //cout<<"Working 2\n";
                              OvlpFlag flag=Init;
                              for(int j=0;j<overlapped_atom_count;j++){ // i != j is unnecessary here. it is ommitted in the previous loop.
                                    jth_atom = AtomArray[overlapped_atom_index[j]];

                                    //double dist         = mid_pt_of_ith->dist(jth_atom);
                                    //double radius_tol   = jth_atom->get_radius_with_tolerance();
                                    //double spread       = ith_pts->data->max_spread;
                                    //cout<<"This of 3\n";
                                    if(jth_atom->isFullyNonTouchingWith(&mid_sphere)){
                                          flag = FullSkip;
                                    }else if(jth_atom->FullyContains(&mid_sphere)){
                                          flag = FullOverlap;
                                          break;
                                    }else{
                                          single_point_index[single_point_count] = overlapped_atom_index[j];
                                          single_point_count ++;
                                    }
                              }
                              if(flag == FullOverlap){
                                    overlap_pts = overlap_pts + ith_pts->data->get_num_points();
                                    //cout<<"Overlap  "<<ith_pts->data->get_num_points()<<endl;

                              }else if(flag == FullSkip && single_point_count == 0){
                                    overlap_pts = overlap_pts + 0;
                                    //cout<<"Skip  "<<endl;
                              }else{
                                    OverlapLinkedList<OvlpPoint3D *>* tmp = ith_pts->data->point_cluster;
                                    OverlapAtom * current_atom;

                                    while(tmp!=NULL){
                                          OvlpPoint3D surface_point = AtomArray[i]->scale_point(tmp->data);
                                          //cout<<"Point "<<surface_point<<endl;
                                          OvlpFlag final_flag=Init;
                                          //cout<<"This of 2222222222\n";

                                          for(int j=0;j<single_point_count;j++){
                                                //cout<<"Single  "<<current_atom->get_serial_no()<<endl;
                                                current_atom = AtomArray[single_point_index[j]];
                                                if( current_atom->dist(&surface_point) <= current_atom->get_radius_with_tolerance()){
                                                      final_flag = SingleInclude;
                                                      break;
                                                }
                                          }
                                          if(final_flag == SingleInclude){
                                                overlap_pts = overlap_pts + 1;
                                          }
                                          tmp = tmp->next;
                                    }
                              }
                              ith_pts = ith_pts->next;
                        }
                  }

                  //cout<<"Working total overlap points for a single atlom = "<<overlap_pts<<endl;

            }
            //cout<<"end of single base pair\n";
            //cout<<"for tresidue no :"<<residue_no<<" total overlap points: "<<overlap_pts<<", total surf_pts= "<<surf_pts<<endl;
            this->total_overlap_points  = overlap_pts;
            this->total_points          = surf_pts;
            this->actual_surface_points = surf_pts - overlap_pts;
      }
      int get_num_atom(){
            return this->num_atom;
      }
};









class OvlpPolymerChain {
private:
      char chain_name;
      int  num_residues;
      int  num_atoms;
      OverlapResidueClass ** ResidueArray;
public:
      ~OvlpPolymerChain(){
            for(int i=0; i<num_residues; ++i){
                  delete ResidueArray[i];
            }
            delete[] ResidueArray;
      }
      OvlpPolymerChain(char chain_name, OverlapStack<OverlapResidueClass *> *ResidueStack){
            this->chain_name    = chain_name;
            this->num_residues  = ResidueStack->get_no_elements();
            this->num_atoms     = 0;
            ResidueArray = new OverlapResidueClass *[num_residues];
            for(int i=num_residues-1; i>=0; i--){
                  ResidueArray[i]  = ResidueStack->pop();
                  this->num_atoms += ResidueArray[i]->get_num_atom();
            }
            if(chain_name == 'Y')
                  this->print();
      }
      int get_num_residues(){
            return num_residues;
      }
      int get_num_atoms(){
            return num_atoms;
      }
      void calc_all_surface(){
            for(int i=0;i<this->num_residues;i++){
                  ResidueArray[i]->compute_surface_area();
            }
      }
      void calc_all_overlap(OvlpPolymerChain * chain){
            for(int i=0;i<this->num_residues;i++){
                  for(int j=0;j<this->num_residues;j++){
                        this->ResidueArray[i]->compute_overlap(chain->ResidueArray[j]);
                  }
            }
      }
      void print(){
            cout<<"CHAIN NAME "<<this->chain_name<<"  Residue numbers "<<this->num_residues<<endl;
            for(int i=0;i<this->num_residues;i++){
                  ResidueArray[i]->print();
            }
      }
};









class OvlpBasePair {
      OverlapResidueClass * base;
      OverlapResidueClass * opposite_base;
      int other[10];
      int primary_base_count;
      int other_base_count;
      int surface_points;
      std::string pair_name;
      char pair_type;
      int total_points;
      int act_surface_pts;
      int total_overlap_pts;
      OverlapAtom ** atom_array;
      int num_atoms;
public:
      ~OvlpBasePair(){
            delete [] atom_array;
            //cout<<"atom_array deleted from OvlpBasePair"<<endl;
            //cout<<"OvlpBasePair deleted"<<endl;
      }
      void set_surface_values(OvlpBasePair * bp){
            this->total_points = bp->total_points;
            this->act_surface_pts = bp->act_surface_pts;
            this->total_overlap_pts = bp->total_overlap_pts;
      }
      int get_original_base_number(){
            return base->get_residue_number();
      }
      bool has_opposite_base(){
            if(this->opposite_base != NULL) return true;
            return false;
      }
      int get_opposite_base_number(){
            assert(this->opposite_base != NULL);
            return opposite_base->get_residue_number();
      }
      bool is_twin(OvlpBasePair * bp){
            if(this->other_base_count == 0 && bp->other_base_count == 0 ){
                  if(this->base->get_residue_number() == bp->opposite_base->get_residue_number()
                     && this->opposite_base->get_residue_number() == bp->base->get_residue_number()){
                        return true;
                  }else{
                        return false;
                  }
            }else{
                  return false;
            }
      }
      int get_num_atoms(){
            return num_atoms;
      }
      OverlapAtom * get_atom(int index){
            return atom_array[index];
      }
      int get_surface_point(){
            return act_surface_pts;
      }
      void compute_surface_area(){
            //cout<<" Trace 1 "<<endl;
            /*if(opposite_base != NULL &&  base->is_enough_close_to(opposite_base, 10) == false){
                cerr<<"Error... bases are not close enough"<<endl;
                exit(11);
            }*/
            int total;
            int surf;
            int ovlp;
            //cout<<" Trace 3 "<<endl;
            OvlpNUtil::compute_surface(atom_array, num_atoms, &total, &ovlp,&surf);
            //cout<<" Trace 4 "<<endl;
            total_points        = total;
            act_surface_pts     = surf;
            total_overlap_pts   = ovlp;
      }
/*    void fprint(FILE* fp){
        char dummy[] = "          ]   ";
        fprintf(fp,"[%6d %c:",base->get_residue_number(),base->get_nucleobase());
        if(opposite_base != NULL){
            fprintf(fp,"%c  %6d ]   ",opposite_base->get_nucleobase(),opposite_base->get_residue_number());
        }else{
            fprintf(fp,"%s",dummy);
        }
    }*/

      void fprint(FILE* fp){
            char dummy[] = "                       ]   ";
            fprintf(fp,"[%6d %3s:",base->get_residue_number(),base->get_nucleobase().c_str());
            if(opposite_base != NULL){
                  fprintf(fp,"%3s  %6d   %6s %c ]   ",opposite_base->get_nucleobase().c_str(),opposite_base->get_residue_number(),this->pair_name.c_str(), this->pair_type);
            }else{
                  fprintf(fp,"%s",dummy);
            }
      }
      OvlpBasePair(OverlapResidueClass * base, OverlapResidueClass * opposite_base,std::string pair_name,char pair_type){
            this->base          = base;
            this->opposite_base = opposite_base;
            this->pair_name     = pair_name;
            this->pair_type     = pair_type;
            this->other_base_count = 0;
            this->total_points = -1; // to check whether the value has been computed or not.
            int num = base->get_num_atom();
            if(opposite_base != NULL){
                  num = num + opposite_base->get_num_atom();
            }
            atom_array = new OverlapAtom *[num];
            this->num_atoms = num;
            int count = 0;
            for(int i=0;i<base->get_num_atom();i++){
                  atom_array[count] = base->get_atom(i);
                  count ++;
            }
            if(opposite_base == NULL){
                  this->primary_base_count = 0;
            }else{
                  for(int i=0;i<opposite_base->get_num_atom();i++){
                        atom_array[count] = opposite_base->get_atom(i);
                        count ++;
                  }
                  this->primary_base_count = 1;
            }
      }

};

class OvlpBasePairArray {
      OvlpBasePair ** base_pair_array;


      int max_residue_number;
      int* pair_index;
      double* pair_result;
      //int* overlap_index;
      //double* bp_bp_overlap_result;
      int num_pairs;
public:
      ~OvlpBasePairArray(){
            delete[] pair_index;
            delete[] pair_result;
      }
      void write_overlap_to_file(string filename){
            string tmp = OvlpUtil::get_base_file_name(filename);
            string bpofile = tmp + ".bpo";
            FILE* fp = fopen(bpofile.c_str(),"w");
            for(int i=0;i<num_pairs-1;i++){
                  int bp = base_pair_array[i]->get_original_base_number();
                  fprintf(fp,"BP ");
                  base_pair_array[i]->fprint(fp);
                  base_pair_array[i+1]->fprint(fp);
                  fprintf(fp,"  : %6.2lf\n",pair_result[bp]);
            }
            fclose(fp);
      }
      OvlpBasePairArray(OverlapStack<OvlpBasePair *>* stack,int max_residue_number){
            this->num_pairs = stack->get_no_elements();
            this->max_residue_number = (max_residue_number>num_pairs)?max_residue_number+1:num_pairs+1;
            pair_index  = new int[this->max_residue_number];
            pair_result = new double[this->max_residue_number];
            //bp_bp_overlap_result = new double[this->max_residue_number];
            for(int i=0;i<this->max_residue_number;i++){
                  pair_index[i] = -1;
                  pair_result[i] = -999.0;
                  //bp_bp_overlap_result[i] = -999.0;
            }
            //pair_result = new double[num_pairs];
            //overlap_result = new double[num_pairs];
            base_pair_array = new OvlpBasePair *[num_pairs];
            for(int i=num_pairs-1;i>=0;i--){
                  base_pair_array[i] = stack->pop();
                  pair_index[base_pair_array[i]->get_original_base_number()] = i;
            }
      }
      void copy_atoms(OverlapAtom ** array,int* num_atoms,int indx1, int indx2){
            int n = base_pair_array[indx1]->get_num_atoms();
            *num_atoms = n;
            int count = 0;
            for(int i=0;i<n;i++){
                  array[count] = base_pair_array[indx1]->get_atom(i);
                  count ++;
            }
            n = base_pair_array[indx2]->get_num_atoms();
            *num_atoms = *num_atoms + n;
            for(int i=0;i<n;i++){
                  array[count] = base_pair_array[indx2]->get_atom(i);
                  count ++;
            }
      }
      void compute_all_bp_surface(){
            //#pragma omp parallel for
            for(int i=0;i<num_pairs;i++){
                  //cout<<"Suspect 2 "<<endl;
                  int flag = 0;
                  //cout<<"BASE ("<<base_pair_array[i]->get_original_base_number()<<","<<base_pair_array[i]->get_opposite_base_number()<<")  :  ("<<base_pair_array[indx]->get_original_base_number()<<","<<base_pair_array[indx]->get_opposite_base_number()<<")"<<endl;
                  if(base_pair_array[i]->has_opposite_base()){
                        if(base_pair_array[i]->get_original_base_number()>base_pair_array[i]->get_opposite_base_number()){
                              int oppo = base_pair_array[i]->get_opposite_base_number();
                              int indx = pair_index[oppo];
                              if(base_pair_array[i]->is_twin(base_pair_array[indx])==true){
                                    //cout<<"BASEBASE    ("<<base_pair_array[i]->get_original_base_number()<<","<<base_pair_array[i]->get_opposite_base_number()<<")  :  ("<<base_pair_array[indx]->get_original_base_number()<<","<<base_pair_array[indx]->get_opposite_base_number()<<")"<<endl;
                                    base_pair_array[i]->set_surface_values(base_pair_array[indx]);
                                    flag = 1;
                              }
                        }
                  }
                  if(flag == 0){
                        base_pair_array[i]->compute_surface_area();
                  }
            }
      }

      void compute_all_bp_bp_overlap(){
            //cout<<"Starts here only\n";

            //cout<<"Anum   "<<num_pairs<<endl;
            //#pragma omp parallel for
            for(int i=0;i<num_pairs-1;i++){  // num_pairs-1 because last but one with last.
                  OverlapAtom * array[200];
                  int num = 0;
                  int total;
                  int surf;
                  int ovlp;
                  //cout<<"Starts here only\n";
                  double overlap = -131.0;
                  int flag = 0;
                  //cout<<"\n\nShound print "<<base_pair_array[i]->get_original_base_number()<<" , "<<endl;
                  if(base_pair_array[i]->has_opposite_base() == true && base_pair_array[i+1]->has_opposite_base() == true){
                        int oppo = base_pair_array[i+1]->get_opposite_base_number();
                        if(oppo < base_pair_array[i]->get_original_base_number()){
                              int indx = pair_index[oppo];
                              if(base_pair_array[indx]->is_twin(base_pair_array[i+1])==true && base_pair_array[indx+1]->is_twin(base_pair_array[i])==true  ){
                                    overlap = this->pair_result[base_pair_array[indx]->get_original_base_number()];
                                    this->pair_result[base_pair_array[i]->get_original_base_number()] = overlap;
                                    //cout<<"From this side??"<<base_pair_array[i]->get_original_base_number()<<endl;
                                    flag = 1;
                              }
                        }

                  }
                  if(flag == 0){
                        copy_atoms(array, &num, i,i+1);
                        //if(base_pair_array[i+1]->is_identical(base_pair_array[i]) == true)
                        OvlpNUtil::compute_surface(array,num, &total,&ovlp, &surf);
                        overlap = (double)(base_pair_array[i]->get_surface_point()+base_pair_array[i+1]->get_surface_point() - surf)/(2.0* OvlpParameters::dots_per_angstrom_sqr);
                        this->pair_result[base_pair_array[i]->get_original_base_number()] = overlap;
                        //cout<<"from normal side"<<endl;

                  }


                  //bp_bp_overlap_result[base_pair_array[i]->get_original_base_number()] = overlap;
                  /*cout<<"BP ";
                  base_pair_array[i]->print();
                  base_pair_array[i+1]->print();
                  cout <<" : " <<overlap<<endl;*/
            }
      }
};





class OvlpRNA_All_Residues {
      int num_residues;
      OverlapResidueClass ** residue_bases;
      //OverlapLinkedList<double> ** ovlp_list;
      //gsl_spmatrix* adj;

      OverlapStack<OvlpFileData *>* file_stack;
      OvlpOutFileRowArray * mOutArray;

public:
      ~OvlpRNA_All_Residues(){

            for(int i=0; i<num_residues; ++i){
                  delete residue_bases[i];
            }
            delete [] residue_bases;
            delete file_stack;
            //cout<<"File stack deleted from OvlpRNA_All_Residues"<<endl;
            delete mOutArray;

            //cout<<"OvlpRNA_All_Residues deleted"<<endl;
      }
      OvlpRNA_All_Residues(OverlapStack<OverlapResidueClass *>* stack,std::string out_file){
            mOutArray = populate_array_from_out_file(out_file);
            this->num_residues  = stack->get_no_elements();
            residue_bases       = new OverlapResidueClass *[this->num_residues];
            //cout<<"This is suspect..."<<endl;
            file_stack = new OverlapStack<OvlpFileData *>(20*this->num_residues); // Assuming max 20 pairs possible for each residue.

            //adj = gsl_spmatrix_alloc(this->num_residues, this->num_residues);
            //gsl_spmatrix_set_zero(adj);
            for(int i=this->num_residues-1; i>=0;i--){
                  //overlap[i] = -999.0;

                  residue_bases[i] = stack->pop();
                  //cout<<i<<endl;
                  //if(i==num_residues-1 || i == 0){
                  //    residue_bases[i]->print();
                  //;cout<<<endl;
                  //}
            }
      }
      void calc_all_surface(){
            for(int i=0;i<this->num_residues;i++){
                  residue_bases[i]->compute_surface_area();
            }
      }
      void calc_all_overlap(){

            double result= 0.0;

            //cout<<"NUM RESIDUE "<<this->num_residues<<endl;
            for(int i=0;i<this->num_residues;i++){
                  //cout<<i<<" of "<<num_residues<<" Started"<<endl;
//cout<<"A2 i="<<i<<endl;//", j="<<j<<endl;
                  for(int j=i+1;j<this->num_residues;j++){
                        result = residue_bases[i]->compute_overlap(residue_bases[j]);
                        //cout<<"RESULT "<<result<<endl;


                        if(result > OvlpParameters::ovlp_cutoff){
                              std::string base1_atom_names[30];
                              std::string base2_atom_names[30];
                              int num_atom_base1 = 0;
                              int num_atom_base2 = 0;
                              bool check = residue_bases[i]->get_atom_overlap_mapping_atoms(residue_bases[j], base1_atom_names, &num_atom_base1, base2_atom_names, &num_atom_base2);
                              if(check == false){
                                    cout<<"No overlapping in atoms... please contact developers..."<<endl;
                                    exit(1);
                              }
                              //cout<<"Working... "<<mOutArray->get_num_rows()<<" "<<endl;
                              OvlpOutFileRow * row_i        = mOutArray->get_outrow(i);
                              OvlpOutFileRow * row_j        = mOutArray->get_outrow(j);

                              std::string chain_name_i = row_i->get_chain_name();
                              std::string chain_name_j = row_j->get_chain_name();
                              int row_i_cor_serial     = row_i->get_cor_serial();
                              int row_j_cor_serial     = row_j->get_cor_serial();
                              std::string base_name    = row_i->get_pair_name(row_i_cor_serial, row_j_cor_serial);
                              //cout<<"A2\n";

                              if(base_name == "NOT_FOUND"){
                                    std::string forward_b_name = row_i->get_pair_type(row_i_cor_serial, row_j_cor_serial-1);
                                    std::string previous_b_name = row_i->get_pair_type(row_i_cor_serial, row_j_cor_serial+1);
                                    //cout<<"("<<row_i_cor_serial<<","<<row_j_cor_serial<<") B_NAME"<<previous_b_name<<endl;
                                    if( forward_b_name == "BP" || previous_b_name == "BP"){
                                          base_name = "BCRS";
                                    }else if( forward_b_name == "TP" || previous_b_name == "TP"){
                                          base_name = "TCRS";
                                    }else if(row_i->get_chain_name() == row_j->get_chain_name()){
                                          int tmpx = row_i->get_cor_serial();
                                          int tmpy = row_j->get_cor_serial();

                                          int val  = (tmpx>tmpy)?tmpx-tmpy:tmpy-tmpx;  // find the absolute difference.
                                          if(val == 1){
                                                base_name = "PSTK"; // probably stacking;
                                          }else{
                                                base_name = "----";
                                          }
                                    }else{
                                          base_name = "----";
                                    }
                              }

                              std::string tmp_base_type = row_i->get_pair_type(row_i_cor_serial, row_j_cor_serial);
                              if(tmp_base_type == "NOT_FOUND"){
                                    tmp_base_type = "--";
                              }

                              OvlpFileData * fdata = new OvlpFileData( this->residue_bases[i]->get_residue_number(),
                                                              row_i->get_pdb_serial(),
                                                              this->residue_bases[i]->get_nucleobase().c_str(),
                                                              chain_name_i,
                                                              this->residue_bases[j]->get_residue_number(),
                                                              row_j->get_pdb_serial(),
                                                              this->residue_bases[j]->get_nucleobase(),
                                                              chain_name_j,
                                                              base_name,
                                                              tmp_base_type,
                                                              result
                              );
                              file_stack->push(fdata);
                              //gsl_spmatrix_set(adj, i, j, result);
                              //gsl_spmatrix_set(adj, j, i, result);
                        }

                  }
                  //cout<<i<<" of "<<num_residues<<" completed"<<endl;
            }

      }

      OverlapResidueClass * heuristic_bin_search(int res_no, int start,int stop){
            if(start > stop) return NULL;
            int heu_mid = start + (stop - start)/2;
            cout<<"res "<<res_no<<" , start "<<start<<" , stop ="<<stop<<endl;
            this->residue_bases[heu_mid]->print();
            int data = this->residue_bases[heu_mid]->get_residue_number();
            if( data == res_no){
                  return this->residue_bases[heu_mid];
            }else if(data < res_no){
                  return heuristic_bin_search(res_no, heu_mid+1, stop);
            }else{
                  return heuristic_bin_search(res_no, start, heu_mid-1);
            }
      }

      OverlapResidueClass * get_residue(int residue_number){
            int no = this->residue_bases[residue_number-1]->get_residue_number();
            if( no == residue_number){
                  return this->residue_bases[residue_number-1];
            }else if(no > residue_number){
                  OverlapResidueClass * tmp = heuristic_bin_search(residue_number, 0, residue_number-2);
                  if(tmp == NULL){
                        cout<<"No valid search result found... "<<no<<" , target "<<residue_number<<endl;
                        exit(1);
                  }
                  assert(tmp != NULL);

                  return tmp;
            }else{
                  OverlapResidueClass * tmp = heuristic_bin_search(residue_number, residue_number, num_residues-1);
                  if(tmp == NULL){
                        cout<<"No valid search result found... "<<no<<" , target "<<residue_number<<endl;
                        exit(1);
                  }
                  //assert(tmp != NULL);
                  return tmp;
            }
      }

      void write_overlap_to_file(FILE * file){
            //gsl_spmatrix_fprintf(file, adj, "%4.5lf");

            fprintf(file,"HEADER    Generated by overlap function\n");
            fprintf(file,"ACCN      %s\n", OvlpParameters::pdb_accn.c_str());
            fprintf(file,"VERSION   %s (Overlap Version)\n", OvlpParameters::file_version.c_str());
            fprintf(file,"TIME      %s\n", OvlpGen::today().c_str());
            fprintf(file,"DVLPR     Parthajit Roy, roy.paryhajit@gmail.com\n");
            fprintf(file,"DVLPR     Dhananjay Bhattacharyya, dhananjay.bhattacharyya@saha.ac.in\n");
            fprintf(file,"COLINFO   (09-14 : 16-21)    .cor file cor serial\n");
            fprintf(file,"COLINFO   (24-29 : 31-36)    .pdb file serial\n");
            fprintf(file,"COLINFO   (40-42 : 44-46)    Nucleobase name:Nucleobase name\n");
            fprintf(file,"COLINFO   (52 : 54)          Chain name:Chain name\n");
            fprintf(file,"COLINFO   (59-62)            Base pair name.\n");
            fprintf(file,"COLINFO   (65-66)            Base pair type. if -- then no type detected.\n");
            fprintf(file,"COLINFO   (74-78)            Overlap Value.\n");
            fprintf(file,"ABVR      PSTK means probably Stacking.\n");
            fprintf(file,"ABVR      BCRS means Probably Cross overlap between BP and the opposite Diagonal.\n");
            fprintf(file,"ABVR      TCRS means Probably Cross overlap between TP and the opposite Diagonal.\n");
            fprintf(file,"NUMROW    %6d\n",file_stack->get_no_elements());
            OverlapStack<OvlpFileData *> ReverseStack = OverlapStack<OvlpFileData *>(this->file_stack->get_no_elements());
            while(this->file_stack->isEmpty() != true){
                  ReverseStack.push(this->file_stack->pop());
            }
            while(ReverseStack.isEmpty() != true){
                  OvlpFileData * fd = ReverseStack.pop();
                  fd->fprint(file);
            }
            //gsl_spmatrix_fwrite(FILE *stream, const gsl_spmatrix *m)
      }
};










class OvlpRNA_NucVatiants {
      std::string rnafile_name;
      std::string rna_fbase;
      std::string ext;
      ifstream    rna_file;

      OvlpStringArray * mAdenineName;
      OvlpStringArray * mGuanineName;
      OvlpStringArray * mCytosineName;
      OvlpStringArray * mUracilName;

      OvlpStringArray * mAdenineAtoms;
      OvlpStringArray * mGuanineAtoms;
      OvlpStringArray * mCytosineAtoms;
      OvlpStringArray * mUracilAtoms;




      void setup_init_values(int* curr_res_no, char curr_residue_name[]){
            //cout<<"Test11234\n";
            std::string line;
            long int start = 0;
            while(true){
                  std::getline(rna_file, line);
                  //cout<<"FIRST SEM "<<line<<endl;
                  if(line.substr(0,4) == "ATOM"){
                        long int curr_len = line.length();
                        start = rna_file.tellg() - curr_len - 1L;
                        rna_file.seekg(start,ios_base::beg);
                        //std::getline(rna_file,line);
                        //cout<<"TESTING  "<<line<<endl;
                        break;
                  }
            }
            rna_file.seekg(18, ios_base::cur); // Go to the chain position in the text file.
            rna_file>>curr_residue_name;
            rna_file>>*curr_res_no;



            rna_file.seekg(start, ios_base::beg);
            //printf("FIRST SEM %s  %d\n",curr_residue_name,*curr_res_no);
      }


/*    void setup_init_values(char* curr_chain_name,int* curr_res_no, char curr_residue_name[]){
        rna_file.seekg(19, ios_base::beg); // Go to the chain position in the text file.
        rna_file>>curr_residue_name;
        rna_file>>*curr_chain_name;
        rna_file>>*curr_res_no;


        rna_file.seekg(0, ios_base::beg);
        printf("%s %s  %d\n",curr_residue_name,curr_chain_name,*curr_res_no);
    }*/
public:
      bool is_valid_atom_for_base(std::string line){
            stringstream ss;
            ss<<line;
            std::string tag;
            int id;
            std::string atom;
            std::string base;
            ss>>tag;
            ss>>id;
            ss>>atom;
            ss>>base;

            if(mAdenineName->search(base) != -1){
                  if(mAdenineAtoms->search(atom) != -1){
                        return true;
                  }else{
                        return false;
                  }
            }
            if(mGuanineName->search(base) != -1){
                  if(mGuanineAtoms->search(atom) != -1){
                        return true;
                  }else{
                        return false;
                  }
            }
            if(mCytosineName->search(base) != -1){
                  if(mCytosineAtoms->search(atom) != -1){
                        return true;
                  }else{
                        return false;
                  }
            }
            if(mUracilName->search(base) != -1){
                  if(mUracilAtoms->search(atom) != -1){
                        return true;
                  }else{
                        return false;
                  }
            }
	    if(base == "PSU") return true;
            cerr<<"Error.... "<<base<<" is not a valid basename in syscon .name files... Please contact developpers\n";
            exit(1);
            return false;
      }
      ~OvlpRNA_NucVatiants(){
            delete mAdenineName;
            delete mGuanineName;
            delete mCytosineName;
            delete mUracilName;

            delete mAdenineAtoms;
            delete mGuanineAtoms;
            delete mCytosineAtoms;
            delete mUracilAtoms;
      }
      void set_rna(std::string rnafile){
            this->rnafile_name  = rnafile;
            this->rna_fbase     = rnafile.substr(0, rnafile.find_last_of("."));

            this->ext = rnafile.substr(rnafile.find_last_of(".") + 1, rnafile.length());
            if(strncasecmp(ext.c_str(), "cor",6)) {
                  cerr<<"Error    ."<<ext<<" file has been supplied insted of .cor file"<<endl;
                  exit(1);
            }
            //cout<<"corfile:"<<this->rnafile_name<<endl;
      }
      OvlpRNA_NucVatiants(){
            char* nucdir = getenv("NUCLEIC_ACID_DIR");
		if(nucdir == NULL){
			fprintf(stderr, "Error... NUCLEIC_ACID_DIR not Defined.\n");
			exit(EXIT_FAILURE);
		}

		string nucfiledir;
		nucfiledir=nucdir;
            mAdenineName = OvlpNUtil::read_words_from_file(nucfiledir+"/AdeVariants.name");
            mGuanineName = OvlpNUtil::read_words_from_file(nucfiledir+"/GuaVariants.name");
            mCytosineName = OvlpNUtil::read_words_from_file(nucfiledir+"/CytVariants.name");
            mUracilName = OvlpNUtil::read_words_from_file(nucfiledir+"/UraVariants.name");


            mAdenineAtoms = OvlpNUtil::read_words_from_file(nucfiledir+"/adenine_atoms.name");
            mGuanineAtoms = OvlpNUtil::read_words_from_file(nucfiledir+"/guanine_atoms.name");
            mCytosineAtoms = OvlpNUtil::read_words_from_file(nucfiledir+"/cytosine_atoms.name");
            mUracilAtoms = OvlpNUtil::read_words_from_file(nucfiledir+"/uracil_atoms.name");
            //cout<<"HERE..."<<endl;


            //cout << rna_fbase<<endl<<ext<<endl;
      }


      OvlpRNA_All_Residues * get_rna_data_structures(OvlpAllSurfacePoints * all_surf,std::string out_filename){
            //all_surf->display();
            rna_file.open(this->rnafile_name.c_str(),std::ios_base::in);
            if(rna_file.is_open() == false){
                  cerr<<"Error... unable to open file named "<<rnafile_name<<endl;
                  exit(1);
            }
            std::string line;
            int curr_residue_no;
            //char curr_chain_name;
            char curr_residue_name[5];


            char tag[10];
            int atom_serial;
            char atom_type[10];
            char residue_name[5];
            int residue_no;
            double x;
            double y;
            double z;
            double occupancy;

            std::string ter = "END";


            // This type of OverlapStack declarations are strictly against parallel programming. This is done for speeding up the sequential processing.
            OverlapStack<OvlpPolymerChain *> StackOfPolymerChain = OverlapStack<OvlpPolymerChain *>(1000);
            OverlapStack<OverlapResidueClass *> StackOfResidue      = OverlapStack<OverlapResidueClass *>(500000);
            OverlapStack<OverlapAtom *>         StackOfAtom         = OverlapStack<OverlapAtom *>(200);

            setup_init_values(&curr_residue_no,curr_residue_name);  // Values have been setup through function call is just because to avoid clumsyness in code and for better readability.
            //residue_no = curr_residue_no;

            //assert(StackOfPolymerChain.isEmpty());
            assert(StackOfResidue.isEmpty());
            //int count = 0;
            std::getline(rna_file,line);
	    line[16] = ' ';
            sscanf(line.c_str(),"%s%d%s%s%d%lf%lf%lf%lf",tag,&atom_serial,atom_type,
                   residue_name,&residue_no,&x,&y,&z,&occupancy);
            OvlpAtomType atomtype = get_atom_base_type(atom_type);
            while(line != ter){
                  assert(StackOfAtom.isEmpty());
                  //bool flag = false;
                  while(curr_residue_no == residue_no){
                        //cout<<"CRUU RES"<<curr_residue_no<<", "<<residue_no<<endl;
                        if(is_valid_atom_for_base(line) == true){
                              //get_atom_base_type()

                              if(atomtype != Unknown_Atom){
                                    //cout<<"LINE "<<line<<endl;
                                    OvlpAtomSurfPoints * surf_pts = all_surf->get_surface_points(atomtype);
                                    OverlapAtom * atom = new OverlapAtom(atom_serial,get_atom_type(atom_type), atom_type,x,y,z, ovlp_get_atom_radius(atom_type[0]),occupancy,-999.0,surf_pts);
                                    StackOfAtom.push(atom);
                                    //if(residue_no)
                                    //flag = true;
                              }
                        }

                        std::getline(rna_file,line);
			line[16] = ' ';
                        sscanf(line.c_str(),"%s%d%s%s%d%lf%lf%lf%lf",tag,&atom_serial,atom_type,
                               residue_name,&residue_no,&x,&y,&z,&occupancy);
                        atomtype = get_atom_base_type(atom_type);

                        if(line == ter){
                              //cout<<"TTTT"<<line<<"BBBB"<<endl;
                              residue_no = -1;
                              break;
                        }
                  }
                  //OverlapResidueClass* residue = new OverlapResidueClass(curr_residue_no,get_nitrogenous_base(curr_residue_name),&StackOfAtom);
                  //cout<<"CURR RES"<<curr_residue_no<<" , "<<curr_residue_name<<", RES NO"<<residue_no<<endl;
                  OverlapResidueClass * residue = new OverlapResidueClass(curr_residue_no,curr_residue_name,&StackOfAtom);


                  StackOfResidue.push(residue);
                  //residue->print();
                  curr_residue_no   = residue_no;

                  strcpy(curr_residue_name , residue_name);

                  //OverlapResidueClass(int residue_number, int num_atom, OvlpNitrogenousBase nucleo_base, OverlapStack<OverlapAtom*> *AtomStack){ // This constructor is non-conventional. This is created because atoms are first placed in a stack when created.
            }
            //cout<<"STK OF RES "<<StackOfResidue.get_no_elements()<<endl;
            OvlpRNA_All_Residues * rna_all_residues = new OvlpRNA_All_Residues(&StackOfResidue,out_filename);
            return rna_all_residues;








/*        while(rna_file.eof() != true){
            assert(StackOfResidue.isEmpty());
            cout<<"Residue no"<<curr_residue_no<<endl;
            int residue_start = curr_residue_no;

            while(curr_chain_name == chain_name){
                assert(StackOfAtom.isEmpty());
                while(curr_residue_no == residue_no){
                    OvlpAtomType atomtype = get_atom_base_type(atom_type[0]);
                    //cout<<"atom type "<<atomtype<<endl;
                    OvlpAtomSurfPoints* surf_pts = all_surf->get_surface_points(atomtype);
                    //surf_pts->test();
                    OverlapAtom* atom = new OverlapAtom(atom_serial,get_atom_type(atom_type),x,y,z,ovlp_get_atom_radius(atom_type[0]),occupancy,temp_factor,surf_pts);
                    //cout<<"atom...  "<<atom_serial<<",  "<<chain_name<<",  "<<residue_no<<endl;
                    StackOfAtom.push(atom);
                    std::getline(rna_file,line);
                    sscanf(line.c_str(),    "%s%d%s%s%c%c%d%lf%lf%lf%lf%lf%s",tag,&atom_serial,atom_type,
                                            residue_name,&pass, &chain_name,&residue_no,&x,&y,&z,&occupancy,
                                            &temp_factor,symbol);
                    if(rna_file.eof() == true){
                        chain_name = '\0';
                        residue_no = -1;
                        //cout<<"ATOM... "<<tag<<", "<<atom_serial<<" "<<residue_no;
                    }
                }
                OverlapResidueClass* residue = new OverlapResidueClass(curr_residue_no,get_nitrogenous_base(curr_residue_name),&StackOfAtom);
                StackOfResidue.push(residue);
                curr_residue_no   = residue_no;

                strcpy(curr_residue_name , residue_name);

                //OverlapResidueClass(int residue_number, int num_atom, OvlpNitrogenousBase nucleo_base, OverlapStack<OverlapAtom*> *AtomStack){ // This constructor is non-conventional. This is created because atoms are first placed in a stack when created.
            }


            OvlpPolymerChain* polymer_chain = new OvlpPolymerChain(curr_chain_name,&StackOfResidue);
            //OvlpPolymerChain(char chain_name[],int num_residue_class,int num_rna_atoms, OverlapStack<OverlapResidueClass*> *ResidueStack){
            StackOfPolymerChain.push(polymer_chain);
            curr_chain_name = chain_name;
            ///cout<<atom_serial<<"  working  "<<residue_no<<endl;
        }
        OvlpRNA* rna = new OvlpRNA(&StackOfPolymerChain);
        cout<<"Upo to this working"<<endl;
        rna_file.close();
        return rna; */
      }
};










class OvlpBasePairNupFile {
      std::string file_name;
      std::string fbase;
      std::string ext;
      ifstream    file;
      void setup_init_values(int* curr_res_no, char curr_residue_name[]){
            file.seekg(19, ios_base::beg); // Go to the chain position in the text file.
            file>>curr_residue_name;
            file>>*curr_res_no;
            file.seekg(0, ios_base::beg);
            //printf("%s  %d\n",curr_residue_name,*curr_res_no);
      }
public:
      OvlpBasePairNupFile(std::string nupfile){
            this->file_name  = nupfile;
            this->fbase      = nupfile.substr(0, nupfile.find_last_of("."));

            this->ext = nupfile.substr(nupfile.find_last_of(".") + 1, nupfile.length());
            if(strncasecmp(ext.c_str(), "nup",6)) {
                  cerr<<"Error    ."<<ext<<" file has been supplied insted of .nup file"<<endl;
                  exit(1);
            }
      }
      OvlpBasePairArray * populate_base_pairs(OvlpRNA_All_Residues * all_bases){
            //all_surf->display();
            file.open(this->file_name.c_str(),std::ios_base::in);
            if(file.is_open() == false){
                  cerr<<"Error... unable to open file "<<file_name<<endl;
                  exit(1);
            }

            int original;
            int primary_pair;
            bool has_pair;
            char pair_name[10];
            char dummy;
            char pair_type;
            int num_secondary_pair;
            int* secondary_pair;
            int len;

            std::string line;
            OverlapStack<OvlpBasePair *> stack = OverlapStack<OvlpBasePair *>(500000);
            //cout<<"trace 1... "<<endl;
            int max_base_no = 0;
            while(file.eof() != true){
                  std::getline(file,line);
                  len = strlen(line.c_str());

                  if(len == 0) break;
                  if(len<=5){
                        sscanf(line.c_str(),    "%d",&original);
                        if(max_base_no < original){
                              max_base_no = original;
                        }
                        has_pair = false;
                        num_secondary_pair = 0;
                        //printf("%d\n",original);
                  }else if(len<=200){
                        sscanf(line.c_str(),    "%d%d%s%c%c",&original,&primary_pair,pair_name,&dummy,&pair_type);
                        if(max_base_no < original){
                              max_base_no = original;
                        }
                        if(max_base_no < primary_pair){
                              max_base_no = primary_pair;
                        }
                        //printf("%d %d %s %c\n",original,primary_pair,pair_name,pair_type);
                        has_pair = true;
                        num_secondary_pair = 0;
                  }else{
                        cout<<"Unknown format in .nup file..."<<endl;
                        exit(1);
                  }
                  OvlpBasePair * bp;
                  if(has_pair == true){
                        bp = new OvlpBasePair(all_bases->get_residue(original), all_bases->get_residue(primary_pair),pair_name,pair_type);
                  }else{
                        bp = new OvlpBasePair(all_bases->get_residue(original), NULL," ",'\0');
                  }
                  stack.push(bp);
            }
            OvlpBasePairArray * base_pair_array = new OvlpBasePairArray(&stack,max_base_no);
            return base_pair_array;
      }
};







void test(){
      /*char x;
      int y;
      char name[10];
      double p,q,r;
      scanf("%c%d%s%lf%lf%lf",&x,&y, name,&p,&q,&r);
      printf("%c --  %d  --  %s --  %lf  --  %lf  --  %lf\n",x,y,name,p,q,r);*/
      /*int sum = 1;
      for(int i=0;i<2000;i++){
          for(int j=0;j<2000;j++){
              sum = sum + 0;
          }
      }*/
      //get_nitrogenous_base("PUANINE");
      int x;
      int y;
      const int* const ptr = &x;
      //*ptr = 5;
      //ptr = &y;

}



void ovlp_base_overlap_comp(string pdb_accn, double ovlp_cutoff,
                            OvlpSurfaceDataFile * surfgen,
                            OvlpAllSurfacePoints * all_surf_points,
                            OvlpRNA_NucVatiants *nucVariants){
      OvlpParameters::ovlp_cutoff = ovlp_cutoff;
      std::string cor_file = pdb_accn+".cor";
      std::string out_file = pdb_accn+".out";


      //OvlpSurfaceDataFile * surfgen = new OvlpSurfaceDataFile("./syscon/surface.xyz",base_type,3);
      //OvlpAllSurfacePoints * all_surf_points = surfgen->generate_surface_points(3);
      //OvlpRNA_NucVatiants *nucVariants = new OvlpRNA_NucVatiants();
      nucVariants->set_rna(cor_file);
      OvlpRNA_All_Residues * rna = nucVariants->get_rna_data_structures(all_surf_points,out_file);

      rna->calc_all_surface();
      rna->calc_all_overlap();
      const string robfile = pdb_accn + ".rob";
      FILE* fp = fopen(robfile.c_str(),"w");
      rna->write_overlap_to_file(fp);
      fclose(fp);
      delete rna;
      //delete all_surf_points;
      //delete surfgen;
}



void ovlp_base_pair_overlap_comp(char **argv, double ovlp_cutoff){
      OvlpParameters::ovlp_cutoff = ovlp_cutoff;
      OvlpAtomType base_type[] = {Carbon,Nitrogen,Oxygen};
      //assert((x==0&&"Working"));

      cout<<"Computation starts..."<<endl;

      OvlpSurfaceDataFile * surfgen = new OvlpSurfaceDataFile("./syscon/surface.xyz",base_type,3);
      OvlpAllSurfacePoints * all_surf_points = surfgen->generate_surface_points(3);
      OvlpRNA_NucVatiants * f1 = new OvlpRNA_NucVatiants();
      f1->set_rna(argv[2]);
      OvlpRNA_All_Residues * rna = f1->get_rna_data_structures(all_surf_points,argv[3]);
      OvlpBasePairNupFile * f2 = new OvlpBasePairNupFile(argv[4]);
      OvlpBasePairArray * all_bp;

      all_bp = f2->populate_base_pairs(rna);
      cout<<"Computation of base pair starts..."<<endl;
      all_bp->compute_all_bp_surface();
      all_bp->compute_all_bp_bp_overlap();
      cout<<"Computation of base pair ends..."<<endl;
      cout<<"Writing to file starts..."<<endl;
      all_bp->write_overlap_to_file(argv[3]);

}





int main1(int argc, char* argv[]){


      OvlpParameters::pdb_accn = argv[2];
      //cout<<"test 1"<<endl;
      //OvlpNUtil::read_words_from_file("./syscon/adenine_atoms.name");

      //all_surf_points->display();

      //OvlpRNADataFile* f  = new OvlpRNADataFile(argv[1]);
      //OvlpRNA* rna        = f->get_rna_data_structures(all_surf_points);


      //cout<<"Starts working..."<<endl;
      if(argc == 1){
            cout<<"Please issue proper command line arguments or run with -h for help."<<endl;
            //populate_array_from_out_file("1asy.out");
            exit(1);
      }

      if(strncasecmp(argv[1], "-ab" , 10) == 0){
            std::string pdb_accn=argv[2];
	    cout<<"test1"<<endl;
            OvlpAtomType base_type[] = {Carbon,Nitrogen,Oxygen};
            OvlpSurfaceDataFile * surfgen = new OvlpSurfaceDataFile("./syscon/surface.xyz",base_type,3);
	    cout<<"test2"<<endl;
            OvlpAllSurfacePoints * all_surf_points = surfgen->generate_surface_points(3);
	    cout<<"test3"<<endl;
            OvlpRNA_NucVatiants *nucVariants = new OvlpRNA_NucVatiants();
	    cout<<"test4"<<endl;
            ovlp_base_overlap_comp(pdb_accn, 25.0, surfgen, all_surf_points, nucVariants);
	    cout<<"test5"<<endl;
            delete all_surf_points;
            delete surfgen;
      }/*else if(strncasecmp(argv[1], "-absummary" , 15){
        std::string pdb_accn=argv[2];
        base_overlap_comp(pdb_accn,2);
    }else if(strncasecmp(argv[1], "-abdetail" , 15){
        std::string pdb_accn=argv[2];
        ovlp_base_overlap_comp(pdb_accn,3);
    }*/
      else if(strncasecmp(argv[1], "-bp" , 10) == 0){
            ovlp_base_pair_overlap_comp(argv, 20.0);
      }else if(strncasecmp(argv[1], "--version" , 10) == 0){
            cout<<"Overlap version "<< OvlpParameters::file_version<<"\n";
      }else if(strncasecmp(argv[1], "-h" , 10) == 0){
            cout<<endl<<"-------------------compiling overlap file---------------------"<<endl<<endl;
            cout<<"g++ overlap.cpp -fopenmp -lgsl -lcblas -lm"<<endl;
            cout<<"NOTE: openmp, gsl needs to be installed in your ubuntu machine"<<endl<<endl;
            cout<<"-------------------running base overlap----------------------"<<endl<<endl;
            cout<<"If you want to compute all base overlaps then do the following"<<endl;
            cout<<"./a.out -ab  <abc.base> "<<endl;
            cout<<"The output file will be stored in .rob with the same accession number"<<endl<<endl;
            cout<<"-------------------running base-pair overlap------------------"<<endl<<endl;
            cout<<"If you want to compute base pair overlap then do the following"<<endl;
            cout<<"./a.out -bp  <abc.base>  <abc.nup>"<<endl;
            cout<<"The output will be stored in .bpo file with same accession number."<<endl<<endl;
            cout<<"-------------------generation of .base file---------------------"<<endl<<endl;
            cout<<"for generation of .base file do the following"<<endl;
            cout<<"gawk -f genbase.awk <abc.cor>"<<endl;
            cout<<"For cor file generation, run bpfind program... see Saha Institute following detail"<<endl;
            cout<< "http://www.saha.ac.in/biop/bioinformatics.html"<<endl<<endl;
            cout<<"---------------------Developpers---------------------------"<<endl<<endl;
            cout<<"Developper: Parthajit Roy and Dr. D. Bhattacharyya"<<endl<<endl;
            cout<<"---------------------Reporting Bugs------------------------"<<endl<<endl;
            cout<<"P. Roy,              e-mail: roy.parthajit@gmail.com or"<<endl;
            cout<<"D. Bhattacharyya     e-mail: dhananjay.bhattacharyya@saha.ac.in"<<endl<<endl;

      }else{
            cout<<"Invalid switch "<<argv[1]<<", Please run the program with -h option."<<endl;
            exit(0);
      }




      return(0);
}

#endif//OVERLAPLIB_OVERLAP_C
