/*
 * geometry.h
 *
 *  Created on: Apr 25, 2019
 *      Author: farri
 */

#ifndef geometry_H_
#define geometry_H_
#include "tri_tri_intersect.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <vector>
using namespace std;

struct HE
{
	int vin;
	int vout;
	int id;
	int type;
	//struct HE *next;
	//struct HE *prev;
	//struct HE *op;

	int nextid;
	int previd;
	int opid;
	double hevec[3];
	double hecent[3];
	double n[3];
	double l;
	bool dout;
	bool din;
	//double nextangle;
};

struct VTX
{
	int vid;
	double co[3];
	vector<int> hein;
	int fusion_vid;
	//int hesurfoutid;
	int hesurfinid;
	vector<int> vneigh;
};
typedef struct HE HE;
typedef struct VTX VTX;

class geometry
{
public:
	bool Test_assembly;
	int Nvlast;
	int Nhelast;
	int Ntype;
	int Nv;
	int Nv5;
	int Nhe;
	int Nsurf;
	int Nd;

	int all_neigh;
	//double gb;
	double dg;
	double mudimer;
	double *mu;
	double mudrug;
	
	//double muAB;
	double drugProb;
	double T;

	double xi;

	double *epsilon;
	double *kappa;
	double *phi0;
	double *kappaPhi;
	double *theta0;
	double *l0;
	double **gb;
	double **gdrug;
	vector<VTX> v; /**< Current positions of the vertices. */
	vector<HE> he; /**< Pairs i,j of vertex indices that are connected. */
	//vector<VTX> tricent; // center of triangles
	int *vidtoindex;
	int *heidtoindex;
	vector<int> surfheid;
	vector<int> surfv;
	vector<int> surfvbond;
	vector<int> fusionv;

	geometry();
	~geometry();

	void initialize(int Ntyp0);

	void update_index();

	void update_neigh();

	void update_neigh_vertex(int vid0);

	void update_surface();

	int is_surface(int heid0);

	int is_bond_in_surface(int heid0);

	int is_bond_out_surface(int heid0);

	int no_bond_surface(int heid0);

	int is_vsurface(int vid);

	int is_bond_vsurface(int vid);

	int open_wedge(int heid0, int *flag);

	int pre_open_wedge(int heid0);

	int next_open_wedge(int heid0);

	int connected(int vid0, int vid1);

	int connectedH(int heid0, int heid1);

	int next_connected(int vid0, int vid1);

	int not_cross_edge(int heid0, int heid1);

	void add_vertex(double *xyz);

	int add_edge_type(int vin, int vout, int etype);

	void add_half_edge_type(int vin0, int vout0, int etype);

	int check_overlap_centerh(double *tempcenter);

	int check_overlap_centerv(double *tempcenter);

	int do_intersect(int heid1, int heid2);

	//int shared_vertices(int heid1, int faceid2);

	int remove_neigh(int vid0, int removevid);
	
	int check_overlap_g(int vid0);

	int find_overlap_g(int vid0);

	int check_overlap_he(int heid0);

	/*int check_overlap_hesurf(int heid0);

	int check_overlap_vsurf(int vid0);

	int check_overlap_vtx(int vid0);*/

	int is_same_triangle(int heid0, int heid1);

	//int get_vindex(int vid);

	//int get_heindex(int heid0);

	void make_triangle();

	int opposite_edge(int heid0);

	void set_prev_next(int heid0, int previd0, int nextid);

	int get_prev_next(int heid0, int *previd0, int *nextid0);

	int add_dimer(int heid0, gsl_rng *r, int typenext, int typeprev);

	int force_add_dimer(int heid0, double *newv, int tyepnext, int typeprev);

	int remove_dimer(int heindex0, int heindexhext0);

	void make_hexamer();

	void he_initialize(int heindex, int heid0, int vin0, int vout0, int etype);

	int add_monomer(int nextofnewid, int prevofnewid, int etype);

	double find_gbb(int etypenew, int etypenextofnew, int etypeprevifnew);

	double find_dg(int type, int typenext, bool drug);

	int force_add_monomer(int nextofnewid, int prevofnewid, int etype);

	int add_monomer_dimer(int heid0);

	int remove_monomer_dimer(int heid0, gsl_rng *r);

	void new_vertex(int heindex0, double *newv);

	void move_p(double *pi, double *pf, gsl_rng *r);

	void move_v(double *pi, double *pf, gsl_rng *r);

	void hecenter(int heindex, double *vcenter);

	void helen(int heindex, double *helen);

	int delete_vertex(int vid0);

	int delete_edge(int heid0);

	int get_normal(int heid0);

	void update_normals();

	void update_edge(int heid0);

	void update_half_edge(int heid0);

	double stretch_energy(int heindex0);

	int check_bind_wedge(int heid0);

	double bend_energy(int heindex0);

	double dimer_bend_energy(int heindex0);

	double dimer_energy(int heid0, int heid1);

	double monomer_energy(int heid0);

	double compute_bind_energy();

	double compute_energy();

	double vertex_energy(int vid0);

	void make_pentamer();

	void dump_parameters();

	int get_fusion_vid(int vid0);

	void update_fusion_pairs();

	void save_vtx(int vid0 , VTX *tempvtx);

	void check_odd_neigh();

};	
void subvec(double *vinit, double *vfin, double *vec);

void addvec(double *vinit, double *vfin, double *vec);

void multvec(double *vinit, double scalar, double *vec);

void centvec(double *vinit, double *vfin, double *vec);

double veclen(double *vinit, double *vfin);

double norm(double *v);

double dot(double *v1, double *v2);

void cross(double *v1, double *v2, double *res);

void randvec(double *v, gsl_rng *r);

void dump_lammps_data_file(geometry &g, int time0);

void move_vertex(geometry &g, gsl_rng *r);

void rotatevec(double *vec, double *axis, double angle, double *vec2);

void read_lammps_data(geometry &g, char filename[]);

void dump_data_frame(geometry &g, FILE *f, int time);

void recenter(geometry &g);

int surfclosev(geometry &g);

//int valid(geometry &g);

//int add_dimer(geometry &g,int hei);

//void add_edge_type(geometry &g,VTX *vin0,VTX *vout0, int etype);

//void make_hexamer(geometry &g);

//void add_vertex(geometry &g, double *xyz);

//void make_triangle(geometry &g);
//int get_vindex(geometry &g, int vid0 ) ;
#endif /* geometry_H_ */
