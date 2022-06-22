/*
 * geometry.cpp
 *
 *  Created on: Apr 25, 2019
 *      Author: farr6
 */

#include "geometry.hpp"
#include "tri_tri_intersect.hpp"
#include <iostream>
#include <vector>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <queue>
#define PI 3.14159265
using namespace std;

geometry::geometry()
{
	// TODO Auto-generated constructor stub
	Test_assembly=0;
	Nvlast = 0;
	Nhelast = 0;
	Ntype = 0;
	Nv = 0;
	Nv5 = 0;
	Nhe = 0;
	epsilon = nullptr;
	kappa = nullptr;
	kappaPhi = nullptr;
	l0 = nullptr;
	theta0 = nullptr;
	phi0 = nullptr;
}

geometry::~geometry()
{
	// TODO Auto-generated destructor stub
	//delete[] v;
	//delete[] he;
	delete[] epsilon;
	delete[] kappa;
	delete[] kappaPhi;
	delete[] theta0;
	delete[] phi0;
	delete[] l0;
	delete[] mu;
	delete[] vidtoindex;
	delete[] heidtoindex;
	delete[] gb;
	delete[] gdrug;
	if (v.size() > 0)
	{
		v.clear();
	}
	if (he.size() > 0)
	{
		he.clear();
	}
	if (surfheid.size() > 0)
	{
		surfheid.clear();
	}
	if (surfv.size() > 0)
	{
		surfv.clear();
	}

	if (fusionv.size() > 0)
	{
		surfv.clear();
	}
}

void geometry::initialize(int Ntype0)
{
	Nvlast = 0;
	Nhelast = 0;
	Ntype = Ntype0;
	Nv = 0;
	Nhe = 0;
	Nv5 = 0;
	Nsurf = 0;
	vidtoindex = new int[1000000];
	heidtoindex = new int[1000000];
	epsilon = new double[Ntype];
	kappa = new double[Ntype];
	kappaPhi = new double[Ntype];
	theta0 = new double[Ntype];
	l0 = new double[Ntype];
	phi0 = new double[Ntype];
	mu = new double[Ntype];
	gb = new double *[Ntype];
	gdrug = new double *[Ntype];

	for (int i = 0; i < Ntype; i++)
	{
		epsilon[i] = -1;
		kappa[i] = -1;
		kappaPhi[i] = -1;
		phi0[i] = -1;
		theta0[i] = -1;
		l0[i] = -1;
		gb[i] = new double[Ntype];
		gdrug[i]= new double[Ntype];
		for (int j = 0; j < Ntype; j++)
		{
			gb[i][j] = 0;
			gdrug[i][j]=0;
		}
	}

	for (int i = 0; i < 1000000; i++)
	{
		vidtoindex[i] = -1;
		heidtoindex[i] = -1;
	}
}
void geometry::dump_parameters()
{

	FILE *pfile;
	pfile = fopen("geo_param.dat", "w");
	fprintf(pfile, "1");
}

void geometry::update_index()
{
	for (vector<VTX>::iterator it = v.begin(); it != v.end(); ++it)
	{
		//it->hesurfinid=-1;
		if (it->hein.size() > 0)
		{
			it->hein.clear();
		}
		vidtoindex[it->vid] = -1;
	}
	//cout <<"T1" <<endl;

	for (vector<HE>::iterator it = he.begin(); it != he.end(); ++it)
	{
		heidtoindex[it->id] = -1;
	}
	//cout <<"T2" <<endl;
	for (vector<VTX>::iterator it = v.begin(); it != v.end(); ++it)
	{

		vidtoindex[it->vid] = distance(v.begin(), it);
		//cout << "vid :" << it->vid << "vindex" <<  distance(v.begin(),it) <<endl;
	}
	//cout <<"T3" <<endl;
	for (vector<HE>::iterator it = he.begin(); it != he.end(); ++it)
	{
		heidtoindex[it->id] = distance(he.begin(), it);
		//cout << "heid :" << it->id << "heindex" <<  distance(he.begin(),it) <<endl;
	}
	//cout <<"T4" <<endl;
	for (vector<HE>::iterator it = he.begin(); it != he.end(); ++it)
	{
		v[vidtoindex[it->vout]].hein.push_back(it->id);
		//if ((is_surface(it->id))>0 && (is_vsurface(it->vout)>0)) { v[vidtoindex[it->vout]].hesurfinid=it->id; }
		
		//v[vidtoindex[it->vin]].hein.push_back(*it);
	}
	//cout <<"T5" <<endl;

}


void geometry::check_odd_neigh()
{
	int alln=0;
	for (vector<VTX>::iterator it = v.begin(); it != v.end(); ++it)
	{
		alln+=it->vneigh.size();
		//cout<< "vindex is " << distance(v.begin(),it) <<  " vid is" << it->vid ;
		//for (vector<int>::iterator itv = it->vneigh.begin(); itv != it->vneigh.end(); itv++)
		//{
			
		//	cout <<"     neighbors are " << *itv << " vindex is " << vidtoindex[*itv];  
		//	if (vidtoindex[*itv]==-1) { cout << "wrong neighbor" <<endl; exit(-1);}
		//}
		//cout <<endl;
	}
	if (alln%2!=0) { update_neigh(); 
		alln=0;
		for (vector<VTX>::iterator it = v.begin(); it != v.end(); ++it)
		{
			alln+=it->vneigh.size();
			//cout<< "vindex is " << distance(v.begin(),it) <<  " vid is" << it->vid ;
			//for (vector<int>::iterator itv = it->vneigh.begin(); itv != it->vneigh.end(); itv++)
			//{
				
			//	cout <<"     neighbors are " << *itv << " vindex is " << vidtoindex[*itv];  
			//	if (vidtoindex[*itv]==-1) { cout << "wrong neighbor" <<endl; exit(-1);}
			//}
			//cout <<endl;
		}
		if (alln%2!=0) { cout <<" check_neigh odd neighbors" <<endl; exit(-1);}
	}
}

void geometry::update_neigh_vertex(int vid0)
{
	int vindex0=vidtoindex[vid0];
	if (v[vindex0].vneigh.size() > 0)
	{
		v[vindex0].vneigh.clear();
	}
	for (vector<VTX>::iterator it = v.begin(); it != v.end(); ++it)
	{
		if (vid0!= it->vid && connected(vid0,it->vid)<0 ){ //} && next_connected(vid0,it->vid)<0){
			if (veclen(v[vindex0].co,it->co)<2*xi){
				v[vindex0].vneigh.push_back(it->vid);
			}	
		}

	}
}

void geometry::update_neigh()
{
	//cout << "update_neigh()"<<endl;
	for (vector<VTX>::iterator it = v.begin(); it != v.end(); ++it)
	{
		if (it->vneigh.size() > 0)
		{
			it->vneigh.clear();
		}
	}

	for (vector<VTX>::iterator it = v.begin(); it != v.end(); ++it)
	{
		for (vector<VTX>::iterator secondit = v.begin(); secondit != v.end(); ++secondit)
		{
			if (it->vid!= secondit->vid && connected(it->vid,secondit->vid)<0 ) { //&& next_connected(it->vid,secondit->vid)<0){
				if (veclen(it->co,secondit->co)<2*xi){
					it->vneigh.push_back(secondit->vid);
				}	
			}

		}
	}
			

}
void geometry::update_surface()
{ // all heid with no next prev
	//Nsurf=0;
	//cout << "in update surface"<<endl;
	all_neigh=0;
	if (surfheid.size() > 0)
	{
		surfheid.clear();
	}
	if (surfv.size() > 0)
	{
		surfv.clear();
	}
	if (surfvbond.size() > 0)
	{
		surfvbond.clear();
	}

	for (vector<VTX>::iterator it = v.begin(); it != v.end(); ++it)
	{
		it->hesurfinid=-1;
		//it->hesurfoutid=-1;
		if (it->hein.size() > 0)
		{
			it->hein.clear();
		}
		all_neigh+=it->vneigh.size();
	}

	//UPDATE INDICES
	update_index();

	// Update Surfheid , SurfV (he->vin) , Surfvbond
	for (vector<HE>::iterator it = he.begin(); it != he.end(); ++it)
	{

		// update vertex

		if (is_surface(it->id) > 0)
		{
			surfheid.push_back(it->id);
			v[vidtoindex[it->vout]].hesurfinid=it->id;

			surfv.push_back(it->vin);
			//cout << "surface updated with edge id " << it->id << endl<< endl;
		
		}

		if (is_bond_in_surface(it->id) >= 0)
		{
			surfvbond.push_back(it->vin);
		}
		
		if (it->vin == -1 || it->vout == -1 ) 
		{
			cout << " update_surface ! error in vin vout of edge " << it->id << " it->vin " << it->vin <<" it->vout " << it->vout << endl;
			exit(-1);
		}

		if (vidtoindex[it->vin] == -1 || vidtoindex[it->vout] == -1) {
			cout << " update_surface ! error in vin vout of edge " << it->id << " vidtoindex[it->vin] " << vidtoindex[it->vin] <<" vidtoindex[it->vout] " <<vidtoindex[it->vout] << endl;
			exit(-1);

		}

	}
	//cout << "AFTER first set of initialization in update surface"<<endl;
	for (vector<HE>::iterator it = he.begin(); it != he.end(); ++it)
	{
		update_half_edge(it->id);
		
	}
	//cout << "AFTER second set of initialization in update surface"<<endl;
	/****************************************************************************/
	/********************use this only in validation steps **********************/
	/****************************************************************************/
	/*for (vector<HE>::iterator it = he.begin(); it != he.end(); ++it)
	{
		
		int optype = he[heidtoindex[it->opid]].type;
		//cout << "Here 0"<<endl;
		//cout <<"id : "<< it->id <<" type is "<< it->type << " opid :  " <<it->opid<< " opidtype is "  << optype <<endl;
		if ((it->type == 1 && optype != 2) || (it->type == 2 && optype != 1) || (it->type == -1) || (it->type == 3 && optype != 0) || (it->type == 0 && optype != 3))
		{
			cout << "id is" << it->id << "type is" << it->type << " opid id  " << it->opid << " wrong opidtype " << optype << endl;
			exit(-1);
		}
		//cout << it->n[2] << "n[2] it is ---???" <<endl;
		if (it->nextid != -1)
		{ //} && it->previd==-1)  { //this and next
			int nextindex = heidtoindex[it->nextid];

			if (he[nextindex].previd == -1)
			{
				cout << " it->id " << it->id << "it->nextid" << it->nextid << "nextindex" << nextindex << endl;
				cout << "wrong g.he[nextindex].previd in update surface edge" << it->id << " he[nextindex].previd " << he[nextindex].previd << endl;
				exit(-1);
			}
		}
		if (it->previd != -1)
		{ //} && it->previd!=-1)  { //this and next
			int previndex = heidtoindex[it->previd];

			if (he[previndex].nextid == -1)
			{
				cout << "it->previd" << it->previd << "previndex" << previndex << endl;
				cout << "wrong g.he[previndex].nextid in update surfaceedge" << it->id << " he[previndex].nextid " << he[previndex].nextid << endl;
				exit(-1);
			}
		}
	}*/
	/*for (vector<int>::iterator it = surfheid.begin(); it != surfheid.end(); ++it)
	{
		
		//if ((is_surface(it->id))>0 && (is_vsurface(it->vout)>0)) { v[vidtoindex[it->vout]].hesurfinid=it->id; }
		cout << "heid" << *it << " next " << he[heidtoindex[*it]].nextid << " prev " << he[heidtoindex[*it]].previd << endl;
		cout << "vin " << he[heidtoindex[*it]].vin << "vout " << he[heidtoindex[*it]].vout <<endl;
		//v[vidtoindex[it->vin]].hein.push_back(*it);
	}*/

	//XXXXX update nearsurf
	/*for (vector<int>::iterator vt = surfv.begin(); vt != surfv.end(); ++vt)
	{
		vindex0=vidtoindex[vt];
		for (vector<VTX>::iterator it = v.begin(); it != v.end(); ++it)
		{
			if v->nearsurf!=1 and veclen(v[vidtoindex[it->vid]==])
		
		//cout << "id " << *vt << " v[vidtoindex[*vt]].hesurfinid" << v[vidtoindex[*vt]].hesurfinid <<endl;
	}*/
	//TEMP DOUBLE CHECK
	

	update_normals();
	//Nsurf=surfheid.size();
	Nsurf = surfheid.size();
}
int geometry::is_vsurface(int vid)
{ // 1 if
	for (vector<int>::iterator it = surfv.begin(); it != surfv.end(); ++it)
	{
		if (*it == vid)
			return 1;
	}
	return -1;
}

int geometry::is_bond_vsurface(int vid)
{
	for (vector<int>::iterator it = surfvbond.begin(); it != surfvbond.end(); ++it)
	{
		if (*it == vid)
			return 1;
	}
	return -1;
}
int geometry::is_surface(int heid0)
{ // it is on the surface if it has no next prevoius
	int heindex0 = heidtoindex[heid0];
	//int opindex0=heidtoindex[he[heindex0].opid);
	if (he[heindex0].nextid == -1 || he[heindex0].previd == -1)
	{ // on surface
		//cout << " he " << heid0 << " on surface " ;
		return 1;
	}
	//cout << " in is_surface normal is " << he[heindex0].n[0] << " " <<endl;
	return -1; // not on surface
}

int geometry::is_bond_in_surface(int heid0)
{ // it is on the surface if it has no next prevoius
	if (heid0==-1) { cout << "id does not exist !!!!!!! in is_bond_in_surface " <<endl; exit(-1); }
	int heindex0 = heidtoindex[heid0];
	//int opindex0=heidtoindex[he[heindex0].opid);
	//cout<<"in Is_bond_surface) for edge heid0 "<<heid0 << endl;
	//cout<<"he[heindex0].nextid" << he[heindex0].nextid << "and  he[heindex0].previd " << he[heindex0].previd <<endl;
	if (he[heindex0].nextid == -1 && he[heindex0].previd != -1)
	{ // not on surface
		//cout << " he " << heid0 << " is bond_in surface " ;
		//cout << "now return he[heindex0].vin" << he[heindex0].vin <<endl;
		return (he[heindex0].vin);
	}

	//cout << " in is_surface normal is " << he[heindex0].n[0] << " " <<endl;
	return -1; // on surface
}

int geometry::is_bond_out_surface(int heid0)
{ // it is on the surface if it has no next prevoius
	if (heid0==-1) { cout << "id does not exist !!!!!!! in is_bond_out_surface " <<endl; exit(-1); }
	int heindex0 = heidtoindex[heid0];
	//int opindex0=heidtoindex[he[heindex0].opid);
	if (he[heindex0].nextid != -1 && he[heindex0].previd == -1)
	{ // not on surface
		//cout << " he " << heid0 << " not on surface " ;
		return he[heindex0].vout;
	}
	//cout << " in is_surface normal is " << he[heindex0].n[0] << " " <<endl;
	return -1; // on surface
}

int geometry::no_bond_surface(int heid0)
{ // it is on the surface if it has no next prevoius
	int heindex0 = heidtoindex[heid0];
	//int opindex0=heidtoindex[he[heindex0].opid);
	if ((he[heindex0].nextid == -1 && he[heindex0].previd == -1))
	{ // not on surface
		//cout << " he " << heid0 << " not on surface " ;
		return 1;
	}
	//cout << " in is_surface normal is " << he[heindex0].n[0] << " " <<endl;
	return -1; // on surface
}

int geometry::is_same_triangle(int heid0, int heid1)
{
	if (heid0==heid1) return 1;
	int heindex0 = heidtoindex[heid0];
	if (heid1==he[heindex0].opid) return 1;
	int heopindex0 = heidtoindex[he[heindex0].opid];
	int opnextid=he[heopindex0].nextid;
	int opprevid=he[heopindex0].previd;
	if (opnextid != -1) {
		if (heid1==opnextid) return 1;
		if (heid1==he[heidtoindex[opnextid]].opid) return 1;
	}
	if (opprevid != -1) {
		if (heid1==opprevid) return 1;
		if (heid1==he[heidtoindex[opprevid]].opid) return 1;
	}
	return -1;
}


int geometry::not_cross_edge(int heid0, int heid1)
{
	int heindex0 = heidtoindex[heid0];
	int heindex1 = heidtoindex[heid1];
	int heopindex0 = heidtoindex[he[heindex0].opid];
	int heopindex1 = heidtoindex[he[heindex1].opid];

	int nextindex0 = heidtoindex[he[heopindex0].nextid];
	int previndex1 = heidtoindex[he[heopindex1].previd];
	int nextindex1 = heidtoindex[he[heopindex1].nextid];
	int previndex0 = heidtoindex[he[heopindex0].previd];

	if (he[nextindex0].opid == he[previndex1].id)
	{ //cout <<"cross1" << endl;
		return -1;
	}
	if (he[nextindex1].opid == he[previndex0].id)
	{ //cout <<"cross2" << endl;
		return -1;
	}
	//cout << "no cross" <<endl;
	return 1;
}

int geometry::next_open_wedge(int heid0)
{
	if (Nhe < 21)
	{
		return -1;
	}
	
	int heindex0 = heidtoindex[heid0]; // this edge
	int vindex0= vidtoindex[he[heindex0].vout];
	if (v[vindex0].hein.size() < 4) { return -1; }
	int opid0=he[heindex0].opid;
	int heidnext=-1;
		
	for (vector<int>::iterator it = v[vindex0].hein.begin(); it != v[vindex0].hein.end(); ++it)
	{

		if (is_surface(he[heidtoindex[*it]].opid)>0) { heidnext=he[heidtoindex[*it]].opid ;}
	}

	if (heidnext==-1) { cout << "wrong geometry in next_open wedge " <<endl; exit(-1); }
	
	int heindex1=heidtoindex[heidnext];
		
	if (he[heindex0].vout != he[heindex1].vin ) { cout << "Wrong geometry vin vout in next_open wedge " <<endl; exit(-1);}
	
	if ( connected(he[heindex0].vin, he[heindex1].vout) < 0)
	{
		double cangle = dot(he[heidtoindex[opid0]].hevec, he[heindex1].hevec) / (he[heidtoindex[he[heindex0].opid]].l * he[heindex1].l);
		if (cangle > 1)
		{
			cout << "in next_open_wedge" << endl;
			exit(-1);
		}
		//cout << "cangle" <<cangle<<endl;
		if (cangle > .2 && not_cross_edge(he[heindex0].id, he[heindex1].id) > 0)
		{
			return heidnext;
			
		}
	
	}
	return -1;
}

int geometry::pre_open_wedge(int heid0)
{
	if (Nhe < 21)
	{
		return -1;
	}
	double *tempv1 = new double[3];
	double *tempv2 = new double[3];
	int heindex0 = heidtoindex[heid0]; // this edge
	int heindex1;
	int opid0 = he[heindex0].opid;
	
	for (vector<int>::iterator it = surfheid.begin(); it != surfheid.end(); ++it)
	{
		if ((*it != heid0) && (*it != opid0))
		{ // other surface edges
			update_half_edge(*it);
			heindex1 = heidtoindex[*it]; // iterating on surface edges
			update_half_edge(he[heindex1].opid);
			if (he[heindex0].vin == he[heindex1].vout && connected(he[heindex0].vout, he[heindex1].vin) < 0)
			{
				if (v[vidtoindex[he[heindex0].vin]].hein.size() > 4)
				{
					double cangle = dot(he[heindex0].hevec, he[heidtoindex[he[heindex1].opid]].hevec) / (he[heindex0].l * he[heidtoindex[he[heindex1].opid]].l);
					if (cangle > 1)
					{
						cout << "in pre_open_wedge" << endl;
						exit(-1);
					}
					//cout << "cangle" <<cangle<<endl;
					if (cangle > 0.2 && not_cross_edge(he[heindex1].id, he[heindex0].id) > 0)
					{
						return *it;
						
					}
				}
			}
		}
	}

	delete[] tempv1;
	delete[] tempv2;
	return -1;
}

int geometry::open_wedge(int heid0, int *flag)
{
	//if (Nhe<21) { return -1;}
	int heindex0 = heidtoindex[heid0]; // this edge
	int heindex1;
	int opid0 = he[heindex0].opid;
	double *tempvec = new double[3];
	//cout << " in open_wedge _ checking heid0 " << heid0  << " heindex0 is "<< heindex0 <<endl;
	for (vector<int>::iterator it = surfheid.begin(); it != surfheid.end(); ++it)
	{
		if ((*it != heid0) && (*it != opid0))
		{								 // other surface edges
			heindex1 = heidtoindex[*it]; // iterating on surface edges
			//cout << " in open_wedge *it is " << *it << endl;
			//cout << "heindex0 is " << heindex0 << endl;
			//cout << "heindex1 is " << heindex1 <<endl;

			if (he[heindex0].vout == he[heindex1].vin && connected(he[heindex0].vin, he[heindex1].vout) < 0)
			{
				subvec(v[vidtoindex[he[heindex0].vin]].co, v[vidtoindex[he[heindex1].vout]].co, tempvec);
				//cout << "norm(tempvec " << norm(tempvec) <<endl;
				if (norm(tempvec) < l0[1] * 1.5 && not_cross_edge(he[heindex0].id, he[heindex1].id) > 0)
				{
					//cout << " not connected first"<<endl;
					//cout << "connection is at edge " << heindex0 << " at " << he[heindex0].vout << endl;
					//cout << "it is shared with edge " << heindex1 << " at " << he[heindex1].vin <<endl;
					//cout << "other end of edge " << heindex0 << " is " << he[heindex0].vin << endl;
					//cout << "other end of edge " << heindex1 << " is " << he[heindex1].vout <<endl;
					*flag = -1;
					delete[] tempvec;
					return *it;
				}
			}
			if (he[heindex0].vin == he[heindex1].vout && connected(he[heindex0].vout, he[heindex1].vin) < 0)
			{
				subvec(v[vidtoindex[he[heindex0].vout]].co, v[vidtoindex[he[heindex1].vin]].co, tempvec);
				if (norm(tempvec) < l0[1] * 1.5 && not_cross_edge(he[heindex1].id, he[heindex0].id) > 0)
				{
					//cout << " not connected second"<<endl;
					//cout << "connection is at " << he[heindex0].vout << endl;
					*flag = 1;
					delete[] tempvec;
					return *it;
				}
			}
		}
	}
	//cout << " no open wedge returning -1" <<endl;
	*flag = 0;
	delete[] tempvec;
	return -1;
}


int geometry::next_connected(int vid0, int vid1)
{
	for (vector<VTX>::iterator it = v.begin() ; it != v.end(); ++it) {
		
		if (vid0!=it->vid && vid1!=it->vid) {
			if (connected(vid0,it->vid)>0 && connected(vid1,it->vid)>0) {
				return it->vid;
			}
		}

	}
	return -1;
}
int geometry::connected(int vid0, int vid1)
{
	for (vector<HE>::iterator it = he.begin(); it != he.end(); ++it)
	{
		if ((it->vin == vid0 && it->vout == vid1) || (it->vout == vid0 && it->vin == vid1))
		{
			return (it->id);
		}
	}
	return -1;
}

int geometry::connectedH(int heid0, int heid1)
{
	if (he[heidtoindex[heid0]].vin==he[heidtoindex[heid1]].vin) return -1;
	if (he[heidtoindex[heid0]].vout==he[heidtoindex[heid1]].vin) return -1;
	if (he[heidtoindex[heid0]].vout==he[heidtoindex[heid1]].vout) return -1;
	if (he[heidtoindex[heid0]].vin==he[heidtoindex[heid1]].vout) return -1;
	return 1;
}

void geometry::add_vertex(double *xyz)
{
	VTX *vtxi;
	vtxi = new VTX;
	vtxi->vid = Nvlast;
	vtxi->co[0] = xyz[0];
	vtxi->co[1] = xyz[1];
	vtxi->co[2] = xyz[2];
	//vtxi->hesurfinid=-1;
	//vtxi->hesurfoutid=-1;
	vtxi->fusion_vid=-1;
	vidtoindex[Nvlast] = Nv;
	v.push_back(*vtxi);
		
	Nv++;
	Nvlast++;
	delete vtxi;
}

int geometry::delete_vertex(int vid0)
{

	int vindex0 = vidtoindex[vid0];

	v.erase(v.begin() + vindex0);
	vidtoindex[vid0] = -1;
	Nv--;

	return 1;
}


int geometry::remove_neigh(int vid0, int removevid) {
	int vindex0=vidtoindex[vid0];
	int counter=-1;
	if (v[vindex0].vneigh.size()==0) { cout << " has no neighbor in remove neigh for " <<vid0 << endl; exit(-1);}
	for (vector<int>::iterator it = v[vindex0].vneigh.begin(); it != v[vindex0].vneigh.end(); ++it)
	{
		if (*it==removevid) { counter=distance(v[vindex0].vneigh.begin(),it);}
	}
	if (counter!=-1) {
		v[vindex0].vneigh.erase(v[vindex0].vneigh.begin()+counter);
	}
	else {
		cout << " not found neigh in remove neigh"<<endl;
		exit(-1);
	}
	return 1;
}

int geometry::delete_edge(int heid0)
{
	//cout << "IN DELETE EDGE " <<endl;

	int heindex0 = heidtoindex[heid0]; //this edge
	if (is_surface(heid0) < 0)
	{
		cout << " edge is not on surface, cannot remove" << endl;
		cout << " edge " << heid0 << " next is " << he[heindex0].nextid << " previous is " << he[heindex0].previd << endl;
		exit(-1);
		return 0;
	}
	//cout << " heidtoindex[he[heindex0].opid); " << heidtoindex[he[heindex0].opid);
	int opid0 = he[heindex0].opid;
	int opindex0 = heidtoindex[opid0]; //opposite edge
	//cout << "deleting opposite edge id "<< opid0 << " with opindex "<<heidtoindex[opid0]<<endl;

	if (opindex0 > heindex0)
	{
		he.erase(he.begin() + opindex0);
		he.erase(he.begin() + heindex0);
		//cout <<" both removed" <<endl;
	}
	else
	{
		he.erase(he.begin() + heindex0);
		he.erase(he.begin() + opindex0);
		//cout <<" both removed" <<endl;
	}
	//cout << "he.size() is " << he.size() <<endl;
	heidtoindex[heid0] = -1;
	heidtoindex[opid0] = -1;
	Nhe--;
	Nhe--;
	//cout << " in delete edge Nhe is " << Nhe <<endl;

	return 1;
}



void geometry::add_half_edge_type(int vin0, int vout0, int etype)
{
	if (vin0 >= Nvlast || vout0 >= Nvlast)
	{
		cout << "ERROR in update_half_type , Nvlast" << endl;
	}
	else if (vin0 < 0 || vout0 < 0)
	{
		cout << "ERROR in update_half_type , vin0 vout0" << endl;
		exit(-1);
	}
	HE *newhei;
	//HE *newheo;
	newhei = new HE;
	he.push_back(*newhei);
	he_initialize(Nhe, Nhelast, vin0, vout0, etype);
	Nhe++;
	Nhelast++;
	delete newhei;
}

int geometry::add_edge_type(int vin0, int vout0, int etype)
{

	if (vin0 >= Nvlast || vout0 >= Nvlast)
	{
		cout << "ERROR in add_edge_type , Nvlast vin0 is" << vin0 << " vout0 is " << vout0 << endl;
		exit(-1);
	}
	else if (vin0 < 0 || vout0 < 0)
	{
		cout << "ERROR in add_edge_type, vin0 vout0" << endl;
		exit(-1);
	}

	HE *newhei;
	HE *newheo;
	newhei = new HE;
	he.push_back(*newhei);
	he_initialize(Nhe, Nhelast, vin0, vout0, etype);
	//heidtoindex[Nhelast]=Nhe;
	Nhe++;
	Nhelast++;
	//now add opposite edge
	int optype = -1;
	if (etype == 0)
	{
		optype = 3;
	}
	else if (etype == 3)
	{
		optype = 0;
	}
	if (etype == 1)
	{
		optype = 2;
	}
	else if (etype == 2)
	{
		optype = 1;
	}
	if (optype == -1)
	{
		cout << "etype is " << etype << "optype in add edge " << optype << endl;
		exit(-1);
	}
	newheo = new HE;
	he.push_back(*newheo);
	he_initialize(Nhe, Nhelast, vout0, vin0, optype);
	//heidtoindex[Nhelast]=Nhe;
	Nhe++;
	Nhelast++;
	he[Nhe - 1].opid = he[Nhe - 2].id;
	he[Nhe - 2].opid = he[Nhe - 1].id;

	delete newhei;
	delete newheo;
	return 1;
}

void geometry::make_triangle()
{
	double xyz0[3];
	xyz0[0] = 0;
	xyz0[1] = 0;
	xyz0[2] = 0;

	for (int i = 0; i < 3; i++)
	{
		add_vertex(xyz0);
		xyz0[0] = cos(i * PI / 3);
		xyz0[1] = sin(i * PI / 3);
		xyz0[2] = 0;
	}
	for (int i = 0; i < 3; i++)
	{

		int j = i + 1;
		if (i == 2)
		{
			j = 0;
		}
		if (i >= Nvlast || j >= Nvlast)
		{
			cout << "ERROR in make triangle" << endl;
			exit(-1);
		}
		add_edge_type(v[i].vid, v[j].vid, 0);
	}

	//HE *nullHe= new HE;
	set_prev_next(he[0].id, he[4].id, he[2].id);
	set_prev_next(he[2].id, he[0].id, he[4].id);
	set_prev_next(he[4].id, he[2].id, he[0].id);
}

void geometry::make_pentamer()
{
	double xyz0[3];
	xyz0[0] = 0;
	xyz0[1] = 0;
	xyz0[2] = 0;

	for (int i = 0; i < 3; i++)
	{
		add_vertex(xyz0);
		xyz0[0] = cos(i * PI / 3);
		xyz0[1] = sin(i * PI / 3);
		xyz0[2] = 0;
	}
	
	add_half_edge_type(v[0].vid, v[1].vid, 2);
	add_half_edge_type(v[1].vid, v[0].vid, 1);
	add_half_edge_type(v[1].vid, v[2].vid, 0);
	add_half_edge_type(v[2].vid, v[1].vid, 0);
	add_half_edge_type(v[2].vid, v[0].vid, 1);
	add_half_edge_type(v[0].vid, v[2].vid, 2);
	

	cout << endl;

	//}

	//HE *nullHe= new HE;
	set_prev_next(he[0].id, he[4].id, he[2].id);
	set_prev_next(he[2].id, he[0].id, he[4].id);
	set_prev_next(he[4].id, he[2].id, he[0].id);
	for (vector<HE>::iterator it = he.begin(); it != he.end(); it++)
	{
		cout << "updating edge" << it->id << endl;
		update_half_edge(it->id);
		cout << it->id << " TYPE " << he[it->id].type << " opid " << it->opid << " OP TYPE " << he[heidtoindex[it->opid]].type << endl;
		//cout<<
		//cout <<
	}

	cout << " HEREH in pentamer 2" << endl;
	double *vco = new double[3];

	for (int x = 0; x < 3; x++)
	{
		update_surface();
		cout << " after update surface" << endl;
		new_vertex(Nhe - 1, vco);
		//cout << vco[0] <<" " << vco[1] << " "<< vco[2] <<endl;

		force_add_dimer(Nhe - 1, vco, 0, 1);
		update_half_edge(Nhelast - 1);
		update_half_edge(Nhelast - 2);
		update_half_edge(Nhelast - 3);
		update_half_edge(Nhelast - 4);
	}

	update_surface();
	force_add_monomer(1, Nhe - 1, 0);
	update_half_edge(Nhelast - 1);
	update_half_edge(Nhelast - 2);

	delete[] vco;
}

int geometry::opposite_edge(int heid0)
{
	int heindex0 = heidtoindex[heid0];
	int hevin = he[heindex0].vin;
	int hevout = he[heindex0].vout;
	for (vector<HE>::iterator it = he.begin(); it != he.end(); ++it)
	{
		if (it->vin == hevout && it->vout == hevin)
		{
			return it->id;
		}
	}
	cout << " No opposite  Invalid geometry" << endl;
	exit(-1);
}

void geometry::set_prev_next(int heid0, int previd0, int nextid0)
{
	int heindex = heidtoindex[heid0];
	//cout << "heindex" <<endl;
	if ((previd0 != -1) && (nextid0 != -1))
	{
		int nextindex = heidtoindex[nextid0];
		int previndex = heidtoindex[previd0];
		//cout << "he[heindex].vout" << he[heindex].vout << " he[nextindex].vin " << he[nextindex].vin << " he[heindex].vin " <<  he[heindex].vin  << " he[previndex].vout "  << he[previndex].vout << endl;
		if ((he[heindex].vout != he[nextindex].vin) || (he[heindex].vin != he[previndex].vout) || (he[previndex].vin != he[nextindex].vout))
		{
			cout << " WRONG VIN VOUT IN SET_PREV_NEXT" << endl;
			exit(-1);
		}
	}
	he[heindex].nextid = nextid0;
	he[heindex].previd = previd0;
}

int geometry::get_prev_next(int heid0, int *previd0, int *nextid0)
{

	int heindex0 = heidtoindex[heid0];
	int hevin = he[heindex0].vin;
	int hevout = he[heindex0].vout;
	//cout << " hevin is " << hevin << " hevout is " << hevout <<endl;
	for (vector<HE>::iterator it = he.begin(); it != he.end(); ++it)
	{
		//cout << " in get_prev_next heid id " << it->id <<endl;
		//cout << " it->vin is " << it->vin << " it->vout is " << it->vout <<endl;
		if (it->id != heid0)
		{
			if (it->vin == hevout && it->nextid != -1 && it->vout != hevin)
			{ // if next has next
				*nextid0 = it->id;
				*previd0 = it->nextid;
				return 1;
			}
			else if (it->vout == hevin && it->previd != -1 && it->vin != hevout)
			{ //if prev has prev
				*previd0 = it->id;
				*nextid0 = it->previd;
				return 1;
			}
			else if (it->vin == hevout)
			{
				for (vector<HE>::iterator secondit = he.begin(); secondit != he.end(); ++secondit)
				{
					if (secondit->id != it->id && secondit->id != heid0)
					{
						if (it->vout == secondit->vin && secondit->vout == hevin)
						{
							*nextid0 = it->id;
							*previd0 = secondit->id;
							return 1;
						}
					}
				}
			}
		}
	}
	return -1;
}

int geometry::add_dimer(int heid0, gsl_rng *r, int et1, int et2)
{ // adds the dimer if it  has no next prev
	if (heid0 >= Nhelast)
	{
		cout << " in add_dimer ERROR\n";
		exit(-1);
	}
	double *newv = new double[3];
	double *tempv1 = new double[3];
	double *tempv2 = new double[3];
	int success = 1;
	int heindex0 = heidtoindex[heid0];
	//cout << "heindex is" << heindex0<<endl;
	if (is_surface(heid0) < 0)
	{
		cout << " he is not on surface " << endl;
		return -1;
	} // return 0;}
	  //else { cout << " add_dimer he is on surface of he with id Nhe is  " << Nhe << " Nv is " << Nv <<endl;  }
	if (he[heindex0].nextid != -1 || he[heindex0].previd != -1)
	{
		cout << "one end ins connected cannot add dimer" << endl;
		exit(-1);
	}
	//cout << "in add_dimer" <<endl;
	new_vertex(heindex0, tempv1);
	move_p(tempv1, newv, r);
	

	if (check_overlap_centerv(newv) < 0)
	{
		//cout << "overalp in adding vertex! \n \n" ;
		delete[] newv;
		delete[] tempv1;
		delete[] tempv2;
		return -1;
	}

	addvec(v[vidtoindex[he[heindex0].vin]].co, newv, tempv1);
	multvec(tempv1, .5, tempv2);

	if (check_overlap_centerh(tempv2) < 0)
	{
		//cout << "overalp in adding edge dimer! \n \n" ;
		delete[] newv;
		delete[] tempv1;
		delete[] tempv2;
		return -1;
	}

	addvec(v[vidtoindex[he[heindex0].vout]].co, newv, tempv1);
	multvec(tempv1, .5, tempv2);
	if (check_overlap_centerh(tempv2) < 0)
	{
		//cout << "overalp in adding edge dimer! \n \n" ;
		delete[] newv;
		delete[] tempv1;
		delete[] tempv2;
		return -1;
	}
	//}
	add_vertex(newv);

	//cout << " add_dimer after adding a vertex Nhe is   " << Nhe << " Nv is " << Nv <<endl;
	delete[] newv;
	delete[] tempv1;
	delete[] tempv2;
	//delete[] tempv;
	//cout << " in add_dimer add_vertex added a new VTX to v \n " ;
	//cout << "in add_dimer "<<endl;
	//cout<< "g.he[hei].vout is "<< he[heindex0].vout <<endl;
	//int hevout= vidtoindex[he[heindex0].vout];
	//int hevin= vidtoindex[he[heindex0].vin];

	//cout <<" vertex added"<< "Nvlast  is " << Nvlast << endl;
	success = add_edge_type(he[heindex0].vout, Nvlast - 1, et1);

	if (success < 1)
	{
		cout << " add_dimer not successfule deleteing vertex  " << Nhe << " Nv is " << Nv << endl;
		delete_vertex(v[Nv - 1].vid);
		//exit(-1);
		return -1;
	}
	success = add_edge_type(Nvlast - 1, he[heindex0].vin, et2);
	if (success < 1)
	{
		cout << " add_dimer not successfule deleteing edge and vertex  " << Nhe << " Nv is " << Nv << endl;
		delete_edge(he[heidtoindex[Nhelast - 1]].id);
		delete_vertex(Nvlast - 1);
		return -1;
	}
	//cout << " add_dimer after adding a vertex Nhe is   " << Nhe << " Nv is " << Nv <<endl;
	update_half_edge(heid0);
	for (int t = 1; t < 4; t++)
	{
		update_half_edge(Nhelast - t);
	}
	set_prev_next(heid0, Nhelast - 2, Nhelast - 4);
	set_prev_next(Nhelast - 4, heid0, Nhelast - 2);
	set_prev_next(Nhelast - 2, Nhelast - 4, heid0);

	//cout << " in add_dime opid of he[Nhe-1]" << he[Nhe-1].id <<" is " << he[Nhe-1].opid<<endl;
	//cout << "he[heidtoindex[he[Nhe-1].opid)].opid " << he[heidtoindex[he[Nhe-1].opid)].opid <<endl;
	//cout << " in add_dime opid of he[Nhe-3]" << he[Nhe-3].id <<" is " << he[Nhe-3].opid<<endl;
	//cout << "he[heidtoindex[he[Nhe-3].opid)].opid " << he[heidtoindex[he[Nhe-3].opid)].opid <<endl;
	//cout << " set prevoius next (heid0) " << heid0 << he[Nhe-2].id <<he[Nhe-4].id <<endl;

	update_surface();
	update_neigh_vertex(Nvlast-1);
	/*int vindex0=vidtoindex[Nvlast-1];
	if (v[vindex0].vneigh.size() > 0)
	{
		for (vector<int>::iterator it = v[vindex0].vneigh.begin(); it != v[vindex0].vneigh.end(); ++it)
		{	
			update_neigh_vertex(*it);
		}
	}*/	
	
	//cout << "dimer added in G.add_dimer" <<endl;
	return success;
}

int geometry::force_add_dimer(int heid0, double *newv, int typenext, int typeprev)
{ // adds the dimer if it  has no next prev
	if (heid0 >= Nhelast)
	{
		cout << " in add_dimer ERROR\n";
		exit(-1);
	}
	//double *newv=new double[3];
	double *tempv1 = new double[3];
	double *tempv2 = new double[3];
	int success = 1;
	int heindex0 = heidtoindex[heid0];
	//cout << newv[0]<<endl;
	//exit(-1);
	cout << "heindex0" << heindex0 << endl;
	if (is_surface(heid0) < 0)
	{
		return -1;
	} //cout << " he is not on surface " <<endl; return 0;}
	
	add_vertex(newv);

	cout << " add_dimer after adding a vertex Nhe is   " << Nhe << " Nv is " << Nv << endl;
	delete[] tempv1;
	delete[] tempv2;
	//delete[] tempv;
	//cout << " in add_dimer add_vertex added a new VTX to v \n " ;
	//cout << "in add_dimer "<<endl;
	//cout<< "g.he[hei].vout is "<< he[heindex0].vout <<endl;
	//int hevout= vidtoindex[he[heindex0].vout];
	//int hevin= vidtoindex[he[heindex0].vin];

	success = add_edge_type(he[heindex0].vout, v[vidtoindex[Nvlast - 1]].vid, typenext);
	if (success < 0)
	{
		cout << "ERROR" << endl;
		exit(-1);
	}

	//if (success<1) {
	//cout << " add_dimer not successfule deleteing vertex  " << Nhe << " Nv is " << Nv <<endl;
	//elete_vertex(v[Nv-1].vid);
	//exit(-1);
	//return -1;
	//}
	success = add_edge_type(v[vidtoindex[Nvlast - 1]].vid, he[heindex0].vin, typeprev);
	if (success < 0)
	{
		cout << "ERROR" << endl;
		exit(-1);
	}
	//if (success<1) {
	//cout << " add_dimer not successfule deleteing edge and vertex  " << Nhe << " Nv is " << Nv <<endl;
	//	delete_edge(he[Nhe-1].id);
	//	delete_vertex(v[Nv-1].vid);
	//	return -1;
	//}

	set_prev_next(heid0, he[Nhe - 2].id, he[Nhe - 4].id);
	set_prev_next(he[Nhe - 4].id, heid0, he[Nhe - 2].id);
	set_prev_next(he[Nhe - 2].id, he[Nhe - 4].id, heid0);
	update_half_edge(heid0);
	update_half_edge(Nhelast - 2);
	update_half_edge(Nhelast - 4);
	update_half_edge(Nhelast - 1);
	update_half_edge(Nhelast - 3);

	//cout << " set prevoius next (heid0) " << heid0 << he[Nhe-2].id <<he[Nhe-4].id <<endl;

	update_surface();
	return success;
}

int geometry::remove_dimer(int heindex0, int heindexnext0)
{
	int vi = he[heindex0].vout;
	if (he[heindexnext0].vin != vi)
	{
		cout << "remove dimer not two consequtive edges!" << endl;
		return -1;
	}
	delete_edge(he[heindex0].id);
	delete_edge(he[heindexnext0].id);
	delete_vertex(v[vi].vid);
	return 1;
}

void geometry::make_hexamer()
{
	/*make_triangle();
	update_surface();
    update_normals();

	//for (int i=0; i<5; i++){ // make pentamer
	double *vco;
	new_vertex(Nhe-1,vco);
	 
	force_add_dimer(Nhe-1,vco,1,1);
	update_surface();
    update_normals();
	cout << "1 dimer added" <<endl;
	new_vertex(Nhe-1,vco);
	force_add_dimer(Nhe-1,vco,1,0);
	update_surface();
	update_normals();
	cout << "2 dimer added" <<endl;
	new_vertex(Nhe-1,vco);
	force_add_dimer(Nhe-1,vco,0,0);
	update_surface();
	update_normals();
	cout << "3 dimer added" <<endl;
	new_vertex(Nhe-1,vco);
	force_add_dimer(Nhe-1,vco,1,1);
	update_surface();
    update_normals();
	cout << "4 dimer added" <<endl;
	force_add_monomer(1,Nhe-1,1);
	update_surface();
    update_normals();
	delete[] vco;
        //for (int rstep=0; rstep<100; rstep++){
        //    move_vertex(g,r);
		//    g.update_surface();
            //dump_lammps_data_file(g, frame++);
        //} 
        //dump_lammps_data_file(g, frame++);
	*/
}

void geometry::he_initialize(int heindex, int heid0, int vin0, int vout0, int etype)
{
	he[heindex].id = heid0;
	he[heindex].opid = -1;
	he[heindex].vin = vin0;
	he[heindex].vout = vout0;
	he[heindex].nextid = -1;
	he[heindex].previd = -1;

	int vindex0 = vidtoindex[vin0];
	int vindex1 = vidtoindex[vout0];
	double *tempv = new double[3];
	subvec(v[vindex0].co, v[vindex1].co, tempv);
	he[heindex].hevec[0] = tempv[0];
	he[heindex].hevec[1] = tempv[1];
	he[heindex].hevec[2] = tempv[2];
	he[heindex].l = norm(tempv);
	he[heindex].type = etype;
	//he[heindex].nextangle=0;
	delete[] tempv;
	he[heindex].din = 0;
	he[heindex].dout = 0;
	heidtoindex[heid0] = heindex;
}

void geometry::update_edge(int heid0)
{ //hevec hecent l n
	int heindex = heidtoindex[heid0];
	int vindex0 = vidtoindex[he[heindex].vin];
	int vindex1 = vidtoindex[he[heindex].vout];
	//cout << "vindex0 " <<vindex0<<endl;
	//cout << "vindex1 " <<vindex1<<endl;

	double *tempv = new double[3];
	subvec(v[vindex0].co, v[vindex1].co, tempv);
	//cout << "tempv " << tempv[0] << " "<< tempv[1] <<" "<< tempv[2] <<endl;
	he[heindex].hevec[0] = tempv[0];
	he[heindex].hevec[1] = tempv[1];
	he[heindex].hevec[2] = tempv[2];
	he[heindex].l = norm(tempv);
	multvec(he[heindex].hevec, .5, tempv);
	he[heindex].hecent[0] = v[vindex0].co[0] + tempv[0];
	he[heindex].hecent[1] = v[vindex0].co[1] + tempv[1];
	he[heindex].hecent[2] = v[vindex0].co[2] + tempv[2];
	if (he[heindex].opid == -1)
	{
		cout << "invalid geometry ";
		exit(-1);
	}

	int opindex = heidtoindex[he[heindex].opid];
	//cout << "opid" <<opindex <<endl;

	if (he[opindex].opid != heid0)
	{
		cout << " opid not set  for op of edge " << heid0 << endl;
		cout << " heindex " << heindex << endl;

		for (vector<HE>::iterator it = he.begin(); it != he.end(); it++)
		{
			cout << "EDGE  " << distance(he.begin(), it) << " id " << it->id << " vin " << vidtoindex[it->vin] << "vout " << vidtoindex[it->vout] << endl;
			cout << "      opid " << it->opid << " nextid " << it->nextid << " previd " << it->previd << endl;
			cout << "      hevec " << it->hevec[0] << " " << it->hevec[1] << " " << it->hevec[2] << " l is " << it->l << endl;
			cout << "      normal " << it->n[0] << " " << it->n[1] << " " << it->n[2] << " " << endl
				 << endl;
		}
		exit(-1);
	}

	subvec(v[vindex0].co, v[vindex1].co, tempv);
	he[opindex].hevec[0] = tempv[0];
	he[opindex].hevec[1] = tempv[1];
	he[opindex].hevec[2] = tempv[2];
	he[opindex].l = norm(tempv);
	multvec(he[opindex].hevec, .5, tempv);
	addvec(v[vindex0].co, tempv, he[opindex].hecent);

	get_normal(heid0);
	get_normal(he[opindex].id);
	//now everything for
	delete[] tempv;
}

void geometry::update_half_edge(int heid0)
{ //hevec hecent l n
	if (heid0 == -1)
	{
		cout << " -1 as heid0 in update" << endl;
		exit(-1);
	}
	int heindex = heidtoindex[heid0];
	//if (heindex==-1)
	int vindex0 = vidtoindex[he[heindex].vin];
	int vindex1 = vidtoindex[he[heindex].vout];
	//cout << "vindex0 " << vindex0 <<endl;
	double *tempv = new double[3];
	subvec(v[vindex0].co, v[vindex1].co, tempv);
	he[heindex].hevec[0] = tempv[0];
	he[heindex].hevec[1] = tempv[1];
	he[heindex].hevec[2] = tempv[2];
	he[heindex].l = norm(tempv);
	multvec(he[heindex].hevec, .5, tempv);
	addvec(v[vindex0].co, tempv, he[heindex].hecent);

	if (he[heindex].opid == -1)
	{
		he[heindex].opid = opposite_edge(heid0);
		cout << " updateing the opposite" << endl;
	}

	delete[] tempv;
}

int geometry::add_monomer(int nextofnewid, int prevofnewid, int et)
{ //nextofnew prevofnew
	int heid0 = prevofnewid;
	//cout << "in add_monomer" << "heid0 is prevofnewid " << heid0 <<endl;
	int heid1 = nextofnewid;
	//cout << "in add_monomer" << "heid1 is nextofnewid " << heid1 <<endl;
	int heprevindex = heidtoindex[heid0]; // this is prev for new edge
	int henextindex = heidtoindex[heid1]; // this is going to be next of new edge
	int newvin = he[heprevindex].vout;
	int newvout = he[henextindex].vin;
	double *tempv1 = new double[3];
	double *tempv2 = new double[3];
	addvec(v[vidtoindex[newvin]].co, v[vidtoindex[newvout]].co, tempv1);
	multvec(tempv1, .5, tempv2);

	if (check_overlap_centerh(tempv2) < 0)
	{
		cout << "overalp in adding edge monomer! \n \n";
		delete[] tempv1;
		delete[] tempv2;
		return -1;
	}

	if (add_edge_type(newvin, newvout, et) < 0)
	{
		cout << "coudnot add the edge " << endl;
		return -1;
	};
	update_half_edge(Nhelast - 2);
	update_half_edge(Nhelast - 1);
	update_half_edge(heid0);
	set_prev_next(Nhelast - 2, heid0, heid1);
	set_prev_next(heid1, Nhelast - 2, heid0);
	set_prev_next(heid0, heid1, Nhelast - 2);
	set_prev_next(Nhelast - 1, -1, -1);

	update_surface();
	//delete[] tempv1; delete[] tempv2;
	return 1;
}

int geometry::force_add_monomer(int nextofnewid, int prevofnewid, int etype)
{ //nextofnew prevofnew
	int heid0 = prevofnewid;
	//cout << "in add_monomer" << "heid0 is prevofnewid " << heid0 <<endl;
	int heid1 = nextofnewid;
	//cout << "in add_monomer" << "heid1 is nextofnewid " << heid1 <<endl;
	int heprevindex = heidtoindex[heid0]; // this is prev for new edge
	int henextindex = heidtoindex[heid1]; // this is going to be next of new edge
	int newvin = he[heprevindex].vout;
	int newvout = he[henextindex].vin;
		
	if (add_edge_type(newvin, newvout, etype) < 0)
	{
		cout << "coudnot add the edge " << endl;
		exit(-1);
	};
	set_prev_next(Nhelast - 2, heid0, heid1);
	set_prev_next(heid1, Nhelast - 2, heid0);
	set_prev_next(heid0, heid1, Nhelast - 2);
	set_prev_next(Nhelast - 1, -1, -1);
	update_half_edge(heid0);
	update_half_edge(Nhelast - 2);
	update_half_edge(Nhelast - 1);
	//update_surface();
	//delete[] tempv1; delete[] tempv2;
	return 1;
}
int geometry::add_monomer_dimer(int heid0)
{

	if (is_surface(heid0) < 0)
	{
		cout << " cannot add not on the surface !" << endl;
		exit(-1);
	}
	update_surface();
	int flagprevnext;
	int xid = open_wedge(heid0, &flagprevnext);

	//cout << "heid is " << heid0 << "xid is" << xid <<endl;
	if (xid >= 0)
	{
		if (flagprevnext > 0)
		{
			//cout << " add-monomer_dimer - monomer with next" <<endl;
			add_monomer(xid, heid0, 2);
			//cout << " monomer added Nhe is " << Nhe <<endl;
			//dump_lammps_data_file(g, frame++);
		}
		else if (flagprevnext < 0)
		{
			//cout << " add-monomer with prev" <<endl;
			add_monomer(heid0, xid, 2);
			//cout << " monomer added Nhe is" << Nhe <<endl;
			//dump_lammps_data_file(g, frame++);
		}
		else
		{
			//cout << "error in main add_monomer_dimer" <<endl;
		}

		//dump_lammps_data_file(g, frame++);
		//exit(-1);
	}
	else
	{
		//cout << "add_dimer Nhe is " << Nhe<<endl;
		//add_dimer(heid0,r);
	}
	return 1;
}


int geometry::remove_monomer_dimer(int heid0, gsl_rng *r)
{
	if (!(is_surface(heid0) > 0))
	{
		cout << "not on surface cannot remove" << endl;
		return -1;
	}
	int heindex0 = heidtoindex[heid0];				 // index of this edge
	int heopindex0 = heidtoindex[he[heindex0].opid]; // indexd of opposite edge
	int nextopid0 = he[heopindex0].nextid;			 // indexd of next of opposite edge
	int prevopid0 = he[heopindex0].previd;			 // indexd of prev of opposite edge
	int heid1 = he[heidtoindex[nextopid0]].opid;	 // now back to this side
	int heid2 = he[heidtoindex[prevopid0]].opid;	 // afetr vertex
	if (is_surface(heid1) < 0 && is_surface(heid2) < 0 && is_vsurface(he[heidtoindex[nextopid0]].vout) < 0)
	{
		int x = delete_edge(heid0);
		if (x < 0)
		{
			cout << "in remove monomer" << endl;
			exit(-1);
		}
		else
		{
			set_prev_next(nextopid0, -1, -1);
			set_prev_next(prevopid0, -1, -1);
			if (true)
			{ //gsl_rng_uniform(r) < crit) {
				return 1;
			}
			else
			{
				add_monomer(nextopid0, prevopid0, 0);
				return -1;
			}
		}
	}
	else
	{

		if (is_surface(heid1) > 0)
		{
			//cout << "\n \n TRYING DIMER REMOVAL - with its next" <<endl;
			/*for (vector<HE>::iterator it =  he.begin() ; it !=  he.end(); it++) {
                cout << "edge " <<distance( he.begin(),it) <<" id " << it->id << " opid " << it-> opid <<"\t";
                cout << "       " <<" vin " << vidtoindex[it->vin] << " vout " << vidtoindex[it->vout] <<"\t";
                cout << "        " << " nextid " << it->nextid << " previd " << it->previd <<endl;
            }*/
			int vi = he[heopindex0].vout;
			//double de=- dimer_energy( he[heidtoindex[heid0)].opid, he[heidtoindex[heid1)].opid);
			//double de=-( stretch_energy(heidtoindex[heid0))+ bend_energy(heidtoindex[heid0)));
			//de-=( stretch_energy(heidtoindex[heid1))+ bend_energy(heidtoindex[heid1)));
			//de-= bend_energy(heidtoindex[ he[heidtoindex[heid1)].nextid));
			//de-= gb*2- mu;
			//cout <<"current triangle with their opp " << endl;
			int success = delete_edge(heid0);
			update_index();
			success *= delete_edge(heid1); // next of op
			update_index();
			//double vp = 4/3.*4*atan(1.)* xi* xi* xi;
			int x = delete_vertex(vi);
			if (x < 0)
			{
				exit(-1);
			}
			else
			{
				//cout  <<"currentlt prev next of prevopid0 = " << prevopid0 << " are " <<  he[heidtoindex[prevopid0)].previd  << " and "  <<  he[heidtoindex[prevopid0)].nextid <<endl;
				//cout << "  set_prev_next(prevopid0,-1,-1); " << prevopid0 << endl;
				set_prev_next(prevopid0, -1, -1); //prev of op
				//cout  <<"after  prev next of prevopid0 = " << prevopid0 << " are " <<  he[heidtoindex[prevopid0)].previd  << " and "  <<  he[heidtoindex[prevopid0)].nextid <<endl;
				// double e2 = compute_energy();

				//cout << " de is " << de <<endl;
				//double crit = exp(-de/ T)/(2*vp);///vp; //* z* K* K);
				//cout << " crit is " << crit << endl;
				if (true)
				{ //gsl_rng_uniform(r) < crit) {
					/*for (vector<HE>::iterator it =  he.begin() ; it !=  he.end(); it++) {
                        cout << "edge " <<distance( he.begin(),it) <<" id " << it->id << " opid " << it-> opid <<"\t";
                        cout << "       " <<" vin " << vidtoindex[it->vin] << " vout " << vidtoindex[it->vout] <<"\t";
                        cout << "        " << " nextid " << it->nextid << " previd " << it->previd <<endl;
                    }*/
					//cout << " crit is " << crit << endl;
					//cout << " DIMER REMOVEd" << endl;
					return 1;
				}
				else
				{
					add_dimer(prevopid0, r, 0, 0);
					return -1;
				}
			}
		}

		else if (is_surface(heid2) > 0)
		{

			int vi = he[heopindex0].vin;
			//double de=- dimer_energy( he[heidtoindex[heid0)].opid, he[heidtoindex[heid1)].opid);

			int success = delete_edge(heid0);
			update_index();
			success *= delete_edge(heid2);
			update_index();

			//double vp = 4/3.*4*atan(1.)* xi* xi* xi;

			int x = delete_vertex(vi);
			if (x < 0 || success < 1)
			{
				exit(-1);
			}
			else
			{
				//cout << "  set_prev_next( he[heidtoindex[heid1)].opid,-1,-1);" <<  he[heidtoindex[heid1)].opid <<endl;
				// cout << " get_index( he[heidtoindex[heid1)].opid)" << heidtoindex[ he[heidtoindex[heid1)].opid) <<endl;
				set_prev_next(nextopid0, -1, -1);
				//double e2 = compute_energy();
				//double de=e2-e1;

				//double de=0;
				//cout << " de is " << de <<endl;
				//double crit = exp(-de/ T)/(2*vp);///vp;
				//cout << " de is " << de <<endl;
				//double crit = exp(-de/ T)/(2*vp* z* K* K);
				//cout << " crit is " << crit << endl;
				if (true)
				{ //gsl_rng_uniform(r) < crit) {
					//cout << " accepted crit is " << crit << endl;

					/*for (vector<HE>::iterator it =  he.begin() ; it !=  he.end(); it++) {
                        cout << "edge " <<distance( he.begin(),it) <<" id " << it->id << " opid " << it-> opid <<endl;
                        cout << "       " <<" vin " << vidtoindex[it->vin] << " vout " << vidtoindex[it->vout] <<endl;
                        cout << "        " << " nextid " << it->nextid << " previd " << it->previd <<endl;
                    }*/
					//cout << " DIMER REMOVED" <<endl;;
					return 1;
				}
				else
				{
					add_dimer(nextopid0, r, 0, 0);
					return -1;
				}
			}
		}
	}
	//return -1;
	return 1;
}

void geometry::move_v(double *pi, double *pf, gsl_rng *r)
{
	double *tempvec, *fvec;
	tempvec = new double[3];
	fvec = new double[3];
	
	//random direction
	randvec(tempvec, r);
	//cout << " tempvec is " << tempvec[0] << " " <<tempvec[1] << " " <<tempvec[2] <<endl;
	double tempveclen=norm(tempvec);
	if (tempveclen>1.001) { cout<< " in move_v veclen > 1 veclen is "<< tempveclen << endl; exit(-1); }
	//double d= .2* gsl_rng_uniform(r);
	double d = gsl_rng_uniform(r) * 2.0 / sqrt(epsilon[0]);
	multvec(tempvec, d, fvec);
	addvec(fvec, pi, pf);

	delete[] tempvec;
	delete[] fvec;
}

void geometry::move_p(double *pi, double *pf, gsl_rng *r)
{
	double *tempvec, *fvec;
	tempvec = new double[3];
	fvec = new double[3];

	//random direction
	randvec(tempvec, r); // random direction with unit length
	//double sigmav=.2;
	double d = xi * (gsl_rng_uniform(r)); // length of change between 0 and xi 
	//double d = xi * gsl_ran_gaussian(r, sigmav);
	multvec(tempvec, d, fvec); // update the vector length
	addvec(fvec, pi, pf); // apply the position change

	delete[] tempvec;
	delete[] fvec;
}
void geometry::new_vertex(int heindex0, double *newv)
{

	double *tempvec, *fvec; // , *p;

	tempvec = new double[3];
	fvec = new double[3];
	//p = new double[3];
	//hecenter(heindex0,vcenter); // find edge center
	int heopindex0 = heidtoindex[he[heindex0].opid];
	//cout << "hecenter is " << he[heindex0].hecent[0] <<" " <<he[heindex0].hecent[1] <<" " <<he[heindex0].hecent[2] <<" " <<endl;

	cross(he[heopindex0].n, he[heindex0].hevec, tempvec); //direction to find the vector
	//cout<< "cross tempvec is "<<tempvec[0] << " " <<  tempvec[1] << " " << tempvec[2] << " " <<endl;
	//cout<< "norm tempvec is" << norm(tempvec)<<endl;
	//cout<< "hevec is " << he[heindex0].hevec[0] <<" " << he[heindex0].hevec[1] <<" "<< he[heindex0].hevec[2] <<endl;
	//cout<< " normal of opedge is " << he[heopindex0].n[0] <<" "<<he[heopindex0].n[1] <<" "<<he[heopindex0].n[2] <<endl;

	multvec(tempvec, 1 * .86 / norm(tempvec), fvec); // length correposonding to one triangle
	//cout<< "fvec is" << fvec[0] <<" " << fvec[1] << " " <<fvec[2] <<endl;
	//cout<< "norm (fvec) is " << norm(fvec)<<endl;

	//multvec(fvec,1,tempvec);

	rotatevec(fvec, he[heopindex0].hevec, (theta0[1] + theta0[0]) / 2, tempvec); // add this later
	//cout<< "new tempvec is "<<tempvec[0] << " " <<  tempvec[1] << " " << tempvec[2] << " " <<endl;
	/*double vx,vy,vz,ex,ey,ez,nm;//,vxr,vyr,vzr;
	vx=fvec[0];
	vy=fvec[1];
	vz=fvec[2];
	//finding normalized vector
	ex = he[heopindex0].hevec[0];
	ey = he[heopindex0].hevec[1];
	ez = he[heopindex0].hevec[2];

	//cout << "hevec is" << he[heindex0].hevec[0] << " " << he[heindex0].hevec[1] << " " << he[heindex0].hevec[2] <<endl;
	nm = norm(he[heopindex0].hevec);
	ex /= nm;
	ey /= nm;
	ez /= nm;
	double ct = cos((theta0[0]+theta0[0])/2);
	double st = sin((theta0[0]+theta0[0])/2);
	double mct = 1-cos((theta0[0]+theta0[0])/2);
	//rotating
	tempvec[0] = vx*(ct+ex*ex*mct) + vy*(ex*ey*mct-ez*st) + vz*(ex*ez*mct+ey*st);
	tempvec[1] = vx*(ey*ex*mct+ez*st) + vy*(ct+ey*ey*mct) + vz*(ey*ez*mct-ex*st);
	tempvec[2] = vx*(ez*ex*mct-ey*st) + vy*(ez*ey*mct+ex*st) + vz*(ct+ez*ez*mct);*/

	//find p
	//cout << "in new_vertex 444 temvec is " <<tempvec[0] <<" "<<tempvec[1] <<" " <<tempvec[2] << endl;
	//cout << "in new_vertex 444bbb hecent is " <<he[heindex0].hecent[0] <<" "<<he[heindex0].hecent[1] <<" " <<he[heindex0].hecent[2] << endl;

	newv[0] = he[heindex0].hecent[0] + tempvec[0];
	newv[1] = he[heindex0].hecent[1] + tempvec[1];
	newv[2] = he[heindex0].hecent[2] + tempvec[2];
	//addvec(he[heindex0].hecent,tempvec,newv);
	//cout<< "in new_vertex 555 newv is " <<newv[0] <<" "<<newv[1] <<" " <<newv[2] << endl;
	delete[] tempvec;
	delete[] fvec;
	//delete[] p;
}

void geometry::hecenter(int heindex, double *vcenter)
{
	//cout << "he[heindex].vin" << he[heindex].vin;
	int vindexi = vidtoindex[he[heindex].vin];
	//cout << " he[heindex].vout " << he[heindex].vout;
	int vindexj = vidtoindex[he[heindex].vout];
	double *tempvec;
	tempvec = new double[3];
	//cout << "in he center" <<endl;
	addvec(v[vindexi].co, v[vindexj].co, tempvec);
	multvec(tempvec, .5, vcenter);

	delete[] tempvec;

	//cout << "in hecenter hei " << vi[0] << endl;
}
void geometry::helen(int heindex, double *helen)
{
	int vindexi = vidtoindex[he[heindex].vin];
	//cout << " he[heindex].vout " << he[heindex].vout;
	int vindexj = vidtoindex[he[heindex].vout];
	double *tempvec;
	tempvec = new double[3];
	//cout << "in he center" <<endl;
	subvec(v[vindexi].co, v[vindexj].co, tempvec);

	*helen = norm(tempvec);
	//cout << " in helen " << *helen<< endl;
	delete[] tempvec;
}


int geometry::get_normal(int heid0)
{
	int heindex0 = heidtoindex[heid0];
	//cout << "in get_normal" <<endl;
	update_half_edge(heid0);
	//int etype=he[heindex0].type;
	if (he[heindex0].nextid == -1 && he[heindex0].previd == -1)
	{
		he[heindex0].n[0] = 0;
		he[heindex0].n[1] = 0;
		he[heindex0].n[2] = 0;
		//cout << " on surface no normal! " <<endl;
		return (-1);
	}
	
	double *vec2 = new double[3];
	if (he[heindex0].nextid != -1) {
		update_half_edge(he[heindex0].nextid);
		int nextindex0 = heidtoindex[he[heindex0].nextid];
		cross(he[heindex0].hevec, he[nextindex0].hevec, vec2);
	}
	else {
		update_half_edge(he[heindex0].previd);
		int previndex0 = heidtoindex[he[heindex0].previd];
		cross(he[previndex0].hevec, he[heindex0].hevec, vec2);
	}
	double nm = norm(vec2);
	

	he[heindex0].n[0] = vec2[0] / nm;
	he[heindex0].n[1] = vec2[1] / nm;
	he[heindex0].n[2] = vec2[2] / nm;

	//cout<< "normal is"  << "he[heindex0].n[0] "<< he[heindex0].n[0] << " he[heindex0].n[1] "<<he[heindex0].n[1] << "he[heindex0].n[2] " <<he[heindex0].n[2] <<endl;
	//cout <<endl;
	//delete[] vec;
	delete[] vec2; //delete[] nNext;
	return 1;
}



void geometry::update_normals()
{
	for (vector<HE>::iterator it = he.begin(); it != he.end(); ++it)
	{
		//if (it->nextid!=-1 && it->previd!=-1) {
		get_normal(it->id); //}
							//cout << "update normal of edge" << it->id << it->n[0] << " "<< it->n[1] << " " <<it->n[2] <<endl;
	}
	//for (vector<HE>::iterator it = he.begin() ; it != he.end(); ++it) {
	//if (it->nextid!=-1 && it->previd!=-1) {
	//get_normal(it->id); //}
	//cout << "After update normal of edge" << it->id << it->n[0] << " "<< it->n[1] << " " <<it->n[2] <<endl;
	//cout << "normal of op edge" << it->opid << he[heidtoindex[it->opid)].n[0] << " " << he[heidtoindex[it->opid)].n[1] <<" " << he[heidtoindex[it->opid)].n[2]<<endl;
	//}
}

/*int geometry::check_overlap_hesurf(int heidsurf0)
{
	if (is_surface(heidsurf0)<0) { cout << " not on surface in check overlap surface"<<endl; exit(-1);}
	if (Nhe <= 6)
		return 1;
	double *vtemp = new double[3];
	////cout<<"check overlap " << endl <<endl;

	//update_half_edge(heidsurf0);

	int heindex0=heidtoindex[heidsurf0];
	for (vector<int>::iterator ithe = surfheid.begin(); ithe != surfheid.end(); ithe++)
	{
		//update_half_edge(*ithe);
		if (*ithe != heidsurf0)
		{ //cout << " he id " << it->id<< endl;
			//cout << "it->id " << it->id <<endl;
			//cout << "heid0 " << heid0 << endl;

			int heindex1=heidtoindex[*ithe];
			subvec(he[heindex0].hecent, he[heindex1].hecent, vtemp);
			//cout <<"it->hecents "<< it->hecent[0] << " " << it->hecent[1] << " " << it->hecent[2] << " " <<endl;
			//cout << "he[heidtoindex[heid0)].hecent " << he[heidtoindex[heid0)].hecent[0] << " " << he[heidtoindex[heid0)].hecent[1] << " "  << he[heidtoindex[heid0)].hecent[2] << " " <<endl;
			//cout << " vtemp[0] " << vtemp[0]<< " vtemp[1] " << vtemp[1]<< " vtemp[2] " << vtemp[2] <<endl;
			//cout << "norm is " <<  norm(vtemp) <<endl;
			if (norm(vtemp) < l0[0] * .2)
			{ //cout << " distance between "<< it->id << " and " << heid0 << " is " << norm(vtemp) <<endl;
				delete[] vtemp;
				return -1;
			}
		}
	}
	

	//cout << "no overlap !\n";	make

	delete[] vtemp;
	return 1;
}*/

int geometry::check_overlap_he(int heid0)
{
	if (Nhe <= 6)
		return 1;
	//double *vtemp = new double[3];
	////cout<<"check overlap " << endl <<endl;
	update_half_edge(heid0);

	int heindex0=heidtoindex[heid0];
	for (vector<HE>::iterator it = he.begin(); it != he.end(); it++)
	{
		update_half_edge(it->id);
		if (!(it->id == heid0) && !(it->opid == heid0))
		{ //cout << " he id " << it->id<< endl;
			//cout << "it->id " << it->id <<endl;
			//cout << "heid0 " << heid0 << endl;
			double temp= veclen(it->hecent, he[heindex0].hecent);
			//cout <<"it->hecents "<< it->hecent[0] << " " << it->hecent[1] << " " << it->hecent[2] << " " <<endl;
			//cout << "he[heidtoindex[heid0)].hecent " << he[heidtoindex[heid0)].hecent[0] << " " << he[heidtoindex[heid0)].hecent[1] << " "  << he[heidtoindex[heid0)].hecent[2] << " " <<endl;
			//cout << " vtemp[0] " << vtemp[0]<< " vtemp[1] " << vtemp[1]<< " vtemp[2] " << vtemp[2] <<endl;
			//cout << "norm is " <<  norm(vtemp) <<endl;
			
			if  (is_same_triangle(it->id,heid0)<0  && temp < l0[0] * .2)
			{ //cout << " distance between "<< it->id << " and " << heid0 << " is " << norm(vtemp) <<endl;
				//delete[] vtemp;
				return -1;
			}

		}
	}
	//cout << " vertices " << v.size();
	//exit(-1);
	int vind = vidtoindex[he[heindex0].vout];
	int vind2 = vidtoindex[he[heindex0].vin];
	for (vector<VTX>::iterator it = v.begin(); it != v.end(); ++it)
	{
		if (it->vid != v[vind].vid && it->vid != v[vind2].vid)
		{
			double temp=veclen(it->co, v[vind].co);
			//cout << " vvvv " << it->vid<< "vtemp[0] " << vtemp[0]<< " vtemp[1] " << vtemp[1]<< " vtemp[2] " << vtemp[2] <<endl;
			//cout << "vvvv " << it->vid<< " norm is " <<  norm(vtemp) <<endl;
			if (temp < l0[0] * .30)
			{
				//delete[] vtemp;
				return -1;
			}
		}
	}

	//cout << "no overlap !\n";	make

	//delete[] vtemp;
	return 1;
}


int geometry::do_intersect(int heid1, int heid2){
	if (is_surface(heid1)>0 or is_surface(heid2)>0){
		return 0;
	}	
	//int n_vertices_shared=shared_vertices( heid1,  heid2);
	//if (n_vertices_shared==0){
	int vi1=he[heidtoindex[heid1]].vin;
	int vi2=he[heidtoindex[heid1]].vout;
	int nextid1=he[heidtoindex[heid1]].nextid;
	int vi3=he[heidtoindex[nextid1]].vout;
	int vj1=he[heidtoindex[heid2]].vin;
	int vj2=he[heidtoindex[heid2]].vout;
	int nextid2=he[heidtoindex[heid2]].nextid;
	int vj3=he[heidtoindex[nextid2]].vout;


	if (vi1==vj1 or vi1==vj2 or vi1==vj3) { return 0;}
	if (vi2==vj1 or vi2==vj2 or vi2==vj3) { return 0;}
	if (vi3==vj1 or vi3==vj2 or vi3==vj3) { return 0;}

	return tri_tri_overlap_test_3d(v[vidtoindex[vi1]].co,v[vidtoindex[vi2]].co ,v[vidtoindex[vi3]].co, 
			v[vidtoindex[vj1]].co,v[vidtoindex[vj2]].co ,v[vidtoindex[vj3]].co);
	
	
	return 0;
}

int geometry::find_overlap_g(int vid0 ) {
	if (Nhe<22) { return 1;}
	//cout << "in _check_overlap_g" <<endl;
	int vindex0=vidtoindex[vid0];   
	double mindis=.1;
	for (vector<int>::iterator itv = v[vindex0].vneigh.begin(); itv != v[vindex0].vneigh.end(); itv++)
	{	
		//cout<< " "
		int vindex1=vidtoindex[*itv];  //neighbor v
		if (vindex1!=-1){
			//if (veclen(v[vindex0].co,v[vindex1].co )<.1) { return -1;} // check vertex vertx overlap
			
			//if (is_vsurface(*itv)>0 && is_vsurface(vid0)>0 ) {mindis=.1; }
			//else { mindis=.2;}
			if (veclen(v[vindex1].co, v[vindex0].co)< l0[0] * mindis) { return -1; }
		

			for (vector<int>::iterator itneigh = v[vindex1].hein.begin(); itneigh != v[vindex1].hein.end(); itneigh++) //neigbor vertex hein s
			{
				int heid1=*itneigh;
				
				if (is_surface(heid1)>0){ heid1=he[heidtoindex[heid1]].opid;}
				int heindex1=heidtoindex[heid1];
				for (vector<int>::iterator it = v[vindex0].hein.begin(); it != v[vindex0].hein.end(); it++) //this vertex hein s
				{
					int heid0=*it;
					int cH=connectedH(heid0,heid1);
					int cvin=connected(vid0,he[heid0].vin);
					int cvout=connected(vid0,he[heid0].vout);
					if (cH==0 && cvin==0 && cvout==0 && heid0!=he[heindex1].id && heid0!=he[heindex1].opid){
						if (is_surface(heid1)>0){ heid0=he[heidtoindex[heid1]].opid;}
						int heindex0=heidtoindex[heid0];
						if (veclen(v[vindex0].co,he[heindex1].hecent )<l0[0] *.4) { cout << "In find_overlap overlap heindex1 vindex0 " << heindex1 <<" " << vindex0 <<endl; exit(-1);}
						if (veclen(v[vindex1].co,he[heindex0].hecent )<l0[0] *.4) { cout << "In find_overlap overlap heindex0 vindex1 " << heindex0 <<" " << vindex1 <<endl; exit(-1);}
						
						if (veclen(he[heindex0].hecent,he[heindex1].hecent)<l0[0] *.4) {  cout << "In find_overlap overlap heindex0 heindex1 " << heindex0 <<" " << heindex1 <<endl; exit(-1);}
						if (do_intersect(heid0,heid1)) {cout << "In find_overlap overlap triangles of  " << heindex0 <<" " << heindex1 <<endl; exit(-1);}
					}
				}

			}

		}
		
	}
	return 1;
}


int geometry::check_overlap_g(int vid0 ) {
	if (Nhe<22) { return 1;}
	//cout << "in _check_overlap_g" <<endl;
	int vindex0=vidtoindex[vid0];   
	double mindis=.1;
	for (vector<int>::iterator itv = v[vindex0].vneigh.begin(); itv != v[vindex0].vneigh.end(); itv++) 
	{	
		//cout<< " "
		int vindex1=vidtoindex[*itv];  //neighbor v
		if (vindex1!=-1){
			//if (veclen(v[vindex0].co,v[vindex1].co )<.1) { return -1;} // check vertex vertx overlap
			
			//if (is_vsurface(*itv)>0 && is_vsurface(vid0)>0 ) {mindis=.1; }
			//else { mindis=.2;}
			if (veclen(v[vindex1].co, v[vindex0].co)< l0[0] * mindis) { return -1; }
		

			for (vector<int>::iterator itneigh = v[vindex1].hein.begin(); itneigh != v[vindex1].hein.end(); itneigh++) //neigbor vertex hein s
			{
				int heid1=*itneigh;
				
				//if (is_surface(heid1)>0){ heid1=he[heidtoindex[heid1]].opid;} //why?
				int heindex1=heidtoindex[heid1];
				for (vector<int>::iterator it = v[vindex0].hein.begin(); it != v[vindex0].hein.end(); it++) //this vertex hein s
				{
					int heid0=*it;
					if (heid0!=he[heindex1].id && heid0!=he[heindex1].opid){
						//if (is_surface(heid0)>0){ heid0=he[heidtoindex[heid1]].opid;} //why?
						int heindex0=heidtoindex[heid0];
						if (veclen(v[vindex0].co,he[heindex1].hecent )<l0[0] *.2) { return -1;} // hecent -vertex overlap
						if (veclen(v[vindex1].co,he[heindex0].hecent )<l0[0] *.2) { return -1;} // hecent -vertex overlap
						float mindishe=.2;
						if (((is_surface(heid1)>0) || (is_surface(he[heidtoindex[heid1]].opid)>0) ) && ((is_surface(heid0)>0) || (is_surface(he[heidtoindex[heid0]].opid)>0))  && connectedH(heid1,heid0)>0){ mindishe=.25;}
						if (veclen(he[heindex0].hecent,he[heindex1].hecent)<l0[0] *mindishe) { return -1;} // hecent - hecent overlap
						if (do_intersect(heid0,heid1)) return -1;
					}
					
				}

			}

		}
		
	}
	return 1;
}


/*int geometry::check_overlap_vsurf(int vid0){
	if (Nhe <= 6)
		return 1;
	//double *vtemp = new double[3];
	//cout <<"check overlap " << endl <<endl;
	for (vector<HE>::iterator it = he.begin(); it != he.end(); it++)
	{
		//cout << " he id " << it->id<< endl;
		if ((it->vin!=vid0 && it->vout!=vid0)){ // prevent self check
			double mindis=.1;
			if (is_surface(it->id)>0 || is_surface(he[heidtoindex[it->id]].opid)>0 ) mindis=.25;
			
			subvec(it->hecent, v[vidtoindex[vid0]].co, vtemp);
			//cout << " vtemp[0] " << vtemp[0]<< " vtemp[1] " << vtemp[1]<< " vtemp[2] " << vtemp[2] <<endl;
			//cout << "norm is " <<  norm(vtemp) <<endl;
			if (norm(vtemp)>0 &&  norm(vtemp) < l0[0] * mindis)
			{ //cout << "OVERLAP HE HE"<<endl;
				delete[] vtemp;
				return -1;
			}
		}	
	}
	for (vector<VTX>::iterator it = v.begin(); it != v.end(); it++)
	{
		if (!(it->vid==vid0)){ // prevent self check
			double mindis=.3;
			if (is_vsurface(it->vid)>0  ) mindis=.1;
			subvec(it->co, v[vidtoindex[vid0]].co, vtemp);
			//cout << " vvvv " << it->vid<< "vtemp[0] " << vtemp[0]<< " vtemp[1] " << vtemp[1]<< " vtemp[2] " << vtemp[2] <<endl;
			//cout << "vvvv " << it->vid<< " norm is " <<  norm(vtemp) <<endl;
			if (norm(vtemp)>0.0 && norm(vtemp) < l0[0] * mindis)
			{ //cout << "OVERLAP HE VERTEX"<<endl;
				delete[] vtemp;
				return -1;
			}
		}
	}		

	//delete[] vtemp;

	//if (check_overlap_centerv(v[vidtoindex[vid0]].co)<0) {
	//	return -1;
	//}
	return 1;
}*/



int geometry::check_overlap_centerv(double *newv) /* checks this v point with all other hecenters and vertices*/
{ 
	if (Nhe <= 6)
		return 1;
	
	for (vector<HE>::iterator it = he.begin(); it != he.end(); it++)
	{
		double mindisvh=.2;
		if (veclen(it->hecent, newv) < l0[0] * mindisvh)	{ return -1;}
	}
	for (vector<VTX>::iterator it = v.begin(); it != v.end(); it++)
	{
		double mindis=.1;
		if (veclen(it->co, newv)< l0[0] * mindis)	{ return -1;}
	}

	return 1;
}

int geometry::check_overlap_centerh(double *newcenter) /* checks this he center point with all other hecenters*/
{ 
	if (Nhe <= 6)	return 1;
	
	for (vector<HE>::iterator it = he.begin(); it != he.end(); it++)
	{
		double mindis=.2;
		if (is_surface(it->id)>0) mindis=.25;
		if (veclen(it->hecent ,newcenter)< l0[0] * mindis) 	{ return -1;}

		double mindisvh=.2;
		if (veclen(v[vidtoindex[it->vin]].co,newcenter ) < l0[0] * mindisvh) { return -1;}

	}
	
	return 1;
}

/*int geometry::check_overlap_vtx(int vid0)
{
	if (Nhe <= 6)
		return 1;
	double *vtemp = new double[3];
	//cout <<"check overlap " << endl <<endl;
	for (vector<HE>::iterator it = he.begin(); it != he.end(); it++)
	{
		//if (!(it->id==heid0) && !(it->opid==heid0)) {//cout << " he id " << it->id<< endl;
		//cout << "it->id " << it->id <<endl;
		//cout << "heid0 " << heid0 << endl;
		subvec(it->hecent, v[vidtoindex[vid0]].co, vtemp);
		//cout <<"it->hecents "<< it->hecent[0] << " " << it->hecent[1] << " " << it->hecent[2] << " " <<endl;
		//cout << "he[heidtoindex[heid0)].hecent " << he[heidtoindex[heid0)].hecent[0] << " " << he[heidtoindex[heid0)].hecent[1] << " "  << he[heidtoindex[heid0)].hecent[2] << " " <<endl;
		//cout << " vtemp[0] " << vtemp[0]<< " vtemp[1] " << vtemp[1]<< " vtemp[2] " << vtemp[2] <<endl;
		//cout << "norm is " <<  norm(vtemp) <<endl;
		if (norm(vtemp) < l0[0] * .1)
		{ //cout << " distance between "<< it->id << " and " << heid0 << " is " << norm(vtemp) <<endl;
			delete[] vtemp;
			return -1;
		}
	}

	//cout << " vertices " << v.size();
	//exit(-1);
	//int vind=get_vindex(he[heidtoindex[heid0)].vout);
	for (vector<VTX>::iterator it = v.begin(); it != v.end(); ++it)
	{
		if (it->vid != vid0)
		{
			subvec(it->co, v[vidtoindex[vid0]].co, vtemp);
			//cout << " vvvv " << it->vid<< "vtemp[0] " << vtemp[0]<< " vtemp[1] " << vtemp[1]<< " vtemp[2] " << vtemp[2] <<endl;
			//cout << "vvvv " << it->vid<< " norm is " <<  norm(vtemp) <<endl;
			if (norm(vtemp) < l0[0] * .1)
			{
				delete[] vtemp;
				return -1;
			}
		}
	}

	//cout << "no overlap !\n";	make

	delete[] vtemp;
	return 1;
}*/

/*ENERGY HELPER FUNCTIONS */

double geometry::find_gbb(int etypenew, int etypenextofnew, int etypeprevifnew)
{
	double e1 = find_dg(etypenew, etypenextofnew, 0);
	double e2 = find_dg(etypenextofnew, etypeprevifnew, 0);
	double e3 = find_dg(etypeprevifnew, etypenew, 0);
	//cout << "e1 is " << e1 << " e2 is " << e2 << " e3 is " <<e3 <<endl;
	return e1 + e2 + e3;
}

double geometry::find_dg(int type, int typenext, bool drug)
{

	/*/if ((typenext == 1 || typenext == 2) and (drug != 0))
	{
		cout << "why drug on 1 2 " << endl;
		exit(-1);
	}*/

	double bindg = gb[type][typenext] + drug*(gdrug[type][typenext]-mudrug);

	/*if (((typenext == 0) || (typenext == 3)) && ((type == 0) || (type == 3))) {
		bindg +=  2*(gdrug-mudrug) * drug;
	}
	else if (((typenext == 0) || (typenext == 3)) && ((type == 1) || (type == 2))) {
		bindg += (gdrug-mudrug) * drug;
	}
	else {
		bindg += (-mudrug) * drug;	
	}*/
	/*if ((type==0 || type==3) &&(typenext==0 || typenext==3)) {
		bindg+= (gdrug-mudrug)*drug;
	}
	else if ((type==2) && (typenext==0 || typenext==3)) {
		bindg+=(gdrug-mudrug)*drug;
	}*/

	/*if (type==3 && typenext==3 ) { //DC-DC //weak GB
		return drug*(gdrug-mudrug)+gb[3][3];
	} 
	else if (type==1 && typenext==2) { //BA-AB  // GB+dg //nodrug
		return gb+dg*gb[1][2];
	}	 
	else if ((type==0 && typenext==1)  )  { //T4 CD-BA  D->B   // GB+dg //nodrug
		return gb +dg*01;
	}  
	else if (type==2 && typenext==0) { //T4  AB-CD   B->C  //weak GB
		return drug*(gdrug-mudrug)+ gb ;
	}

	else if (  (type==3 && typenext==1))  { // T3  DC-BA   C->B  // GB+dg //nodrug
		
		return gb;
	}
	else if ((type==2 && typenext==3) ) { 	// T3 AB-(CC) DC  B->C(D) //weak GB 
		return drug*(gdrug-mudrug)+gb ;
	}*/

	//else if ( (type==0 && typenext==0)){//} || (type==0 && typenext==3) || (type==3 && typenext==0)) { //hexagonal sheet //DC-DC
	//	return 2*drug*(gdrug-mudrug)+ gb-5*dg*gb;
	//}
	//else if ((type=0 || type==0 )){
	//return drug*(gdrug-mudrug)+ gb-2*dg*gb;
	//}
	//} && (typenext==0 || typenext==3) ){  //very weak but dug binds both
	//cout << " should not exist " << "type " << type << "nexttype " << typenext<<endl;
	//exit(-1);
	//return 2*drug*(gdrug-mudrug);//gb-5*dg*gb;
	//return drug*(gdrug-mudrug)+ gb-1*dg*gb;
	//}
	//else if ((typenext==0 || typenext==3) ){
	//cout << " should not exist " << "type " << type << "nexttype " << typenext<<endl;
	//exit(-1);
	//return drug*(gdrug-mudrug);//gb-5*dg*gb;
	//}
	//else {
	//	return 0;
	//}
	return (bindg);
}

double geometry::stretch_energy(int heindex0)
{
	int et = he[heindex0].type;

	//cout<<  "epsilon[et]" << epsilon[et]  <<endl;
	//cout<<  "(he[heindex0].l " << he[heindex0].l <<endl ;
	//cout <<"stretch energy is " <<  0.5 * epsilon[et] * (he[heindex0].l - l0[et]) * (he[heindex0].l - l0[et]) <<endl;
	return 0.5 * epsilon[et] * (he[heindex0].l - l0[et]) * (he[heindex0].l - l0[et]);
}
int geometry::check_bind_wedge(int heid0){
	int heindex0=heidtoindex[heid0];
	get_normal(he[heindex0].id);
	get_normal(he[heindex0].opid);
	int opindex0 = heidtoindex[he[heindex0].opid];
	int nextindex0 = heidtoindex[he[heindex0].nextid];
	//int previndex0 = heidtoindex[he[heindex0].previd];
	int opnextindex0 = heidtoindex[he[opindex0].nextid];

	
	double *tempvec = new double[3];
	subvec(v[vidtoindex[he[nextindex0].vout]].co, v[vidtoindex[he[opnextindex0].vout]].co, tempvec);
	
	//double ndot = dot(he[opindex0].n, he[heindex0].n);
	if (dot(tempvec, he[heindex0].n) > 0) {
		delete[] tempvec;
		return -1;
	}
	delete[] tempvec;
	return 1;
}

double geometry::bend_energy(int heindex0)
{
	//cout << " \n\n     normal " << he[heindex0].n[0] << " " << he[heindex0].n[1] << " " << he[heindex0].n[2] << " " <<endl ;
	//cout << "\n in bend_energy edge index is "  << heindex0<< endl;
	int opindex0 = heidtoindex[he[heindex0].opid];
	if (he[heindex0].nextid == -1 && he[heindex0].previd == -1)
	{
		return 0;
	}
	if (he[opindex0].nextid == -1 && he[opindex0].previd == -1)
	{
		return 0;
	}
	//if (is_surface(he[heindex0].id)>0 || is_surface(he[heindex0].opid)>0)  { return 0;}
	//cout << "in bend_energy 1111111" << endl;
	double bendE = 0;
	double ndot;
	//cout << "first BEND E IS " << bendE <<endl;
	//int et=he[heindex0].type; // type 0 is always no curvature
	//cout << " EDGE TYPE is " << he[heindex0].type;
	//if (he[heindex0].nextid==-1 || he[heindex0].previd==-1) { return 0; }

	get_normal(he[heindex0].id);
	get_normal(he[heindex0].opid);
	int previndex0 = -1;
	int nextindex0 = -1;
	int nexttype = -1;
	int prevtype = -1;
	int etype = he[heindex0].type;
	int vid0=-1;
	int opvid0=-1;
	if (he[heindex0].nextid!=-1) { 
		nextindex0 = heidtoindex[he[heindex0].nextid];
		nexttype = he[nextindex0].type;
		vid0= he[nextindex0].vout;


	}

	if (he[heindex0].previd!=-1) { 
		previndex0 = heidtoindex[he[heindex0].previd];
		prevtype = he[previndex0].type;
		vid0= he[previndex0].vin;
	}
	
	int opnextindex0 = -1;
	int opprevindex0 = -1;
	int opnexttype = -1;
	int opprevtype = -1;
	int opetype = he[opindex0].type;

	if (he[opindex0].nextid!=-1) {
		opnextindex0 = heidtoindex[he[opindex0].nextid];
		opnexttype = he[opnextindex0].type;
		opvid0=he[opnextindex0].vout;
	}
	if (he[opindex0].previd!=-1) {
		opprevindex0 = heidtoindex[he[opindex0].previd];
		opprevtype = he[opprevindex0].type;
		opvid0= he[opprevindex0].vin;
	}
	
	
	
	//for (vector<HE>::iterator it = g.he.begin() ; it != g.he.end(); it++) {
	//cout << "EDGE index " <<heindex0 <<" id " << he[heindex0].id << " vin " << vidtoindex[he[heindex0].vin] << "vout " << vidtoindex[he[heindex0].vout] <<endl;
	//cout << "      opid " << he[heindex0].opid << " nextid " << he[heindex0].nextid << " previd " << he[heindex0].previd <<endl;
	//cout << "      hevec " << he[heindex0].hevec[0] << " " << he[heindex0].hevec[1]  << " "  << he[heindex0].hevec[2] <<" l is " << he[heindex0].l << endl;
	//cout << "      normal " << he[heindex0].n[0] << " " << he[heindex0].n[1] << " " << he[heindex0].n[2] << " " <<endl ;
	//cout << "edge type " << he[heindex0].type << " op type " << he[heidtoindex[he[heindex0].opid)].type <<endl;
	//cout << " next type" << he[heidtoindex[he[heindex0].nextid)].type << " prev type" << he[heidtoindex[he[heindex0].previd)].type <<endl<<endl;
	//}
	
	// test for convex
	double *tempvec = new double[3];
	
	subvec(v[vidtoindex[vid0]].co, v[vidtoindex[opvid0]].co, tempvec);
	//cout << " normal of this edge is " << he[heindex0].n[0] <<" "<< he[heindex0].n[1] << " "<< he[heindex0].n[2] <<endl;
	//cout << " normal of op edge is " << he[opindex0].n[0] <<" "<< he[opindex0].n[1] << " "<< he[opindex0].n[2] <<endl;

	ndot = dot(he[opindex0].n, he[heindex0].n);
	//cout << " after dot normal of this edge is " << he[heindex0].n[0] <<" "<< he[heindex0].n[1] << " "<< he[heindex0].n[2] <<endl;
	//cout << " normal of op edge is " << he[opindex0].n[0] <<" "<< he[opindex0].n[1] << " "<< he[opindex0].n[2] <<endl;
	double theta;
	if (ndot >= 1)
	{
		theta = 0;
	}
	else
	{
		theta = acos(ndot);
	}
	// prevents numerical precision errors
	//double theta = acos(ndot);
	//cout << "ndot is "  << ndot << " theta is " << theta << endl;
	//if (dot(tempvec,he[heindex0].n)<0) {
	//	theta*=-1;
	//}

	int angle0 = -1; //theta0[0]/2;
	//CD-BA-AB :: DC-BA-AB in T3
	if (((etype == 0 || etype == 3) && nexttype == 1 && prevtype == 2) && ((opetype == 3 || opetype == 0) && opnexttype == 1 && opprevtype == 2))
	{
		angle0 = 0;
	}
	// BA-AB-CD :: AB-CD-AB //T3 and T4 5fold

	else if ((etype == 1 && nexttype == 2 && (prevtype == 0 || prevtype == 3)) && (opetype == 2 && (opnexttype == 0 || opnexttype == 3) && opprevtype == 1))
	{
		angle0 = 1;
		//cout <<"/ BA-AB-CD :: AB-CD-BA" <<endl;
	}
	//AB-AB-CD :: CD-AB-AB

	else if ((etype == 2 && (nexttype == 0 || nexttype == 3) && prevtype == 1) && (opetype == 1 && opnexttype == 2 && (opprevtype == 0 || opprevtype == 3)))
	{
		angle0 = 1;
		//cout <<"/ BA-AB-CD :: AB-CD-BA" <<endl;
	}

	else if ((etype == 0 && nexttype == 1 && prevtype == 2) && (opetype == 3 && opnexttype == 3 && opprevtype == 3))
	{ //AB-AB-CD :: CD-CD-CD
		angle0 = 0;
		//cout <<"//CD_BA_AB :: CD-CD-CD"<<endl;
	}
	else if ((etype == 3 && nexttype == 3 && prevtype == 3) && (opetype == 0 && opnexttype == 1 && opprevtype == 2))
	{ // CD-CD-CD :: CD-AB-AB
		angle0 = 0;
		//cout <<"// CD-CD-CD :: CD-BA-AB  "<<endl;
	}

	// are these necessary???->?????????
	else if ((etype == 3 && nexttype == 1 && prevtype == 2) && (opetype == 0 && opnexttype == 0 && opprevtype == 0))
	{ //AB-AB-CD :: CD-CD-CD
		angle0 = 0;
		//cout <<"//CD_BA_AB :: CD-CD-CD"<<endl;
	}
	else if ((etype == 0 && nexttype == 0 && prevtype == 0) && (opetype == 3 && opnexttype == 1 && opprevtype == 2))
	{ // CD-CD-CD :: CD-AB-AB
		angle0 = 0;
		//cout <<"// CD-CD-CD :: CD-BA-AB  "<<endl;
	}

	
	else if (((etype == 0 || etype == 3) && (nexttype == 0 || nexttype == 3) && (prevtype == 0 || prevtype == 3) ) && ( (opetype == 3 || opetype == 0)  &&  (opnexttype == 0 || opnexttype == 3)  && (opprevtype == 0 || opprevtype == 3)))
	{
		angle0 = 2;
	}

	else
	{
		//cout << " etype " <<  etype << " optype " << opetype <<endl;
		angle0 = 3;
		//cout << "//all other   "<<endl;
	}
	//else if ((etype==1 && nexttype==0 && prevtype==1) &&	(opetype==1 && opnexttype==0 && opprevtype==1)) {//  CD-CD-CD :: CD-CD-CD
	//	angle0=-theta0[0];//trying
	//}*/

	//bendE = kappa[et] * (1-ndot);
	//double theta = acos(ndot);
	//cout << " preferred angle is " << angle0 <<endl;
	//cout << "kappa[0] is " << kappa[0] << endl;

	/*if ((etype == 0 && optype == 3)  || ( etype == 3 && optype == 0 )) { angle0=0;}
	else if ((etype==1 && optype==2) || (etype==2 && optype==1)) {angle0=1;}
	else { cout << "types don't match" <<endl;}*/

	//bendE = .5*kappa[angle0] * (theta - theta0[angle0])*(theta - theta0[angle0]);

	/*if (dot(tempvec, he[heindex0].n) > 0)
	{
		theta=PI+theta;
		//cout <<" convex"<<endl;
		
	}
	bendE = .5*kappa[angle0] * (theta - theta0[angle0])*(theta - theta0[angle0]);*/

	if (theta0[angle0] > 0)
	{
		bendE = kappa[angle0] * (1 - cos(theta - theta0[angle0]));
		

		if (dot(tempvec, he[heindex0].n) > 0)
		{
			//cout <<"convex" <<endl;
			//theta=PI+theta;
			//bendE = kappa[angle0] * (theta - theta0[angle0])*(theta - theta0[angle0]);
			bendE = kappa[angle0] * (1 - cos(-theta - theta0[angle0]));

			if (theta<theta0[angle0]) {bendE*=1000; }
			else {bendE*=100;  }
			//bendE*=1000;
		}
		//bendE = kappa[angle0] * (theta - theta0[angle0])*(theta - theta0[angle0]);
	}
	else
	{
		bendE = kappa[angle0] * (1 - ndot);
	}
	 // cout << " convex"<<endl;}
	  /** ???? correct for convex**/
	  //cout << " angle0 is "<< angle0<<endl;
	  //cout << "BEND E IS " << bendE <<endl;
	if (bendE < 0)
	{
		cout << "BEND E IS " << bendE << endl;
		exit(-1);
	}

	delete[] tempvec;
	
	return bendE;
}

double geometry::dimer_bend_energy(int heindex0)
{
	if (he[heindex0].nextid == -1)
	{
		return 0;
	}

	double DbendE;
	int nextindex = heidtoindex[he[heindex0].nextid];
	//cout << "nextindex: " << nextindex<<endl;
	int opindex = heidtoindex[he[heindex0].opid];
	//cout << "opindex: " << opindex<<endl;
	update_half_edge(he[opindex].id);
	update_half_edge(he[nextindex].id);
	double ndot = (dot(he[opindex].hevec, he[nextindex].hevec) / (he[opindex].l * he[nextindex].l));
	//cout << "ndot "<<ndot<<endl;
	//if (ndot<-.9) { cout <<" not accepted triangle edge he[opindex].hevec "<<he[opindex].hevec[0] <<endl; exit(-1);}
	if (ndot > 1)
	{
		ndot = 1;
	}
	double phi = acos(ndot);
	//cout <<phi<< " is phi "<<endl;
	int phitype = -1;
	int etype = he[heindex0].type;
	int nexttype = he[nextindex].type;
	//cout << " type " << it->type << "nexttype" << nexttype <<endl;
	if ((etype == 0 && nexttype == 0) || (etype == 3 && nexttype == 3) || (etype == 3 && nexttype == 0) || (etype == 0 && nexttype == 3))
	{ //(cd-cd)  // T4 and hexamer sheet
		phitype = 0;
		//cout << "PHI 0 //CD -CD" <<endl;
	}
	else if ((etype == 1 && nexttype == 2))
	{ //(ab-ab)
		phitype = 1;
		//cout << "PHI 1 //BA-AB " <<endl;
	}
	else if ((etype == 0 && nexttype == 1) || (etype == 3 && nexttype == 1))
	{ //( CD-BA) (DC-BA) (<60 drug doesnt binds)
		phitype = 2;
		//cout << "PHI 2 //DC-BA" <<endl;
	}
	else if ((etype == 2 && nexttype == 0) || (etype == 2 && nexttype == 3))
	{ //AB-DC  and T3  (=60 drug binds)
		phitype = 2;
	}

	else
	{
		//cout << "PHI 1 //all other" <<endl;
		//cout << "phitype other  etype is " << etype << " nexttype is " << nexttype <<endl;
		phitype = 0;
		//exit(-1);
	}

	if (he[nextindex].din == 1)
		phitype = 3;
	//cout <<phitype<< " is phitype "<<endl;
	//if phi0[phitype]==0
	
	double kPhi = kappaPhi[phitype];
	//if (he[nextindex].din==1) kPhi*=10;
	if (phi<phi0[phitype]) {
		kPhi/=5;
	}
	//DbendE = kPhi * (1 - cos(phi - phi0[phitype])); // }
	DbendE = .5*kPhi * (phi - phi0[phitype])* (phi - phi0[phitype]);
	//if (phi<phi0[phitype]) { DbendE*=.5; }
	//else {DbendE= kPhi*.1 * (1 - cos(phi-phi0[phitype]));}
	if (DbendE < 0)
	{
		cout << "DBEND E IS " << DbendE << endl;
		exit(-1);
	}
	return DbendE;
}
double geometry::compute_bind_energy()
{
	double tot_bind = 0;
	for (vector<HE>::iterator it = he.begin(); it != he.end(); ++it)
	{
		if (it->nextid != -1)
		{
			int nextid = heidtoindex[it->nextid];
			tot_bind += find_dg(it->type, he[nextid].type, he[nextid].din);
		}
	}
	return tot_bind;
}
double geometry::compute_energy()
{
	double tot_eng = 0;
	//update_normals();
	//for (int i=0; i<he.size(); i++) {
	for (vector<HE>::iterator it = he.begin(); it != he.end(); ++it)
	{
		int heindex0 = heidtoindex[it->id];
		// add the contribution from stretching the bonds
		tot_eng += stretch_energy(heindex0) * .5;
		//cout << " edge i "  << i << " stretch E " << stretch_energy(i) <<endl;
		// compute the normal vectors for the triangles along edge i
		double temp = bend_energy(heindex0) * .5;
		//cout << " edge i "  << i << " bend E " << bend_energy(i) <<endl;
		//cout << "i is " << i  << " + stretch total energy is" << tot_eng <<endl;
		tot_eng += temp;
		//cout << "i is " << i  << " + dbend total energy is" << tot_eng <<endl;
		double temp2 = dimer_bend_energy(heindex0);
		tot_eng += temp2;
	}
	return tot_eng; // double count
}

double geometry::monomer_energy(int heid0)
{
	//heid0 is on the surface
	int heindex0 = heidtoindex[heid0];
	double etot = stretch_energy(heindex0);
	int heopindex0 = heidtoindex[he[heindex0].opid];
	int nextindex0 = heidtoindex[he[heopindex0].nextid];
	int previndex0 = heidtoindex[he[heopindex0].previd];
	if (nextindex0 == -1)
	{
		cout << " no other side !!!!!" << endl;
		exit(-1);
		return etot;
	}
	etot += bend_energy(nextindex0);
	if (previndex0 == -1)
	{
		cout << " no other side !!!!!" << endl;
		exit(-1);
		return etot;
	}
	//etot+=bend_energy(nextindex0);
	etot += bend_energy(previndex0);

	return (etot);
}

double geometry::dimer_energy(int heid0, int heid1)
{
	//cout << "in dimer energy" << endl;
	int heindex0 = heidtoindex[heid0];
	int heindex1 = heidtoindex[heid1];
	double etot = stretch_energy(heindex0);
	etot += stretch_energy(heindex1);
	int heid2 = -1;
	//int heopindex0 = heidtoindex[he[heindex0].opid);
	//int heopindex1 = heidtoindex[he[heindex1].opid);
	if (he[heindex0].vin == he[heindex1].vout)
	{
		heid2 = he[heindex0].nextid;
	}
	else if (he[heindex0].vout == he[heindex1].vin)
	{
		heid2 = he[heindex0].previd;
	}
	if (heid2 == -1)
	{
		cout << " in dimer energy !!!!!! no bend " << endl;
		exit(-1);
		return etot;
	}
	etot += bend_energy(heidtoindex[heid2]);
	return (etot);
}

double geometry::vertex_energy(int vid0)
{
	double v_eng = 0;
	int vindex0=vidtoindex[vid0];
	for (vector<int>::iterator ithe = v[vindex0].hein.begin(); ithe != v[vindex0].hein.end(); ++ithe)
    { 
        int heindex = heidtoindex[*ithe];
                
        v_eng += stretch_energy(heindex);
        v_eng += bend_energy(heindex);
        v_eng += dimer_bend_energy(heindex);

        if (he[heindex].previd != -1)
        {
            
            int previndex = heidtoindex[he[heindex].previd];
            v_eng +=   dimer_bend_energy(previndex);

			//if (is_surface(he[heindex].previd)<0) {
			v_eng += bend_energy(previndex);
			//}
        }
        if (he[heindex].nextid!= -1)
        {
            int nextindex= heidtoindex[he[heindex].nextid];
            v_eng +=  dimer_bend_energy(nextindex);
        }
		//if (is_surface(he[heindex].opid)>0) // ?????
        //{
		//	int opindex= heidtoindex[he[heindex].nextid];
          //  v_eng +=  dimer_bend_energy(opindex);
		//}

    }
	
	
	return v_eng;
}


int geometry::get_fusion_vid(int vidsurf0){

	if (is_bond_vsurface(vidsurf0)>0)  { return -1;}
	//cout << "in get_fusion vid for vid0 " << vidsurf0 <<endl;   
    int vsurfindex0=vidtoindex[vidsurf0];
	//cout << "vindex is " << vsurfindex0 <<endl;

	//TESTNEW WAY
	if  (v[vsurfindex0].vneigh.size()==0) { return -1;}
	
	int vhecount=0; /* prevent 7 fold */
	for (vector<int>::iterator it = v[vsurfindex0].hein.begin(); it != v[vsurfindex0].hein.end(); ++it)
	{	
		vhecount++;
	}
	

	
	int hesurfindextemp=-1;
	int vsurfindextemp=vsurfindex0;

	for (int i=0; i<3; i++){
		//cout << "hesurfinid " << v[vsurfindextemp].hesurfinid <<endl;
		if (v[vsurfindextemp].hesurfinid==-1) { cout <<"something is wrong for this ! " << endl; exit(-1);}
		hesurfindextemp=heidtoindex[v[vsurfindextemp].hesurfinid];
		//cout << "hesurfindextemp " << hesurfindextemp << " id  " << he[hesurfindextemp].id <<endl;
		vsurfindextemp=vidtoindex[he[hesurfindextemp].vin];
		//cout << "vsurfindextemp " << vsurfindextemp << " vid " << v[vsurfindextemp].vid  <<endl;
	}
	for (vector<int>::iterator it = v[vsurfindextemp].hein.begin(); it != v[vsurfindextemp].hein.end(); ++it){	
		vhecount++;
	}
	if (vhecount>6) { return -1;}


	if (is_bond_vsurface(v[vsurfindextemp].vid)>0) { return -1;}
	if (connected(v[vsurfindextemp].vid, vidsurf0) >0) {return -1;}
	double d= veclen(v[vsurfindex0].co, v[vsurfindextemp].co);
    if (vsurfindextemp==vsurfindex0) { cout << " something went wrong " <<endl; exit(-1);}
    //cout << "vertices distance is " << d << endl;
    
    /* temp restrict numer of bonds in fusion*/
	if (d < (xi) ){
        return v[vsurfindextemp].vid;
    }

	return -1;
}

void geometry::update_fusion_pairs(){

	//cout <<"in update_fusionP_pairs" <<endl;
	if (fusionv.size()>0) {
		fusionv.clear();
	}
	for (vector<int>::iterator it = surfv.begin(); it != surfv.end(); ++it)
	{
		v[vidtoindex[*it]].fusion_vid=-1;

	}
	//cout <<"in update_fusionP_pairs surfv.size " << surfv.size() <<endl;

	for (vector<int>::iterator it = surfv.begin(); it != surfv.end(); ++it)
	{
		
		if (v[vidtoindex[*it]].fusion_vid==-1) {
			int tempfusionvid=get_fusion_vid(*it);
			
			
			if (tempfusionvid!=-1 && (v[vidtoindex[tempfusionvid]].fusion_vid)==-1) {
				//cout << "vid is" << *it <<endl;
				//cout << "tempfusionvid fusion pair!" << tempfusionvid <<endl;
				v[vidtoindex[*it]].fusion_vid=tempfusionvid;
				v[vidtoindex[tempfusionvid]].fusion_vid=*it;
				fusionv.push_back(*it);
				fusionv.push_back(tempfusionvid);
			}
		}

	}

}

void geometry::save_vtx(int vid0 , VTX *tempvtx)
{

	int vindex0=vidtoindex[vid0];
	
	tempvtx->co[0] = v[vindex0].co[0];
    tempvtx->co[1] = v[vindex0].co[1];
    tempvtx->co[2] = v[vindex0].co[2];
	

	for (vector<int>::iterator ithe = v[vindex0].hein.begin(); ithe != v[vindex0].hein.end(); ++ithe)
    {
		tempvtx->hein.push_back(*ithe); 
	}
	for (vector<int>::iterator itv = v[vindex0].vneigh.begin(); itv != v[vindex0].vneigh.end(); ++itv)
    {
		tempvtx->vneigh.push_back(*itv); 
	}
	tempvtx->fusion_vid=v[vindex0].fusion_vid;
	tempvtx->hein=v[vindex0].hein;
	tempvtx->hesurfinid=v[vindex0].hesurfinid;
	tempvtx->vid=-1; //decide on this later
}

void subvec(double *vinit, double *vfin, double *vec)
{
	vec[0] = vfin[0] - vinit[0];
	vec[1] = vfin[1] - vinit[1];
	vec[2] = vfin[2] - vinit[2];
}

void addvec(double *vinit, double *vfin, double *vec)
{
	vec[0] = vinit[0] + vfin[0];
	vec[1] = vinit[1] + vfin[1];
	vec[2] = vinit[2] + vfin[2];
}

void multvec(double *vinit, double scalar, double *vec)
{
	vec[0] = vinit[0] * scalar;
	vec[1] = vinit[1] * scalar;
	vec[2] = vinit[2] * scalar;
}

void centvec(double *vinit, double *vfin, double *vec)
 {
	vec[0] = (vinit[0] + vfin[0])/2;
	vec[1] = (vinit[1] + vfin[1])/2;
	vec[2] = (vinit[2] + vfin[2])/2;
 }

double veclen(double *vinit, double *vfin){
	double *vec=new double[3];
	vec[0] = vfin[0] - vinit[0];
	vec[1] = vfin[1] - vinit[1];
	vec[2] = vfin[2] - vinit[2];
	double vlen=norm(vec);
	delete[] vec;
	return(vlen);

}

double norm(double *v)
{
	return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

double dot(double *v1, double *v2)
{
	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

void cross(double *v1, double *v2, double *res)
{
	res[0] = v1[1] * v2[2] - v1[2] * v2[1];
	res[1] = v1[2] * v2[0] - v1[0] * v2[2];
	res[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

void randvec(double *v, gsl_rng *r)
{
	double psi1 = gsl_rng_uniform(r);
	double psi2 = gsl_rng_uniform(r);
	double theta = 2 * PI * psi2;
	double phi = acos(1 - 2 * psi1);
	v[0] = sin(phi) * cos(theta);
	v[1] = sin(phi) * sin(theta);
	v[2] = cos(phi);
	//cout <<"v[0] is " << v[0] <<endl;
	//cout << " norm v in randvec " <<norm(v)<<endl;
	
	//double vlen=norm(v);
	//v[0]/=vlen;
	//cout <<"normalized v[0] is  " << v[0] <<endl;
	//v[1]/=vlen;
	//v[2]/=vlen;
}

void dump_lammps_data_file(geometry &g, int time0)
{
	char filename[80];
	float box = 3.0;
	g.Nv5 = 0;
	sprintf(filename, "lammps_%07d.dat", time0);
	FILE *f;
	f = fopen(filename, "w");
	//fprintf(f,"@<TRIPOS>MOLECULE\n");

	fprintf(f, "LAMMPSDescription-Generated by Shelfrich at time_step=%d\n", time0);
	fprintf(f, "\n%d atoms", g.Nv + g.Nd);
	fprintf(f, "\n%d bonds", g.Nhe / 2);
	//fprintf(f,"\n%d bonds",g.Nhe/2+g.Nsurf);
	fprintf(f, "\n");
	fprintf(f, "\n4 atom types");
	fprintf(f, "\n4 bond types");
	fprintf(f, "\n");
	fprintf(f, "\n%8.3f %8.3f xlo xhi", -box, box);
	fprintf(f, "\n%8.3f %8.3f ylo yhi", -box, box);
	fprintf(f, "\n%8.3f %8.3f zlo zhi", -box, box);
	fprintf(f, "\n");
	fprintf(f, "\nAtoms");
	fprintf(f, "\n");
	//cout << "here in dump 000"<<endl;
	for (vector<VTX>::iterator it = g.v.begin(); it != g.v.end(); ++it)
	{
		//cout <<" it->co[0]"<< it->co[0]<< endl;
		//exit(-1);
		if (g.is_bond_vsurface(it->vid) > 0)
		{
			fprintf(f, "\n%li 3 %8.3f %8.3f %8.3f", distance(g.v.begin(), it) + 1, it->co[0], it->co[1], it->co[2]);
		}

		else if (it->hein.size() == 5 && g.is_vsurface(it->vid) < 0)
		{
			fprintf(f, "\n%li 1 %8.3f %8.3f %8.3f", distance(g.v.begin(), it) + 1, it->co[0], it->co[1], it->co[2]);
			g.Nv5++;
		}
		else
		{
			fprintf(f, "\n%li 2 %8.3f %8.3f %8.3f", distance(g.v.begin(), it) + 1, it->co[0], it->co[1], it->co[2]);
		}
		//fprintf(stderr,"\n%li 1 %8.3f %8.3f %8.3f", distance(g.v.begin(),it)+1 ,it->co[0], it->co[1], it->co[2]);
	}
	int counter = g.Nv + 1;
	for (vector<HE>::iterator it = g.he.begin(); it != g.he.end(); ++it)
	{
		if (it->din == 1)
		{
			int vindex = g.vidtoindex[it->vin];
			double x0 = 0;
			double x1 = 0;
			double x2 = 0;
			if (it->previd != -1)
			{
				int preindex = g.heidtoindex[it->previd];
				x0 = -.1 * g.he[preindex].hevec[0];
				x1 = -.1 * g.he[preindex].hevec[1];
				x2 = -.1 * g.he[preindex].hevec[2];
			}
			fprintf(f, "\n%d 4 %8.3f %8.3f %8.3f", counter++, x0 + (g.v[vindex]).co[0] + .15 * (it->hevec[0]), x1 + g.v[vindex].co[1] + .15 * (it->hevec[1]), x2 + g.v[vindex].co[2] + .15 * (it->hevec[2]));
		}
	}
	
	fprintf(f, "\n");
	fprintf(f, "\nBonds");
	fprintf(f, "\n");
	//cout <<" "<<endl;
	//cout << "here in dump 222"<<endl;
	for (vector<HE>::iterator it = g.he.begin(); it != g.he.end(); ++it)
	{
		//cout << "edge " <<distance(g.he.begin(),it)+1 <<" " << it->id << "vin vid" << g.vidtoindex[it->vin] << "vout vid" << g.vidtoindex[it->vout] <<endl;
		if (it->vin == -1 || it->vout == -1 || g.vidtoindex[it->vin] == -1 || g.vidtoindex[it->vout] == -1)
		{
			cout << " dump_data ! error in vin vout of edge " << it->id << endl;
			exit(-1);
		}
		int btype = it->type + 1;
		//if ( it->type==2) {  btype=2 ;}
		//if ( it->type==3) {  btype=1 ;}
		if ((it->type == 1) || (it->type == 0))
		{
			fprintf(f, "\n%li %d %d %d", distance(g.he.begin(), it) + 1, btype, g.vidtoindex[it->vin] + 1, g.vidtoindex[it->vout] + 1);
		}
		//}
		//if ( g.is_surface(it->id)<0 && g.is_surface(it->opid)<0) {
		//fprintf(f, "\n%li 1 %d %d",distance(g.he.begin(),it)+1 , g.vidtoindex[it->vin]+1, g.vidtoindex[it->vout]+1);
		//fprintf(stderr, "\n%li 1 %d %d",distance(g.he.begin(),it)+1 , g.vidtoindex[it->vin]+1, g.vidtoindex[it->vout]+1);
		//}
		//else if ( g.is_surface(it->id)<0 && g.is_surface(it->opid)>0) {
		//fprintf(f, "\n%li 2 %d %d",distance(g.he.begin(),it)+1 , g.vidtoindex[it->vin]+1, g.vidtoindex[it->vout]+1);
		//fprintf(stderr, "\n%li 2 %d %d",distance(g.he.begin(),it)+1 , g.vidtoindex[it->vin]+1, g.vidtoindex[it->vout]+1);
		//}
		//else {
		//fprintf(f, "\n%li 2 %d %d",distance(g.he.begin(),it)+1 , g.vidtoindex[it->vin]+1, g.vidtoindex[it->vout]+1);
		//fprintf(stderr, "\n%li 2 %d %d",distance(g.he.begin(),it) , g.vidtoindex[it->vin], g.vidtoindex[it->vout]);
		//}
	}

	//exit(-1);
	fprintf(f, "\n");
	fclose(f);
}

void rotatevec(double *vec, double *axis, double angle, double *vec2)
{
	double vx, vy, vz, ex, ey, ez, nm; //,vxr,vyr,vzr;
	vx = vec[0];
	vy = vec[1];
	vz = vec[2];
	//finding normalized vector
	ex = axis[0];
	ey = axis[1];
	ez = axis[2];

	//cout << "hevec is" << he[heindex0].hevec[0] << " " << he[heindex0].hevec[1] << " " << he[heindex0].hevec[2] <<endl;
	nm = norm(axis);
	ex /= nm;
	ey /= nm;
	ez /= nm;
	double ct = cos(angle);
	double st = sin(angle);
	double mct = 1 - cos(angle);
	//rotating
	vec2[0] = vx * (ct + ex * ex * mct) + vy * (ex * ey * mct - ez * st) + vz * (ex * ez * mct + ey * st);
	vec2[1] = vx * (ey * ex * mct + ez * st) + vy * (ct + ey * ey * mct) + vz * (ey * ez * mct - ex * st);
	vec2[2] = vx * (ez * ex * mct - ey * st) + vy * (ez * ey * mct + ex * st) + vz * (ct + ez * ez * mct);
}

void read_lammps_data(geometry &g, char filename[])
{

	int fNv = 42;
	int fNe = 240;
	//float v[Nv][3];
	//int edge[Ne][2];
	//int etype[Ne];
	char s[100];
	char index[5], vtype[5], temp1[20], temp2[20], temp0[20];
	char TT[] = "Bonds";
	char AA[] = "Atoms";
	//int hetype=0;
	//int vin=-1;
	//int vout=-1;
	FILE *file;
	file = fopen(filename, "r");
	int x = 0;
	while (fscanf(file, "%s", s) == 1)
	{
		if (strcmp(s, AA) == 0)
		{
			break;
		}
	}
	double *vec = new double[3];
	for (int i = 0; i < fNv; i++)
	{

		x = fscanf(file, "%s %s %s %s %s\n", index, vtype, temp0, temp1, temp2);

		vec[0] = atof(temp0);
		vec[1] = atof(temp1);
		vec[2] = atof(temp2);
		//fprintf(stderr,"%s %s %s %s %s\n" ,index,vtype ,temp0,temp1,temp2);
		g.add_vertex(vec);
		//fprintf(stderr, "%d  %f %f %f \n",i, g.v[i][0],g.v[i][1],g.v[i][2]);
	}
	delete[] vec;
	//fprintf(stdout, "read all Vertices\n graph has %d vertices", g.Nv);

	cout << "now edges" << endl;
	while (fscanf(file, "%s", s) == 1)
	{ //s[0]!='1') {
		if (strcmp(s, TT) == 0)
		{
			break;
		}
	}
	for (int i = 0; i < fNe; i++)
	{
		x = fscanf(file, "%s %s %s %s\n", index, temp0, temp1, temp2);
		//fprintf(stderr,"%s %s %s %s \n" ,index,btype ,temp0,temp1);
		//g.e[i][0]=atoi(temp0)-1;
		//g.e[i][1]=atoi(temp1)-1;
		//g.etype[i]=atoi(btype);
		//if (atoi(btype)==1) {
		//fetype=2;
		//} else {
		//fetype=0;
		//}
		g.add_half_edge_type(atoi(temp1) - 1, atoi(temp2) - 1, atoi(temp0) - 1);
		//fprintf(stderr, "add edge %d  %d %d %d\n",i,atoi(temp1)-1, atoi(temp2)-1,atoi(temp0)-1 );
		//cout << i <<endl;
	}
	//int previd0=-1;
	//int nextid0=-1;
	cout << "unused x" << x << endl;

	for (vector<HE>::iterator it = g.he.begin(); it != g.he.end(); ++it)
	{
		g.update_half_edge(it->id);
	}
	//SEt prevoius next
	int Nt = 0;
	vector<HE>::iterator it = g.he.begin();
	//if (it==g.he.end()) { break;}
	if (!(it->nextid != -1 && it->previd != -1))
	{
		vector<HE>::iterator nextit = g.he.begin();
		while (true)
		{
			//cout << (nextit->id)<<endl;
			if (nextit == g.he.end())
			{
				break;
			}
			if (!(nextit->nextid != -1 && nextit->previd != -1))
			{
				vector<HE>::iterator previt = g.he.begin();
				while (true)
				{
					//cout << (previt->id)<<endl;
					if (previt == g.he.end())
					{
						break;
					}
					if (!(previt->nextid != -1 && previt->previd != -1))
					{
						if (it->vout == nextit->vin && nextit->vout == previt->vin && previt->vout == it->vin)
						{
							g.set_prev_next(it->id, previt->id, nextit->id);
							g.set_prev_next(previt->id, nextit->id, it->id);
							g.set_prev_next(nextit->id, it->id, previt->id);
							//cout << "set prev_next " << it->id << " "<< previt->id<< " " <<nextit->id<<endl;
							Nt++;
							//++it;
							//break;
							//cout << "Nt is" << Nt <<endl;
							//cout << "set prev_next " << it->id << " "<< previt->id<< " " <<nextit->id<<endl;
							//cout << "set prev_next " <<  previt->id<< " " <<nextit->id<<" "<< it->id <<endl;
							//cout << "set prev_next " <<  nextit->id<<" "<< it->id <<" "<< previt->id<<endl;

							++previt;
							++it;
							++nextit;
						}
					}

					++previt;
				}
			}
			++nextit;
		}
	}
	// it is done in two steps to avoid double counting triangles
	while (true)
	{

		//cout << (it->id)<<endl;
		if (it == g.he.end())
		{
			it = g.he.begin();
		}
		if (Nt == 80)
		{
			break;
		}
		//else { cout<< " running Nt is " << Nt<<endl;}
		int opindex = g.heidtoindex[it->opid];
		if ((it->nextid == -1 && it->previd == -1) && (g.he[opindex].nextid != -1 && g.he[opindex].previd != -1))
		{
			vector<HE>::iterator nextit = g.he.begin();
			while (true)
			{
				//cout << (nextit->id)<<endl;
				if (nextit == g.he.end())
				{
					break;
				}
				if (nextit->nextid == -1 && nextit->previd == -1)
				{
					vector<HE>::iterator previt = g.he.begin();
					while (true)
					{
						//cout << (previt->id)<<endl;
						if (previt == g.he.end())
						{
							break;
						}
						if (previt->nextid == -1 && previt->previd == -1)
						{
							if (it->vout == nextit->vin && nextit->vout == previt->vin && previt->vout == it->vin)
							{
								if (it->vin != nextit->vout && nextit->vin != previt->vout && previt->vin != it->vout)
								{
									g.set_prev_next(it->id, previt->id, nextit->id);
									g.set_prev_next(previt->id, nextit->id, it->id);
									g.set_prev_next(nextit->id, it->id, previt->id);
									//cout << "Nt is" << Nt <<endl;
									//cout << "set prev_next " << it->id << " "<< previt->id<< " " <<nextit->id<<endl;
									//cout << "set prev_next " <<  previt->id<< " " <<nextit->id<<" "<< it->id <<endl;
									//cout << "set prev_next " <<  nextit->id<<" "<< it->id <<" "<< previt->id<<endl;
									Nt++;
									//++it;
								}
							}
						}

						++previt;
					}
				}
				++nextit;
			}
		}
		++it;
	}
	/*for (vector<HE>::iterator it = g.he.begin() ; it != g.he.end(); it++) {
		cout << "EDGE  " <<distance(g.he.begin(),it) <<" id " << it->id << " vin " << g.vidtoindex[it->vin] << "vout " << g.vidtoindex[it->vout] <<endl;
		cout << "      opid " << it->opid << " nextid " << it->nextid << " previd " << it->previd <<endl;
		cout << "      hevec " << it->hevec[0] << " " << it->hevec[1]  << " "  << it->hevec[2] <<" l is " << it->l << endl;
		//cout << "      normal " << it->n[0] << " " << it->n[1] << " " << it->n[2] << " " <<endl <<endl;
	}
	cout << "Nt " << Nt <<endl;;*/
	//fprintf(stderr," edges : %d\n",g.Nhe);

	//exit(-1);
}

void dump_data_frame(geometry &g, FILE *f, int time)
{
	double avgL0 = 0, avgL1 = 0, avgTheta0 = 0, avgTheta1 = 0, avgPhi00 = 0, avgPhi11 = 0, avgPhi01 = 0;
	int L0 = 0, L1 = 0, Theta0 = 0, Theta1 = 0, Phi00 = 0, Phi11 = 0, Phi01 = 0;

	fprintf(f, "<configuration time_step=\"%d\">\n", time);
	fprintf(f, "<Edges num=\"%d\">\n", g.Nhe);

	for (vector<HE>::iterator it = g.he.begin(); it != g.he.end(); ++it)
	{

		fprintf(f, "%d %li %.4f\n", it->type, distance(g.he.begin(), it), it->l);
		//fprintf(stderr, "%d %d %.4f\n",  it->type,distance(g.he.begin(),it),it->l);
		if (it->type == 0)
		{
			avgL0 += it->l;
			L0 += 1;
		}
		else if (it->type == 1)
		{
			avgL1 += it->l;
			L1 += 1;
		}
	}
	fprintf(f, "<Theta>\n");
	//fprintf(stderr,"<Theta>\n");
	double theta;
	//for (int edge=0; edge<Ne; edge++)
	for (vector<HE>::iterator it = g.he.begin(); it != g.he.end(); ++it)
	{

		//{
		//if (t[edge][1] != -1 && t[edge][0] != -1) {

		//cout << "      normal " << it->n[0] << " " << it->n[1] << " " << it->n[2] << " " <<endl <<endl;
		//cout << "other  normal" << g.he[g.heidtoindex[it->opid)].n[0] << " " << g.he[g.heidtoindex[it->opid)].n[1] << " " << g.he[g.heidtoindex[it->opid)].n[2] << endl;
		//cout << dot(it->n,g.he[g.heidtoindex[it->opid)].n) <<endl;
		double ndot = dot(it->n, g.he[g.heidtoindex[it->opid]].n);
		if (ndot < -1)
		{
			ndot = -1;
		}
		//cout << "ndot is"  << ndot <<endl;
		theta = acos(ndot);
		//}
		fprintf(f, "%d %li %.4f\n", it->type, distance(g.he.begin(), it), theta);
		//fprintf(stderr, "%d %d %.4f\n", it->type, distance(g.he.begin(),it), theta);
		if (it->type == 0)
		{
			avgTheta0 += theta;
			Theta0 += 1;
			//cout << "Theta0  " << Theta0 <<endl;
			//cout <<  "avgTheta0" << avgTheta0 <<endl;
		}
		else if (it->type == 1)
		{
			avgTheta1 += theta;
			Theta1 += 1;
			//cout << "Theta1  " << Theta1 <<endl;
			//cout << "avgTheta1  " << avgTheta1 <<endl;
		}
		else
		{
			cout << "ERRRRRRRRRRRRRRRRRRRORRRRRRRRRRRRRRR , it->id" << endl;
		}
	}

	//cout << "avgTheta1/Theta1 " << avgTheta1/Theta1 <<endl;
	//cout << "avgTheta0/Theta0 " << avgTheta0/Theta0 <<endl;

	fprintf(f, "<Phi>\n");
	//fprintf(stderr,"<Phi>\n");
	//update_Phi();
	int phitype = -1;
	double phi;
	for (vector<HE>::iterator it = g.he.begin(); it != g.he.end(); ++it)
	{

		int nextindex = g.heidtoindex[it->nextid];
		int opindex = g.heidtoindex[it->opid];
		double ndot = (dot(g.he[opindex].hevec, g.he[nextindex].hevec) / (g.he[opindex].l * g.he[nextindex].l));
		if (ndot < -1)
		{
			ndot = -1;
		}
		if (ndot > 1)
		{
			ndot = 1;
		}
		phi = acos(ndot);
		int nexttype = g.he[nextindex].type;
		//cout << " type " << it->type << "nexttype" << nexttype <<endl;
		if (it->type == 0 && nexttype == 0)
		{
			phitype = 0;
			avgPhi00 += phi;
			Phi00 += 1;
			//cout << ":Phi00" <<Phi00 <<endl;
		}
		else if ((it->type == 0 && nexttype == 1) || (it->type == 1 && nexttype == 0))
		{
			phitype = 2;
			avgPhi01 += phi;
			Phi01 += 1;
		}
		else if ((it->type == 1 && nexttype == 1))
		{
			phitype = 1;
			avgPhi11 += phi;
			Phi11 += 1;
		}
		fprintf(f, "%d %d %d %.4f\n", phitype, it->id, it->nextid, phi);
		//fprintf(stderr, "%d %d %d %.4f\n", phitype, it->id, it->nextid,phi);
	}
	fprintf(stderr, " L0 %.d L1 %.d Theta0 %.d Theta1 %.d Phi00 %.d Phi11 %.d Phi01 %.d \n", L0, L1, Theta0, Theta1, Phi00, Phi11, Phi01);
	fprintf(stderr, " L0 %.3f L1 %.3f Theta0 %.3f Theta1 %.3f Phi00 %.3f Phi11 %.3f Phi01 %.3f \n", avgL0 / L0, avgL1 / L1, avgTheta0 / Theta0, avgTheta1 / Theta1, avgPhi00 / Phi00, avgPhi11 / Phi11, avgPhi01 / Phi01);
}

void recenter(geometry &g) {

	double XCM=0;
	double YCM=0;
	double ZCM=0;
	for (vector<VTX>::iterator it = g.v.begin(); it != g.v.end(); ++it)
	{    
		XCM +=it->co[0];
		YCM +=it->co[1];
		ZCM +=it->co[2];

	}
	XCM/=g.Nv;
	YCM/=g.Nv;
	ZCM/=g.Nv;

	//cout << "HERE2"<<endl;
	double vx = XCM;//g.v[0].co[0];
	double vy = YCM;//g.v[0].co[1];
	double vz = ZCM;//g.v[0].co[2];
	for (vector<VTX>::iterator it = g.v.begin(); it != g.v.end(); ++it)
	{
		it->co[0] -= vx;
		it->co[1] -= vy;
		it->co[2] -= vz;
	}
	g.update_surface();
}

int surfclosev(geometry &g){
	int alln=0;
	for (vector<VTX>::iterator it = g.v.begin(); it != g.v.end(); ++it){
		alln+=it->vneigh.size();
		if (it->vneigh.size()>0) {
			cout<< "vindex is " << distance(g.v.begin(),it) <<  " vid is" << it->vid ;
			for (vector<int>::iterator itv = it->vneigh.begin(); itv != it->vneigh.end(); itv++)
			{
				
				cout <<"     neighbors are " << *itv << " vindex is " << g.vidtoindex[*itv] ;  
				if (g.vidtoindex[*itv]==-1) { cout << "wrong neighbor" <<endl; exit(-1);}
			}
			cout<<endl;
		}
	}
		
	if (alln%2!=0) { cout <<" odd neighbors" <<endl; exit(-1);}
	
	
	int surfclosevCount=0;
	for (vector<int>::iterator it = g.surfv.begin(); it != g.surfv.end(); ++it)
	{
		if (g.v[g.vidtoindex[*it]].vneigh.size()>0) { 
			for (vector<int>::iterator itv = g.v[g.vidtoindex[*it]].vneigh.begin(); itv != g.v[g.vidtoindex[*it]].vneigh.end(); itv++)
			{    
				if (veclen(g.v[g.vidtoindex[*itv]].co,g.v[g.vidtoindex[*it]].co)<.5*g.l0[0]) { 
					surfclosevCount+=1; 
					break;
				}
			}
		}
	}
	return(surfclosevCount);
}
//void show_status(geometry &g, int frame, int sweep, int seconds ){




/*int valid(geometry &g){
	for (vector<HE>::iterator it = g.he.begin() ; it != g.he.end(); ++it) {
		if ((it->nextid==-1 && it->previd!=-1) || (it->nextid==-1 && it->previd!=-1)) {
			cout << " WRONG geometry " << endl;
			return -1;
		}
		if (it->opid==-1) { 
			cout<< " NO OPPOSIT EDGE" <<endl;
			return -1;
		}

		if (it->type!=g.he[g.heidtoindex[it->opid)].type) {
			cout << " TYPE MISMATCH";
			return -1;
		}


	}
	return 1;

}*/