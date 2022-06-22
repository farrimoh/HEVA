/*
 * Geometry.cpp
 *
 *  Created on: Apr 25, 2019
 *      Author: farri
 */

#include "Geometry.hpp"
#include "tri_tri_intersect.hpp"
#include <iostream>
#include <vector>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <queue>
#define PI 3.14159265
using namespace std;

Geometry::Geometry()
{
	// TODO Auto-generated constructor stub
	Test_assembly = 0;
	Nvlast = 0;
	Nhelast = 0;
	Ntype = 0;
	Nv = 0;
	Nv5 = 0;
	Nv6 = 0;
	Nhe = 0;
	epsilon = nullptr;
	kappa = nullptr;
	kappaPhi = nullptr;
	l0 = nullptr;
	theta0 = nullptr;
	phi0 = nullptr;
}

Geometry::~Geometry()
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
	if (boundary.size() > 0)
	{
		boundary.clear();
	}
	if (boundaryv.size() > 0)
	{
		boundaryv.clear();
	}

	if (fusionv.size() > 0)
	{
		boundaryv.clear();
	}
}

void Geometry::initialize(int Ntype0)
{
	Test_assembly=0; // for trial runs =1
	Nvlast=0;
	Nhelast=0;
	Ntype=Ntype0;
	Nv=0;
	Nv5=0;
	Nv6=0;
	Nhe=0;
	Nsurf=0;
	Nd=0;
	NAB=0;
	NAB_in=0;
	NCD=0;
	NCD_T4=0;
	NCD_T3=0;
	NCD_T4_in=0;
	NCD_T3_in=0;
	NCD_Hex=0;
	NCD_other=0;
	Nv_in=0;
	Nhe_in=0;
	Nboundary=1;
	Nboundarylast=0;
	accepted_vmove=0;
	rejected_vmove=0;

	all_neigh=0;
	gb0=0;
	dg=0;
	mudimer=0;
	mudrug=0;

	gaussian_sigma=0;
	l_thermal_sigma=0;
	l_thermal_kappa=0;
	theta_thermal_kappa=0;

	drugProb=0;
	T=1;

	xi=0;

	lenpoints=200000;
	dist_points = new double *[lenpoints];

	for (int i = 0; i < lenpoints; i++)
	{
		dist_points[i] = new double[3];
		for (int j = 0; j < 3; j++)
		{
			dist_points[i][j]=0.0;
		}
	}

	vidtoindex = new int[2000000000];
	heidtoindex = new int[2000000000];
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
		gdrug[i] = new double[Ntype];
		for (int j = 0; j < Ntype; j++)
		{
			gb[i][j] = 0;
			gdrug[i][j] = 0;
		}
	}

	for (int i = 0; i < 1000000; i++)
	{
		vidtoindex[i] = -1;
		heidtoindex[i] = -1;
	}
}

void Geometry::dump_parameters()
{

	FILE *pfile;
	pfile = fopen("geo_param.dat", "w");
	fprintf(pfile, "1");
}

void Geometry::update_index()
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
		//if ((is_boundary(it->id))>0 && (is_vboundary(it->vout)>0)) { v[vidtoindex[it->vout]].hesurfinid=it->id; }

		//v[vidtoindex[it->vin]].hein.push_back(*it);
	}
	//cout <<"T5" <<endl;
}

void Geometry::check_odd_neigh()
{
	int alln = 0;
	for (vector<VTX>::iterator it = v.begin(); it != v.end(); ++it)
	{
		alln += it->vneigh.size();
		//cout<< "vindex is " << distance(v.begin(),it) <<  " vid is" << it->vid ;
		//for (vector<int>::iterator itv = it->vneigh.begin(); itv != it->vneigh.end(); itv++)
		//{

		//	cout <<"     neighbors are " << *itv << " vindex is " << vidtoindex[*itv];
		//	if (vidtoindex[*itv]==-1) { cout << "wrong neighbor" <<endl; exit(-1);}
		//}
		//cout <<endl;
	}
	if (alln % 2 != 0)
	{
		update_neigh();
		alln = 0;
		for (vector<VTX>::iterator it = v.begin(); it != v.end(); ++it)
		{
			alln += it->vneigh.size();
			//cout<< "vindex is " << distance(v.begin(),it) <<  " vid is" << it->vid ;
			//for (vector<int>::iterator itv = it->vneigh.begin(); itv != it->vneigh.end(); itv++)
			//{

			//	cout <<"     neighbors are " << *itv << " vindex is " << vidtoindex[*itv];
			//	if (vidtoindex[*itv]==-1) { cout << "wrong neighbor" <<endl; exit(-1);}
			//}
			//cout <<endl;
		}
		if (alln % 2 != 0)
		{
			cout << " check_neigh odd neighbors" << endl;
			exit(-1);
		}
	}
}

void Geometry::update_neigh_vertex(int vid0)
{
	int vindex0 = vidtoindex[vid0];
	if (v[vindex0].vneigh.size() > 0)
	{
		v[vindex0].vneigh.clear();
	}
	for (vector<VTX>::iterator it = v.begin(); it != v.end(); ++it)
	{
		if (vid0 != it->vid && connected(vid0, it->vid) < 0)
		{ //} && next_connected(vid0,it->vid)<0){
			if (veclen(v[vindex0].co, it->co) < 3* xi)
			{
				v[vindex0].vneigh.push_back(it->vid);
			}
		}
	}
}

void Geometry::update_neigh_vertex_and_neigh(int vid0)
{
	/* update current neighbors */
	for (vector<int>::iterator vit = v[vidtoindex[vid0]].vneigh.begin(); vit != v[vidtoindex[vid0]].vneigh.end(); ++vit)
	{
		update_neigh_vertex(*vit);
	}

	/* update this vertex */
	update_neigh_vertex(vid0);

	/* update new neighbors */
	for (vector<int>::iterator vit = v[vidtoindex[vid0]].vneigh.begin(); vit != v[vidtoindex[vid0]].vneigh.end(); ++vit)
	{
		update_neigh_vertex(*vit);
	}
}

void Geometry::update_neigh()
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
		update_neigh_vertex(it->vid);
		/*for (vector<VTX>::iterator secondit = v.begin(); secondit != v.end(); ++secondit)
		{
			if (it->vid != secondit->vid && connected(it->vid, secondit->vid) < 0)
			{ //&& next_connected(it->vid,secondit->vid)<0){
				if (veclen(it->co, secondit->co) < 2 * xi)
				{
					it->vneigh.push_back(secondit->vid);
				}
			}
		}*/
	}
}

/****************** update boundary **********************/
/* This function updates the boundary related parameters */
/* Since I am keeping track of everything, (ToDo) it is possible to just keep it as validation step */

void Geometry::update_boundary()
{
	Nsurf = 0;
	NAB = 0;
	NCD = 0;

	//cout << "in update boundary"<<endl;

	//clear the boundary
	if (boundary.size() > 0)
		boundary.clear();
	if (boundaryv.size() > 0)
		boundaryv.clear();
	if (boundaryvbond.size() > 0)
		boundaryvbond.clear();

	//clear hein, hebundaryoutid hebundaryoutid2 (for double boundary) , update all_neigh
	all_neigh = 0;
	for (vector<VTX>::iterator it = v.begin(); it != v.end(); ++it)
	{
		it->heboundaryoutid = -1;
		it->heboundaryoutid2 = -1;
		if (it->hein.size() > 0)
			it->hein.clear();
		all_neigh += it->vneigh.size();
	}

	//UPDATE INDICES
	update_index(); // this also update hein for vertices // why not here?

	// Update boundaryheid , boundaryv (he->vin) , boundaryvbond
	// if !double surf ->

	for (vector<HE>::iterator it = he.begin(); it != he.end(); ++it)
	{
		////cout << " heid" << it->id <<" boundary_index " << it->boundary_index <<endl;
		//it->nextid_boundary=-1;
		//it->previd_boundary=-1;

		if (it->type == 0)
			NCD++;
		if (it->type == 1)
			NAB++;

		// update vertex hein and double check boundary_index

		if (is_boundary(it->id) > 0)
		{
			boundary.push_back(it->id);

			//all vertices have heboundaryoutid , doubleboundaries have 2
			v[vidtoindex[it->vin]].heboundaryoutid = it->id;

			boundaryv.push_back(it->vin);
			//cout << "boundary updated with edge id " << it->id << endl<< endl;
			//if (it->boundary_index==-1) { it->boundary_index=0;}
			//cout << " in update boundary - wrong boundary boundary index"<<endl; cout << " heid" << it->id <<" boundary_index " << it->boundary_index <<endl; exit(-1);}
		}

		//else {
		//if (it->boundary_index!=-1) { it->boundary_index=-1;}
		//cout << " in update boundary - wrong inside boundary_index"<<endl; cout << " heid" << it->id <<" boundary_index " << it->boundary_index <<endl;exit(-1);}
		//}

		if (is_bond_in_boundary(it->id) >= 0)
		{
			boundaryvbond.push_back(it->vin);
		}

		if (it->vin == -1 || it->vout == -1)
		{
			cout << " update_boundary ! error in vin vout of edge " << it->id << " it->vin " << it->vin << " it->vout " << it->vout << endl;
			exit(-1);
		}

		if (vidtoindex[it->vin] == -1 || vidtoindex[it->vout] == -1)
		{
			cout << " update_boundary ! error in vin vout of edge " << it->id << " vidtoindex[it->vin] " << vidtoindex[it->vin] << " vidtoindex[it->vout] " << vidtoindex[it->vout] << endl;
			exit(-1);
		}
	}
	//cout << "in update boundary 222"<<endl;
	int N_doubleboundary = 0;
	// update doubleboundary vertices
	for (vector<int>::iterator vt = boundaryv.begin(); vt != boundaryv.end(); ++vt)
	{
		int vindex0 = vidtoindex[*vt];

		v[vindex0].doubleboundary = -1;

		int nboundary = 0;

		for (vector<int>::iterator it = v[vindex0].hein.begin(); it != v[vindex0].hein.end(); ++it)
		{
			if (is_boundary(*it) > 0)
				nboundary++;
		}

		if (nboundary == 2)
		{
			//cout << "*vt " << *vt << "doubleboundary" << v[vindex0].doubleboundary << endl;
			v[vindex0].doubleboundary = 1;
			N_doubleboundary++;
		}
	}
	//g.Nboundary=N_doubleboundary+1;
	// now go over the boundary again and add the heboundaryoutid2
	for (vector<int>::iterator it = boundary.begin(); it != boundary.end(); ++it)
	{

		int vin0 = he[heidtoindex[*it]].vin;
		if (v[vidtoindex[vin0]].doubleboundary == 1)
		{
			if ((v[vidtoindex[vin0]].heboundaryoutid != *it) && (v[vidtoindex[vin0]].heboundaryoutid2 == -1))
			{
				v[vidtoindex[vin0]].heboundaryoutid2 = *it;
			}
		}
	}
	//cout << "in update boundary 333"<<endl;
	//if (Nboundary != 1) cout << "AFTER first set of initialization in update boundary"<<endl;
	for (vector<HE>::iterator it = he.begin(); it != he.end(); ++it)
	{
		update_half_edge(it->id);

		if (Nhe == 6) //(Nboundary == 1)
		{
			if (is_boundary(it->id) > 0)
			{
				if (it->boundary_index == -1)
				{
					it->boundary_index = 0;
				}
				//cout << " in update boundary - wrong boundary boundary index"<<endl; cout << " heid" << it->id <<" boundary_index " << it->boundary_index <<endl; exit(-1);}
				//update boundary_nextindex
				it->nextid_boundary = v[vidtoindex[it->vout]].heboundaryoutid;
				he[heidtoindex[v[vidtoindex[it->vout]].heboundaryoutid]].previd_boundary = it->id;
			}
			else
			{
				if (it->boundary_index != -1)
				{
					it->boundary_index = -1;
				}
				//cout << " in update boundary - wrong inside boundary_index"<<endl; cout << " heid" << it->id <<" boundary_index " << it->boundary_index <<endl;exit(-1);}
			}
		}
		else
		{

			if (is_boundary(it->id) > 0)
			{
				if (it->boundary_index == -1)
				{
					cout << " in update boundary - wrong boundary boundary index" << endl;
					cout << " heid " << it->id << " boundary_index " << it->boundary_index << endl;
					exit(-1);
				}

				/*************** Temp Test the nextid_boundary ******/
				if (it->boundary_index != he[heidtoindex[it->nextid_boundary]].boundary_index || it->boundary_index != he[heidtoindex[it->previd_boundary]].boundary_index)
				{
					cout << "error in boundary index" << endl;
					cout << " heid " << it->id << " boundary_index " << it->boundary_index;
					cout << " nextid_boundary " << it->nextid_boundary << " boundary index of nextid_boundary " << he[heidtoindex[it->nextid_boundary]].boundary_index;
					cout << " previd_boundary " << it->previd_boundary << " boundary index of previd_boundary " << he[heidtoindex[it->previd_boundary]].boundary_index << endl;
					exit(-1);
				}

				/*************** Temp Test the nextid_boundary ******/
				if (it->nextid != -1)
				{
					if (it->nextid_boundary != it->nextid)
					{
						cout << "next previous dont match nextid_boundary previd_boundary" << endl;
						cout << " heid " << it->id;
						cout << " nextid " << it->nextid << " nextid_boundary " << it->nextid_boundary;
						cout << " previd " << it->previd << " previd_boundary " << it->previd_boundary << endl;
						exit(-1);
					}
				}
				if (it->previd != -1)
				{
					if (it->previd_boundary != it->previd)
					{
						cout << "next previous dont match nextid_boundary previd_boundary" << endl;
						cout << " heid " << it->id;
						cout << " nextid " << it->nextid << " nextid_boundary " << it->nextid_boundary;
						cout << " previd " << it->previd << " previd_boundary " << it->previd_boundary << endl;
						exit(-1);
					}
				}
			}
			else
			{
				if (it->boundary_index != -1)
				{
					cout << " error in update boundary - wrong inside boundary_index" << endl;
					cout << " heid" << it->id << " boundary_index " << it->boundary_index << endl;
					exit(-1);
				}
			}
		}
	}

	//cout << "AFTER second set of initialization in update boundary"<<endl;
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
				cout << "wrong g.he[nextindex].previd in update boundary edge" << it->id << " he[nextindex].previd " << he[nextindex].previd << endl;
				exit(-1);
			}
		}
		if (it->previd != -1)
		{ //} && it->previd!=-1)  { //this and next
			int previndex = heidtoindex[it->previd];

			if (he[previndex].nextid == -1)
			{
				cout << "it->previd" << it->previd << "previndex" << previndex << endl;
				cout << "wrong g.he[previndex].nextid in update boundaryedge" << it->id << " he[previndex].nextid " << he[previndex].nextid << endl;
				exit(-1);
			}
		}
	}*/
	/*for (vector<int>::iterator it = boundary.begin(); it != boundary.end(); ++it)
	{
		
		//if ((is_boundary(it->id))>0 && (is_vboundary(it->vout)>0)) { v[vidtoindex[it->vout]].heboundaryinid=it->id; }
		cout << "heid" << *it << " next " << he[heidtoindex[*it]].nextid << " prev " << he[heidtoindex[*it]].previd << endl;
		cout << "vin " << he[heidtoindex[*it]].vin << "vout " << he[heidtoindex[*it]].vout <<endl;
		//v[vidtoindex[it->vin]].hein.push_back(*it);
	}*/

	//TEMP DOUBLE CHECK
	/*
	for (vector<int>::iterator ht = boundary.begin(); ht != boundary.end(); ++ht)
	{
		//cout << "on boundary" <<endl;
		//cout << "id " << *ht << " next id " << he[heidtoindex[*ht]].nextid << " prev id " << he[hei:update_boundary
		he[heidtoindex[*ht]].previd_boundary=he[heidtoindex[*ht]].previd;
	}*/
	update_normals();
	update_excluder_top();
	Nsurf = boundary.size();
}

int Geometry::is_vboundary(int vid)
{ // 1 if
	for (vector<int>::iterator it = boundaryv.begin(); it != boundaryv.end(); ++it)
	{
		if (*it == vid)
			return 1;
	}
	return -1;
}

int Geometry::is_bond_vboundary(int vid)
{
	for (vector<int>::iterator it = boundaryvbond.begin(); it != boundaryvbond.end(); ++it)
	{
		if (*it == vid)
			return 1;
	}
	return -1;
}
int Geometry::is_boundary(int heid0)
{ // it is on the boundary if it has no next prevoius
	int heindex0 = heidtoindex[heid0];
	//int opindex0=heidtoindex[he[heindex0].opid);
	if (he[heindex0].nextid == -1 || he[heindex0].previd == -1)
	{ // on boundary
		//cout << " he " << heid0 << " on boundary " ;
		return 1;
	}
	//cout << " in is_boundary normal is " << he[heindex0].n[0] << " " <<endl;
	return -1; // not on boundary
}

int Geometry::is_bond_in_boundary(int heid0)
{ // it is on the boundary if it has no next prevoius
	if (heid0 == -1)
	{
		cout << "id does not exist !!!!!!! in is_bond_in_boundary " << endl;
		exit(-1);
	}
	int heindex0 = heidtoindex[heid0];
	//int opindex0=heidtoindex[he[heindex0].opid);
	//cout<<"in Is_bond_boundary) for edge heid0 "<<heid0 << endl;
	//cout<<"he[heindex0].nextid" << he[heindex0].nextid << "and  he[heindex0].previd " << he[heindex0].previd <<endl

	if (he[heindex0].nextid == -1 && he[heindex0].previd != -1)
	{ // not on boundary
		//cout << " he " << heid0 << " is bond_in boundary " ;
		//cout << "now return he[heindex0].vin" << he[heindex0].vin <<endl;
		return (he[heindex0].vin);
	}

	//cout << " in is_boundary normal is " << he[heindex0].n[0] << " " <<endl;
	return -1; // on boundary
}

int Geometry::is_bond_out_boundary(int heid0)
{ // it is on the boundary if it has no next prevoius
	if (heid0 == -1)
	{
		cout << "id does not exist !!!!!!! in is_bond_out_boundary " << endl;
		exit(-1);
	}
	int heindex0 = heidtoindex[heid0];
	//int opindex0=heidtoindex[he[heindex0].opid);
	if (he[heindex0].nextid != -1 && he[heindex0].previd == -1)
	{ // not on boundary
		//cout << " he " << heid0 << " not on boundary " ;
		return he[heindex0].vout;
	}
	//cout << " in is_boundary normal is " << he[heindex0].n[0] << " " <<endl;
	return -1; // on boundary
}

int Geometry::no_bond_boundary(int heid0)
{ // it is on the boundary if it has no next prevoius
	int heindex0 = heidtoindex[heid0];
	//int opindex0=heidtoindex[he[heindex0].opid);
	if ((he[heindex0].nextid == -1 && he[heindex0].previd == -1))
	{ // not on boundary
		//cout << " he " << heid0 << " not on boundary " ;
		return 1;
	}
	//cout << " in is_boundary normal is " << he[heindex0].n[0] << " " <<endl;
	return -1; // on boundary
}

int Geometry::is_same_triangle(int heid0, int heid1)
{
	if (heid0 == heid1)
		return 1;
	int heindex0 = heidtoindex[heid0];
	if (heid1 == he[heindex0].opid)
		return 1;
	int heopindex0 = heidtoindex[he[heindex0].opid];
	int opnextid = he[heopindex0].nextid;
	int opprevid = he[heopindex0].previd;
	if (opnextid != -1)
	{
		if (heid1 == opnextid)
			return 1;
		if (heid1 == he[heidtoindex[opnextid]].opid)
			return 1;
	}
	if (opprevid != -1)
	{
		if (heid1 == opprevid)
			return 1;
		if (heid1 == he[heidtoindex[opprevid]].opid)
			return 1;
	}
	return -1;
}

int Geometry::not_cross_edge(int heid0, int heid1)
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

int Geometry::next_open_wedge(int heid0)
{
	//cout <<" in next_open_wedge"<<endl;

	if (Nhe < 10)
	{
		return -1;
	}

	int heindex_prev = heidtoindex[heid0]; // this edge
	if (he[heindex_prev].nextid != -1)
	{
		return -1;
	}
	if (he[heindex_prev].previd != -1)
	{
		return -1;
	}
	int vindex0 = vidtoindex[he[heindex_prev].vout];
	//if (v[vindex0].doubleboundary!=-1) { return -1;}  // sould not be double surf
	if (v[vindex0].hein.size() < 4)
	{
		return -1;
	} //  not the spike vertices
	int opid_prev = he[heindex_prev].opid;
	int heidnext = -1;
	heidnext = he[heindex_prev].nextid_boundary;

	/*for (vector<int>::iterator it = v[vindex0].hein.begin(); it != v[vindex0].hein.end(); ++it)
	{
		int opid0 = he[heidtoindex[*it]].opid; // becaus it looks at hein
		//cout << "in next open wedge, opid is "<< opid0<<endl;
		if (is_boundary(opid0) > 0 && he[heidtoindex[opid0]].previd == -1)
		{
			heidnext = opid0;
		} // find the one that doesn't have previd !!! there should be just one!
	}*/
	//cout << "in next open wedge, heidnext "<< heidnext<<endl;
	if (heidnext == -1)
	{
		cout << " in get_next_wedge! error no nextid_boundary" << endl;
		exit(-1);
	}

	int heindex_next = heidtoindex[heidnext];
	if (he[heindex_next].nextid != -1)
	{
		return -1;
	} // next should not be bound
	// temp double check
	if (he[heindex_prev].vout != he[heindex_next].vin)
	{
		cout << "Wrong geometry vin vout in next_open wedge " << endl;
		exit(-1);
	}

	double cangle = dot(he[heidtoindex[opid_prev]].hevec, he[heindex_next].hevec) / (he[heidtoindex[opid_prev]].l * he[heindex_next].l);
	//temp check
	if (cangle > 1)
	{
		cout << "in next_open_wedge" << endl;
		exit(-1);
	}
	//cout << "cangle" <<cangle<<endl;
	if (cangle > .2 && not_cross_edge(he[heindex_prev].id, he[heindex_next].id) > 0) // do we need not cross wedge?
	{
		return heidnext;
	}

	return -1;
}

int Geometry::pre_open_wedge(int heid0)
{
	if (Nhe < 10)
	{
		return -1;
	}

	int heindex_next = heidtoindex[heid0]; // this edge
	if (he[heindex_next].previd != -1)
	{
		return -1;
	}
	//temp double check
	if (he[heindex_next].nextid != -1)
	{
		return -1;
	}
	int vindex0 = vidtoindex[he[heindex_next].vin];
	//if (v[vindex0].doubleboundary!=-1) { return -1;}  // sould not be double surf
	if (v[vindex0].hein.size() < 4)
	{
		return -1;
	} //  not the spike vertices

	int heidprev = -1;
	heidprev = he[heindex_next].previd_boundary;
	/* 
	for (vector<int>::iterator it = v[vindex0].hein.begin(); it != v[vindex0].hein.end(); ++it)
	{
		//find prev
		if (is_boundary(*it) > 0 && he[heidtoindex[*it]].nextid == -1)
		{
			heidprev = *it;
		} // find the one that doesn't have previd !!! there should be just one!
	}
	*/
	//cout << "in pre open wedge, heidprev "<< heidprev<<endl;
	if (heidprev == -1)
	{
		cout << " in get_pre_wedge! error no previd_boundary" << endl;
		exit(-1);
	}

	int heindex_prev = heidtoindex[heidprev];
	int opid_prev = he[heindex_prev].opid;
	if (he[heindex_prev].previd != -1)
	{
		return -1;
	} // prev should not be bound to prev
	// temp double check
	if (he[heindex_prev].vout != he[heindex_next].vin)
	{
		cout << "Wrong geometry vin vout in next_open wedge " << endl;
		exit(-1);
	}

	double cangle = dot(he[heidtoindex[opid_prev]].hevec, he[heindex_next].hevec) / (he[heidtoindex[opid_prev]].l * he[heindex_next].l);
	//temp check
	if (cangle > 1)
	{
		cout << "in pre_open_wedge" << endl;
		exit(-1);
	}
	//cout << "cangle" <<cangle<<endl;
	if (cangle > .2 && not_cross_edge(he[heindex_prev].id, he[heindex_next].id) > 0) // do we need not cross wedge?
	{
		return heidprev;
	}

	return -1;
}

int Geometry::open_wedge(int heid0, int *flag) // ised in add_monomer_dimer
{
	//if (Nhe<21) { return -1;}
	int heindex0 = heidtoindex[heid0]; // this edge
	int heindex1;
	int opid0 = he[heindex0].opid;
	double *tempvec = new double[3];
	//cout << " in open_wedge _ checking heid0 " << heid0  << " heindex0 is "<< heindex0 <<endl;
	//int nb=he[heindex0].boundary_index;
	for (vector<int>::iterator it = boundary.begin(); it != boundary.end(); ++it)
	{
		if ((*it != heid0) && (*it != opid0))
		{								 // other boundary edges
			heindex1 = heidtoindex[*it]; // iterating on boundary edges
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

int Geometry::next_connected_boundary(int vid0, int vid1)
{
	/*for (vector<VTX>::iterator it = v.begin() ; it != v.end(); ++it) {
		
		if (vid0!=it->vid && vid1!=it->vid) {
			if (connected(vid0,it->vid)>0 && connected(vid1,it->vid)>0) {
				return it->vid;
			}
		}

	}*/

	for (vector<int>::iterator it = boundary.begin(); it != boundary.end(); ++it)
	{
		int heindex0 = heidtoindex[*it];
		if (he[heidtoindex[he[heindex0].nextid_boundary]].vout == vid0 && he[heidtoindex[he[heindex0].previd_boundary]].vin == vid1)
			return 1;
		if (he[heidtoindex[he[heindex0].nextid_boundary]].vout == vid1 && he[heidtoindex[he[heindex0].previd_boundary]].vin == vid0)
			return 1;
	}
	return -1;
}
int Geometry::connected(int vid0, int vid1)
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

int Geometry::connectedH(int heid0, int heid1)
{
	if (he[heidtoindex[heid0]].vin == he[heidtoindex[heid1]].vin)
		return -1;
	if (he[heidtoindex[heid0]].vout == he[heidtoindex[heid1]].vin)
		return -1;
	if (he[heidtoindex[heid0]].vout == he[heidtoindex[heid1]].vout)
		return -1;
	if (he[heidtoindex[heid0]].vin == he[heidtoindex[heid1]].vout)
		return -1;
	return 1;
}

void Geometry::add_vertex(double *xyz)
{
	VTX *vtxi;
	vtxi = new VTX;
	vtxi->vid = Nvlast;
	vtxi->co[0] = xyz[0];
	vtxi->co[1] = xyz[1];
	vtxi->co[2] = xyz[2];
	//vtxi->hesurfinid=-1;
	vtxi->heboundaryoutid = -1;
	vtxi->heboundaryoutid2 = -1;
	//vtxi->fusion_vid=-1;
	vtxi->doubleboundary = -1;
	vidtoindex[Nvlast] = Nv;
	v.push_back(*vtxi);

	Nv++;
	Nvlast++;
	delete vtxi;
}

int Geometry::delete_vertex(int vid0)
{

	int vindex0 = vidtoindex[vid0];

	v.erase(v.begin() + vindex0);
	vidtoindex[vid0] = -1;
	Nv--;

	return 1;
}

int Geometry::remove_neigh(int vid0, int removevid)
{
	int vindex0 = vidtoindex[vid0];
	int counter = -1;
	if (v[vindex0].vneigh.size() == 0)
	{
		cout << " has no neighbor in remove neigh for " << vid0 << endl;
		exit(-1);
	}
	for (vector<int>::iterator it = v[vindex0].vneigh.begin(); it != v[vindex0].vneigh.end(); ++it)
	{
		if (*it == removevid)
		{
			counter = distance(v[vindex0].vneigh.begin(), it);
		}
	}
	if (counter != -1)
	{
		v[vindex0].vneigh.erase(v[vindex0].vneigh.begin() + counter);
	}
	else
	{
		cout << " not found neigh in remove neigh" << endl;
		exit(-1);
	}
	return 1;
}

int Geometry::delete_edge(int heid0)
{
	//cout << "IN DELETE EDGE " <<endl;

	int heindex0 = heidtoindex[heid0]; //this edge
	if (is_boundary(heid0) < 0)
	{
		cout << " edge is not on boundary, cannot remove" << endl;
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

void Geometry::add_half_edge_type(int vin0, int vout0, int etype, int b_index)
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
	he_initialize(Nhe, Nhelast, vin0, vout0, etype, b_index);
	Nhe++;
	Nhelast++;
	delete newhei;
}

int Geometry::add_edge_type(int vin0, int vout0, int etype)
{
	int b_index = -1;
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
	he_initialize(Nhe, Nhelast, vin0, vout0, etype, b_index);
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
	he_initialize(Nhe, Nhelast, vout0, vin0, optype, b_index);
	//heidtoindex[Nhelast]=Nhe;
	Nhe++;
	Nhelast++;
	he[Nhe - 1].opid = he[Nhe - 2].id;
	he[Nhe - 2].opid = he[Nhe - 1].id;

	delete newhei;
	delete newheo;
	return 1;
}

int Geometry::opposite_edge(int heid0)
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
	cout << " No opposite  Invalid GEometry" << endl;
	exit(-1);
}

void Geometry::set_prev_next(int heid0, int previd0, int nextid0)
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

	he[heindex].boundary_index = -1;
	he[heindex].nextid_boundary = -1;
	he[heindex].previd_boundary = -1;
}

void Geometry::set_prev_next_boundary(int previd0, int nextid0)
{

	if ((is_boundary(previd0) < 0) || (is_boundary(nextid0) < 0))
	{
		cout << " wrong in set prev next boundary" << endl;
		int previndex0 = heidtoindex[previd0];
		cout << "he[previndex0].previd" << he[previndex0].previd << " he[previndex0].nextid " << he[previndex0].nextid << endl;
		int nextindex0 = heidtoindex[nextid0];
		cout << "he[nextindex0].previd" << he[nextindex0].previd << " he[nextindex0].nextid " << he[nextindex0].nextid << endl;
		exit(-1);
	}

	int nextindex = heidtoindex[nextid0];
	int previndex = heidtoindex[previd0];
	//cout << "he[heindex].vout" << he[heindex].vout << " he[nextindex].vin " << he[nextindex].vin << " he[heindex].vin " <<  he[heindex].vin  << " he[previndex].vout "  << he[previndex].vout << endl;
	if ((he[previndex].vout != he[nextindex].vin))
	{
		cout << " WRONG VIN VOUT IN SET_PREV_NEXT" << endl;
		exit(-1);
	}
	if ((he[previndex].boundary_index != he[nextindex].boundary_index))
	{
		cout << " WRONG boundary index in set prev next boundary" << endl;
		exit(-1);
	}

	he[previndex].nextid_boundary = nextid0;
	he[nextindex].previd_boundary = previd0;
}

int Geometry::get_prev_next(int heid0, int *previd0, int *nextid0)
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

double Geometry::add_dimer(int heid0, gsl_rng *r, int et1, int et2 , double *distance_vector)
{ // adds the dimer if it  has no next prev
	if (heid0 >= Nhelast)
	{
		cout << " in add_dimer ERROR\n";
		exit(-1);
	}
	double *newv = new double[3];
	double *tempv1 = new double[3];
	double *tempv2 = new double[3];
	
	
	double success = 1;
	int heindex0 = heidtoindex[heid0];

	int nextboundary = he[heindex0].nextid_boundary;
	int prevboundary = he[heindex0].previd_boundary;

	
	if (is_boundary(heid0) < 0)
	{
		cout << " he is not on boundary " << endl;
		delete[] newv;
		delete[] tempv1;
		delete[] tempv2;
		return -1;
	} 
	 
	if (he[heindex0].nextid != -1 || he[heindex0].previd != -1)
	{
		cout << "one end ins connected cannot add dimer" << endl;
		exit(-1);
	}
	//read boundary_index and update it for heid and added halfedges
	int bi = he[heindex0].boundary_index;

	//find the eqy=uilibrium point of new dimer XXXX
	//cout << " in add_dimer" <<endl;
	new_vertex_edge(heindex0, tempv1, et1);
	//double d=new_vertex_edge_and_move(heindex0, newv, et1,r);
	
	
	//double d=move_p(xi, tempv1, newv, r);
	double d=move_p_gaussian(xi,tempv1, newv, r);
	//subvec(tempv1,newv,distance_vector);
	//cout << "now move it" <<endl;
	//double d=move_p_gaussian_axes(heindex0,xi,tempv1, newv, r, distance_vector); //  also return the projections in each direction  for calculation of vp in detailed balance

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

	he[heidtoindex[Nhelast - 1]].boundary_index = bi;
	he[heidtoindex[Nhelast - 3]].boundary_index = bi;
	set_prev_next_boundary(prevboundary, Nhelast - 1);
	set_prev_next_boundary(Nhelast - 3, nextboundary);
	set_prev_next_boundary(Nhelast - 1, Nhelast - 3);
	//cout << " in add_dimer update boundary " << endl;
	update_boundary();
	update_neigh_vertex(Nvlast - 1);
	//update
	/*int vindex0=vidtoindex[Nvlast-1];
	if (v[vindex0].vneigh.size() > 0)
	{
		for (vector<int>::iterator it = v[vindex0].vneigh.begin(); it != v[vindex0].vneigh.end(); ++it)
		{	
			update_neigh_vertex(*it);
		}
	}*/
	if (success>0) success=d;
	//cout << "dimer added in G.add_dimer" <<endl;
	return success;
}

int Geometry::force_add_dimer(int heid0, double *newv, int typenext, int typeprev)
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
	if (is_boundary(heid0) < 0)
	{
		delete[] tempv1;
		delete[] tempv2;
		return -1;
	} //cout << " he is not on boundary " <<endl; return 0;}

	add_vertex(newv);

	//cout << " add_dimer after adding a vertex Nhe is   " << Nhe << " Nv is " << Nv << endl;
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

	update_boundary();
	return success;
}

int Geometry::remove_dimer(int heindex0, int heindexnext0)
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

void Geometry::make_hexamer()
{
	/*make_triangle();
	update_boundary();
    update_normals();

	//for (int i=0; i<5; i++){ // make pentamer
	double *vco;
	new_vertex(Nhe-1,vco);
	 
	force_add_dimer(Nhe-1,vco,1,1);
	update_boundary();
    update_normals();
	cout << "1 dimer added" <<endl;
	new_vertex(Nhe-1,vco);
	force_add_dimer(Nhe-1,vco,1,0);
	update_boundary();
	update_normals();
	cout << "2 dimer added" <<endl;
	new_vertex(Nhe-1,vco);
	force_add_dimer(Nhe-1,vco,0,0);
	update_boundary();
	update_normals();
	cout << "3 dimer added" <<endl;
	new_vertex(Nhe-1,vco);
	force_add_dimer(Nhe-1,vco,1,1);
	update_boundary();
    update_normals();
	cout << "4 dimer added" <<endl;
	force_add_monomer(1,Nhe-1,1);
	update_boundary();
    update_normals();
	delete[] vco;
        //for (int rstep=0; rstep<100; rstep++){
        //    move_vertex(g,r);
		//    g.update_boundary();
            //dump_lammps_data_file(g, frame++);
        //} 
        //dump_lammps_data_file(g, frame++);
	*/
}

void Geometry::he_initialize(int heindex, int heid0, int vin0, int vout0, int etype, int b_index)
{
	he[heindex].id = heid0;
	he[heindex].opid = -1;
	he[heindex].vin = vin0;
	he[heindex].vout = vout0;
	he[heindex].nextid = -1;
	he[heindex].previd = -1;
	he[heindex].nextid_boundary = -1;
	he[heindex].previd_boundary = -1;
	he[heindex].next_fusion_heid = -1;
	he[heindex].prev_fusion_heid = -1;
	he[heindex].next_wedge_fusion_heid = -1;
	he[heindex].prev_wedge_fusion_heid = -1;
	int vindex0 = vidtoindex[vin0];
	int vindex1 = vidtoindex[vout0];
	double *tempv = new double[3];
	subvec(v[vindex0].co, v[vindex1].co, tempv);
	he[heindex].hevec[0] = tempv[0];
	he[heindex].hevec[1] = tempv[1];
	he[heindex].hevec[2] = tempv[2];
	he[heindex].l = norm(tempv);
	he[heindex].type = etype;
	he[heindex].boundary_index = b_index;
	delete[] tempv;
	he[heindex].din = 0;
	//he[heindex].dout = 0;
	heidtoindex[heid0] = heindex;
}

void Geometry::update_edge(int heid0)
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

void Geometry::update_half_edge(int heid0)
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

int Geometry::add_monomer(int nextofnewid, int prevofnewid, int et)
{ //nextofnew prevofnew

	int heid0 = prevofnewid;
	//cout << "in add_monomer" << "heid0 is prevofnewid " << heid0 <<endl;
	int heid1 = nextofnewid;

	//cout << "in add_monomer" << "heid1 is nextofnewid " << heid1 <<endl;
	int heprevindex = heidtoindex[heid0]; // this is prev for new edge
	int henextindex = heidtoindex[heid1]; // this is going to be next of new edge

	//read boundary index and edit it for added halfedges
	int biprev = he[heprevindex].boundary_index;
	int binext = he[henextindex].boundary_index;
	if (biprev != binext)
	{
		cout << "boundary index not matching. " << endl;
		exit(-1);
	}
	//read and keep nextid_boundary and previd_boundary
	int bound_nextid = he[heprevindex].nextid_boundary;
	int bound_previd = he[henextindex].previd_boundary;

	//double check
	if ((he[heprevindex].previd != nextofnewid) || (he[henextindex].nextid != prevofnewid))
	{
		cout << "wrong chosen half edges for add monomer" << endl;
		exit(-1);
	}
	int newvin = he[heprevindex].vout;
	int newvout = he[henextindex].vin;
	double *tempv1 = new double[3];
	double *tempv2 = new double[3];
	addvec(v[vidtoindex[newvin]].co, v[vidtoindex[newvout]].co, tempv1);
	multvec(tempv1, .5, tempv2);

	if (check_overlap_centerh(tempv2) < 0)
	{
		//cout << "overalp in adding edge monomer! \n \n";
		delete[] tempv1;
		delete[] tempv2;
		return -1;
	}

	if (add_edge_type(newvin, newvout, et) < 0)
	{
		//cout << "coudnot add the edge " << endl;
		delete[] tempv1;
		delete[] tempv2;
		return -1;
	}
	update_half_edge(Nhelast - 2);
	update_half_edge(Nhelast - 1);
	update_half_edge(heid0);
	set_prev_next(Nhelast - 2, heid0, heid1);
	set_prev_next(heid1, Nhelast - 2, heid0);
	set_prev_next(heid0, heid1, Nhelast - 2);
	set_prev_next(Nhelast - 1, -1, -1);
	he[heprevindex].boundary_index = -1;
	he[henextindex].boundary_index = -1;
	he[heidtoindex[Nhelast - 1]].boundary_index = biprev;
	//cout << "in add_monomer" <<endl;

	set_prev_next_boundary(Nhelast - 1, bound_nextid);

	set_prev_next_boundary(bound_previd, Nhelast - 1);

	//cout << "in g.add_monomer , monomer added"<<endl;
	update_boundary();
	update_normals();
	delete[] tempv1;
	delete[] tempv2;
	return 1;
}

int Geometry::force_add_monomer(int nextofnewid, int prevofnewid, int etype)
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
	//update_boundary();
	//delete[] tempv1; delete[] tempv2;
	return 1;
}
int Geometry::add_monomer_dimer(int heid0)
{

	if (is_boundary(heid0) < 0)
	{
		cout << " cannot add not on the boundary !" << endl;
		exit(-1);
	}
	update_boundary();
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
	//else
	//{
	//cout << "add_dimer Nhe is " << Nhe<<endl;
	//add_dimer(heid0,r,3,3);
	//}
	return 1;
}

int Geometry::remove_monomer_dimer(int heid0, gsl_rng *r) /* NOT USED */
{
	/*if (!(is_boundary(heid0) > 0))
	{
		cout << "not on boundary cannot remove" << endl;
		return -1;
	}
	int heindex0 = heidtoindex[heid0];																		   // index of this edge
	int heopindex0 = heidtoindex[he[heindex0].opid];														   // indexd of opposite edge
	int nextopid0 = he[heopindex0].nextid;																	   // indexd of next of opposite edge
	int prevopid0 = he[heopindex0].previd;																	   // indexd of prev of opposite edge
	int heid1 = he[heidtoindex[nextopid0]].opid;															   // now back to this side
	int heid2 = he[heidtoindex[prevopid0]].opid;															   // afetr vertex
	if (is_boundary(heid1) < 0 && is_boundary(heid2) < 0 && is_vboundary(he[heidtoindex[nextopid0]].vout) < 0) // this prevents dissociation
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

		if (is_boundary(heid1) > 0)
		{
			//cout << "\n \n TRYING DIMER REMOVAL - with its next" <<endl;
			
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
				//double crit = exp(-de/ T)/(2*vp);///vp; // z* K* K);
				//cout << " crit is " << crit << endl;
				if (true)
				{ //gsl_rng_uniform(r) < crit) {
					
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

		else if (is_boundary(heid2) > 0)
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
	//return -1;*/
	return 1;
}

void Geometry::move_v_epsilon(double len_move, double *pi, double *pf, gsl_rng *r)
{
	double *tempvec, *fvec;
	tempvec = new double[3];
	fvec = new double[3];

	//random direction
	randvec(tempvec, r); // random direction with unit length
	//double sigmav=.2;
	//double d = pow(xi,(1/3)) * gsl_rng_uniform(r); // length of change between 0 and xi
	double d = len_move* cbrt(gsl_rng_uniform(r)) * 2.0 / sqrt(epsilon[0]);
	//double d = xi * gsl_ran_gaussian(r, sigmav);
	//cout << " len d in move_v_epsilon " << d << endl;
	multvec(tempvec, d, fvec); // update the vector length
	addvec(fvec, pi, pf);	   // apply the position change

	delete[] tempvec;
	delete[] fvec;
}

double Geometry::move_p(double *pi, double *pf, gsl_rng *r)
{
	double *x;//tempvec, *fvec;
	//tempvec = new double[3];
	x = new double[3];

	//random direction
	//randvec(tempvec, r); // random direction with unit length
	//double sigmav=.2;
	//double x = gsl_rng_uniform(r);
	//cout << " x is " << x << endl;
	//cout << " x^1/3 is" << cbrt(x)<<endl<<endl;
	//double d = len_v * cbrt(x); // length of change between 0 and xi
	//cout << " len d in move_p " << d << endl;
	//double d = xi * gsl_ran_gaussian(r, sigmav);
	//multvec(tempvec, d, fvec); // update the vector length

	//x[0] = l_thermal_sigma* (gsl_rng_uniform(r)-.5);
	//x[1] = l_thermal_sigma* (gsl_rng_uniform(r)-.5);

	x[0] =  l_thermal_kappa* (gsl_rng_uniform(r)-.5);
	x[1] =  l_thermal_kappa* (gsl_rng_uniform(r)-.5);


	//this is theta

	x[2] =  l_thermal_kappa* (gsl_rng_uniform(r)-.5);

	double d=norm(x);

	addvec(x, pi, pf);	   // apply the position change

	//delete[] tempvec;
	delete[] x;
	
	return d;
}
double Geometry::move_p_gaussian(double len_v, double *pi, double *pf, gsl_rng *r)
{
	double *x;
	x = new double[3];

	double sig = gaussian_sigma;
	//cout << " sig is " <<sig<<endl;

	x[0] = gsl_ran_gaussian(r, sig);
	x[1] = gsl_ran_gaussian(r, sig);
	x[2] = gsl_ran_gaussian(r, sig);

	double d = norm(x);
	addvec(pi, x, pf); // apply the position change
	delete[] x;

	return d;
}

//find the projection of the vector in direction of 3 rotating axes
// use this for gaussian correction
double Geometry::find_project_dist_axes(int heindex0, double *newv, double *oldv , double *dis_vector_project )
{

	double *fvecx, *fvecy, *fvecz, *tempvec_dis , *tempvec;
	fvecx = new double[3];
	fvecy = new double[3];
	fvecz = new double[3];
	tempvec_dis = new double[3];
	tempvec = new double[3];

	subvec(oldv,newv,tempvec_dis);
	cout << "tempvec " <<  tempvec_dis[0] <<" "<< tempvec_dis[1] << " "<<tempvec_dis[2]<<endl;
	//direction of x hevec 
	multvec(he[heindex0].hevec, 1/norm(he[heindex0].hevec),fvecx);
	dis_vector_project[0]=dot(tempvec_dis,fvecx);

	//direction of y pi -hecent
	subvec(he[heindex0].hecent,oldv,tempvec);
	multvec(tempvec, 1/norm(tempvec),fvecy);
	dis_vector_project[1]=dot(tempvec_dis,fvecy);

	//derection of z cross(hevec, y)
	cross(he[heindex0].hevec, fvecy, tempvec); 
	multvec(tempvec, 1/norm(tempvec),fvecz);
	dis_vector_project[2]=dot(tempvec_dis,fvecy);

	double d= norm(tempvec_dis);
	delete[] fvecx;
	delete[] fvecy;
	delete[] fvecz;
	delete[] tempvec_dis;
	delete[] tempvec;
	

	return(d);
}


double Geometry::move_p_rotate_axes(int heindex0,double len_v, double *pi, double *pf, gsl_rng *r, double *dis_vector)
{
	// this is incomplete XXX
	double *tempvec,*fvecx, *fvecy, *fvecz; // , *p;

	double *x;
	x = new double[3];

	x[0] = l_thermal_sigma* (gsl_rng_uniform(r)-.5);
	x[1] = l_thermal_sigma* (gsl_rng_uniform(r)-.5);

	//this is theta

	x[2] = theta_thermal_kappa* (gsl_rng_uniform(r)-.5);


	multvec(x,1,dis_vector);
	double d = norm(x);


	tempvec = new double[3];
	fvecx = new double[3];
	fvecy = new double[3];
	fvecz = new double[3];
		
	//direction of x hevec 
	multvec(he[heindex0].hevec, x[0]/norm(he[heindex0].hevec),fvecx);
	

	//direction of y pi -hecent
	subvec(he[heindex0].hecent,pi,tempvec);
	multvec(tempvec, x[1]/norm(tempvec),fvecy);

	//derection of z cross(hevec, y)
	cross(he[heindex0].hevec, fvecy, tempvec); 
	multvec(tempvec, x[2]/norm(tempvec),fvecz);

	// add to hecent

	tempvec[0]=fvecx[0]+fvecy[0];
	tempvec[1]=fvecx[1]+fvecy[1];
	tempvec[2]=fvecx[2]+fvecy[2];

	// apply the position change
	addvec(pi, tempvec, pf); 

	// now rotate it 

	//cout<< "in new_vertex newv is " <<newv[0] <<" "<<newv[1] <<" " <<newv[2] << endl;
	
	delete[] tempvec;
	delete[] fvecx;
	delete[] fvecy;
	delete[] fvecz;
	delete[] x;

	return d;

}



double Geometry::move_p_gaussian_axes(int heindex0,double len_v, double *pi, double *pf, gsl_rng *r, double *dis_vector)
{

	double *tempvec,*fvecx, *fvecy, *fvecz; // , *p;

	double *x;
	x = new double[3];

	//double sig = gaussian_sigma;
	//cout << " sig is " <<sig<<endl;

	double sig = gaussian_sigma;
	x[0] = gsl_ran_gaussian(r, sig);
	x[1] = gsl_ran_gaussian(r, sig);

	x[2] = gsl_ran_gaussian(r, sig);


	 // we send this back for claculation of vp

	//for now keep it a box ! XXX

	/*x[0] = l_thermal_sigma* (gsl_rng_uniform(r)-.5);
	x[1] = l_thermal_sigma* (gsl_rng_uniform(r)-.5);
	x[2] = l_thermal_kappa* (gsl_rng_uniform(r)-.5);*/


	multvec(x,1,dis_vector);
	double d = norm(x);


	tempvec = new double[3];
	fvecx = new double[3];
	fvecy = new double[3];
	fvecz = new double[3];
		
	//direction of x hevec 
	multvec(he[heindex0].hevec, x[0]/norm(he[heindex0].hevec),fvecx);
	

	//direction of y pi -hecent
	subvec(he[heindex0].hecent,pi,tempvec);
	multvec(tempvec, x[1]/norm(tempvec),fvecy);

	//derection of z cross(hevec, y)
	cross(he[heindex0].hevec, fvecy, tempvec); 
	multvec(tempvec, x[2]/norm(tempvec),fvecz);

	// add to hecent

	tempvec[0]=fvecx[0]+fvecy[0]+fvecz[0];
	tempvec[1]=fvecx[1]+fvecy[1]+fvecz[1];
	tempvec[2]=fvecx[2]+fvecy[2]+fvecz[2];

	// apply the position change
	addvec(pi, tempvec, pf); 
	//cout<< "in new_vertex newv is " <<newv[0] <<" "<<newv[1] <<" " <<newv[2] << endl;
	
	delete[] tempvec;
	delete[] fvecx;
	delete[] fvecy;
	delete[] fvecz;
	delete[] x;

	return d;

}



double Geometry::new_vertex_edge_and_move(int heindex0, double *newv, int etnew,gsl_rng *r)
{
	//cout << " in new_vertex_edge_and_move heindex0 is " << heindex0<<endl;
	double *tempvec, *fvec , *fvecx , *fvecy; 

	tempvec = new double[3];
	fvec = new double[3];
	fvecx = new double[3];
	fvecy = new double[3];


	double *x;
	x = new double[3];

	x[0] = 2.0 * l_thermal_sigma* (gsl_rng_uniform(r)-0.5);
	x[1] = 2.0 * l_thermal_sigma* (gsl_rng_uniform(r)-0.5);

	//cout << "gsl_rng_uniform(r)-0.5 " << gsl_rng_uniform(r)-0.5 <<endl;
	//cout << "l_thermal_sigma " << l_thermal_sigma << endl;
	

	//this is theta

	x[2] = 2.0 * theta_thermal_kappa* (gsl_rng_uniform(r)-0.5);

	//cout << "X " << x[0] << " " << x[1] << " " << x[2] <<endl;

	int heopindex0 = heidtoindex[he[heindex0].opid];
	//cout << "heopindex0 is  "<< heopindex0 <<endl;
	get_normal(he[heopindex0].id);
	int et = he[heindex0].type;
	
	
	//cout <<"heopindex0 " <<heopindex0 <<endl;
	// find the angle

	/*g.phi0[0] = 1.05;
    g.phi0[1] = 1.17;
    g.phi0[2] = .98;
    g.phi0[3] = 1.05;*/
	double angle=1.05;
	if (et==1 && etnew==2) angle=phi0[1]; // 1 -> 2 T3 and T4 (5fold vertex)
	else if (et==2 && etnew==0 ) angle=phi0[2]; // 2 -> 0 T3 and T4
	else if (et==0 && etnew==1 ) angle=phi0[2]; // 0 -> 1 T3 and T4
	else if (et==3 && etnew==1 ) angle=phi0[2]; // 3 -> 1 T3 
	else if (et==2 && etnew==3 ) angle=phi0[2]; // 2 -> 3 T3 

	// XXX should consider 0 -> 1 ? T3 ?
	// the rest are considered 60 degrees

	// hevec is x direction to find the vector 
	// - sign is because of considering this edge and next edge
	double ll=l0[et];

	multvec(he[heindex0].hevec,(ll*(1+x[0])*(-cos(angle)))/norm(he[heindex0].hevec),fvecx);
	//cout << " he[heindex0].hevec " << he[heindex0].hevec[0] <<" "<< he[heindex0].hevec[1] << he[heindex0].hevec[2] <<endl;
	//cout << "ll " << ll <<endl;
	//cout << "x[0] " << x[0] <<endl;
	//cout << " -cos(angle) " << -cos(angle) <<endl; 
	//cout << "ll*(1+x[0])*(-cos(angle)))/norm(he[heindex0].hevec) " << ll*(1+x[0])*(-cos(angle))/norm(he[heindex0].hevec) <<endl;
	//cout<< " fvecx" << fvecx[0] <<" " << fvecx[1] <<" " << fvecx[2] <<" " << endl;
	
	// tempvec is y direction to find the vector
	cross(he[heopindex0].n, he[heindex0].hevec, tempvec); 
	
	//cout<< " tempvec" << tempvec[0] <<" " << tempvec[1] <<" " << tempvec[2] <<" " << endl;
	//cout<< " he[heopindex0].n" << he[heopindex0].n[0] <<" " << he[heopindex0].n[1] <<" " << he[heopindex0].n[2] <<" " << endl;
	if (he[heopindex0].n[0]==0 && he[heopindex0].n[1]==0 && he[heopindex0].n[2]==0) exit(-1);
	multvec(tempvec,ll*sin(angle)*(1+x[1])/norm(tempvec),fvecy);
	//cout<< " fvecy" << fvecy[0] <<" " << fvecy[1] <<" " << fvecy[2] <<" " << endl;

	addvec(fvecx,fvecy,fvec);

	//cout << "len new edge is " << norm(fvec) <<endl;
	//exit(-1);
 
	// find the rotate angle

	/*gg.theta0[0] = atof(argv[5]);
    g.theta0[1] = atof(argv[6]);
    g.theta0[2] = g.theta0[0];   // CD-CD
    g.theta0[3] = g.theta0[0];  */


	double theta = theta0[0]; // .24 for all angles
	if (et==1 || et==2) theta=theta0[1];
	rotatevec(fvec, he[heopindex0].hevec, theta*(1+x[2]), tempvec); 
	
	addvec(v[vidtoindex[he[heindex0].vout]].co,tempvec,newv);

	
	//cout << "theta" << theta <<endl;

	delete[] tempvec;
	delete[] fvec;
	delete[] fvecx;
	delete[] fvecy;
	delete[] x;

	double d= sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]*.86*.86);
	return(d);
	
}


void Geometry::new_vertex_edge(int heindex0, double *newv, int etnew)
{
	
	double *tempvec, *fvec , *fvecx , *fvecy; 

	tempvec = new double[3];
	fvec = new double[3];
	fvecx = new double[3];
	fvecy = new double[3];

	int heopindex0 = heidtoindex[he[heindex0].opid];
	get_normal(he[heindex0].id);
	int et = he[heindex0].type;
	
	//cout <<"heopindex0 " <<heopindex0 <<endl;
	// find the angle

	/*g.phi0[0] = 1.05;
    g.phi0[1] = 1.17;
    g.phi0[2] = .98;
    g.phi0[3] = 1.05;*/
	double angle=1.05;
	if (et==1 && etnew==2) angle=phi0[1]; // 1 -> 2 T3 and T4 (5fold vertex)
	else if (et==2 && etnew==0 ) angle=phi0[2]; // 2 -> 0 T3 and T4
	else if (et==0 && etnew==1 ) angle=phi0[2]; // 0 -> 1 T3 and T4
	else if (et==3 && etnew==1 ) angle=phi0[2]; // 3 -> 1 T3 
	else if (et==2 && etnew==3 ) angle=phi0[2]; // 2 -> 3 T3 

	// XXX should consider 0 -> 1 ? T3 ?
	// the rest are considered 60 degrees

	// hevec is x direction to find the vector 
	// - sign is because of considering this edge and next edge

	multvec(he[heindex0].hevec,-cos(angle)*l0[et]/norm(he[heindex0].hevec),fvecx);
	//cout<< " fvecy" << fvecx[0] <<" " << fvecx[1] <<" " << fvecx[2] <<" " << endl;
	
	// tempvec is y direction to find the vector
	cross(he[heopindex0].n, he[heindex0].hevec, tempvec); 
	
	//cout<< " tempvec" << tempvec[0] <<" " << tempvec[1] <<" " << tempvec[2] <<" " << endl;
	//cout<< " he[heopindex0].n" << he[heopindex0].n[0] <<" " << he[heopindex0].n[1] <<" " << he[heopindex0].n[2] <<" " << endl;
	if (he[heopindex0].n[0]==0 && he[heopindex0].n[1]==0 && he[heopindex0].n[2]==0) exit(-1);
	multvec(tempvec,sin(angle)*l0[et]/norm(tempvec),fvecy);
	//cout<< " fvecy" << fvecy[0] <<" " << fvecy[1] <<" " << fvecy[2] <<" " << endl;

	addvec(fvecx,fvecy,fvec);

	// find the rotate angle

	/*gg.theta0[0] = atof(argv[5]);
    g.theta0[1] = atof(argv[6]);
    g.theta0[2] = g.theta0[0];   // CD-CD
    g.theta0[3] = g.theta0[0];  */


	double theta = theta0[0]; // .24 for all angles
	if (et==1 || et==2) theta=theta0[1];
	rotatevec(fvec, he[heopindex0].hevec, theta, tempvec); 
	
	addvec(v[vidtoindex[he[heindex0].vout]].co,tempvec,newv);

	

	delete[] tempvec;
	delete[] fvec;
	delete[] fvecx;
	delete[] fvecy;
	
}

void Geometry::new_vertex(int heindex0, double *newv)
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
	//cout<< "in new_vertex newv is " <<newv[0] <<" "<<newv[1] <<" " <<newv[2] << endl;
	delete[] tempvec;
	delete[] fvec;
	//delete[] p;
}


void Geometry::new_vertex_points(int heindex0, int ind ,double *newv)
{

	double *tempvec,*fvecx, *fvecy, *fvecz; // , *p;

	tempvec = new double[3];
	fvecx = new double[3];
	fvecy = new double[3];
	fvecz = new double[3];
	
	int heopindex0 = heidtoindex[he[heindex0].opid];
	
	//direction of x hevec 
	multvec(he[heindex0].hevec, dist_points[ind][0],fvecx);
	

	//direction of y tempvec

	cross(he[heopindex0].n, he[heindex0].hevec, tempvec); 
	multvec(tempvec, dist_points[ind][1],fvecy);

	//derection of z n
	multvec(he[heindex0].n, dist_points[ind][2],fvecz);

	// add to hecent

	tempvec[0]=fvecx[0]+fvecy[0]+fvecz[0];
	tempvec[1]=fvecx[1]+fvecy[1]+fvecz[1];
	tempvec[2]=fvecx[2]+fvecy[2]+fvecz[2];

	addvec( tempvec , he[heopindex0].hecent, newv); 

	
	cout<< "in new_vertex newv is " <<newv[0] <<" "<<newv[1] <<" " <<newv[2] << endl;
	delete[] tempvec;
	delete[] fvecx;
	delete[] fvecy;
	delete[] fvecz;

}


void Geometry::hecenter(int heindex, double *vcenter)
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
void Geometry::helen(int heindex, double *helen)
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

int Geometry::get_normal(int heid0)
{
	int heindex0 = heidtoindex[heid0];
	//cout << "in get_normal" <<endl;
	update_half_edge(heid0);
	//int etype=he[heindex0].type;
	if (he[heindex0].nextid == -1 || he[heindex0].previd == -1)
	{
		he[heindex0].n[0] = 0;
		he[heindex0].n[1] = 0;
		he[heindex0].n[2] = 0;
		//cout << " on boundary no normal! " <<endl;
		return (-1);
	}

	double *vec2 = new double[3];
	if (he[heindex0].nextid != -1)
	{
		update_half_edge(he[heindex0].nextid);
		int nextindex0 = heidtoindex[he[heindex0].nextid];
		cross(he[heindex0].hevec, he[nextindex0].hevec, vec2);
	}
	else
	{
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

void Geometry::update_excluder_top_he(int heid)
{
	double *vtemp = new double[3];
	double *vtemp2 = new double[3];

	int heindex0 = heidtoindex[heid];
	int heindexop = heidtoindex[he[heindex0].opid];

	if (is_boundary(heid) > 0)
	{
		multvec(he[heindexop].n, .425, vtemp);
		addvec(he[heindex0].hecent, vtemp, he[heindex0].hetop);
	}
	else if (is_boundary(he[heindex0].opid) > 0)
	{
		multvec(he[heindex0].n, .425, vtemp);
		addvec(he[heindex0].hecent, vtemp, he[heindex0].hetop);
	}
	else
	{
		addvec(he[heindex0].n, he[heindexop].n, vtemp2); //averaging the two normals
		multvec(vtemp2, .2125, vtemp);
		addvec(he[heindex0].hecent, vtemp, he[heindex0].hetop);
	}

	delete[] vtemp;
	delete[] vtemp2;
}

void Geometry::update_excluder_top()
{

	for (vector<HE>::iterator it = he.begin(); it != he.end(); ++it)
	{
		update_excluder_top_he(it->id);
	}
}

void Geometry::update_normals()
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

void Geometry::update_normals_vertex(int vindex0)
{
	for (vector<int>::iterator ithe = v[vindex0].hein.begin(); ithe != v[vindex0].hein.end(); ithe++)
	{
		get_normal(*ithe);
		int heindex = heidtoindex[*ithe];
		get_normal(he[heindex].opid);
		if (he[heindex].previd != -1)
		{
			get_normal(he[heindex].previd);
			get_normal(he[heidtoindex[he[heindex].previd]].opid);
		}
		if (he[heindex].nextid != -1)
		{
			get_normal(he[heindex].nextid);
			get_normal(he[heidtoindex[he[heindex].nextid]].opid);
		}

		if (is_boundary(he[heindex].opid) > 0)
		{
			int opindex = heidtoindex[he[heindex].opid];
			if (he[opindex].previd != -1)
			{
				get_normal(he[opindex].previd);
				get_normal(he[heidtoindex[he[opindex].previd]].opid);
			}
			if (he[opindex].nextid != -1)
			{
				get_normal(he[opindex].nextid);
				get_normal(he[heidtoindex[he[opindex].nextid]].opid);
			}
		}
	}
}

void Geometry::update_geometry_vertex(int vindex0)
{
	//cout << "in update geometry vertex vindex0 " <<vindex0 <<endl;
	for (vector<int>::iterator ithe = v[vindex0].hein.begin(); ithe != v[vindex0].hein.end(); ithe++)
	{

		int heindex = heidtoindex[*ithe];

		update_half_edge(*ithe);
		update_half_edge(he[heindex].opid);
		if (he[heindex].previd != -1)
		{
			update_half_edge(he[heindex].previd);
			update_half_edge(he[heidtoindex[he[heindex].previd]].opid);
		}

		if (he[heindex].nextid != -1)
		{
			update_half_edge(he[heindex].nextid);
			update_half_edge(he[heidtoindex[he[heindex].nextid]].opid);
		}

		if (is_boundary(he[heindex].opid) > 0)
		{
			int opindex = heidtoindex[he[heindex].opid];
			if (he[opindex].previd != -1)
			{
				update_half_edge(he[opindex].previd);
				update_half_edge(he[heidtoindex[he[opindex].previd]].opid);
			}
			if (he[opindex].nextid != -1)
			{
				update_half_edge(he[opindex].nextid);
				update_half_edge(he[heidtoindex[he[opindex].nextid]].opid);
			}
		}
	}
}

void Geometry::update_excluder_top_vertex(int vindex0)
{

	for (vector<int>::iterator ithe = v[vindex0].hein.begin(); ithe != v[vindex0].hein.end(); ithe++)
	{
		update_excluder_top_he(*ithe);
	}
}

int Geometry::check_inside_overlap(int heid0)
{
	if (Nhe <= 180)
		return 1;
	int heindex0 = heidtoindex[heid0];
	int opindex = heidtoindex[he[heindex0].opid];
	int vi1 = he[opindex].vin;
	int vi2 = he[opindex].vout;
	int nextopid = he[opindex].nextid;
	int vi3 = he[heidtoindex[nextopid]].vout;
	if ((vi3 == vi1) || (vi3 == vi2))
	{
		cout << " wrong in inside !" << endl;
		exit(-1);
	}
	double xyz0[3];

	double XCM = 0;
	double YCM = 0;
	double ZCM = 0;
	// find CM
	/*for (vector<VTX>::iterator it = v.begin(); it != v.end(); ++it)
	{    
		XCM +=it->co[0];
		YCM +=it->co[1];
		ZCM +=it->co[2];

	}*/

	xyz0[0] = XCM / Nv;
	xyz0[1] = YCM / Nv;
	xyz0[2] = ZCM / Nv;
	//int nooverlap=0;
	for (vector<HE>::iterator it = he.begin(); it != he.end(); it++)
	{

		if ((it->id != heid0) && (it->opid != heid0))
		{
			//cout << "checkoverlap "<< it->id << " and " << heid0 <<endl;
			int vj1 = it->vin;
			int vj2 = it->vout;

			if (vi1 != vj1 && vi1 != vj2 && vi2 != vj1 && vi2 != vj2 && vi3 != vj1 && vi3 != vj2)
			{

				if (tri_tri_overlap_test_3d(v[vidtoindex[vi1]].co, v[vidtoindex[vi2]].co, v[vidtoindex[vi3]].co,
											v[vidtoindex[vj1]].co, v[vidtoindex[vj2]].co, xyz0) > 0)
				{
					return -1;
				}

				// double check ovelap ?? TEMP 1022

				int tempnextid = it->nextid;
				if (tempnextid == -1)
				{
					tempnextid = he[heidtoindex[it->opid]].nextid;
				}
				if (tempnextid == -1)
				{
					cout << " wrong geometry! inside overlap" << endl;
				}
				if (tri_tri_overlap_test_3d(v[vidtoindex[vi1]].co, v[vidtoindex[vi2]].co, v[vidtoindex[vi3]].co,
											it->hecent, he[heidtoindex[tempnextid]].hecent, xyz0) > 0)
				{
					return -1;
				}
				//else
				//{
				//	nooverlap++;
				//}
			}
			//else{
			//	cout << "vi1 " << vi1 <<" vi2 " << vi2  <<" vi3 " << vi3 << " vj1 " << vj1<< " vj2 " << vj2<< endl;
			//}
		}
	}

	//cout << "nooverlapInsideTest " << nooverlap << " - Nhe " << Nhe << " -- " <<endl;
	return 1;
}

int Geometry::do_intersect(int heid1, int heid2)
{
	if (is_boundary(heid1) > 0 || is_boundary(heid2) > 0)
	{
		return 0;
	}
	//int n_vertices_shared=shared_vertices( heid1,  heid2);
	//if (n_vertices_shared==0){
	int vi1 = he[heidtoindex[heid1]].vin;
	int vi2 = he[heidtoindex[heid1]].vout;
	int nextid1 = he[heidtoindex[heid1]].nextid;
	int vi3 = he[heidtoindex[nextid1]].vout;
	int vj1 = he[heidtoindex[heid2]].vin;
	int vj2 = he[heidtoindex[heid2]].vout;
	int nextid2 = he[heidtoindex[heid2]].nextid;
	int vj3 = he[heidtoindex[nextid2]].vout;

	if (vi1 == vj1 || vi1 == vj2 || vi1 == vj3)
	{
		return 0;
	}
	if (vi2 == vj1 || vi2 == vj2 || vi2 == vj3)
	{
		return 0;
	}
	if (vi3 == vj1 || vi3 == vj2 || vi3 == vj3)
	{
		return 0;
	}

	int overlap = 0;

	overlap = tri_tri_overlap_test_3d(v[vidtoindex[vi1]].co, v[vidtoindex[vi2]].co, v[vidtoindex[vi3]].co,
									  v[vidtoindex[vj1]].co, v[vidtoindex[vj2]].co, v[vidtoindex[vj3]].co);

	if (overlap == 1)
		return overlap;

	/* check intersect with top triangle */

	overlap = tri_tri_overlap_test_3d(v[vidtoindex[vi1]].co, v[vidtoindex[vi2]].co, he[heidtoindex[heid1]].hetop,
									  v[vidtoindex[vj1]].co, v[vidtoindex[vj2]].co, he[heidtoindex[heid2]].hetop);

	if (overlap == 1)
		return overlap;

	overlap = tri_tri_overlap_test_3d(v[vidtoindex[vi1]].co, v[vidtoindex[vi2]].co, v[vidtoindex[vi3]].co,
									  v[vidtoindex[vj1]].co, v[vidtoindex[vj2]].co, he[heidtoindex[heid2]].hetop);

	if (overlap == 1)
		return overlap;

	overlap = tri_tri_overlap_test_3d(v[vidtoindex[vi1]].co, v[vidtoindex[vi2]].co, he[heidtoindex[heid1]].hetop,
									  v[vidtoindex[vj1]].co, v[vidtoindex[vj2]].co, v[vidtoindex[vj3]].co);

	return overlap;
}

int Geometry::find_overlap_all()
{
	for (vector<VTX>::iterator it = v.begin(); it != v.end(); it++)
	{
		check_overlap_g(it->vid);
	}
	return 0;
}

int Geometry::find_overlap_g(int vid0)
{
	if (Nhe < 22)
	{
		return 1;
	}
	//cout << "in _find_overlap_g" <<endl;
	int vindex0 = vidtoindex[vid0];
	double mindis = .2;
	for (vector<int>::iterator itv = v[vindex0].vneigh.begin(); itv != v[vindex0].vneigh.end(); itv++)
	{

		int vindex1 = vidtoindex[*itv]; //neighbor v
		if (vindex1 != -1)
		{
			//if (veclen(v[vindex0].co,v[vindex1].co )<.1) { return -1;} // check vertex vertx overlap

			//if (is_vboundary(*itv)>0 && is_vboundary(vid0)>0 ) {mindis=.1; }
			//else { mindis=.2;}
			if (veclen(v[vindex1].co, v[vindex0].co) < l0[0] * mindis)
			{
				cout << "In find_overlap overlap vindex0 vindex1 " << vindex0 << " " << vindex1 << endl;
				exit(-1);
			}

			for (vector<int>::iterator itneigh = v[vindex1].hein.begin(); itneigh != v[vindex1].hein.end(); itneigh++) //neigbor vertex hein s
			{
				int heid1 = *itneigh;

				if (is_boundary(heid1) > 0)
				{
					heid1 = he[heidtoindex[heid1]].opid;
				}
				int heindex1 = heidtoindex[heid1];
				for (vector<int>::iterator it = v[vindex0].hein.begin(); it != v[vindex0].hein.end(); it++) //this vertex hein s
				{
					int heid0 = *it;
					int cH = connectedH(heid0, heid1);
					int cvin = connected(vid0, he[heid0].vin);
					int cvout = connected(vid0, he[heid0].vout);
					if (cH == 0 && cvin == 0 && cvout == 0 && heid0 != he[heindex1].id && heid0 != he[heindex1].opid)
					{
						if (is_boundary(heid1) > 0)
						{
							heid0 = he[heidtoindex[heid1]].opid;
						}
						int heindex0 = heidtoindex[heid0];
						if (veclen(v[vindex0].co, he[heindex1].hecent) < l0[0] * .2)
						{
							cout << "In find_overlap overlap heindex1 vindex0 " << heindex1 << " " << vindex0 << endl;
							exit(-1);
						}
						if (veclen(v[vindex1].co, he[heindex0].hecent) < l0[0] * .2)
						{
							cout << "In find_overlap overlap heindex0 vindex1 " << heindex0 << " " << vindex1 << endl;
							exit(-1);
						}

						if (veclen(he[heindex0].hecent, he[heindex1].hecent) < l0[0] * .2)
						{
							cout << "In find_overlap overlap heindex0 heindex1 " << heindex0 << " " << heindex1 << endl;
							exit(-1);
						}
						if (do_intersect(heid0, heid1) > 0)
						{
							cout << "In find_overlap overlap triangles of  " << heindex0 << " " << heindex1 << endl;
							exit(-1);
						}
					}
				}
			}
		}
	}
	return 1;
}

int Geometry::check_overlap_g(int vid0)
{
	if (Nhe < 22)
	{
		return 1;
	}
	//cout << "in _check_overlap_g" <<endl;
	int vindex0 = vidtoindex[vid0];
	double mindis = .2;

	//check overlap with neighbors
	for (vector<int>::iterator itv = v[vindex0].vneigh.begin(); itv != v[vindex0].vneigh.end(); itv++)
	{
		//cout<< " "
		int vindex1 = vidtoindex[*itv]; //neighbor v
		
		if (vindex1==-1) {
			
			cout << "in checkoverlap wrong neigh updating neigh!"<<endl;
			//exit(-1);
			//update_index();

			update_neigh_vertex(vid0);
			update_neigh_vertex_and_neigh(vid0);
			//update_index();
			//update_neigh();
			break;
			//update_neigh_vertex_and_neigh(*itv);
		}
	}
	for (vector<int>::iterator itv = v[vindex0].vneigh.begin(); itv != v[vindex0].vneigh.end(); itv++)
	{
		//cout<< " "
		int vindex1 = vidtoindex[*itv]; //neighbor v
		
		if (vindex1==-1) {
			
			cout << "in checkoverlap wrong neigh !"<<endl;
			exit(-1);
			
		}
		else
		{
			//if (veclen(v[vindex0].co,v[vindex1].co )<.1) { return -1;} // check vertex vertx overlap

			//if (is_vboundary(*itv)>0 && is_vboundary(vid0)>0 ) {mindis=.1; }
			//else { mindis=.2;}
			if (veclen(v[vindex1].co, v[vindex0].co) < l0[0] * mindis) // vertex vertex overlap
			{
				return -1;
			}

			for (vector<int>::iterator itneigh = v[vindex1].hein.begin(); itneigh != v[vindex1].hein.end(); itneigh++) //neigbor vertex hein s
			{
				int heid1 = *itneigh;

				//if (is_boundary(heid1)>0){ heid1=he[heidtoindex[heid1]].opid;} //why?
				int heindex1 = heidtoindex[heid1];
				/* 10/29 add spike ecluded volume */
				for (vector<int>::iterator it = v[vindex0].hein.begin(); it != v[vindex0].hein.end(); it++) //this vertex hein s
				{
					int heid0 = *it;
					if (heid0 != he[heindex1].id && heid0 != he[heindex1].opid)
					{
						//if (is_boundary(heid0)>0){ heid0=he[heidtoindex[heid1]].opid;} //why?
						int heindex0 = heidtoindex[heid0];

						if (veclen(v[vindex0].co, he[heindex1].hetop) < l0[0] * .2) 
						{
							return -1;
						} // hecent -vertex overlap
						if (veclen(v[vindex1].co, he[heindex0].hetop) < l0[0] * .2)
						{
							return -1;
						} // hecent -vertex overlap

						if (veclen(v[vindex0].co, he[heindex1].hecent) < l0[0] * .2)
						{
							return -1;
						} // hecent -vertex overlap
						if (veclen(v[vindex1].co, he[heindex0].hecent) < l0[0] * .2)
						{
							return -1;
						} // hecent -vertex overlap
						float mindishe = .2;
						//if (((is_boundary(heid1)>0) || (is_boundary(he[heidtoindex[heid1]].opid)>0) ) && ((is_boundary(heid0)>0) || (is_boundary(he[heidtoindex[heid0]].opid)>0))  && connectedH(heid1,heid0)>0){ mindishe=.25;}
						if (veclen(he[heindex0].hecent, he[heindex1].hecent) < l0[0] * mindishe)
						{
							return -1;
						} // hecent - hecent overlap
						if (veclen(he[heindex0].hetop, he[heindex1].hecent) < l0[0] * .2)
						{
							return -1;
						} // hetop - hecent overlap
						if (veclen(he[heindex0].hecent, he[heindex1].hetop) < l0[0] * .2)
						{
							return -1;
						} // hecent - hetop overlap

						if (do_intersect(heid0, heid1) > 0)
							return -1;
						
					}
				}
			}
		}
		
	}
	return 1;
}

int Geometry::check_overlap_centerv(double *newv) /* checks this v point with all other hecenters and vertices*/
{
	if (Nhe <= 6)
		return 1;

	for (vector<HE>::iterator it = he.begin(); it != he.end(); it++)
	{
		double mindisvh = .2;
		if (veclen(it->hecent, newv) < l0[0] * mindisvh)
		{
			return -1;
		}
		if (veclen(it->hetop, newv) < l0[0] * mindisvh)
		{
			return -1;
		}
	}
	for (vector<VTX>::iterator it = v.begin(); it != v.end(); it++)
	{
		double mindis = .2;
		if (veclen(it->co, newv) < l0[0] * mindis)
		{
			return -1;
		}
	}

	return 1;
}

int Geometry::check_overlap_centerh(double *newcenter) /* checks this he center point with all other hecenters*/
{
	if (Nhe <= 6)
		return 1;

	for (vector<HE>::iterator it = he.begin(); it != he.end(); it++)
	{
		double mindis = .25;
		//if (is_boundary(it->id)>0) mindis=.25;
		if (veclen(it->hecent, newcenter) < l0[0] * mindis)
		{
			return -1;
		}
		if (veclen(it->hetop, newcenter) < l0[0] * .2)
		{
			return -1;
		}
		double mindisvh = .2;
		if (veclen(v[vidtoindex[it->vin]].co, newcenter) < l0[0] * mindisvh)
		{
			return -1;
		}
	}

	return 1;
}

/*ENERGY HELPER FUNCTIONS */

double Geometry::find_gbb(int etypenew, int etypenextofnew, int etypeprevifnew)
{
	double e1 = find_dg(etypenew, etypenextofnew, 0);
	double e2 = find_dg(etypenextofnew, etypeprevifnew, 0);
	double e3 = find_dg(etypeprevifnew, etypenew, 0);
	//cout << "e1 is " << e1 << " e2 is " << e2 << " e3 is " <<e3 <<endl;
	return e1 + e2 + e3;
}

double Geometry::find_dg(int type, int typenext, bool drug)
{

	/*/if ((typenext == 1 || typenext == 2) and (drug != 0))
	{
		cout << "why drug on 1 2 " << endl;
		exit(-1);
	}*/

	double bindg = gb[type][typenext] + drug * (gdrug[type][typenext] - mudrug);

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

double Geometry::stretch_energy(int heindex0)
{
	int et = he[heindex0].type;

	//cout<<  "epsilon[et]" << epsilon[et]  <<endl;
	//cout<<  "(he[heindex0].l " << he[heindex0].l <<endl ;
	//cout <<"stretch energy is " <<  0.5 * epsilon[et] * (he[heindex0].l - l0[et]) * (he[heindex0].l - l0[et]) <<endl;
	return 0.5 * epsilon[et] * (he[heindex0].l - l0[et]) * (he[heindex0].l - l0[et]);
}

int Geometry::check_bind_wedge(int heid0)
{
	int heindex0 = heidtoindex[heid0];
	get_normal(he[heindex0].id);
	get_normal(he[heindex0].opid);
	int opindex0 = heidtoindex[he[heindex0].opid];
	int nextindex0 = heidtoindex[he[heindex0].nextid];
	//int previndex0 = heidtoindex[he[heindex0].previd];
	int opnextindex0 = heidtoindex[he[opindex0].nextid];

	double *tempvec = new double[3];
	subvec(v[vidtoindex[he[nextindex0].vout]].co, v[vidtoindex[he[opnextindex0].vout]].co, tempvec);

	//double ndot = dot(he[opindex0].n, he[heindex0].n);
	if (dot(tempvec, he[heindex0].n) > 0)
	{
		delete[] tempvec;
		return -1;
	}
	delete[] tempvec;
	return 1;
}

double Geometry::bend_energy(int heindex0)
{
	//cout << " \n\n     normal " << he[heindex0].n[0] << " " << he[heindex0].n[1] << " " << he[heindex0].n[2] << " " <<endl ;
	//cout << "\n in bend_energy edge index is "  << heindex0<< endl;
	int opindex0 = heidtoindex[he[heindex0].opid];
	if (he[heindex0].nextid == -1 || he[heindex0].previd == -1)
	{
		return 0;
	}
	if (he[opindex0].nextid == -1 || he[opindex0].previd == -1)
	{
		return 0;
	}
	//if (is_boundary(he[heindex0].id)>0 || is_boundary(he[heindex0].opid)>0)  { return 0;}
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
	int vid0 = -1;
	int opvid0 = -1;
	if (he[heindex0].nextid != -1)
	{
		nextindex0 = heidtoindex[he[heindex0].nextid];
		nexttype = he[nextindex0].type;
		vid0 = he[nextindex0].vout;
	}

	if (he[heindex0].previd != -1)
	{
		previndex0 = heidtoindex[he[heindex0].previd];
		prevtype = he[previndex0].type;
		vid0 = he[previndex0].vin;
	}

	int opnextindex0 = -1;
	int opprevindex0 = -1;
	int opnexttype = -1;
	int opprevtype = -1;
	int opetype = he[opindex0].type;

	if (he[opindex0].nextid != -1)
	{
		opnextindex0 = heidtoindex[he[opindex0].nextid];
		opnexttype = he[opnextindex0].type;
		opvid0 = he[opnextindex0].vout;
	}
	if (he[opindex0].previd != -1)
	{
		opprevindex0 = heidtoindex[he[opindex0].previd];
		opprevtype = he[opprevindex0].type;
		opvid0 = he[opprevindex0].vin;
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

	else if (((etype == 0 || etype == 3) && (nexttype == 0 || nexttype == 3) && (prevtype == 0 || prevtype == 3)) && ((opetype == 3 || opetype == 0) && (opnexttype == 0 || opnexttype == 3) && (opprevtype == 0 || opprevtype == 3)))
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

			//if (theta<theta0[angle0]) {bendE*=1000; }
			//else {bendE*=100;  }
			bendE *= 1000;
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

double Geometry::dimer_bend_energy(int heindex0)
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
	if (phi < phi0[phitype])
	{
		kPhi /= 5;
	}
	//DbendE = kPhi * (1 - cos(phi - phi0[phitype])); // }
	DbendE = .5 * kPhi * (phi - phi0[phitype]) * (phi - phi0[phitype]);
	//if (phi<phi0[phitype]) { DbendE*=.5; }
	//else {DbendE= kPhi*.1 * (1 - cos(phi-phi0[phitype]));}
	if (DbendE < 0)
	{
		cout << "DBEND E IS " << DbendE << endl;
		exit(-1);
	}
	return DbendE;
}

double Geometry::compute_bind_energy()
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

double Geometry::compute_energy()
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

double Geometry::monomer_energy(int heid0)
{
	//heid0 is on the boundary
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

double Geometry::dimer_energy(int heid0, int heid1)
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

double Geometry::vertex_energy(int vid0)
{
	double v_eng = 0;
	int vindex0 = vidtoindex[vid0];
	//cout << "vertex_energy vindex0" << vindex0 <<endl;
	for (vector<int>::iterator ithe = v[vindex0].hein.begin(); ithe != v[vindex0].hein.end(); ++ithe)
	{
		int heindex = heidtoindex[*ithe];

		v_eng += stretch_energy(heindex);
		//cout << "vertex_energy heindex0 " << heindex << "stretch energy        " << stretch_energy(heindex) << endl;
		v_eng += bend_energy(heindex);
		//cout << "vertex_energy heindex0 " << heindex << "bend energy           " << bend_energy(heindex) << endl;
		v_eng += dimer_bend_energy(heindex);
		//cout << "vertex_energy heindex0 " << heindex << "dimer bend energy     " << dimer_bend_energy(heindex) << endl;

		if (he[heindex].previd != -1)
		{

			int previndex = heidtoindex[he[heindex].previd];
			v_eng += dimer_bend_energy(previndex);
			//cout << "vertex_energy previndex " << previndex << "dimer bend energy     " <<dimer_bend_energy(previndex) <<endl;
			//if (is_boundary(he[heindex].previd)<0) {
			v_eng += bend_energy(previndex);
			//cout << "vertex_energy previndex " << previndex << "bend energy          " <<bend_energy(previndex) <<endl;
			//}
		}
		if (he[heindex].nextid != -1)
		{
			int nextindex = heidtoindex[he[heindex].nextid];
			v_eng += dimer_bend_energy(nextindex);
			//cout << "vertex_energy nextindex " << nextindex << "dimer bend energy      " << dimer_bend_energy(nextindex) <<endl;
		}
		/* in case of existance of bound surface add it */
		if (is_boundary(he[heindex].opid) > 0)
		{
			int opindex = heidtoindex[he[heindex].opid];
			v_eng += dimer_bend_energy(opindex);
			//cout << "vertex_energy opindex " << opindex << "dimer bend energy           " <<dimer_bend_energy(opindex) <<endl;
			if (he[opindex].nextid != -1)
			{
				int nextopindex = heidtoindex[he[opindex].nextid];
				v_eng += bend_energy(nextopindex);
				//cout << "vertex_energy nextopindex " << nextopindex << "bend energy              " <<bend_energy(opindex) <<endl;
			}
		}
	}

	//cout <<"v_eng " << v_eng <<endl;
	return v_eng;
}

/*void Geometry::get_prev_fusion_heid(int heidsurf0){ //not used
	//prev_fusion_heid
	int heindexsurf=heidtoindex[heidsurf0];
	int vid_in = he[heindexsurf].vin; //if this is next
	int vindex_in=vidtoindex[vid_in];
	

	if (v[vindex_in].vneigh.size()>0) {
		for (vector<int>::iterator itv = v[vindex_in].vneigh.begin(); itv != v[vindex_in].vneigh.end(); ++itv){ //loop over neighbors
			int itvindex=vidtoindex[*itv];
			if (veclen(v[vindex_in].co,v[itvindex].co)<xi) { //find close neighbor
				if (v[itvindex].doubleboundary==-1){  // not a doubleboundary vertex
					int heindex_prev=heidtoindex[v[itvindex].hesurfinid]; // this he could be he_prev
					
					int opid_prev=he[heindex_prev].opid;
					//check the angle:
					double cangle = dot(he[heidtoindex[opid_prev]].hevec, he[heindexsurf].hevec) / (he[heidtoindex[opid_prev]].l, he[heindexsurf].l);
					
					if (cangle > .2 ){ //&& not_cross_edge(he[heindexsurf].id, he[heindex_prev].id) > 0){
						if ((connected(he[heindexsurf].vout, he[heindex_prev].vin) >0) && (is_bond_vboundary(he[heindexsurf].vout) || is_bond_vboundary(he[heindex_prev].vin))) {
							he[heindexsurf].prev_wedge_fusion_heid=he[heindex_prev].id;
						}
						else
						{
							he[heindexsurf].prev_fusion_heid=he[heindex_prev].id;
						}
						//he[heindex_prev].next_fusion_heid=he[heindexsurf].id;
						//if (connected(he[heindexsurf].vout, he[heindex_prev].vin) >0) { //now update fusion_vid //maybe not?
						//	v[vidtoindex[vid_in]].fusion_vid=*itv;
						//	v[vidtoindex[*itv]].fusion_vid=vid_in;
					}
					
				}
			}
		}
	}
}*/

/****************get_next_fusion_heid ****************/
/* find close vertex from neighborlist of vout of the edge */
/* if the outgoing edge of that vertx has correct angle that consider half_edge  */
/* if one of the half edges is bond and vertices on the other ends are not connected -> wedge_fusion_heid */
/* otherwise if none of the two are bond on the other end,  fusion_heid*/

void Geometry::get_next_fusion_heid(int heidsurf0)
{

	//next_fusion_heid
	int heindexsurf = heidtoindex[heidsurf0];
	int opindexsurf = heidtoindex[he[heindexsurf].opid];
	//cout << "00 01" <<endl;
	if (heindexsurf == -1)
	{
		cout << "00 erroe in get_next_fusion_heid" << endl;
		exit(-1);
	}

	int vid_out = he[heindexsurf].vout;
	int vindex_out = vidtoindex[vid_out];
	//cout << "00 in get_next_fusion_heid for heid "<< heidsurf0 <<endl;
	//cout << "00 in get_next_fusion_heid vid_out is  "<< vid_out <<endl;
	//if vout has neighbors, check each neighbor vertex
	/* heindexsurf -> vout   */

	if (v[vindex_out].vneigh.size() > 0 && v[vindex_out].doubleboundary == -1)
	{
		//cout << "00 in get_next_fusion_heid for heid "<< heidsurf0 <<endl;
		//cout << "00 in get_next_fusion_heid vid_out is  "<< vid_out <<endl;
		//cout<< "00 in get_next_fusion_heid number of neighbors are "<< v[vindex_out].vneigh.size()<<endl;
		//loop over neighbors
		/* heindexsurf -> vout -> vneigh (close but not connected) */
		for (vector<int>::iterator itv = v[vindex_out].vneigh.begin(); itv != v[vindex_out].vneigh.end(); ++itv)
		{
			//cout <<endl<< "00 03 in get_next_fusion_heid vid_neigh is  "<< *itv <<endl;
			int itvindex = vidtoindex[*itv];

			//find close neighbor , size of hein should be <6
			//next_connected_boundary(vid_out, *itv)<0 &&
			if ((veclen(v[vindex_out].co, v[itvindex].co) < xi) && (v[vindex_out].hein.size() + v[itvindex].hein.size()) < 7)
			{
				//cout << "00 04 veclen(v[vindex_out].co,v[itvindex].co"<< veclen(v[vindex_out].co,v[itvindex].co) <<endl;
				// not a doubleboundary vertex
				if (v[itvindex].doubleboundary == -1)
				{
					// choose the out going halfedge that can be the heid_next for fusion_heid
					/* heindexsurf -> vout -> vneigh (close but not connected) -> heboundaryoutid */

					int heindex_next = heidtoindex[v[itvindex].heboundaryoutid];
					//cout << "00 05in get_next_fusion_heid v[itvindex].heboundaryoutid is "<< v[itvindex].heboundaryoutid <<endl;
					//check the angle between opposit of two half edges:
					double cangle = dot(he[opindexsurf].hevec, he[heindex_next].hevec) / (he[opindexsurf].l, he[heindex_next].l);
					//cout << "00 06 check the angle between opposit of two half edges" << cangle << endl;
					//cout << "00 07 check boundary index of two half edges, he[heindex_next].boundary_index " << he[heindex_next].boundary_index << endl;
					//cout << "00 07 check boundary index of two half edges, he[heindexsurf].boundary_index " << he[heindexsurf].boundary_index << endl;
					if (he[heindex_next].boundary_index == he[heindexsurf].boundary_index && cangle > .2)
					{ // && not_cross_edge(he[heindexsurf].id, he[heindex_next].id) > 0){
						//cout << "00 08 he[heindexsurf].nextid " << he[heindexsurf].nextid << endl;
						//cout << "00 08 he[heindex_next].previd " << he[heindex_next].previd << endl;
						// if one of the half edges is bond and vertices on the other ends are not connected -> wedge_fusion_heid
						if ((he[heindexsurf].previd != -1 || he[heindex_next].nextid != -1) && connected(he[heindexsurf].vin, he[heindex_next].vout) > 0)
						{
							//cout << "00 connected " << he[heindexsurf].vin << " and " <<  he[heindex_next].vout << " is " << connected(he[heindexsurf].vin, he[heindex_next].vout) <<endl;
							he[heindexsurf].next_wedge_fusion_heid = he[heindex_next].id;
							//cout << "00 07" <<" heidsurf0 "<<heidsurf0<< " he[heindexsurf].next_wedge_fusion_heid "<< he[heindexsurf].next_wedge_fusion_heid << endl;
						}
						//if none of the two are bond on the other end,  fusion_heid
						else if (he[heindexsurf].previd == -1 && he[heindex_next].nextid == -1 && connected(he[heindexsurf].vin, he[heindex_next].vout) < 0)
						{
							//cout << "00 08" <<endl;
							he[heindexsurf].next_fusion_heid = he[heindex_next].id;
							//cout << "00 07" <<" heidsurf0 "<<heidsurf0<< " he[heindexsurf].next_fusion_heid "<< he[heindexsurf].next_fusion_heid << endl;
						}
						//he[heindex_next].prev_fusion_heid=he[heindexsurf].id;
						//if (connected(he[heindexsurf].vin, he[heindex_next].vout) >0) { //maybe not!
						//	v[vidtoindex[vid_out]].fusion_vid=*itv;
						//	v[vidtoindex[*itv]].fusion_vid=vid_out;
						//}
					}
				}
			}
		}
	}
}

/*int Geometry::get_fusion_vid(int vidsurf0){

	if (is_bond_vboundary(vidsurf0)>0)  { return -1;} 
	//cout << "in get_fusion vid for vid0 " << vidsurf0 <<endl;   
    int vsurfindex0=vidtoindex[vidsurf0];
	//cout << "vindex is " << vsurfindex0 <<endl;

	//TESTNEW WAY
	if  (v[vsurfindex0].vneigh.size()==0) { return -1;}
	
	int vhecount=0; // prevent 7 fold 
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


	if (is_bond_vboundary(v[vsurfindextemp].vid)>0) { return -1;}
	if (connected(v[vsurfindextemp].vid, vidsurf0) >0) {return -1;}
	double d= veclen(v[vsurfindex0].co, v[vsurfindextemp].co);
    if (vsurfindextemp==vsurfindex0) { cout << " something went wrong " <<endl; exit(-1);}
    //cout << "vertices distance is " << d << endl;
    
    // temp restrict number of bonds in fusion
	if (d < (xi) ){
        return v[vsurfindextemp].vid;
    }

	return -1;
}*/

void Geometry::update_fusion_pairs_he()
{

	//clear fusionhe (vector of all fusion edges)
	//cout <<"in update fusionpair "<<endl;
	if (fusionhe.size() > 0)
	{
		fusionhe.clear();
	}
	if (fusionwedgehe.size() > 0)
	{
		fusionwedgehe.clear();
	}
	// set all fusion_heid_nex_prev=-1

	for (vector<int>::iterator it = boundary.begin(); it != boundary.end(); ++it)
	{
		//cout << "setting pairs for " << *it <<endl;
		he[heidtoindex[*it]].next_fusion_heid = -1;
		he[heidtoindex[*it]].prev_fusion_heid = -1;
		he[heidtoindex[*it]].next_wedge_fusion_heid = -1;
		he[heidtoindex[*it]].prev_wedge_fusion_heid = -1;
	}
	// go over all boundary edges
	//update next_fusion_heid for each vertex

	for (vector<int>::iterator it = boundary.begin(); it != boundary.end(); ++it)
	{
		//cout << "setting next_fusion_heid for " << *it <<endl;
		// only update if fusion_vid==-1
		int heinddex0 = heidtoindex[*it];
		if (he[heinddex0].next_fusion_heid == -1 && he[heinddex0].next_wedge_fusion_heid == -1)
		{
			get_next_fusion_heid(*it);

			// regular fusion
			int temp_next_fusion_heid = he[heinddex0].next_fusion_heid;
			//cout << "temp_next_fusion_heid " << temp_next_fusion_heid <<endl;
			if (temp_next_fusion_heid != -1)
			{
				if (he[heidtoindex[temp_next_fusion_heid]].prev_fusion_heid == -1)
				{
					//cout << "heid is" << *it <<endl;
					//cout << "temp_next_fusion_heid fusion pair!" << temp_next_fusion_heid <<endl;

					he[heidtoindex[temp_next_fusion_heid]].prev_fusion_heid = *it;
					fusionhe.push_back(*it);
					fusionhe.push_back(temp_next_fusion_heid);
				}
				else
				{
					he[heinddex0].next_fusion_heid = -1;
				}
			}
			// wedeg fusion
			int temp_next_wedge_fusion_heid = he[heinddex0].next_wedge_fusion_heid;
			//cout << "temp_next_wedge_fusion_heid" << temp_next_wedge_fusion_heid <<endl;
			if (temp_next_wedge_fusion_heid != -1)
			{
				if (he[heidtoindex[temp_next_wedge_fusion_heid]].prev_wedge_fusion_heid == -1)
				{
					//cout << "heid is" << *it <<endl;
					//cout << "temp_next_fusion_heid fusion pair!" << temp_next_fusion_heid <<endl;

					he[heidtoindex[temp_next_wedge_fusion_heid]].prev_wedge_fusion_heid = *it;

					fusionwedgehe.push_back(*it);
					fusionwedgehe.push_back(temp_next_wedge_fusion_heid);
				}
				else
				{
					he[heinddex0].next_wedge_fusion_heid = -1;
				}
			}
		}
	}

	//temporary double check  this needs an update in function "get_prev_fusion_heid(*it);""
	/*for (vector<int>::iterator it = boundary.begin(); it != boundary.end(); ++it)
	{
		// only update if fusion_vid==-1
		if (he[heidtoindex[*it]].prev_fusion_heid!=-1) {
			int temp_prev_fusion_heid=get_prev_fusion_heid(*it);
			if (he[heidtoindex[temp_prev_fusion_heid]].next_fusion_heid!=*it) {
				cout << "what is wrong with fusion heid" <<endl;
			}
			
			
		}

	}*/
}

/*
void Geometry::update_fusion_pairs(){


	//clear fusionv (vector of all fusion vertices)
	if (fusionv.size()>0) {
		fusionv.clear();
	}

	// set all fusion_vid=-1
	for (vector<int>::iterator it = boundaryv.begin(); it != boundaryv.end(); ++it)
	{
		v[vidtoindex[*it]].fusion_vid=-1;

	}
	// go over all boundary vertices
	//update fusion_vid for each vertex 
	//update fusionv vector of all 
	

	for (vector<int>::iterator it = boundaryv.begin(); it != boundaryv.end(); ++it)
	{
		// only update if fusion_vid==-1
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

}*/

void Geometry::save_vtx(int vid0, VTX *tempvtx)
{

	int vindex0 = vidtoindex[vid0];

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
	//tempvtx->fusion_vid=v[vindex0].fusion_vid;
	tempvtx->hein = v[vindex0].hein;
	tempvtx->doubleboundary = v[vindex0].doubleboundary;
	//tempvtx->hesurfinid=v[vindex0].hesurfinid;
	tempvtx->vid = -1; //decide on this later
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
	vec[0] = (vinit[0] + vfin[0]) / 2;
	vec[1] = (vinit[1] + vfin[1]) / 2;
	vec[2] = (vinit[2] + vfin[2]) / 2;
}

double veclen(double *vinit, double *vfin)
{
	double *vec = new double[3];
	vec[0] = vfin[0] - vinit[0];
	vec[1] = vfin[1] - vinit[1];
	vec[2] = vfin[2] - vinit[2];
	double vlen = norm(vec);
	delete[] vec;
	return (vlen);
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

void dump_lammps_traj(Geometry &g, int time0)
{
	char filename[80];
	float box = 3.0;
	sprintf(filename, "trajlammps_bonds.dat");
	FILE *f;
	f = fopen(filename, "a");
	//fprintf(f,"@<TRIPOS>MOLECULE\n");

	fprintf(f, "LAMMPSDescription-Generated by HEVA at time_step=%d\n", time0);
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
		if (g.is_bond_vboundary(it->vid) > 0)
		{
			fprintf(f, "\n%li 3 %10.6f %10.6f %10.6f", distance(g.v.begin(), it) + 1, it->co[0], it->co[1], it->co[2]);
		}

		else if (it->hein.size() == 5 && g.is_vboundary(it->vid) < 0)
		{
			fprintf(f, "\n%li 1 %10.6f %10.6f %10.6f", distance(g.v.begin(), it) + 1, it->co[0], it->co[1], it->co[2]);
			
		}
		else
		{
			fprintf(f, "\n%li 2 %10.6f %10.6f %10.6f", distance(g.v.begin(), it) + 1, it->co[0], it->co[1], it->co[2]);
		}
		//fprintf(stderr,"\n%li 1 %10.6f %10.6f %10.6f", distance(g.v.begin(),it)+1 ,it->co[0], it->co[1], it->co[2]);
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
			fprintf(f, "\n%d 4 %10.6f %10.6f %10.6f", counter++, x0 + (g.v[vindex]).co[0] + .15 * (it->hevec[0]), x1 + g.v[vindex].co[1] + .15 * (it->hevec[1]), x2 + g.v[vindex].co[2] + .15 * (it->hevec[2]));
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
		//if ( g.is_boundary(it->id)<0 && g.is_boundary(it->opid)<0) {
		//fprintf(f, "\n%li 1 %d %d",distance(g.he.begin(),it)+1 , g.vidtoindex[it->vin]+1, g.vidtoindex[it->vout]+1);
		//fprintf(stderr, "\n%li 1 %d %d",distance(g.he.begin(),it)+1 , g.vidtoindex[it->vin]+1, g.vidtoindex[it->vout]+1);
		//}
		//else if ( g.is_boundary(it->id)<0 && g.is_boundary(it->opid)>0) {
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

void dump_lammps_traj_restart(Geometry &g, int time0)
{ //currently no drug
	char filename[80];
	float box = 3.0;
	sprintf(filename, "trajlammps_restart.dat");
	FILE *f;
	f = fopen(filename, "a");
	//fprintf(f,"@<TRIPOS>MOLECULE\n");

	fprintf(f, "LAMMPSDescription-Generated by HEVA at time_step=%d\n", time0);
	fprintf(f, "\n%d atoms", g.Nv);
	fprintf(f, "\n%d bonds", g.Nhe);
	fprintf(f, "\n%d angles", g.Nhe);				  //next _ prev
	fprintf(f, "\n%li impropers", g.boundary.size()); // prev_boundary this next_boundary
	fprintf(f, "\n");
	fprintf(f, "\n1 atom types"); //vertex
	fprintf(f, "\n4 bond types");
	fprintf(f, "\n1 angle types");
	fprintf(f, "\n%d improper types", g.Nboundary);
	fprintf(f, "\n");
	fprintf(f, "\n%8.3f %8.3f xlo xhi", -box, box);
	fprintf(f, "\n%8.3f %8.3f ylo yhi", -box, box);
	fprintf(f, "\n%8.3f %8.3f zlo zhi", -box, box);
	fprintf(f, "\n");
	fprintf(f, "\nAtoms");
	fprintf(f, "\n");

	for (vector<VTX>::iterator it = g.v.begin(); it != g.v.end(); ++it)
	{
		fprintf(f, "\n%li 1 %10.6f %10.6f %10.6f", distance(g.v.begin(), it) + 1, it->co[0], it->co[1], it->co[2]);
	}

	fprintf(f, "\n");
	fprintf(f, "\nBonds");
	fprintf(f, "\n");

	for (vector<HE>::iterator it = g.he.begin(); it != g.he.end(); ++it)
	{

		if (it->vin == -1 || it->vout == -1 || g.vidtoindex[it->vin] == -1 || g.vidtoindex[it->vout] == -1)
		{
			cout << " dump_data ! error in vin vout of edge " << it->id << endl;
			exit(-1);
		}
		int btype = it->type + 1;
		fprintf(f, "\n%li %d %d %d", distance(g.he.begin(), it) + 1, btype, g.vidtoindex[it->vin] + 1, g.vidtoindex[it->vout] + 1);
	}
	fprintf(f, "\n");
	fprintf(f, "\nAngles"); // this is he - next -prev
	fprintf(f, "\n");
	for (vector<HE>::iterator it = g.he.begin(); it != g.he.end(); ++it)
	{

		int atype = 1;
		int henext = -1;
		int heprev = -1;
		if (it->nextid != -1)
		{
			henext = g.heidtoindex[it->nextid];
		}
		if (it->previd != -1)
		{
			heprev = g.heidtoindex[it->previd];
		}

		fprintf(f, "\n%li %d %d %d %d", distance(g.he.begin(), it) + 1, atype, g.heidtoindex[it->id] + 1, henext + 1, heprev + 1);
	}

	fprintf(f, "\n");
	fprintf(f, "\nImpropers"); // this is he - prev_boundary this next_boundary
	fprintf(f, "\n");
	for (vector<int>::iterator it = g.boundary.begin(); it != g.boundary.end(); ++it)
	{
		int heindex0 = g.heidtoindex[*it];
		int btype = 0; //ToDo should be updated!
		fprintf(f, "\n%li %d %d %d %d", distance(g.boundary.begin(), it) + 1, btype, g.heidtoindex[g.he[heindex0].previd_boundary] + 1, heindex0 + 1, g.he[heindex0].boundary_index);
	}
	fprintf(f, "\n");
	fclose(f);
}

void dump_lammps_data_file(Geometry &g, int time0)
{
	char filename[80];
	float box = 3.0;
	sprintf(filename, "snap_%07d.dat", time0);
	FILE *f;
	f = fopen(filename, "w");
	//fprintf(f,"@<TRIPOS>MOLECULE\n");

	fprintf(f, "LAMMPSDescription-Generated by HEVA at time_step=%d\n", time0);
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
		if (g.is_bond_vboundary(it->vid) > 0)
		{
			fprintf(f, "\n%li 3 %10.6f %10.6f %10.6f", distance(g.v.begin(), it) + 1, it->co[0], it->co[1], it->co[2]);
		}

		else if (it->hein.size() == 5 && g.is_vboundary(it->vid) < 0)
		{
			fprintf(f, "\n%li 1 %10.6f %10.6f %10.6f", distance(g.v.begin(), it) + 1, it->co[0], it->co[1], it->co[2]);
			
		}
		else
		{
			fprintf(f, "\n%li 2 %10.6f %10.6f %10.6f", distance(g.v.begin(), it) + 1, it->co[0], it->co[1], it->co[2]);
		}
		//fprintf(stderr,"\n%li 1 %10.6f %10.6f %10.6f", distance(g.v.begin(),it)+1 ,it->co[0], it->co[1], it->co[2]);
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
			fprintf(f, "\n%d 4 %10.6f %10.6f %10.6f", counter++, x0 + (g.v[vindex]).co[0] + .15 * (it->hevec[0]), x1 + g.v[vindex].co[1] + .15 * (it->hevec[1]), x2 + g.v[vindex].co[2] + .15 * (it->hevec[2]));
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
		//if ( g.is_boundary(it->id)<0 && g.is_boundary(it->opid)<0) {
		//fprintf(f, "\n%li 1 %d %d",distance(g.he.begin(),it)+1 , g.vidtoindex[it->vin]+1, g.vidtoindex[it->vout]+1);
		//fprintf(stderr, "\n%li 1 %d %d",distance(g.he.begin(),it)+1 , g.vidtoindex[it->vin]+1, g.vidtoindex[it->vout]+1);
		//}
		//else if ( g.is_boundary(it->id)<0 && g.is_boundary(it->opid)>0) {
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


void update_geometry_parameters(Geometry &g)
{
	g.NAB=0;
	g.NAB_in = 0;
	g.NCD_Hex = 0;
	g.NCD_other = 0;
	g.NCD_T3 = 0;
	g.NCD_T4 = 0;
	g.NCD_T3_in = 0;
	g.NCD_T4_in = 0;
	g.Nv_in = 0;
	g.Nhe_in = 0;
	g.Nv5 =0;
	g.Nv6=0;
	//cout <<"in dump  initialization"<<endl;
	for (vector<HE>::iterator it = g.he.begin(); it != g.he.end(); ++it)
	{
		//if (it->type==0) g.NCD++;
		int heindex0 = g.heidtoindex[it->id];
		int etype = g.he[heindex0].type;

		int nexttype = -1;
		int prevtype = -1;

		if (it->nextid != -1)
			nexttype = g.he[g.heidtoindex[it->nextid]].type;
		if (it->previd != -1)
			prevtype = g.he[g.heidtoindex[it->previd]].type;

		int opindex0 = g.heidtoindex[it->opid];
		int opetype = g.he[opindex0].type;

		int opnexttype = -1;
		int opprevtype = -1;
		//cout <<"in dump  here 000"<<endl;
		if (g.he[opindex0].nextid != -1)
			opnexttype = g.he[g.heidtoindex[g.he[opindex0].nextid]].type;
		if (g.he[opindex0].previd != -1)
			opprevtype = g.he[g.heidtoindex[g.he[opindex0].previd]].type;
		//cout <<"in dump  here 011"<<endl;
		//consider only the internal structure

		if ((etype == 1) || (etype == 2))
			g.NAB++;


		if (((etype == 0 || etype == 3) && nexttype == 1 && prevtype == 2) && ((opetype == 3 || opetype == 0) && opnexttype == 1 && opprevtype == 2))
				//cout <<"in dump  here 013"<<endl;
				g.NCD_T3 += 1; //CD-BA-AB :: DC-BA-AB in T3
			else if (((etype == 0 || etype == 3) && nexttype == 1 && prevtype == 2) && (opetype == 3 || opetype == 0) && (opnexttype == 0 || opnexttype == 3) && (opprevtype == 0 || opprevtype == 3))
				//cout <<"in dump  here 014"<<endl;
				g.NCD_T4 += 1; //CD_BA_AB :: CD-CD-CD"<<endl;
			else if ((etype == 0 || etype == 3) && (nexttype == 0 || nexttype == 3) && (prevtype == 0 || prevtype == 3) && ((opetype == 3 || opetype == 0) && opnexttype == 1 && opprevtype == 2))
				//cout <<"in dump  here 015"<<endl;
				g.NCD_T4 += 1; // CD-CD-CD :: CD-AB-AB
			else if (((etype == 0 || etype == 3) && (nexttype == 0 || nexttype == 3) && (prevtype == 0 || prevtype == 3)) && ((opetype == 3 || opetype == 0) && (opnexttype == 0 || opnexttype == 3) && (opprevtype == 0 || opprevtype == 3)))
				//cout <<"in dump  here 016"<<endl;
				g.NCD_Hex += 1;
			else if ((etype == 0 || etype == 3))
				//cout <<"in dump  here 017"<<endl;
				g.NCD_other += 1;

		if (g.he[heindex0].nextid != -1 && g.he[opindex0].nextid != -1 && g.v[g.vidtoindex[g.he[g.heidtoindex[g.he[heindex0].nextid]].vout]].hein.size() > 2 && g.v[g.vidtoindex[g.he[g.heidtoindex[g.he[opindex0].nextid]].vout]].hein.size() > 2)
		{
			//cout <<"in dump  here 012"<<endl;
			
			if (((etype == 0 || etype == 3) && nexttype == 1 && prevtype == 2) && ((opetype == 3 || opetype == 0) && opnexttype == 1 && opprevtype == 2))
				//cout <<"in dump  here 013"<<endl;
				g.NCD_T3_in += 1; //CD-BA-AB :: DC-BA-AB in T3
			else if (((etype == 0 || etype == 3) && nexttype == 1 && prevtype == 2) && (opetype == 3 || opetype == 0) && (opnexttype == 0 || opnexttype == 3) && (opprevtype == 0 || opprevtype == 3))
				//cout <<"in dump  here 014"<<endl;
				g.NCD_T4_in += 1; //CD_BA_AB :: CD-CD-CD"<<endl;
			else if ((etype == 0 || etype == 3) && (nexttype == 0 || nexttype == 3) && (prevtype == 0 || prevtype == 3) && ((opetype == 3 || opetype == 0) && opnexttype == 1 && opprevtype == 2))
				//cout <<"in dump  here 015"<<endl;
				g.NCD_T4_in += 1; // CD-CD-CD :: CD-AB-AB
		}

		if (g.v[g.vidtoindex[it->vin]].hein.size() > 2 && g.v[g.vidtoindex[it->vout]].hein.size() > 2)
		{
			g.Nhe_in++;
			if ((etype == 1) || (etype == 2))
				g.NAB_in++;

		}
	}

	g.NAB /= 2;
	g.NCD_T4 /= 2;
	g.NCD_T3 /= 2;
	g.NCD_T4_in /= 2;
	g.NCD_T3_in /= 2;
	g.NCD_Hex /= 2;
	g.NCD_other /= 2;
	g.Nhe_in /= 2;
	g.NAB_in /=2;


	for (vector<VTX>::iterator it = g.v.begin(); it != g.v.end(); ++it)
	{
		if (it->hein.size()>2) g.Nv_in++;
		if (g.is_vboundary(it->vid)<0) {
			if (it->hein.size()==6) g.Nv6++;
			else if (it->hein.size()==5) g.Nv5++;
		}
	}
	
}


void dump_lammps_traj_dimers(Geometry &g, int time0)
{

	char filename[80];
	float box = 3.0;
	
	sprintf(filename, "trajlammps.dat");
	FILE *f;
	f = fopen(filename, "a");
	//fprintf(f,"@<TRIPOS>MOLECULE\n");

	fprintf(f, "LAMMPSDescription-Generated by HEVA at time_step=%d\n", time0);
	fprintf(f, "\n%d atoms", g.Nhe + g.Nhe + g.Nhe + g.Nd + 8);
	fprintf(f, "\n%d bonds", g.Nhe);
	//fprintf(f,"\n%d bonds",g.Nhe/2+g.Nsurf);
	fprintf(f, "\n");
	fprintf(f, "\n10 atom types");
	fprintf(f, "\n4 bond types");
	fprintf(f, "\n");
	fprintf(f, "\n%8.3f %8.3f xlo xhi", -box, box);
	fprintf(f, "\n%8.3f %8.3f ylo yhi", -box, box);
	fprintf(f, "\n%8.3f %8.3f zlo zhi", -box, box);
	fprintf(f, "\n");
	fprintf(f, "\nAtoms");
	fprintf(f, "\n");
	//vin of each half edge
	int counter = 1;

	//cout <<"in dump  here 222"<<endl;
	for (vector<HE>::iterator it = g.he.begin(); it != g.he.end(); ++it)
	{

		//double x0 = g.v[g.vidtoindex[it->vin]].co[0];
		//double x1 = g.v[g.vidtoindex[it->vin]].co[1];
		//double x2 = g.v[g.vidtoindex[it->vin]].co[2];

		double x0 = .9 * (g.v[g.vidtoindex[it->vin]].co[0]) + .1 * (g.v[g.vidtoindex[it->vout]].co[0]);
		double x1 = .9 * (g.v[g.vidtoindex[it->vin]].co[1]) + .1 * (g.v[g.vidtoindex[it->vout]].co[1]);
		double x2 = .9 * (g.v[g.vidtoindex[it->vin]].co[2]) + .1 * (g.v[g.vidtoindex[it->vout]].co[2]);

		int atype = it->type + 1;

		if (it->type == 1)
		{ // AB A
			fprintf(f, "\n%d %d %10.6f %10.6f %10.6f", counter++, atype, x0, x1, x2);
		}
		else if (it->type == 2)
		{ //AB B
			fprintf(f, "\n%d %d %10.6f %10.6f %10.6f", counter++, atype, x0, x1, x2);
		}
		else if (it->type == 0)
		{ //CD C
			fprintf(f, "\n%d %d %10.6f %10.6f %10.6f", counter++, atype, x0, x1, x2);
		}
		else if (it->type == 3)
		{ //CD D
			if (g.is_boundary(it->id) < 0)
			{
				int nexttype = g.he[g.heidtoindex[it->nextid]].type;
				int prevtype = g.he[g.heidtoindex[it->previd]].type;

				if ((nexttype == 1) && (prevtype == 2))
					atype = 1;
			}
			fprintf(f, "\n%d %d %10.6f %10.6f %10.6f", counter++, atype, x0, x1, x2);
		}

		if (g.v[g.vidtoindex[it->vin]].hein.size() > 2 && g.v[g.vidtoindex[it->vout]].hein.size() > 2)
		{
			g.Nhe_in++;
		}
	}
	//cout <<"in dump   here 333"<<endl;
	g.Nhe_in /= 2;
	//int counter = g.Nhe + 1;

	//he center beads
	for (vector<HE>::iterator it = g.he.begin(); it != g.he.end(); ++it)
	{

		double x0 = it->hecent[0];
		double x1 = it->hecent[1];
		double x2 = it->hecent[2];
		int atype = it->type + 1;
		if (it->type == 1)
		{ // AB A
			fprintf(f, "\n%d %d %10.6f %10.6f %10.6f", counter++, atype, x0, x1, x2);
		}
		else if (it->type == 2)
		{ //AB B
			fprintf(f, "\n%d %d %10.6f %10.6f %10.6f", counter++, atype, x0, x1, x2);
		}
		else if (it->type == 0)
		{ //CD C
			fprintf(f, "\n%d %d %10.6f %10.6f %10.6f", counter++, atype, x0, x1, x2);
		}
		else if (it->type == 3)
		{ //CD D
			if (g.is_boundary(it->id) < 0)
			{
				int nexttype = g.he[g.heidtoindex[it->nextid]].type;
				int prevtype = g.he[g.heidtoindex[it->previd]].type;

				if ((nexttype == 1) && (prevtype == 2))
					atype = 1;
			}
			fprintf(f, "\n%d %d %10.6f %10.6f %10.6f", counter++, atype, x0, x1, x2);
		}
	}
	//cout <<"in dump   here 444"<<endl;
	//drug beads
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
			fprintf(f, "\n%d 5 %10.6f %10.6f %10.6f", counter++, x0 + (g.v[vindex]).co[0] + .15 * (it->hevec[0]), x1 + g.v[vindex].co[1] + .15 * (it->hevec[1]), x2 + g.v[vindex].co[2] + .15 * (it->hevec[2]));
		}
	}
	//cout <<"in dump   here 555"<<endl;

	fprintf(f, "\n%d 6 %10.6f %10.6f %10.6f", 2 * g.Nhe + g.Nd + 1, box, box, box);
	fprintf(f, "\n%d 6 %10.6f %10.6f %10.6f", 2 * g.Nhe + g.Nd + 2, -box, box, box);
	fprintf(f, "\n%d 6 %10.6f %10.6f %10.6f", 2 * g.Nhe + g.Nd + 3, box, -box, box);
	fprintf(f, "\n%d 6 %10.6f %10.6f %10.6f", 2 * g.Nhe + g.Nd + 4, box, box, -box);
	fprintf(f, "\n%d 6 %10.6f %10.6f %10.6f", 2 * g.Nhe + g.Nd + 5, -box, -box, box);
	fprintf(f, "\n%d 6 %10.6f %10.6f %10.6f", 2 * g.Nhe + g.Nd + 6, -box, box, -box);
	fprintf(f, "\n%d 6 %10.6f %10.6f %10.6f", 2 * g.Nhe + g.Nd + 7, box, -box, -box);
	fprintf(f, "\n%d 6 %10.6f %10.6f %10.6f", 2 * g.Nhe + g.Nd + 8, -box, -box, -box);

	counter += 8;
	//temp
	for (vector<HE>::iterator it = g.he.begin(); it != g.he.end(); ++it)
	{
		double x0 = it->hetop[0];
		double x1 = it->hetop[1];
		double x2 = it->hetop[2];
		//int atype=-1;
		//if ((it->type == 1) | (it->type == 2)) atype=2;

		int atype = 7 + it->type;

		fprintf(f, "\n%d %d %10.6f %10.6f %10.6f", counter++, atype, x0, x1, x2);
	}
	//cout <<"in dump   here 666"<<endl;
	// now bonds
	fprintf(f, "\n");
	fprintf(f, "\nBonds");
	fprintf(f, "\n");

	//edge_counter=1;
	for (vector<HE>::iterator it = g.he.begin(); it != g.he.end(); ++it)
	{
		// 1 A 2 B 3 C 4 D
		int btype = it->type + 1;
		//vin and the middle
		int index = distance(g.he.begin(), it);
		fprintf(f, "\n%d %d %d %d", index + 1, btype, index + 1, g.Nhe + index + 1);
	}
	//cout <<"in dump   here 777"<<endl;
	fprintf(f, "\n");
	fclose(f);
}

void dump_lammps_data_dimers(Geometry &g, int time0)
{

	char filename[80];
	float box = 3.0;
	sprintf(filename, "snap_%07d.dat", time0);
	FILE *f;
	f = fopen(filename, "w");
	//fprintf(f,"@<TRIPOS>MOLECULE\n");

	fprintf(f, "LAMMPSDescription-Generated by HEVA at time_step=%d\n", time0);
	fprintf(f, "\n%d atoms", g.Nhe + g.Nhe + g.Nd + 8);
	fprintf(f, "\n%d bonds", g.Nhe);
	//fprintf(f,"\n%d bonds",g.Nhe/2+g.Nsurf);
	fprintf(f, "\n");
	fprintf(f, "\n6 atom types");
	fprintf(f, "\n4 bond types");
	fprintf(f, "\n");
	fprintf(f, "\n%8.3f %8.3f xlo xhi", -box, box);
	fprintf(f, "\n%8.3f %8.3f ylo yhi", -box, box);
	fprintf(f, "\n%8.3f %8.3f zlo zhi", -box, box);
	fprintf(f, "\n");
	fprintf(f, "\nAtoms");
	fprintf(f, "\n");
	//cout << "here in dump 000"<<endl;
	int counter = 1;
	for (vector<HE>::iterator it = g.he.begin(); it != g.he.end(); ++it)
	{

		//double x0 = g.v[g.vidtoindex[it->vin]].co[0];
		//double x1 = g.v[g.vidtoindex[it->vin]].co[1];
		//double x2 = g.v[g.vidtoindex[it->vin]].co[2];

		double x0 = .9 * (g.v[g.vidtoindex[it->vin]].co[0]) + .1 * (g.v[g.vidtoindex[it->vout]].co[0]);
		double x1 = .9 * (g.v[g.vidtoindex[it->vin]].co[1]) + .1 * (g.v[g.vidtoindex[it->vout]].co[1]);
		double x2 = .9 * (g.v[g.vidtoindex[it->vin]].co[2]) + .1 * (g.v[g.vidtoindex[it->vout]].co[2]);

		int atype = it->type + 1;

		if (it->type == 1)
		{ // AB A
			fprintf(f, "\n%d %d %10.6f %10.6f %10.6f", counter++, atype, x0, x1, x2);
		}
		else if (it->type == 2)
		{ //AB B
			fprintf(f, "\n%d %d %10.6f %10.6f %10.6f", counter++, atype, x0, x1, x2);
		}
		else if (it->type == 0)
		{ //CD C
			fprintf(f, "\n%d %d %10.6f %10.6f %10.6f", counter++, atype, x0, x1, x2);
		}
		else if (it->type == 3)
		{ //CD D
			if (g.is_boundary(it->id) < 0)
			{
				int nexttype = g.he[g.heidtoindex[it->nextid]].type;
				int prevtype = g.he[g.heidtoindex[it->previd]].type;

				if ((nexttype == 1) && (prevtype == 2))
					atype = 1;
			}
			fprintf(f, "\n%d %d %10.6f %10.6f %10.6f", counter++, atype, x0, x1, x2);
		}
	}
	//int counter = g.Nhe + 1;

	//he center beads
	for (vector<HE>::iterator it = g.he.begin(); it != g.he.end(); ++it)
	{

		double x0 = it->hecent[0];
		double x1 = it->hecent[1];
		double x2 = it->hecent[2];
		int atype = it->type + 1;
		if (it->type == 1)
		{ // AB A
			fprintf(f, "\n%d %d %10.6f %10.6f %10.6f", counter++, atype, x0, x1, x2);
		}
		else if (it->type == 2)
		{ //AB B
			fprintf(f, "\n%d %d %10.6f %10.6f %10.6f", counter++, atype, x0, x1, x2);
		}
		else if (it->type == 0)
		{ //CD C
			fprintf(f, "\n%d %d %10.6f %10.6f %10.6f", counter++, atype, x0, x1, x2);
		}
		else if (it->type == 3)
		{ //CD D
			if (g.is_boundary(it->id) < 0)
			{
				int nexttype = g.he[g.heidtoindex[it->nextid]].type;
				int prevtype = g.he[g.heidtoindex[it->previd]].type;

				if ((nexttype == 1) && (prevtype == 2))
					atype = 1;
			}
			fprintf(f, "\n%d %d %10.6f %10.6f %10.6f", counter++, atype, x0, x1, x2);
		}
	}

	//drug beads
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
			fprintf(f, "\n%d 5 %10.6f %10.6f %10.6f", counter++, x0 + (g.v[vindex]).co[0] + .15 * (it->hevec[0]), x1 + g.v[vindex].co[1] + .15 * (it->hevec[1]), x2 + g.v[vindex].co[2] + .15 * (it->hevec[2]));
		}
	}

	fprintf(f, "\n%d 6 %10.6f %10.6f %10.6f", 2 * g.Nhe + g.Nd + 1, box, box, box);
	fprintf(f, "\n%d 6 %10.6f %10.6f %10.6f", 2 * g.Nhe + g.Nd + 2, -box, box, box);
	fprintf(f, "\n%d 6 %10.6f %10.6f %10.6f", 2 * g.Nhe + g.Nd + 3, box, -box, box);
	fprintf(f, "\n%d 6 %10.6f %10.6f %10.6f", 2 * g.Nhe + g.Nd + 4, box, box, -box);
	fprintf(f, "\n%d 6 %10.6f %10.6f %10.6f", 2 * g.Nhe + g.Nd + 5, -box, -box, box);
	fprintf(f, "\n%d 6 %10.6f %10.6f %10.6f", 2 * g.Nhe + g.Nd + 6, -box, box, -box);
	fprintf(f, "\n%d 6 %10.6f %10.6f %10.6f", 2 * g.Nhe + g.Nd + 7, box, -box, -box);
	fprintf(f, "\n%d 6 %10.6f %10.6f %10.6f", 2 * g.Nhe + g.Nd + 8, -box, -box, -box);
	// now bonds
	fprintf(f, "\n");
	fprintf(f, "\nBonds");
	fprintf(f, "\n");

	//edge_counter=1;
	for (vector<HE>::iterator it = g.he.begin(); it != g.he.end(); ++it)
	{
		// 1 A 2 B 3 C 4 D
		int btype = it->type + 1;
		//vin and the middle
		int index = distance(g.he.begin(), it);
		fprintf(f, "\n%d %d %d %d", index + 1, btype, index + 1, g.Nhe + index + 1);
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


void read_points(Geometry &g)
{

	char filename[80];
	sprintf(filename, "points.txt");
	char temp1[20], temp2[20], temp0[20];
	//g.lenpoints=200000;
	FILE *file;
	file = fopen(filename, "r");
	
	int x=-1;
	for (int i = 0; i < g.lenpoints; i++)
	{

		x = fscanf(file, "%s %s %s\n", temp0, temp1, temp2);
		if (x==-1){
			cout <<"error in read points" <<endl;
			exit(-1);
		}
		g.dist_points[i][0] = atof(temp0);
		g.dist_points[i][1] = atof(temp1);
		g.dist_points[i][2] = atof(temp2);
		
	}
	
	fclose(file);
}



void read_lammps_data(Geometry &g, char filename[])
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
		g.add_half_edge_type(atoi(temp1) - 1, atoi(temp2) - 1, atoi(temp0) - 1, -1);
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
	fclose(file);
	//exit(-1);
}

/* this fuction should be updated with boundary index */
int read_restart_lammps_data_file(Geometry &g, char filename[])
{

	int fNv = 0;
	int fNhe = 0;

	char s[100];
	char temp1[20], temp2[20], temp0[20];
	char TT[] = "Bonds";
	//char AA[] = "Atoms";
	char GG[] = "Angles";

	char II[] = "Impropers";

	FILE *file;
	if ((file = fopen(filename, "r")))
	{
		cout << "reading restart file" << endl;
		//int x = 0;
		int x = fscanf(file, "%*s %*s %s\n", temp0);
		int boundarysize = -1;
		if (x == 0)
		{
			exit(-1);
		}
		int step = atoi(temp0);
		x = fscanf(file, "\n%d %*s", &fNv);			 //atoms
		x = fscanf(file, "\n%d %*s", &fNhe);		 //bonds
		x = fscanf(file, "\n%*s %*s");				 //angles next_ prev
		x = fscanf(file, "\n%d %*s", &boundarysize); //impropers prev_boundary this next_boundary
		x = fscanf(file, "\n");
		x = fscanf(file, "\n%*s %*s %*s"); //vertex
		x = fscanf(file, "\n%*s %*s %*s");
		x = fscanf(file, "\n%*s %*s %*s");
		x = fscanf(file, "\n%*s %*s %*s"); //impropers
		x = fscanf(file, "\n");
		x = fscanf(file, "\n%*s %*s %*s %*s");
		x = fscanf(file, "\n%*s %*s %*s %*s");
		x = fscanf(file, "\n%*s %*s %*s %*s");
		x = fscanf(file, "\n");
		x = fscanf(file, "\n%*s");
		x = fscanf(file, "\n");

		double *vec = new double[3];
		for (int i = 0; i < fNv; i++)
		{

			x = fscanf(file, "%*s %*s %s %s %s\n", temp0, temp1, temp2);

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
		/* this part needs update for boundary edges with boundary index*/
		for (int i = 0; i < fNhe; i++)
		{
			x = fscanf(file, "%*s %s %s %s\n", temp0, temp1, temp2);
			g.add_half_edge_type(atoi(temp1) - 1, atoi(temp2) - 1, atoi(temp0) - 1, -1);
			//fprintf(stderr, "add edge %d  %d %d %d\n",i,atoi(temp1)-1, atoi(temp2)-1,atoi(temp0)-1 );
			//cout << i <<endl;
		}

		cout << "now next previus from angles" << endl;
		while (fscanf(file, "%s", s) == 1)
		{ //s[0]!='1') {
			if (strcmp(s, GG) == 0)
			{
				break;
			}
		}
		for (int i = 0; i < fNhe; i++)
		{
			x = fscanf(file, "%*s %*s %s %s %s\n", temp0, temp1, temp2);
			g.set_prev_next(atoi(temp0) - 1, atoi(temp2) - 1, atoi(temp1) - 1);
		}

		cout << "now next previus boundary from impropers" << endl;
		while (fscanf(file, "%s", s) == 1)
		{ //s[0]!='1') {
			if (strcmp(s, II) == 0)
			{
				break;
			}
		}
		for (int i = 0; i < boundarysize; i++)
		{
			x = fscanf(file, "%*s %*s %s %s %s \n", temp0, temp1, temp2);
			g.he[atoi(temp0) - 1].boundary_index = atoi(temp2);
			g.he[atoi(temp1) - 1].boundary_index = atoi(temp2);
			g.set_prev_next_boundary(atoi(temp0) - 1, atoi(temp1) - 1);
			// this should get updated for double boundary
		}

		g.update_boundary();

		/*for (vector<HE>::iterator it = g.he.begin() ; it != g.he.end(); it++) {
			cout << "EDGE  " <<distance(g.he.begin(),it) <<" id " << it->id << " vin " << g.vidtoindex[it->vin] << "vout " << g.vidtoindex[it->vout] <<endl;
			cout << "      opid " << it->opid << " nextid " << it->nextid << " previd " << it->previd <<endl;
			cout << "      hevec " << it->hevec[0] << " " << it->hevec[1]  << " "  << it->hevec[2] <<" l is " << it->l << endl;
			//cout << "      normal " << it->n[0] << " " << it->n[1] << " " << it->n[2] << " " <<endl <<endl;
		}
		cout << "Nt " << Nt <<endl;;*/
		//fprintf(stderr," edges : %d\n",g.Nhe);
		fclose(file);
		return (step);
	}
	else
	{
		make_initial_triangle(g);

	}
	return (0);
}

/* this fuction should be updated with boundary index */
int read_restart_lammps_data_traj(Geometry &g, FILE *trajfile, int step = -1)
{

	char s[100];
	char temp1[20], temp2[20], temp0[20];
	char TT[] = "Bonds";
	char GG[] = "Angles";
	char II[] = "Impropers";
	int initialstep = step;

	//initialstep=0 read trajectory and do analysis
	//initialstep=-1 only read trajectory
	//initailstep>0 read one fram and do analysis

	FILE *analysis_file;
	analysis_file = fopen("analysis.dat", "w");
	cout << "analysis file openned" << endl;
	
	if (trajfile)
	{
		cout << "reading restart file" << endl;
		;
		while (fscanf(trajfile, "%*s %*s %s\n", temp0) == 1)
		{
			int x = 0;
			int fNv = -1;
			int fNhe = -1;
			int boundarysize = -1;
			step = atoi(temp0); //timestep -- sweep
			cout << "step " << step << endl;
			x = fscanf(trajfile, "\n%*s");					 //atoms
			if (x==0) break;

			x = fscanf(trajfile, "\n%*s");					 //bonds
			x = fscanf(trajfile, "\n%d %*s", &fNv);			 //atoms
			x = fscanf(trajfile, "\n%d %*s", &fNhe);		 //bonds
			x = fscanf(trajfile, "\n%*s %*s");				 //angles next_ prev
			x = fscanf(trajfile, "\n%d %*s", &boundarysize); //impropers prev_boundary this next_boundary
			x = fscanf(trajfile, "\n");
			x = fscanf(trajfile, "\n%*s %*s %*s"); //vertex
			x = fscanf(trajfile, "\n%*s %*s %*s");
			x = fscanf(trajfile, "\n%*s %*s %*s");
			x = fscanf(trajfile, "\n%*s %*s %*s"); //impropers
			x = fscanf(trajfile, "\n");
			x = fscanf(trajfile, "\n%*s %*s %*s %*s");
			x = fscanf(trajfile, "\n%*s %*s %*s %*s");
			x = fscanf(trajfile, "\n%*s %*s %*s %*s");
			x = fscanf(trajfile, "\n");
			x = fscanf(trajfile, "\n%*s");
			x = fscanf(trajfile, "\n");

			cout << " fNv is " << fNv << endl;
			cout << " fNhe is " << fNhe << endl;
			cout << "boundaysize is " << boundarysize << endl;
			double *vec = new double[3];
			for (int i = 0; i < fNv; i++)
			{

				x = fscanf(trajfile, "%*s %*s %s %s %s\n", temp0, temp1, temp2);

				vec[0] = atof(temp0);
				vec[1] = atof(temp1);
				vec[2] = atof(temp2);
				//fprintf(stderr,"%s %s %s %s %s\n" ,index,vtype ,temp0,temp1,temp2);
				g.add_vertex(vec);
				//fprintf(stderr, "%d  %f %f %f \n",i, g.v[i][0],g.v[i][1],g.v[i][2]);
			}
			delete[] vec;
			fprintf(stderr, "read all Vertices\n graph has %d vertices", g.Nv);

			//cout << "now edges" << endl;
			while (fscanf(trajfile, "%s", s) == 1)
			{ //s[0]!='1') {
				if (strcmp(s, TT) == 0)
				{
					break;
				}
			}
			/* this part needs update for boundary edges with boundary index*/
			cout << "fNhe is " << fNhe << endl;
			for (int i = 0; i < fNhe; i++)
			{
				x = fscanf(trajfile, "%*s %s %s %s\n", temp0, temp1, temp2);
				g.add_half_edge_type(atoi(temp1) - 1, atoi(temp2) - 1, atoi(temp0) - 1, -1);
				fprintf(stderr, "add edge %d  %d %d %d\n", i, atoi(temp1) - 1, atoi(temp2) - 1, atoi(temp0) - 1);
				//cout << i <<endl;
			}

			cout << "now next previus from angles" << endl;
			while (fscanf(trajfile, "%s", s) == 1)
			{ //s[0]!='1') {
				if (strcmp(s, GG) == 0)
				{
					break;
				}
			}
			for (int i = 0; i < fNhe; i++)
			{
				x = fscanf(trajfile, "%*s %*s %s %s %s\n", temp0, temp1, temp2);
				g.set_prev_next(atoi(temp0) - 1, atoi(temp2) - 1, atoi(temp1) - 1);
			}

			cout << "now next previus boundary from impropers" << endl;
			while (fscanf(trajfile, "%s", s) == 1)
			{ //s[0]!='1') {
				if (strcmp(s, II) == 0)
				{
					break;
				}
			}
			cout << " here boundary size is" << boundarysize << endl;
			cout << g.he.size() << endl;
			for (int i = 0; i < boundarysize; i++)
			{
				x = fscanf(trajfile, "%*s %*s %s %s %s \n", temp0, temp1, temp2);
				cout << "read trajfile boundary index nest prev boundary" << endl;
				cout << "g.he[atoi(temp0) - 1].boundary_index=atoi(temp2);" << g.he[atoi(temp0) - 1].boundary_index << " = " << atoi(temp2) << endl;
				g.he[atoi(temp0) - 1].boundary_index = atoi(temp2);
				g.he[atoi(temp1) - 1].boundary_index = atoi(temp2);
				cout << "set next prev boundary" << endl;
				g.set_prev_next_boundary(atoi(temp0) - 1, atoi(temp1) - 1);
				// this should get updated for double boundary
			}

			cout << "update neigh" << endl;
			g.update_neigh();
			cout << "update normals" << endl;
			g.update_normals();

			if (initialstep == -1)
				break;

			cout << "now do analysis" << endl;

			dump_analysis(g, analysis_file, step, -1, -1);

			if (initialstep > 0)
				break;

			cout << "step" << step << endl;
			step++;
		}
		return (step);
	}
	
	fclose(analysis_file);
	return (0);
}

void dump_restart_lammps_data_file(Geometry &g, int time0)
{ //currently no drug
	char filename[80];
	float box = 3.0;
	sprintf(filename, "restart_lammps.dat");
	FILE *f;
	f = fopen(filename, "w");

	fprintf(f, "HEVA-LAMMPSDescription-Generated  time_step= %d\n", time0);
	fprintf(f, "\n%d atoms", g.Nv);
	fprintf(f, "\n%d bonds", g.Nhe);
	fprintf(f, "\n%d angles", g.Nhe);				  //next _ prev
	fprintf(f, "\n%li impropers", g.boundary.size()); // prev_boundary this next_boundary
	fprintf(f, "\n");
	fprintf(f, "\n1 atom types"); //vertex
	fprintf(f, "\n4 bond types");
	fprintf(f, "\n1 angle types");
	fprintf(f, "\n%d improper types", g.Nboundary);
	fprintf(f, "\n");
	fprintf(f, "\n%8.3f %8.3f xlo xhi", -box, box);
	fprintf(f, "\n%8.3f %8.3f ylo yhi", -box, box);
	fprintf(f, "\n%8.3f %8.3f zlo zhi", -box, box);
	fprintf(f, "\n");
	fprintf(f, "\nAtoms");
	fprintf(f, "\n");

	for (vector<VTX>::iterator it = g.v.begin(); it != g.v.end(); ++it)
	{
		fprintf(f, "\n%li 1 %10.6f %10.6f %10.6f", distance(g.v.begin(), it) + 1, it->co[0], it->co[1], it->co[2]);
	}

	fprintf(f, "\n");
	fprintf(f, "\nBonds");
	fprintf(f, "\n");

	for (vector<HE>::iterator it = g.he.begin(); it != g.he.end(); ++it)
	{

		if (it->vin == -1 || it->vout == -1 || g.vidtoindex[it->vin] == -1 || g.vidtoindex[it->vout] == -1)
		{
			cout << " dump_data ! error in vin vout of edge " << it->id << endl;
			exit(-1);
		}
		int btype = it->type + 1;
		fprintf(f, "\n%li %d %d %d", distance(g.he.begin(), it) + 1, btype, g.vidtoindex[it->vin] + 1, g.vidtoindex[it->vout] + 1);
	}
	fprintf(f, "\n");
	fprintf(f, "\nAngles"); // this is he - next -prev
	fprintf(f, "\n");
	for (vector<HE>::iterator it = g.he.begin(); it != g.he.end(); ++it)
	{

		int atype = 1;
		int henext = -1;
		int heprev = -1;
		if (it->nextid != -1)
		{
			henext = g.heidtoindex[it->nextid];
		}
		if (it->previd != -1)
		{
			heprev = g.heidtoindex[it->previd];
		}

		fprintf(f, "\n%li %d %d %d %d", distance(g.he.begin(), it) + 1, atype, g.heidtoindex[it->id] + 1, henext + 1, heprev + 1);
	}

	fprintf(f, "\n");
	fprintf(f, "\nImpropers"); // this is he - prev_boundary this next_boundary
	fprintf(f, "\n");
	for (vector<int>::iterator it = g.boundary.begin(); it != g.boundary.end(); ++it)
	{
		int heindex0 = g.heidtoindex[*it];
		int btype = 0; //ToDo should be updated!
		fprintf(f, "\n%li %d %d %d %d", distance(g.boundary.begin(), it) + 1, btype, g.heidtoindex[g.he[heindex0].previd_boundary] + 1, heindex0 + 1, g.he[heindex0].boundary_index);
	}
	fprintf(f, "\n");
	fclose(f);
}

void dump_data_frame(Geometry &g, FILE *f, int time)
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

void dump_analysis(Geometry &g, FILE *ofile, int sweep = -1, int seed = -1, int seconds = -1)
{

	if (sweep == 0)
		fprintf(ofile, "sweep,seed,seconds,epsilon,kappa,kappaPhi,theta0,theta1,gb0,mu,dmu,dg,theta2,energy,binding_energy,Nv5,Nv6,NAB,NAB_in,NCD_Hex,NCD_other,NVin,Nhein,NCD_T4_in, NCD_T3_in,NCD_T4,NCD_T3,Nv,NE,Nsurf,Nboundary\n");
																			//Nv5,	Nv6,	NAB,	NAB_in,	NCD_Hex, 	NCD_other, 	NVin,	Nhein,NCD_T4,	NCD_T3,	Nv,	NE,	Nsurf, Nboundary\n");
	update_geometry_parameters(g);

	fprintf(ofile, "%d,%d,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.5f,%.5f,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",
			sweep, seed, seconds, g.epsilon[0], g.kappa[0], g.kappaPhi[0], g.theta0[0], g.theta0[1], g.gb0, g.mu[0], g.mu[1] - g.mu[0], g.dg, g.theta0[2],
			g.compute_energy(),g.compute_bind_energy(), g.Nv5, g.Nv6, g.NAB, g.NAB_in,  g.NCD_Hex, g.NCD_other, g.Nv_in,g.Nhe_in,g.NCD_T4_in, g.NCD_T3_in,g.NCD_T4, g.NCD_T3,  g.Nv, g.Nhe / 2,g.Nsurf, g.Nboundary);
			                    //Nv5,	Nv6,	NAB,	NAB_in,	NCD_Hex, 	NCD_other, 		NVin,	Nhein,NCD_T4,	NCD_T3,		Nv,		NE,		Nsurf, 	Nboundary\n");
	fflush(ofile);
}

void recenter(Geometry &g)
{

	double XCM = 0;
	double YCM = 0;
	double ZCM = 0;
	for (vector<VTX>::iterator it = g.v.begin(); it != g.v.end(); ++it)
	{
		XCM += it->co[0];
		YCM += it->co[1];
		ZCM += it->co[2];
	}
	XCM /= g.Nv;
	YCM /= g.Nv;
	ZCM /= g.Nv;

	//cout << "HERE2"<<endl;
	double vx = XCM; //g.v[0].co[0];
	double vy = YCM; //g.v[0].co[1];
	double vz = ZCM; //g.v[0].co[2];
	for (vector<VTX>::iterator it = g.v.begin(); it != g.v.end(); ++it)
	{
		it->co[0] -= vx;
		it->co[1] -= vy;
		it->co[2] -= vz;
	}
	g.update_boundary();
}

int surfclosev(Geometry &g)
{
	int alln = 0;
	for (vector<VTX>::iterator it = g.v.begin(); it != g.v.end(); ++it)
	{
		alln += it->vneigh.size();
		if (it->vneigh.size() > 0)
		{
			//cout<< "vindex is " << distance(g.v.begin(),it) <<  " vid is" << it->vid ;
			for (vector<int>::iterator itv = it->vneigh.begin(); itv != it->vneigh.end(); itv++)
			{

				//cout <<"     neighbors are " << *itv << " vindex is " << g.vidtoindex[*itv] ;
				if (g.vidtoindex[*itv] == -1)
				{
					cout << "wrong neighbor" << endl;
					exit(-1);
				}
			}
			//cout<<endl;
		}
	}

	if (alln % 2 != 0)
	{
		cout << " odd neighbors" << endl;
		exit(-1);
	}

	int surfclosevCount = 0;
	for (vector<int>::iterator it = g.boundaryv.begin(); it != g.boundaryv.end(); ++it)
	{
		if (g.v[g.vidtoindex[*it]].vneigh.size() > 0)
		{
			for (vector<int>::iterator itv = g.v[g.vidtoindex[*it]].vneigh.begin(); itv != g.v[g.vidtoindex[*it]].vneigh.end(); itv++)
			{
				if (veclen(g.v[g.vidtoindex[*itv]].co, g.v[g.vidtoindex[*it]].co) < .5 * g.l0[0])
				{
					surfclosevCount += 1;
					break;
				}
			}
		}
	}
	return (surfclosevCount);
}
//void show_status(Geometry &g, int frame, int sweep, int seconds ){

void make_initial_triangle(Geometry &g)
{
	double xyz0[3];
	xyz0[0] = 0;
	xyz0[1] = 0;
	xyz0[2] = 0;

	for (int i = 0; i < 3; i++)
	{
		g.add_vertex(xyz0);
		xyz0[0] = cos(i * PI / 3);
		xyz0[1] = sin(i * PI / 3);
		xyz0[2] = 0;
	}
	//for (int i=0; i<3; i++) {

	//	int j=i+1;
	//	if (i==2) { j=0;}
	//	if (i>=Nvlast || j>=Nvlast) { cout << "ERROR in make triangle" << endl; exit(-1); }
	g.add_half_edge_type(g.v[0].vid, g.v[1].vid, 2, -1);
	g.add_half_edge_type(g.v[1].vid, g.v[0].vid, 1, 0); //boundary
	g.add_half_edge_type(g.v[1].vid, g.v[2].vid, 0, -1);
	g.add_half_edge_type(g.v[2].vid, g.v[1].vid, 3, 0); //boundary
	g.add_half_edge_type(g.v[2].vid, g.v[0].vid, 1, -1);
	g.add_half_edge_type(g.v[0].vid, g.v[2].vid, 2, 0); //boundary

	g.set_prev_next(g.he[0].id, g.he[4].id, g.he[2].id);
	g.set_prev_next(g.he[2].id, g.he[0].id, g.he[4].id);
	g.set_prev_next(g.he[4].id, g.he[2].id, g.he[0].id);

	for (vector<HE>::iterator it = g.he.begin(); it != g.he.end(); it++)
	{
		cout << "in make triangle updating edge" << it->id << endl;

		g.update_half_edge(it->id);
		cout << it->id << "in make triangle  TYPE " << g.he[it->id].type << " opid " << it->opid << " OP TYPE " << g.he[g.heidtoindex[it->opid]].type << endl;
		cout << it->id << "in make triangle  ID " << it->id << " nextid " << it->nextid << " previd " << it->previd << endl;
		cout << it->id << "in make triangle  boundary Index " << it->boundary_index << " next_boundary " << it->nextid_boundary << " previd boundary " << it->previd_boundary << endl;
	}

	g.update_boundary();
}

void make_initial_pentamer(Geometry &g)
{
	double xyz0[3];
	xyz0[0] = 0;
	xyz0[1] = 0;
	xyz0[2] = 0;

	for (int i = 0; i < 3; i++)
	{
		g.add_vertex(xyz0);
		xyz0[0] = cos(i * PI / 3);
		xyz0[1] = sin(i * PI / 3);
		xyz0[2] = 0;
	}
	//for (int i=0; i<3; i++) {

	//	int j=i+1;
	//	if (i==2) { j=0;}
	//	if (i>=Nvlast || j>=Nvlast) { cout << "ERROR in make pentamer" << endl; exit(-1); }
	g.add_half_edge_type(g.v[0].vid, g.v[1].vid, 2, -1);
	g.add_half_edge_type(g.v[1].vid, g.v[0].vid, 1, 0); //boundary
	g.add_half_edge_type(g.v[1].vid, g.v[2].vid, 0, -1);
	g.add_half_edge_type(g.v[2].vid, g.v[1].vid, 3, 0); //boundary
	g.add_half_edge_type(g.v[2].vid, g.v[0].vid, 1, -1);
	g.add_half_edge_type(g.v[0].vid, g.v[2].vid, 2, 0); //boundary

	g.set_prev_next(g.he[0].id, g.he[4].id, g.he[2].id);
	g.set_prev_next(g.he[2].id, g.he[0].id, g.he[4].id);
	g.set_prev_next(g.he[4].id, g.he[2].id, g.he[0].id);

	for (vector<HE>::iterator it = g.he.begin(); it != g.he.end(); it++)
	{
		cout << "in make pentame updating edge" << it->id << endl;
		;
		g.update_half_edge(it->id);
		cout << it->id << "in make pentamer  TYPE " << g.he[it->id].type << " opid " << it->opid << " OP TYPE " << g.he[g.heidtoindex[it->opid]].type << endl;
	}

	g.update_boundary();
	//exit(-1);

	cout << " HEREH in pentamer 2" << endl;
	double *vco = new double[3];

	for (int x = 0; x < 3; x++)
	{
		g.update_boundary();
		//dump_lammps_data_file(g, frame++);

		cout << " after g.update boundary" << endl;
		g.new_vertex(g.Nhe - 1, vco);
		cout << vco[0] << " " << vco[1] << " " << vco[2] << endl;
		//exit(-1);
		g.force_add_dimer(g.Nhe - 1, vco, 0, 1);
		g.update_half_edge(g.Nhelast - 1);
		g.update_half_edge(g.Nhelast - 2);
		g.update_half_edge(g.Nhelast - 3);
		g.update_half_edge(g.Nhelast - 4);
	}
	//exit(-1);
	g.update_boundary();
	g.force_add_monomer(1, g.Nhe - 1, 0);
	g.update_half_edge(g.Nhelast - 1);
	g.update_half_edge(g.Nhelast - 2);

	delete[] vco;
	//exit(-1);
	//int vind=-1;
	g.update_boundary();
}

int check_bind_triangle(Geometry &g) //
{
	//cout << "in attempt_bind_triangle heid0 " << heid0 << endl;
	for (vector<int>::iterator it = g.boundary.begin(); it != g.boundary.end(); ++it)
	{
		int heid0 = *it;
		int heindex0 = g.heidtoindex[heid0]; // this edge on boundary
		/* if triangle */
		//int bi=g.he[heindex0].boundary_index;

		//cout << "in attempt_bind_triangle heindex0 " << heindex0 <<  " boundary_index " <<bi << endl;

		int nextboundaryid0 = g.he[heindex0].nextid_boundary;
		int prevboundaryid0 = g.he[heindex0].previd_boundary;

		if (nextboundaryid0 == -1 || prevboundaryid0 == -1)
		{
			cout << "error in attempt_bind_triangle heindex0 " << endl;
			exit(-1);
		}

		//int nextboundaryindex0=g.heidtoindex[nextboundaryid0]; // next of heid0
		int prevboundaryindex0 = g.heidtoindex[prevboundaryid0]; // prev of heid0
		if (g.he[prevboundaryindex0].previd == g.he[heindex0].nextid_boundary)
		{

			g.set_prev_next(heid0, prevboundaryid0, nextboundaryid0);
			g.set_prev_next(prevboundaryid0, nextboundaryid0, heid0);
			g.set_prev_next(nextboundaryid0, heid0, prevboundaryid0);
			g.Nboundary--;
			return 1;
		}
	}
	return 0;
}