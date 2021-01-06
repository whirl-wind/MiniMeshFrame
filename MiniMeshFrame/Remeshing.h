#pragma once
#ifndef REMESHING_H_
#define REMESHING_H_

#include <vector>
#include <cmath>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Handles.hh>
#include "DifferentialGeometry.h"

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

class Remeshing
{
public:
	Remeshing();
	~Remeshing();
	void LoadDDG(DifferentialGeometry* ddg);
	virtual void remesh();

	MyMesh* ptr_mesh_;
	DifferentialGeometry* DDG;
};

class IncrementalRemeshing :
	public Remeshing 
{
public:
	IncrementalRemeshing(MyMesh* ptr, DifferentialGeometry* ddg);
	~IncrementalRemeshing();

	void remesh(double target_edge_length);
	void split_long_edges(double high);
	void collapse_short_edges(double low, double high);
	void equalize_valences();
	void tangential_relaxation();
	void project_to_surface();
};


#endif // !REMESHING_H_
