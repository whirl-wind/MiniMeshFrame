#pragma once
#ifndef MINIMALSURFACE_H
#define MINIMALSURFACE_H

#include <OpenMesh/Core/IO/MeshIO.hh>
//#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Handles.hh>
#include <vector>

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;
using std::vector;

class MinimalSurface
{
public:
	MinimalSurface();
	MinimalSurface(MyMesh *ptr_mesh_);
	~MinimalSurface();

	void local_methods();
	void global_methods();
private:
	MyMesh *ptr_object_;
};

#endif