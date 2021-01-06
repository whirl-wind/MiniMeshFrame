#pragma once
#ifndef TRIANGULATION_H_
#define TRIANGULATION_H_

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Handles.hh>
#include <vector>

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

class Triangulation
{
public:
	Triangulation();
	~Triangulation();

	virtual void triangulate();

	MyMesh* ptr_mesh_;
protected:
	bool is_flatten_;

};

class DelaunayTriangulation :
	public Triangulation
{
public:
	DelaunayTriangulation();
	DelaunayTriangulation(MyMesh* ptr);
	~DelaunayTriangulation();

	void triangulate();
private:
	bool is_inside_tri(MyMesh::VertexHandle p, MyMesh::VertexHandle v1, MyMesh::VertexHandle v2, MyMesh::VertexHandle v3, bool strict);
	void init_convexhull();
	void init_tri();
	OpenMesh::VPropHandleT< bool > is_added;
	std::vector<MyMesh::VertexHandle> convexhull;
};

#endif // !TRIANGULATION_H_