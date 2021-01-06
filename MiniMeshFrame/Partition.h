#pragma once
#ifndef PARTITION_H_
#define PARTITION_H_

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Handles.hh>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#include <vector>
#include <queue>
#include <iostream>
#include "DifferentialGeometry.h"

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

class Partition
{
public:
	Partition();
	~Partition();
	void LoadDDG(DifferentialGeometry* ddg);
	int get_classnum();
	virtual void partition();
	virtual void Draw();
protected:
	int classnum;
	MyMesh* ptr_mesh_;
	DifferentialGeometry* DDG;
};

class Pa_VSA ://Variational Shape Approximation
	public Partition
{
public:
	Pa_VSA(MyMesh *ptr, DifferentialGeometry* ddg, int n);
	~Pa_VSA();

	void init();
	float proxy();
	void partition();
	void sort_into_oneclass(std::queue<MyMesh::FaceHandle> &oneclass, int classidx);
	void Draw();
private:
	std::vector<Eigen::Vector3f> P;
	std::vector< std::vector<int> > R;
	std::vector<bool> conqueue;

	OpenMesh::FPropHandleT< int > label;
};

#endif // !PARTITION_H_
