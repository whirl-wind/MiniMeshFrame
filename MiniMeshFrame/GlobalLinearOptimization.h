#pragma once
#ifndef GOLOBAL_LO_H
#define GOLOBAL_LO_H

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Handles.hh>
#include<Eigen/Sparse>
#include <vector>
#include "Cluster.h"

//using Eigen::MatrixXf;
using namespace Eigen;
using namespace Eigen::internal;
using namespace Eigen::Architecture;
using std::vector;

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;
typedef Eigen::SparseMatrix<float> SparseMatrixFType;
typedef Eigen::Triplet<float> TF;

class GlobalLinearOptimization
{
private:
	float LM_NonLinearSearch();
	float LinearOptimization();
	void Set_Vertices();

	Cluster *mesh_cluster;
	MyMesh *ptr_mesh_;
	MyMesh *ptr_mesh_base_;
	SparseMatrixFType *Eg;
	VectorXf *Cg;
	float *Ag;
	SparseMatrixFType *Ec;
	VectorXf *Cc;
	float *Ac;
	SparseMatrixFType *Eb;
	VectorXf *Cb;
	float *Ab;

	bool NonLinear;
public:
	GlobalLinearOptimization();
	GlobalLinearOptimization(Cluster *c,float opt_err);
	~GlobalLinearOptimization();

	void Run();
	void Opt_Cluster();
	float One_Step();
	float get_loss();

	Cluster *get_cluster();

	float Get_Eg(float *p, int m);
	float Get_Ec(float *p, int m);
	float Get_Eb(float *p, int m);

	vector<Vector3f> new_vectors;
	float loss;
	float t_err;
};

#endif