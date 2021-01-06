#pragma once
#ifndef DEFORMATION_H
#define DEFORMATION_H

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Handles.hh>
#include "DifferentialGeometry.h"
#include "Interaction.h"
#include <vector>

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;
typedef Eigen::SparseMatrix<double> SparseMatrixType;
typedef Eigen::Triplet<double> T;
typedef Eigen::SimplicialLLT<SparseMatrixType> Solve_LLT;
typedef Eigen::SparseLU<SparseMatrixType> Solve_LU;

class Deformation
{
public:
	Deformation();
	Deformation(MyMesh* ptr_mesh_, vector<Interaction*>	IXs, Interaction* select_ix);
	~Deformation();
	void LoadDDG(DifferentialGeometry* ddg);
	virtual void deform();
	virtual void Init();
	virtual void update_IXs(vector<Interaction*> IXs, Interaction* select_ix);

	MyMesh* ptr_mesh_;
	DifferentialGeometry* DDG;
	OpenMesh::VPropHandleT<bool> fixed;
protected:
	
};

class Df_ARAP :
	public Deformation
{
public:
	Df_ARAP(MyMesh* ptr_mesh_, vector<Interaction*>	IXs, Interaction* select_ix, DifferentialGeometry* ddg);
	~Df_ARAP();
	void deform();
	void Init();
	void update_IXs(vector<Interaction*>	IXs, Interaction* select_ix);
	
	void setTimes(int time);
	void LocalSetL();
	void GlobalSolveU(Solve_LLT& solver);
	void GlobalMatrixA(SparseMatrixType& A);

private:
	std::vector < Eigen::Matrix3d> L;
	int times;
	Solve_LLT solver;

};


#endif //DEFORMATION_H