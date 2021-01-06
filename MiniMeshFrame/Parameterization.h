#pragma once
#ifndef PARAMETERIZATION_H
#define PARAMETERIZATION_H

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Handles.hh>
#include "DifferentialGeometry.h"
#include <vector>
#include <QPainter>

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;
typedef Eigen::SparseMatrix<double> SparseMatrixType;
typedef Eigen::Triplet<double> T;
typedef Eigen::SimplicialLLT<SparseMatrixType> Solve_LLT;
typedef Eigen::SparseLU<SparseMatrixType> Solve_LU;

class Parameterization
{
public:
	Parameterization();
	Parameterization(MyMesh* ptr_mesh_, int dw, int dh);
	~Parameterization();

	void LoadDDG(DifferentialGeometry* ddg);
	void Draw(QPainter &painter, int dw = 0, int dh = 0);
	void SetTexcoord();
	virtual void parameterization();

	bool is_good_;
	MyMesh* ptr_object_;
	DifferentialGeometry* DDG;
	OpenMesh::VPropHandleT<MyMesh::Point> cogs;
protected:
	int w, h;

	void SetBoundary(int dw, int dh);
	void FitIn(float x_min, float x_max, float y_min, float y_max, bool same_scale = false);
	bool EqualVector(OpenMesh::Vec3f v1, OpenMesh::Vec3f v2);
};

class Pt_uniform :
	public Parameterization
{
public:
	Pt_uniform(MyMesh* ptr_mesh_, int dw, int dh);
	Pt_uniform(MyMesh* ptr_mesh_, int dw, int dh, DifferentialGeometry* ddg);
	~Pt_uniform();
	void parameterization();
};

class Pt_weightedleastsquares :
	public Parameterization
{
public:
	Pt_weightedleastsquares(MyMesh* ptr_mesh_, int dw, int dh);
	~Pt_weightedleastsquares();
	void parameterization();
};

class Pt_shapepreserving :
	public Parameterization
{
public:
	Pt_shapepreserving(MyMesh* ptr_mesh_, int dw, int dh);
	~Pt_shapepreserving();
	void parameterization();

private:
	OpenMesh::FPropHandleT<float> ang;
};

class Pt_ASAP :
	public Parameterization
{
public:
	Pt_ASAP(MyMesh* ptr_mesh_, int dw, int dh);
	Pt_ASAP(MyMesh* ptr_mesh_, int dw, int dh, DifferentialGeometry* ddg);
	~Pt_ASAP();
	void parameterization();
	void flatten();

private:
	std::vector<std::vector<double>> cot; 
	std::vector<std::vector<double>> diff_x;
	std::vector<std::vector<double>> diff_y;
};

class Pt_ARAP :
	public Parameterization
{
public:
	Pt_ARAP(MyMesh* ptr_mesh_, int dw, int dh, DifferentialGeometry* ddg);
	~Pt_ARAP();
	void parameterization();
	void flatten();
	void setTimes(int time);
	void LocalSetL();
	void GlobalSolveU(Solve_LLT& solver);
	void GlobalMatrixA(SparseMatrixType& A);

private:
	std::vector<std::vector<double>> cot;
	std::vector< std::vector<Eigen::Vector2f> > x_plane;
	std::vector < Eigen::Matrix2d> L;
	int times;
	std::vector<size_t> fixed_vertices_;
	std::vector<Eigen::Vector2f> fixed_coords_;
};

class Pt_SLIM :
	public Parameterization
{
public:
	Pt_SLIM(MyMesh* ptr_mesh_, int dw, int dh, DifferentialGeometry* ddg);
	~Pt_SLIM();
	void parameterization();
	void flatten();
	void setTimes(int time);
	void LocalSetL();
	void GlobalSolveU(Solve_LLT& solver);
	void UpdateGlobalMatrixA(SparseMatrixType& A);
	void LineSearch(SparseMatrixType& A);

private:
	std::vector<std::vector<double>> cot;
	std::vector< std::vector<Eigen::Vector2f> > x_plane;
	std::vector < Eigen::Matrix2d> L;
	std::vector < Eigen::Matrix2d> W;
	int times;
	std::vector<size_t> fixed_vertices_;
	std::vector<Eigen::Vector2f> fixed_coords_;
	std::vector<MyMesh::Point> p;
};

#endif