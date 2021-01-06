#ifndef SMOOTHING_H_
#define SMOOTHING_H_

#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Handles.hh>
#include "DifferentialGeometry.h"

typedef Eigen::SparseMatrix<double> SparseMatrixType;
typedef Eigen::Triplet<double> T;
typedef Eigen::SparseLU<SparseMatrixType> Solve_LU;
typedef Eigen::SimplicialLLT<SparseMatrixType> Solve_LLT;

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

using std::vector;

class Smoothing
{
public:
	Smoothing();
	virtual ~Smoothing();

	void LoadMesh(MyMesh* ptr_);
	void LoadDDG(DifferentialGeometry *ddg);
	void AddNoise(double scale);
	virtual void SmoothMesh();

	MyMesh* ptr_mesh_;
protected:
	DifferentialGeometry *DDG;
};

class BilateralMeshDenoising: public Smoothing 
{
public:
	BilateralMeshDenoising();
	~BilateralMeshDenoising();

	void Set_LocalParameter(int num=20, double sig=0.5);
	void Set_GlobalParameter(double n_lam=0.5, double p_lam=0.5, double sig = 0.5);
	void SmoothMesh();
	void ScaleMesh(double rate);
	void SetGlobal(bool is_global_);

	vector<Eigen::Vector3d> fnormal;
private:
	void update_NormalField();
	void update_VertexPosition();
	void Local_Iterative();
	void Global_EnergyMinimization();

	int iter_num;

	bool is_global_;
	double sig_s;
	double sig_c;
	double n_lam;
	double p_lam;
};






#endif // !SMOOTHING_H_

