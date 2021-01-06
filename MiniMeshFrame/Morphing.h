#ifndef MORPHING_H_
#define MORPHING_H_

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Handles.hh>
#include "DifferentialGeometry.h"
#include <vector>

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;
typedef Eigen::SparseMatrix<double> SparseMatrixType;
typedef Eigen::Triplet<double> T;
typedef Eigen::SimplicialLLT<SparseMatrixType> Solve_LLT;
typedef Eigen::SparseLU<SparseMatrixType> Solve_LU;


class Morphing
{
public:
	Morphing();
	Morphing(MyMesh* ptr_mesh_, MyMesh* ptr_target_, float t = 0.01);
	~Morphing();
	void LoadDDG(DifferentialGeometry* ddg);
	virtual void morphing(float t_now);
	virtual void Init();

	float t; //step
protected:
	MyMesh* ptr_mesh_;
	MyMesh* ptr_target_mesh_;
	DifferentialGeometry* DDG;
};

class Mp_ARAP :
	public Morphing
{
public:
	Mp_ARAP(MyMesh* ptr_mesh_, MyMesh* ptr_target_, DifferentialGeometry* ddg, float t = 0.01);
	~Mp_ARAP();

	virtual void morphing(float t_now, int i);
	virtual void Init();
	void LocalSetL();
	void GlobalSolveU(Solve_LLT& solver, float t_now, int i);
private:
	Solve_LLT solver;
	std::vector < Eigen::Matrix3d> L;
	std::vector < Eigen::Matrix3d> S;
	std::vector<size_t> fixed_vertices_;
	std::vector<Eigen::Vector3f> fixed_position_;

	//MyMesh* ptr_backup_;
};

#endif // !MORPHING_H_
