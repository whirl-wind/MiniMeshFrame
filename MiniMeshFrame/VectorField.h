#pragma once
#ifndef VECTORFIELD_H_
#define VECTORFIELD_H_

#include <vector>
#include <cmath>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Handles.hh>

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

using namespace std;

class VectorField
{
public:
	VectorField();
	VectorField(MyMesh* ptr, int mm = 1, int nn = 1); // m: degree  n: free
	~VectorField();

	void BuildLocal();
	virtual void SetFix(vector < size_t > fidx, vector < vector< double > > fvectorfield);
	virtual void UpdateMesh(MyMesh* ptr);
	virtual void BuildVectorField();
	void SetColor();
	void SetColorIter(MyMesh::FaceHandle f_it, MyMesh::HalfedgeHandle fh_it);

	virtual void Draw(float size);

	int m, n;
	MyMesh* ptr_mesh_;
	OpenMesh::FPropHandleT< vector< double > > vectorfield_theta;
	OpenMesh::FPropHandleT< vector< size_t > > vectorfield_idx;
	OpenMesh::HPropHandleT< std::complex<double> > localcoord;

protected:
	vector<size_t> fix_idx;
	vector< vector<double> > fix_vectorfield;
};

class NpolyVectorFields :
	public VectorField
{
public:
	NpolyVectorFields(MyMesh* ptr, int nn = 1, int mm = 1);
	~NpolyVectorFields();

	void BuildLaplace();
	void BuildVectorField();
	void UpdateMesh(MyMesh* ptr);

private:
	vector< Eigen::SparseMatrix< std::complex<double> > > Lm;
	//vector< Eigen::SparseLU< Eigen::SparseMatrix< std::complex<double> > >* > solver;
};

#endif // !VECTORFIELD_H_