#pragma once
#ifndef CLUSTER_H
#define CLUSTER_H

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Handles.hh>
#include "colorbar.h"
#include <Eigen/Dense>
#include <vector>

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

//using Eigen::MatrixXf;
using namespace Eigen;
using namespace Eigen::internal;
using namespace Eigen::Architecture;
using std::vector;

enum NonlinearType
{
	kDefault = 0,
	kT = 1,
	kLM = 2,
};

class Cluster
{
private:
	void Set_Canonical_Triangle();
	bool addCentroids();
	void init();
	void updateCluster();
	void LM_NonLinearSearch();

	float distTriangles(Matrix3f A, Matrix3f B);
	float distTriangles(Vector3f A, Vector3f B);

	vector<int> centroids_index;

public:
	Cluster();
	Cluster(MyMesh *ptr_mesh_,int cla);
	~Cluster();

	int get_minIdx(int idx_f);
	void kmeans();
	void updateCentroids();
	void updateTriangles();
	void Set_Rigid_Transformation();

	float distTriangles(Vector3f A, Matrix3f B); 
	Vector3f Cluster::Into_Canonical_Triangle(Matrix3f t);
	Matrix3f Cluster::Into_Triangle(Vector3f t);
	void Set_Faces_Colors();
	void show_info();
	void Draw();

	typedef struct Node
	{
		int min_Index; //the index of each node
		float min_Dist;
		Node(int idx, float dist) :min_Index(idx), min_Dist(dist) {}
	}tNode;

	typedef struct transNode
	{
		Vector3f T;
		Matrix3f R;
		Matrix3f R_;
		transNode(Vector3f t, Matrix3f r, Matrix3f r_) :T(t), R(r), R_(r_) {}
	}transNode;

	vector<tNode>  clusterAssment;
	vector<Matrix3f> triangles;
	vector<Vector3f> canonical_triangles;
	vector<Vector3f> canonical_centroids_triangles;
	vector<transNode> rigid_transformations;
	NonlinearType nonlinear_type_;
	MyMesh *ptr_object_;
	vector<int> clusterNumber;
	int class_num;
};

#endif