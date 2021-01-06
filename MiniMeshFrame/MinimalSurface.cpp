#include "MinimalSurface.h"
#include<Eigen/Sparse>
#include<Eigen/SparseLU>
#include<iostream>

typedef Eigen::SparseMatrix<double> SparseMatrixType; 
typedef Eigen::Triplet<double> T; 
typedef Eigen::SparseLU<SparseMatrixType> Solve_LU; //����ԭ���õ�Eigen::SimplicialCholesky���Ͳ��У���֪��Ϊʲô��

MinimalSurface::MinimalSurface()
{
	ptr_object_ = new MyMesh();
}

MinimalSurface::MinimalSurface(MyMesh * ptr_mesh_)
{
	ptr_object_ = ptr_mesh_;
}


MinimalSurface::~MinimalSurface()
{
}

void MinimalSurface::global_methods()
{
	if (ptr_object_->n_faces() == 0)return;

	int row, col;
	int num = ptr_object_->n_vertices();
	row = col = num;
	
	//����������ı���;	
	SparseMatrixType A(row, col);
	Eigen::VectorXd x;
	Eigen::VectorXd y;
	Eigen::VectorXd z;
	Eigen::VectorXd X;
	Eigen::VectorXd Y;
	Eigen::VectorXd Z;
	std::vector<T> tripletList;
	
	//������b��ֵ;	
	X.resize(row);
	Y.resize(row);
	Z.resize(row);
	int i = 0;
	for (auto v_it = ptr_object_->vertices_begin(); v_it != ptr_object_->vertices_end(); v_it++) {
		if (ptr_object_->is_boundary(*v_it)) {
			auto vert = ptr_object_->point(v_it.handle());
			X(i) = vert[0];
			Y(i) = vert[1];
			Z(i) = vert[2];
		}
		else {
			X(i) = 0;
			Y(i) = 0;
			Z(i) = 0;
		}
		i++;
	}

	//��ϡ�����A��ֵ;
	tripletList.reserve(num);
	i = 0;
	for (auto v_it = ptr_object_->vertices_begin(); v_it != ptr_object_->vertices_end(); v_it++) {
		if (ptr_object_->is_boundary(*v_it)) tripletList.push_back(T(i, i, 1));
		else {
			int nb_num = 0;
			for (auto vv_it = ptr_object_->vv_iter(*v_it); vv_it.is_valid(); ++vv_it) {
				tripletList.push_back(T(i, vv_it.handle().idx(), -1)); ;
				nb_num++;
			}
			tripletList.push_back(T(i, i, nb_num));
		}
		i++;
	}
	A.setFromTriplets(tripletList.begin(), tripletList.end()); //������Ԫ���б�����ϡ�����
	A.makeCompressed();

	//���;	
	Solve_LU solver;
	solver.compute(A);
	if (solver.info() != Eigen::Success) {
		std::cout << "decomposition failed!" << std::endl;
		return;
	}
	x = solver.solve(X);
	if (solver.info() != Eigen::Success) {
		std::cout << "solving failed!" << std::endl;
		return;
	}
	y = solver.solve(Y);
	if (solver.info() != Eigen::Success) {
		std::cout << "solving failed!" << std::endl;
		return;
	}
	z = solver.solve(Z);
	if (solver.info() != Eigen::Success) {
		std::cout << "solving failed!" << std::endl;
		return;
	}

	//��ֵ
	i = 0;
	for (auto v_it = ptr_object_->vertices_begin(); v_it != ptr_object_->vertices_end(); v_it++) {
		MyMesh::Point new_vert(x[i],y[i],z[i]);
		ptr_object_->set_point(*v_it, new_vert);
		i++;
	}

	ptr_object_->update_normals();
}

void MinimalSurface::local_methods()
{
	if (ptr_object_->n_faces() == 0)return;

	int N = 3000; //��������
	float r = 0.3; //ѧϰ��

	OpenMesh::VPropHandleT<MyMesh::Point> cogs;
	ptr_object_->add_property(cogs);
	MyMesh::Scalar valence;

	for (int i = 0; i < N; i++) {
		for (auto v_it = ptr_object_->vertices_begin(); v_it != ptr_object_->vertices_end(); v_it++) {
			ptr_object_->property(cogs, *v_it).vectorize(0.0f);
			valence = 0.0;
			for (auto vv_it = ptr_object_->vv_iter(*v_it);  vv_it.is_valid(); vv_it++) {
				ptr_object_->property(cogs, *v_it) += ptr_object_->point(*vv_it);
				++valence;
			}
			ptr_object_->property(cogs, *v_it) /= valence;
		}
		for (auto v_it = ptr_object_->vertices_begin(); v_it != ptr_object_->vertices_end(); v_it++)
			if (!ptr_object_->is_boundary(*v_it)) {
				MyMesh::Point new_point = ptr_object_->point(*v_it) + r * (ptr_object_->property(cogs, *v_it) - ptr_object_->point(*v_it));
				ptr_object_->set_point(*v_it, new_point);
			}
	}
	ptr_object_->remove_property(cogs);
	ptr_object_->update_normals();
}