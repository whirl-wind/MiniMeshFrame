#include "GlobalLinearOptimization.h"
#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Handles.hh>
#include <limits.h> //for INT_MIN INT_MAX
#include <math.h>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include "levmar.h"

#ifdef _cplusplus
extern "C" {
#include "f2c.h"
#include "clapack.h"
}
#endif

#ifndef LM_DBL_PREC
#error Example program assumes that levmarhas been compiled with double precision, see LM_DBL_PREC!
#endif

#define M_PI 3.14
#undef REPEATABLE_RANDOM
#define DBL_RAND_MAX (double)(RAND_MAX)

#ifdef _MSC_VER // MSVC
#include <process.h>
#define GETPID  _getpid
#elif defined(__GNUC__) // GCC
#include <sys/types.h>
#include <unistd.h>
#define GETPID  getpid
#else
#warning Do not know the name of thefunction returning the process id for your OS / compiler combination
#define GETPID  0
#endif

#ifdef REPEATABLE_RANDOM
#define INIT_RANDOM(seed) srandom(seed)
#else
#define INIT_RANDOM(seed)srandom((int)GETPID()) // seed unused
#endif

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;
typedef Eigen::SparseMatrix<float> SparseMatrixFType;
typedef Eigen::Triplet<float> TF;

GlobalLinearOptimization::GlobalLinearOptimization()
{

}

GlobalLinearOptimization::GlobalLinearOptimization(Cluster * c)
{
	NonLinear = false;

	mesh_cluster = c;
	ptr_mesh_ = mesh_cluster->ptr_object_;
	ptr_mesh_base_ = new MyMesh();
	ptr_mesh_base_->assign(*ptr_mesh_, true);
//	ptr_mesh_->assign(*ptr_mesh_update_, true);
}

GlobalLinearOptimization::~GlobalLinearOptimization()
{
}

void GlobalLinearOptimization::Run()
{
	std::cout << "Modifying... " << std::endl;
	int MAX_C = 2000;
	float error = 0.075;

	vector<float> lo_mem;
	float lo = INT_MAX;
	int count = 0;
	lo_mem.clear();
	while (lo > error && count < MAX_C) {

		lo = One_Step();

		count++;
		std::cout << count << ") loss: " << lo << std::endl;

		lo_mem.push_back(lo);
		if (lo_mem.size() >= 10) {
			int count_lo = 0;
			for (int k = 0; k < lo_mem.size() - 1; k++) {
				if (lo_mem[k] <= lo) count_lo++;
				lo_mem[k] = lo_mem[k + 1];
			}
			if (count_lo > 5) break;
			lo_mem.pop_back();
		}
	}

	loss = lo;

	mesh_cluster->updateTriangles();
	mesh_cluster->kmeans();

	std::cout << "------------------------------------------------- " << std::endl;
	std::cout << "Cluster Num: " << mesh_cluster->class_num << std::endl;
	std::cout << "Total Loss: " << loss << std::endl;
	mesh_cluster->show_info();
	std::cout << "------------------------------------------------- \n" << std::endl;
}

float GlobalLinearOptimization::One_Step()
{
	//mesh_cluster->kmeans();
	mesh_cluster->updateCentroids();
	//mesh_cluster->show_info();
	std::cout << "Loss_cl: " << get_loss() << std::endl;
	if (NonLinear) {
		loss = LM_NonLinearSearch();
	}
	else {
		loss = LinearOptimization();
	}
	loss = get_loss();
	std::cout << "Loss_op: " << loss << std::endl;
	Set_Vertices();
	mesh_cluster->updateTriangles();
	loss = get_loss();
	std::cout << "Loss_sv: " << loss << std::endl;
	return loss;
}

float GlobalLinearOptimization::get_loss()
{
	int num = ptr_mesh_->n_vertices();
	int m = 3 * num;
	float a, b, c;
	a = 1.000; b = 0.001; c = 0.010;
	//给稀疏矩阵G,C,A赋值;
	float* p = new float[m];
	Eigen::VectorXf q(m, 1);
	int v_idx = 0;
	for (MyMesh::VertexIter v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); ++v_it) {
		p[0 + v_idx * 3] = ptr_mesh_->point(*v_it).data()[0];
		p[1 + v_idx * 3] = ptr_mesh_->point(*v_it).data()[1];
		p[2 + v_idx * 3] = ptr_mesh_->point(*v_it).data()[2];
		q(0 + v_idx * 3) = p[0 + v_idx * 3];
		q(1 + v_idx * 3) = p[1 + v_idx * 3];
		q(2 + v_idx * 3) = p[2 + v_idx * 3];
		v_idx++;
	}

	float eg, ec, eb;
	eg = Get_Eg(p, m);
	ec = Get_Ec(p, m);
	eb = Get_Eb(p, m);

	return a * eg + b * ec + c * eb;
}

Cluster * GlobalLinearOptimization::get_cluster()
{
	return mesh_cluster;
}

float GlobalLinearOptimization::Get_Eg(float *p, int m)
{
	float rotate_rate = 0.0;

	float eg = 0;
	Eg = new SparseMatrixFType(m, m);
	Eg->setZero();
	Cg = new VectorXf(m, 1);
	Cg->setZero();
	Ag = new float(0);
	SparseMatrixFType E(9, m);
	MatrixXf F(9, 1);
	std::vector<TF> tripletList;

	int *idx_;
	idx_ = new int[3];
	float *x_, *y_, *z_;
	x_ = new float[3];
	y_ = new float[3];
	z_ = new float[3];
	float *fx_, *fy_, *fz_;
	fx_ = new float[3];
	fy_ = new float[3];
	fz_ = new float[3];
	Vector3f J = Vector3f::Ones(3);

	ptr_mesh_->update_normals();
	int idx = 0;
	VectorXf pp(m, 1);
	for (int i = 0; i < m; i++) {
		pp(i, 0) = p[i];
	}
	for (MyMesh::FaceIter f_it = ptr_mesh_->faces_begin(); f_it != ptr_mesh_->faces_end(); ++f_it) {
		int i = 0;
		Vector3f n_f = Vector3f(ptr_mesh_->normal(f_it.handle()).data());
		n_f.normalize();
		for (MyMesh::FaceVertexIter fv_it = ptr_mesh_->fv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
			idx_[i] = fv_it.handle().idx();
			x_[i] = ptr_mesh_->point(ptr_mesh_->vertex_handle(idx_[i])).data()[0];
			y_[i] = ptr_mesh_->point(ptr_mesh_->vertex_handle(idx_[i])).data()[1];
			z_[i] = ptr_mesh_->point(ptr_mesh_->vertex_handle(idx_[i])).data()[2];
			fx_[i] = p[idx_[i] * 3 + 0];
			fy_[i] = p[idx_[i] * 3 + 1];
			fz_[i] = p[idx_[i] * 3 + 2];
			i++;
		}

		float S0 = 0.5 * mesh_cluster->canonical_triangles[idx](0) * mesh_cluster->canonical_triangles[idx](2);
		Vector3f a = Vector3f(x_[0] - x_[1], y_[0] - y_[1], z_[0] - z_[1]);
		Vector3f b = Vector3f(x_[2] - x_[0], y_[2] - y_[0], z_[2] - z_[0]);
		Vector3f c_na = n_f.cross(a);
		Vector3f c_nb = n_f.cross(b);
		Vector3f fx, fy, fz;
		fx = ((fx_[2] - fx_[0])*c_na + (fx_[0] - fx_[1])*c_nb) / (2 * S0);
		fy = ((fy_[2] - fy_[0])*c_na + (fy_[0] - fy_[1])*c_nb) / (2 * S0);
		fz = ((fz_[2] - fz_[0])*c_na + (fz_[0] - fz_[1])*c_nb) / (2 * S0);
		Matrix3f P;
		P.row(0) = fx.transpose();
		P.row(1) = fy.transpose();
		P.row(2) = fz.transpose();

		Matrix3f et;
		et = mesh_cluster->rigid_transformations[idx].R*mesh_cluster->Into_Triangle(mesh_cluster->canonical_centroids_triangles[mesh_cluster->clusterAssment[idx].min_Index])*mesh_cluster->rigid_transformations[idx].R_ + mesh_cluster->rigid_transformations[idx].T*J.transpose();
		
		//Rotate
		float dist = INT_MAX;
		int chose_i = 0;
		Vector3f mi = Vector3f(ptr_mesh_base_->normal(f_it.handle()).data());
		for (int k = 0; k < 3; k++) {
			float chose_d = (x_[k] - P(0))*(x_[k] - P(0)) + (y_[k] - P(1))*(y_[k] - P(1)) + (z_[k] - P(2))*(z_[k] - P(2));
			if (chose_d < dist) {
				chose_i = idx_[k];
				dist = chose_d;
			}
		}
		Vector3f ni = Vector3f((ptr_mesh_base_->vertex_normals()[chose_i]).data());
		mi.normalize();
		ni.normalize();
		Vector3f ei = mi.cross(ni);
		float thetai = asin((ei).norm());
		ei.normalize();	
		Matrix3f ei_;
		Matrix3f I;
		I << 1, 0, 0,
			 0, 1, 0,
			 0, 0, 1;
		ei_ <<	0,			-ei(2, 0), 	ei(1, 0),
				ei(2, 0),	0,			-ei(0, 0),
				-ei(1, 0),	ei(0, 0),	0;
		thetai = rotate_rate * thetai;
		Matrix3f Rotate = ei * ei.transpose() + cos(thetai) * (I - ei * ei.transpose()) + sin(thetai) * ei_;
		et = Rotate * et;

		Vector3f cx, cy, cz;
		cx = ((et(0, 2) - et(0, 0))*c_na + (et(0, 0) - et(0, 1))*c_nb) / (2 * S0);
		cy = ((et(1, 2) - et(1, 0))*c_na + (et(1, 0) - et(1, 1))*c_nb) / (2 * S0);
		cz = ((et(2, 2) - et(2, 0))*c_na + (et(2, 0) - et(2, 1))*c_nb) / (2 * S0);
		Matrix3f Q;
		Q.row(0) = cx.transpose();
		Q.row(1) = cy.transpose();
		Q.row(2) = cz.transpose();

		tripletList.clear();
		for (int k = 0; k < 3; k++) {
			tripletList.push_back(TF(k * 3 + 0, idx_[0] * 3 + k, (c_nb(0) / (2 * S0)) - (c_na(0) / (2 * S0))));
			tripletList.push_back(TF(k * 3 + 0, idx_[1] * 3 + k, -(c_nb(0) / (2 * S0))));
			tripletList.push_back(TF(k * 3 + 0, idx_[2] * 3 + k, (c_na(0) / (2 * S0))));

			tripletList.push_back(TF(k * 3 + 1, idx_[0] * 3 + k, (c_nb(1) / (2 * S0)) - (c_na(1) / (2 * S0))));
			tripletList.push_back(TF(k * 3 + 1, idx_[1] * 3 + k, -(c_nb(1) / (2 * S0))));
			tripletList.push_back(TF(k * 3 + 1, idx_[2] * 3 + k, (c_na(1) / (2 * S0))));

			tripletList.push_back(TF(k * 3 + 2, idx_[0] * 3 + k, (c_nb(2) / (2 * S0)) - (c_na(2) / (2 * S0))));
			tripletList.push_back(TF(k * 3 + 2, idx_[1] * 3 + k, -(c_nb(2) / (2 * S0))));
			tripletList.push_back(TF(k * 3 + 2, idx_[2] * 3 + k, (c_na(2) / (2 * S0))));

			F(k * 3 + 0, 0) = Q(k, 0);
			F(k * 3 + 1, 0) = Q(k, 1);
			F(k * 3 + 2, 0) = Q(k, 2);
		}
		E.setFromTriplets(tripletList.begin(), tripletList.end()); //根据三元组列表生成稀疏矩阵
		E.makeCompressed();

		*Eg += E.transpose() * E * S0;
		*Cg += -2 * E.transpose() * F * S0;
		*Ag += (F.transpose() * F)(0, 0) * S0;
		eg += (P - Q).squaredNorm() * S0;
		//std::cout << (pp.transpose()*E.transpose() * E *pp)(0, 0) * S0 - (2 * F.transpose() * E *pp)(0, 0)* S0 + (F.transpose() * F)(0, 0) * S0 << std::endl;
		//std::cout << (P - Q).squaredNorm() * S0 << std::endl;
		idx++;
	}

	*Eg *= 2.0;
	return eg;
}

float GlobalLinearOptimization::Get_Ec(float * p, int m)
{
	float ec = 0;
	Ec = new SparseMatrixFType(m, m);
	Ec->setZero();
	Cc = new VectorXf(m, 1);
	Cc->setZero();
	Ac = new float(0);
	SparseMatrixFType J(3, m);
	std::vector<TF> tripletList;

	int *idx_;
	idx_ = new int[3];
	float *x_, *y_, *z_;
	x_ = new float[3];
	y_ = new float[3];
	z_ = new float[3];
	float *fx_, *fy_, *fz_;
	fx_ = new float[3];
	fy_ = new float[3];
	fz_ = new float[3];

	int idx = 0;
	for (MyMesh::FaceIter f_it = ptr_mesh_base_->faces_begin(); f_it != ptr_mesh_base_->faces_end(); ++f_it) {
		int i = 0;
		Vector3f n_f;
		for (MyMesh::FaceVertexIter fv_it = ptr_mesh_base_->fv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
			idx_[i] = fv_it.handle().idx();
			x_[i] = ptr_mesh_base_->point(ptr_mesh_base_->vertex_handle(idx_[i])).data()[0];
			y_[i] = ptr_mesh_base_->point(ptr_mesh_base_->vertex_handle(idx_[i])).data()[1];
			z_[i] = ptr_mesh_base_->point(ptr_mesh_base_->vertex_handle(idx_[i])).data()[2];
			fx_[i] = p[idx_[i] * 3 + 0];
			fy_[i] = p[idx_[i] * 3 + 1];
			fz_[i] = p[idx_[i] * 3 + 2];
			i++;
		}

		Vector3f P = Vector3f((fx_[0] + fx_[1] + fx_[2]) / 3.0, (fy_[0] + fy_[1] + fy_[2]) / 3.0, (fz_[0] + fz_[1] + fz_[2]) / 3.0);
		
		float dist = INT_MAX;
		int chose_i = 0;
		for (int k = 0; k < 3; k++) {
			float chose_d = (x_[k] - P(0))*(x_[k] - P(0)) + (y_[k] - P(1))*(y_[k] - P(1)) + (z_[k] - P(2))*(z_[k] - P(2));
			if (chose_d < dist) {
				chose_i = idx_[k];
				//n_f = Vector3f((ptr_mesh_base_->vertex_normals()[chose_i]).data());
				dist = chose_d;
			}
		}
/*		for (MyMesh::FaceIter f_it_2 = ptr_mesh_base_->faces_begin(); f_it_2 != ptr_mesh_base_->faces_end(); ++f_it_2) {
			for (MyMesh::FaceVertexIter fv_it_2 = ptr_mesh_base_->fv_iter(*f_it_2); fv_it_2.is_valid(); ++fv_it_2) {
				float chose_d = (Vector3f(ptr_mesh_base_->point(fv_it_2.handle()).data()) - P).squaredNorm();;
				if (chose_d < dist) {
					chose_i = fv_it_2.handle().idx();
					n_f = Vector3f(ptr_mesh_base_->normal(fv_it_2.handle()).data());
					dist = chose_d;
				}
			}
		}
		for (MyMesh::VertexIter v_it = ptr_mesh_base_->vertices_begin(); v_it != ptr_mesh_base_->vertices_end(); ++v_it) {
			float chose_d = (Vector3f(ptr_mesh_base_->point(v_it.handle()).data()) - P).squaredNorm();;
			if (chose_d < dist) {
				chose_i = v_it.handle().idx();
				//n_f = Vector3f(ptr_mesh_base_->normal(v_it.handle()).data());
				dist = chose_d;
			}
			//std::cout << Vector3f(ptr_mesh_base_->point(v_it.handle()).data()).transpose() << std::endl;
		}
*/
		n_f = Vector3f(ptr_mesh_base_->normal(f_it.handle()).data());
		Vector3f xi = Vector3f(ptr_mesh_base_->point(ptr_mesh_base_->vertex_handle(chose_i)).data());
		float ec_i = (n_f.transpose() * (xi - P))(0, 0);

		tripletList.clear();
		for (int k = 0; k < 3; k++) {
			tripletList.push_back(TF(k, idx_[0] * 3 + k, 1 / 3.0));
			tripletList.push_back(TF(k, idx_[1] * 3 + k, 1 / 3.0));
			tripletList.push_back(TF(k, idx_[2] * 3 + k, 1 / 3.0));
		}
		J.setFromTriplets(tripletList.begin(), tripletList.end()); //根据三元组列表生成稀疏矩阵
		J.makeCompressed();

		*Ec += ((J.transpose()) * n_f.sparseView()) * (n_f.transpose().sparseView() * (J));
		*Cc += -2 * J.transpose() * n_f * n_f.transpose() * xi;
		*Ac += (xi.transpose() * n_f * n_f.transpose() * xi)(0,0);

		ec += ec_i * ec_i;
		idx++;
	}

	*Ec *= 2.0;
	return ec;
}

float GlobalLinearOptimization::Get_Eb(float * p, int m)
{
	float eb = 0;
	Eb = new SparseMatrixFType(m, m);
	Eb->setZero();
	Cb = new VectorXf(m, 1);
	Cb->setZero();
	Ab = new float(0);
	SparseMatrixFType J(3, m);
	std::vector<TF> tripletList;

	int b_num = 0;
	MyMesh::VertexHandle b_vert;
	for (MyMesh::VertexIter v_it = ptr_mesh_base_->vertices_begin(); v_it != ptr_mesh_base_->vertices_end(); v_it++) {
		if (ptr_mesh_base_->is_boundary(*v_it)) {
			if (b_num == 0) b_vert = v_it.handle();
			b_num++;
		}
	}

	if (b_num == 0) {
		return eb;
	}

	MyMesh::HalfedgeHandle b_edge = ptr_mesh_base_->halfedge_handle(b_vert);
	b_vert = ptr_mesh_base_->to_vertex_handle(b_edge);
	while (ptr_mesh_base_->face_handle(b_edge).is_valid()) {
		b_edge = ptr_mesh_base_->next_halfedge_handle(ptr_mesh_base_->opposite_halfedge_handle(b_edge));
		b_vert = ptr_mesh_base_->to_vertex_handle(b_edge);
	}
	b_vert = ptr_mesh_base_->from_vertex_handle(b_edge);

	int idx;
	float fx, fy, fz;
	float y1_x, y1_y, y1_z;
	float y2_x, y2_y, y2_z;
	float y3_x, y3_y, y3_z;
	MyMesh::VertexHandle begin_vert = b_vert;
	do {
		b_vert = ptr_mesh_base_->from_vertex_handle(b_edge);

		MyMesh::VertexHandle x0_vert = ptr_mesh_base_->to_vertex_handle(b_edge);
		MyMesh::VertexHandle x2_vert = ptr_mesh_base_->to_vertex_handle(ptr_mesh_base_->next_halfedge_handle(b_edge));
		idx = x0_vert.idx();
		fx = p[idx * 3 + 0];
		fy = p[idx * 3 + 1];
		fz = p[idx * 3 + 2];
		y1_x = ptr_mesh_base_->point(x0_vert).data()[0];
		y1_y = ptr_mesh_base_->point(x0_vert).data()[1];
		y1_z = ptr_mesh_base_->point(x0_vert).data()[2];
		y2_x = ptr_mesh_base_->point(x2_vert).data()[0];
		y2_y = ptr_mesh_base_->point(x2_vert).data()[1];
		y2_z = ptr_mesh_base_->point(x2_vert).data()[2];
		y3_x = ptr_mesh_base_->point(b_vert).data()[0];
		y3_y = ptr_mesh_base_->point(b_vert).data()[1];
		y3_z = ptr_mesh_base_->point(b_vert).data()[2];
		if ((y2_x - fx)*(y2_x - fx) + (y2_y - fy)*(y2_y - fy) + (y2_z - fz)*(y2_z - fz) > (y3_x - fx)*(y3_x - fx) + (y3_y - fy)*(y3_y - fy) + (y3_z - fz)*(y3_z - fz)) {
			y2_x = y3_x;
			y2_y = y3_y;
			y2_z = y3_z;
		}

		Vector3f pl = Vector3f(fx, fy, fz);
		Vector3f y1 = Vector3f(y1_x, y1_y, y1_z);
		Vector3f y2 = Vector3f(y2_x, y2_y, y2_z);
		float pr = ((((pl - y1).transpose()*(y2 - y1))(0, 0)) / ((y2 - y1).squaredNorm()));

		Vector3f yy = -y1 + ((y1.transpose()*(y2 - y1))(0, 0) / ((y2 - y1).squaredNorm()))*(y2 - y1);
		tripletList.clear();
		for (int k = 0; k < 3; k++) {
			tripletList.push_back(TF(k, idx * 3 + k, 1));
		}
		J.setFromTriplets(tripletList.begin(), tripletList.end()); //根据三元组列表生成稀疏矩阵
		J.makeCompressed();

		eb += (pl - y1 - pr * (y2 - y1)).squaredNorm();
		*Eb += (J - ((y2 - y1)*(y2 - y1).transpose()).sparseView() / ((y2 - y1).squaredNorm()) * J).transpose() * (J - ((y2 - y1)*(y2 - y1).transpose()).sparseView() / ((y2 - y1).squaredNorm()) * J);
		*Cb += 2 * (J - ((y2 - y1)*(y2 - y1).transpose()) / ((y2 - y1).squaredNorm()) * J).transpose() * yy;
		*Ab += (yy.transpose() * yy)(0, 0);

		b_edge = ptr_mesh_base_->next_halfedge_handle(b_edge);
	} while (b_vert != begin_vert);
	
	*Eb *= 2.0;
	return eb;
}

struct glodata {
	GlobalLinearOptimization *g;
};

void Loss(float *p, float *x, int m, int n, void *data)
{
	register int i;
	struct glodata *dat;
	dat = (struct glodata *)data;
	float eg = dat->g->Get_Eg(p, m);
	float ec = dat->g->Get_Ec(p, m);
	float eb = dat->g->Get_Eb(p, m);
	float a, b, c;
	a = 1.000; b = 0.001; c = 0.010;
	for (i = 0; i < n; ++i)
	{
		x[i] = a * eg + b * ec + c * eb;
	}
}

float GlobalLinearOptimization::LM_NonLinearSearch()
{
	int m = 3 * mesh_cluster->ptr_object_->n_vertices(), n = 3 * mesh_cluster->ptr_object_->n_vertices();     //m表示待求参数的维度，n表示测量值的维度

	float *p = new float[m];
	float *x = new float[n];
	float opts[LM_OPTS_SZ], info[LM_INFO_SZ];
	opts[0] = LM_INIT_MU;
	struct glodata data;
	data.g = this;
	//待求参数的初值
	int v_idx = 0;
	for (MyMesh::VertexIter v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); ++v_it) {
		p[0 + v_idx * 3] = ptr_mesh_->point(*v_it).data()[0];
		p[1 + v_idx * 3] = ptr_mesh_->point(*v_it).data()[1];
		p[2 + v_idx * 3] = ptr_mesh_->point(*v_it).data()[2];
		v_idx++;
	}
	//模拟n次测量的结果，由于是最小值，设为零
	for (int i = 0; i < n; i++) x[i] = 0.0;
	// 调用迭代入口函数
	int ret = slevmar_dif(
		Loss,      //描述测量值之间关系的函数指针
		p,          //初始化的待求参数，结果一并保存在其中
		x,          //测量值
		m,          //参数维度
		n,          //测量值维度
		1000,       //最大迭代次数
		opts,       //迭代的一些参数
		info,       //关于最小化结果的一些参数，不需要设为NULL
		NULL, NULL,//一些内存的指针，暂时不需要，以后再学习这个具体由什么用
		(void *)&data
	);
	
	printf("Levenberg-Marquardtreturned in %g iter, reason %g, sumsq %g [%g]\n", info[5], info[6], info[1], info[0]);

	new_vectors.clear();
	for (int i = 0; i < m / 3; i++) {
		new_vectors.push_back(Vector3f(p[0 + 3 * i], p[1 + 3 * i], p[2 + 3 * i]));
	}

	return info[1];
}

float GlobalLinearOptimization::LinearOptimization()
{
	float lam = 1;
	float a, b, c;
	a = 1.000; b = 0.001; c = 0.010;
	int row, col;
	int num = ptr_mesh_->n_vertices();
	int m = 3 * num;
	row = col = m;

	//声明方程组的变量;
	//0.5*Xt*G*X-Ct*X+A
	SparseMatrixFType G(row, col);
	VectorXf C(row, 1);
	float A;

	//给稀疏矩阵G,C,A赋值;
	float *p = new float[m];
	Eigen::VectorXf q(m, 1);
	int v_idx = 0;
	for (MyMesh::VertexIter v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); ++v_it) {
		p[0 + v_idx * 3] = ptr_mesh_->point(*v_it).data()[0];
		p[1 + v_idx * 3] = ptr_mesh_->point(*v_it).data()[1];
		p[2 + v_idx * 3] = ptr_mesh_->point(*v_it).data()[2];
		q(0 + v_idx * 3) = p[0 + v_idx * 3];
		q(1 + v_idx * 3) = p[1 + v_idx * 3];
		q(2 + v_idx * 3) = p[2 + v_idx * 3];
		v_idx++;
	}

	float eg, ec, eb;
	eg = Get_Eg(p, m);
	ec = Get_Ec(p, m);
	eb = Get_Eb(p, m);

	G =  a * *Eg + b * *Ec + c * *Eb;
	C = -a * *Cg - b * *Cc - c * *Cb;
	A =  a * *Ag + b * *Ac + c * *Ab;

	//求解
	G.makeCompressed();
	//loss = (0.5 * q.transpose() * G * q - C.transpose() * q)(0, 0) + A;
	//std::cout << "Loss0: " << loss << std::endl;
	
//	Eigen::SparseQR<SparseMatrixFType, AMDOrdering < int > > solver;
	Eigen::SparseLU<SparseMatrixFType> solver;
	solver.compute(G);
	Eigen::VectorXf min_q(m, 1);
	min_q = solver.solve(C);

	if (solver.info() != Success) std::cout << "Oh, Very bad!" << std::endl;
	//else std::cout << "Succeed!" << std::endl;

	q = (1 - lam) * q + lam * min_q;

	new_vectors.clear();
	for (int i = 0; i < m / 3; i++) {
		new_vectors.push_back(Vector3f(q(0 + 3 * i), q(1 + 3 * i), q(2 + 3 * i)));
	}
	//loss = a * eg + b * ec + c * eb;
	//std::cout << "Loss1: " << loss << std::endl;
	loss = (0.5 * q.transpose()*G*q - C.transpose()*q)(0, 0)+A;
	//std::cout << "Loss: "<< loss <<std::endl;

	return loss;
	
}

void GlobalLinearOptimization::Set_Vertices()
{
	//赋值
	int i = 0;
	for (auto v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); v_it++) {
		MyMesh::Point new_vert(new_vectors[i](0), new_vectors[i](1), new_vectors[i](2));
		ptr_mesh_->set_point(*v_it, new_vert);
		i++;
	}

	ptr_mesh_->update_normals();
}

