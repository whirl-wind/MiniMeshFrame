#include "Parameterization.h"
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#include <cmath>
#include <iostream>

Parameterization::Parameterization()
{
	ptr_object_ = new MyMesh();
	ptr_object_->add_property(cogs);
	is_good_ = true;
}

Parameterization::Parameterization(MyMesh* ptr_mesh_, int dw, int dh)
{
	ptr_object_ = ptr_mesh_;
	ptr_object_->add_property(cogs);
	w = dw;
	h = dh;
	is_good_ = true;
	SetBoundary(dw, dh);
}

Parameterization::~Parameterization()
{
	ptr_object_->remove_property(cogs);
}

void Parameterization::SetBoundary(int dw, int dh)
{
	int b_num = 0;
	MyMesh::VertexHandle b_vert;
	for (MyMesh::VertexIter v_it = ptr_object_->vertices_begin(); v_it != ptr_object_->vertices_end(); v_it++) {
		ptr_object_->property(cogs, *v_it).vectorize(0.0f);
		if (ptr_object_->is_boundary(*v_it)) {
			if (b_num == 0) b_vert = v_it.handle();
			b_num++;
		}
	}
	if (b_num == 0) {
		is_good_ = false;
		return;
	}
	MyMesh::HalfedgeHandle b_edge = ptr_object_->halfedge_handle(b_vert);
	b_vert = ptr_object_->to_vertex_handle(b_edge);
	while(ptr_object_->face_handle(b_edge).is_valid()) {
		b_edge = ptr_object_->next_halfedge_handle(ptr_object_->opposite_halfedge_handle(b_edge));
		b_vert = ptr_object_->to_vertex_handle(b_edge);
	} 
	//std::cout << "Boundary verts: " << b_num << std::endl;

	int d = 2 * (dw + dh) / b_num;
	if (dw % d < d / 2 && dw % d > 0 && dh % d < d / 2 && dh % d > 0) d = 2 * (dw + dh) / (b_num + 4);

	int pattern = 0;
	//int count = 0;
	int x, y;
	x = y = 0;

	b_vert = ptr_object_->from_vertex_handle(b_edge);
	MyMesh::VertexHandle begin_vert = b_vert;
	do{
		bool s = true;
		if (x == 0 && y == 0) s = false;
		MyMesh::Point pt_vert(x, y, 0);
		ptr_object_->property(cogs, b_vert) = pt_vert;

		switch (pattern) {
			case 0:
				x += d;
				if (x >= dw && x - dw < d / 2) {
					x = dw;
					pattern++;
				}
				else if (x - dw > d / 2) {
					x = dw;
					b_vert = ptr_object_->from_vertex_handle(b_edge);
					//count--;
					s = false;
					MyMesh::Point pt_vert_(x, y, 0);
					ptr_object_->property(cogs, b_vert) = pt_vert_;
					y += d;
					pattern++;
				}
				break;
			case 1:
				y += d;
				if (y >= dh && y - dh < d / 2) {
					y = dh;
					pattern++;
				}
				else if (y - dh > d / 2) {
					y = dh;
					b_vert = ptr_object_->from_vertex_handle(b_edge);
					//count--;
					s = false;
					MyMesh::Point pt_vert_(x, y, 0);
					ptr_object_->property(cogs, b_vert) = pt_vert_;
					x -= d;
					pattern++;
				}
				break;
			case 2:
				x -= d;
				if (x <= 0 && -x < d / 2) {
					x = 0;
					pattern++;
				}
				else if (-x > d / 2) {
					x = 0;
					b_vert = ptr_object_->from_vertex_handle(b_edge);
					//count--;
					s = false;
					MyMesh::Point pt_vert_(x, y, 0);
					ptr_object_->property(cogs, b_vert) = pt_vert_;
					y -= d;
					pattern++;
				}
				break;
			case 3:
				y -= d;
				if (y <= 0 && -y < d / 2) {
					y = 0;
					pattern++;
				}
				else if (-y > d / 2) {
					y = 0;
					b_vert = ptr_object_->from_vertex_handle(b_edge);
					//count--;
					s = false;
					MyMesh::Point pt_vert_(x, y, 0);
					ptr_object_->property(cogs, b_vert) = pt_vert_;
					pattern++;
				}
				break;
			default:
				break;
		}
		if(s) b_edge = ptr_object_->next_halfedge_handle(b_edge);
		b_vert = ptr_object_->to_vertex_handle(b_edge);
		//count++;
		//if (count >= b_num) break;
	}while (b_vert != begin_vert);
}

void Parameterization::FitIn(float x_min, float x_max, float y_min, float y_max, bool same_sacle)
{
	float scale_x = w / (x_max - x_min);
	float scale_y = h / (y_max - y_min);
	if (same_sacle) {
		if (scale_x > scale_y) scale_x = scale_y;
		else scale_y = scale_x;
	}
	for (auto v_it = ptr_object_->vertices_begin(); v_it != ptr_object_->vertices_end(); v_it++) {
		MyMesh::Point new_vert = ptr_object_->property(cogs, *v_it);
		new_vert[0] -= x_min;
		new_vert[0] *= scale_x;
		new_vert[1] -= y_min;
		new_vert[1] *= scale_y;
		ptr_object_->property(cogs, *v_it) = new_vert;
	}
}

bool Parameterization::EqualVector(OpenMesh::Vec3f v1, OpenMesh::Vec3f v2)
{
	if (abs(v1[0] - v2[0]) < 1e-6 && abs(v1[2] - v2[2]) < 1e-6 && abs(v1[1] - v2[1]) < 1e-6)
		return true;
	else
		return false;
}

void Parameterization::parameterization()
{
}

void Parameterization::LoadDDG(DifferentialGeometry* ddg)
{
	DDG = ddg;
}

void Parameterization::Draw(QPainter & painter, int dw, int dh)
{
	QColor  greenblue(0, 150, 200, 191);
	QColor  lightblue(137, 207, 240, 100);
	QColor  blue(0, 0, 255, 100);
	QPen pen1(blue);
	QPen pen2(greenblue);
	pen1.setWidth(1);
	pen2.setWidth(2);
	//Edges
	for (MyMesh::HalfedgeIter he_it = ptr_object_->halfedges_begin(); he_it != ptr_object_->halfedges_end(); ++he_it) {
		auto startv = ptr_object_->from_vertex_handle(*he_it);
		auto endv = ptr_object_->to_vertex_handle(*he_it);
		MyMesh::Point startvert = ptr_object_->property(cogs, startv);
		MyMesh::Point endvert = ptr_object_->property(cogs, endv);
		QPoint startp(startvert[0] + dw, startvert[1] + dh);
		QPoint endp(endvert[0] + dw, endvert[1] + dh);
		
		painter.setPen(pen1);
		painter.drawLine(startp.x(), startp.y(), endp.x(), endp.y());
		painter.setPen(pen2);
		painter.drawPoint(startp);
		painter.drawPoint(endp);
	}
}

void Parameterization::SetTexcoord()
{
	ptr_object_->request_vertex_texcoords2D();
	for (MyMesh::VertexIter v_it = ptr_object_->vertices_begin(); v_it != ptr_object_->vertices_end(); v_it++) {
		OpenMesh::Vec2f coord(ptr_object_->property(cogs, v_it.handle())[0] / w, ptr_object_->property(cogs, v_it.handle())[1] / h);
		ptr_object_->set_texcoord2D(v_it.handle(), coord);
	}
}

Pt_uniform::Pt_uniform(MyMesh* ptr_mesh_, int dw, int dh) {
	ptr_object_ = ptr_mesh_;
	ptr_object_->add_property(cogs);
	w = dw;
	h = dh;
	SetBoundary(dw, dh);
}

Pt_uniform::Pt_uniform(MyMesh* ptr_mesh_, int dw, int dh, DifferentialGeometry* ddg)
{
	ptr_object_ = ptr_mesh_;
	ptr_object_->add_property(cogs);
	w = dw;
	h = dh;
	SetBoundary(dw, dh);
	LoadDDG(ddg);
}

Pt_uniform::~Pt_uniform()
{
}

void Pt_uniform::parameterization()
{
	int row, col;
	int num = ptr_object_->n_vertices();
	row = col = num;

	//声明方程组的变量;	
	SparseMatrixType A(row, col);
	Eigen::VectorXd x;
	Eigen::VectorXd y;
	Eigen::VectorXd z;
	Eigen::VectorXd X;
	Eigen::VectorXd Y;
	Eigen::VectorXd Z;
	std::vector<T> tripletList;

	//给向量b赋值;	
	X.resize(row);
	Y.resize(row);
	Z.resize(row);
	int i = 0;
	for (auto v_it = ptr_object_->vertices_begin(); v_it != ptr_object_->vertices_end(); v_it++) {
		if (ptr_object_->is_boundary(*v_it)) {
			MyMesh::Point vert = ptr_object_->property(cogs, *v_it);
			X(i) = vert[0];
			Y(i) = vert[1];
		}
		else {
			X(i) = 0;
			Y(i) = 0;
		}
		i++;
	}

	//给稀疏矩阵A赋值;
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
	A.setFromTriplets(tripletList.begin(), tripletList.end()); //根据三元组列表生成稀疏矩阵
	A.makeCompressed();

	//求解;	
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
	
	//赋值
	i = 0;
	for (auto v_it = ptr_object_->vertices_begin(); v_it != ptr_object_->vertices_end(); v_it++) {
		MyMesh::Point new_vert(x[i], y[i], 0);
		ptr_object_->property(cogs, *v_it) = new_vert;
		i++;
	}
}

Pt_weightedleastsquares::Pt_weightedleastsquares(MyMesh * ptr_mesh_, int dw, int dh)
{
	ptr_object_ = ptr_mesh_;
	ptr_object_->add_property(cogs);
	w = dw;
	h = dh;
	SetBoundary(dw, dh);
}

Pt_weightedleastsquares::~Pt_weightedleastsquares()
{
}

void Pt_weightedleastsquares::parameterization()
{
	int q = 2;

	int row, col;
	int num = ptr_object_->n_vertices();
	row = col = num;

	//声明方程组的变量;	
	SparseMatrixType A(row, col);
	Eigen::VectorXd x;
	Eigen::VectorXd y;
	Eigen::VectorXd z;
	Eigen::VectorXd X;
	Eigen::VectorXd Y;
	Eigen::VectorXd Z;
	std::vector<T> tripletList;

	//给向量b赋值;	
	X.resize(row);
	Y.resize(row);
	Z.resize(row);
	int i = 0;
	for (auto v_it = ptr_object_->vertices_begin(); v_it != ptr_object_->vertices_end(); v_it++) {
		if (ptr_object_->is_boundary(*v_it)) {
			MyMesh::Point vert = ptr_object_->property(cogs, *v_it);
			X(i) = vert[0];
			Y(i) = vert[1];
		}
		else {
			X(i) = 0;
			Y(i) = 0;
		}
		i++;
	}

	//给稀疏矩阵A赋值;
	tripletList.reserve(num);
	i = 0;
	for (auto v_it = ptr_object_->vertices_begin(); v_it != ptr_object_->vertices_end(); v_it++) {
		if (ptr_object_->is_boundary(*v_it)) tripletList.push_back(T(i, i, 1));
		else {
			float sig_w= 0;
			for (auto vv_it = ptr_object_->vv_iter(*v_it); vv_it.is_valid(); ++vv_it) {
				float pow = (ptr_object_->point(*v_it) - ptr_object_->point(*vv_it)).length();
				for (int j = 1; j < q; j++) pow *= pow;
				float w_ij =  1 / pow;
				tripletList.push_back(T(i, vv_it.handle().idx(), -w_ij)); ;
				sig_w += w_ij;
			}
			tripletList.push_back(T(i, i, sig_w));
		}
		i++;
	}
	A.setFromTriplets(tripletList.begin(), tripletList.end()); //根据三元组列表生成稀疏矩阵
	A.makeCompressed();

	//求解;	
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

	//赋值
	i = 0;
	for (auto v_it = ptr_object_->vertices_begin(); v_it != ptr_object_->vertices_end(); v_it++) {
		MyMesh::Point new_vert(x[i], y[i], 0);
		ptr_object_->property(cogs, *v_it) = new_vert;
		i++;
	}
}

Pt_shapepreserving::Pt_shapepreserving(MyMesh * ptr_mesh_, int dw, int dh)
{
	ptr_object_ = ptr_mesh_;
	ptr_object_->add_property(cogs);
	w = dw;
	h = dh;
	SetBoundary(dw, dh);
}

Pt_shapepreserving::~Pt_shapepreserving()
{
}

void Pt_shapepreserving::parameterization()
{
	int row, col;
	int num = ptr_object_->n_vertices();
	row = col = num;

	ptr_object_->add_property(ang);
	for (auto f_it = ptr_object_->faces_begin(); f_it != ptr_object_->faces_end(); f_it++) {
		ptr_object_->property(ang, *f_it) = 0;
	}

	//声明方程组的变量;	
	SparseMatrixType A(row, col);
	Eigen::VectorXd x;
	Eigen::VectorXd y;
	Eigen::VectorXd z;
	Eigen::VectorXd X;
	Eigen::VectorXd Y;
	Eigen::VectorXd Z;
	std::vector<T> tripletList;

	//给向量b赋值;	
	X.resize(row);
	Y.resize(row);
	Z.resize(row);
	int i = 0;
	for (auto v_it = ptr_object_->vertices_begin(); v_it != ptr_object_->vertices_end(); v_it++) {
		if (ptr_object_->is_boundary(*v_it)) {
			MyMesh::Point vert = ptr_object_->property(cogs, *v_it);
			X(i) = vert[0];
			Y(i) = vert[1];
		}
		else {
			X(i) = 0;
			Y(i) = 0;
		}
		i++;
	}

	//给稀疏矩阵A赋值;
	tripletList.reserve(num);
	float **mu = new float*[num];
	for (int j = 0; j < num; j++) {
		mu[j] = new float[num];
		for (int i = 0; i < num; i++) {
			mu[j][i] = 0;
		}
	}
	i = 0;
	for (auto v_it = ptr_object_->vertices_begin(); v_it != ptr_object_->vertices_end(); v_it++) {
		if (ptr_object_->is_boundary(*v_it)) tripletList.push_back(T(i, i, 1));
		else {
			int di = 0;	
			float theta = 0;
			MyMesh::HalfedgeHandle e_next = ptr_object_->halfedge_handle(*v_it);
			do {
				MyMesh::HalfedgeHandle e_oppo = ptr_object_->opposite_halfedge_handle(e_next);
				e_next = ptr_object_->next_halfedge_handle(e_oppo);
				MyMesh::FaceHandle f = ptr_object_->face_handle(e_next);
				MyMesh::Point next = ptr_object_->point(ptr_object_->to_vertex_handle(e_next)) - ptr_object_->point(ptr_object_->from_vertex_handle(e_next));
				MyMesh::Point oppo = ptr_object_->point(ptr_object_->from_vertex_handle(e_oppo)) - ptr_object_->point(ptr_object_->to_vertex_handle(e_oppo));
				float anger = acos((next[0] * oppo[1] + next[1] * oppo[0]) / (next.length() * oppo.length()));
				ptr_object_->property(ang, f) = anger;
				theta += anger;
				di++;
			} while (e_next != ptr_object_->halfedge_handle(*v_it));

			do {
				e_next = ptr_object_->next_halfedge_handle(ptr_object_->opposite_halfedge_handle(e_next));
				MyMesh::VertexHandle v = ptr_object_->to_vertex_handle(e_next);
				MyMesh::HalfedgeHandle ve_next = e_next;
				MyMesh::HalfedgeHandle ve_last = e_next;
				MyMesh::FaceHandle f = ptr_object_->face_handle(ve_next);
				float judge_anger = 0;
				do {
					ve_last = ve_next;
					ve_next = ptr_object_->next_halfedge_handle(ptr_object_->opposite_halfedge_handle(ve_next));
					f = ptr_object_->face_handle(ve_next);
					judge_anger += ptr_object_->property(ang, f) * 2 * M_PI / theta;
					if (judge_anger > M_PI) break;
				} while (ve_next != e_next);
				MyMesh::VertexHandle v1 = ptr_object_->to_vertex_handle(ve_last);
				MyMesh::VertexHandle v2 = ptr_object_->to_vertex_handle(ve_next);

				
				if (abs(M_PI - judge_anger) < 1e-6) {
					float tri_v = (ptr_object_->point(v2) - ptr_object_->point(v_it.handle())).length();
					float tri_v1 = 0;
					float tri_v2 = (ptr_object_->point(v) - ptr_object_->point(v_it.handle())).length();
					float tri_sum = tri_v + tri_v1 + tri_v2;
					mu[v.idx()][v.idx()] = tri_v / tri_sum;
					mu[v1.idx()][v.idx()] = tri_v1 / tri_sum;
					mu[v2.idx()][v.idx()] = tri_v2 / tri_sum;
				}
				else if (abs(M_PI - judge_anger + ptr_object_->property(ang, f) * 2 * M_PI / theta) < 1e-6) {
					float tri_v = (ptr_object_->point(v1) - ptr_object_->point(v_it.handle())).length();
					float tri_v1 = (ptr_object_->point(v) - ptr_object_->point(v_it.handle())).length();
					float tri_v2 = 0;
					float tri_sum = tri_v + tri_v1 + tri_v2;
					mu[v.idx()][v.idx()] = tri_v / tri_sum;
					mu[v1.idx()][v.idx()] = tri_v1 / tri_sum;
					mu[v2.idx()][v.idx()] = tri_v2 / tri_sum;
				}
				else {
					float tri_v = (ptr_object_->point(v1) - ptr_object_->point(v_it.handle())).length() * (ptr_object_->point(v2) - ptr_object_->point(v_it.handle())).length() * sin(ptr_object_->property(ang, f) * 2 * M_PI / theta) / 2;
					float tri_v1 = (ptr_object_->point(v2) - ptr_object_->point(v_it.handle())).length() * (ptr_object_->point(v) - ptr_object_->point(v_it.handle())).length() * sin(2 * M_PI - judge_anger) / 2;
					float tri_v2 = (ptr_object_->point(v) - ptr_object_->point(v_it.handle())).length() * (ptr_object_->point(v1) - ptr_object_->point(v_it.handle())).length() * sin(judge_anger - ptr_object_->property(ang, f) * 2 * M_PI / theta) / 2;
					float tri_sum = tri_v + tri_v1 + tri_v2;
					mu[v.idx()][v.idx()] = tri_v / tri_sum;
					mu[v1.idx()][v.idx()] = tri_v1 / tri_sum;
					mu[v2.idx()][v.idx()] = tri_v2 / tri_sum;
				}

			} while (e_next != ptr_object_->halfedge_handle(*v_it));
			
			tripletList.push_back(T(i, i, -di));
			for (auto vv_it = ptr_object_->vv_iter(*v_it); vv_it.is_valid(); ++vv_it) {
				float lmd_iv = 0;
				for (int j = 0; j < num; j++) lmd_iv += mu[vv_it.handle().idx()][j];
				tripletList.push_back(T(i, vv_it.handle().idx(), lmd_iv));
			}
			for (int j = 0; j < num; j++) {
				for (int k = 0; k < num; k++) {
					mu[j][k] = 0;
				}
			}
		}
		i++;
	}
	A.setFromTriplets(tripletList.begin(), tripletList.end()); //根据三元组列表生成稀疏矩阵
	A.makeCompressed();


	//求解;	
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

	//赋值
	i = 0;
	for (auto v_it = ptr_object_->vertices_begin(); v_it != ptr_object_->vertices_end(); v_it++) {
		MyMesh::Point new_vert(x[i], y[i], 0);
		ptr_object_->property(cogs, *v_it) = new_vert;
		i++;
	}
}

Pt_ASAP::Pt_ASAP(MyMesh* ptr_mesh_, int dw, int dh)
{
	ptr_object_ = ptr_mesh_;
	ptr_object_->add_property(cogs);
	w = dw;
	h = dh;

	flatten();
}

Pt_ASAP::Pt_ASAP(MyMesh* ptr_mesh_, int dw, int dh, DifferentialGeometry* ddg)
{
	ptr_object_ = ptr_mesh_;
	ptr_object_->add_property(cogs);
	w = dw;
	h = dh;
	LoadDDG(ddg);

	flatten();
}

Pt_ASAP::~Pt_ASAP()
{
}

void Pt_ASAP::parameterization()
{
	int nT = static_cast<int>(ptr_object_->n_faces());
	int nV = static_cast<int>(ptr_object_->n_vertices());

	if (nV < 2)
		return;

	int index1 = static_cast<int>(ptr_object_->from_vertex_handle(DDG->Boundaries[0][0]).idx());
	size_t a = DDG->Boundaries[0].size() / 2;
	int index2 = static_cast<int>(ptr_object_->to_vertex_handle(DDG->Boundaries[0][a]).idx());

	std::vector<double> two_points = { 0, 0, 1, 1 };

	// construct the equations
	std::vector<T> A_Triplet;
	Eigen::VectorXd b = Eigen::VectorXd::Zero(2 * (size_t)nV + 2 * (size_t)nT);

	for (MyMesh::FaceIter f_it = ptr_object_->faces_begin(); f_it != ptr_object_->faces_end(); ++f_it) {
		int t = (*f_it).idx();
		double tmp = 0;
		int i = 0;
		for (MyMesh::FaceVertexIter fv_it = ptr_object_->fv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
			int v_index = (*fv_it).idx();
			tmp += cot[t][i] * (diff_x[t][i] * diff_x[t][i] + diff_y[t][i] * diff_y[t][i]);

			if (v_index != index1 && v_index != index2)
			{
				A_Triplet.push_back(T(2 * nV + 2 * t, 2 * v_index,
					-cot[t][i] * diff_x[t][i] + cot[t][i == 0 ? 2 : i - 1] * diff_x[t][i == 0 ? 2 : i - 1]));
				A_Triplet.push_back(T(2 * nV + 2 * t, 2 * v_index + 1,
					-cot[t][i] * diff_y[t][i] + cot[t][i == 0 ? 2 : i - 1] * diff_y[t][i == 0 ? 2 : i - 1]));
				A_Triplet.push_back(T(2 * nV + 2 * t + 1, 2 * v_index,
					-cot[t][i] * diff_y[t][i] + cot[t][i == 0 ? 2 : i - 1] * diff_y[t][i == 0 ? 2 : i - 1]));
				A_Triplet.push_back(T(2 * nV + 2 * t + 1, 2 * v_index + 1,
					cot[t][i] * diff_x[t][i] - cot[t][i == 0 ? 2 : i - 1] * diff_x[t][i == 0 ? 2 : i - 1]));

				A_Triplet.push_back(T(2 * v_index, 2 * nV + 2 * t,
					-cot[t][i] * diff_x[t][i] + cot[t][i == 0 ? 2 : i - 1] * diff_x[t][i == 0 ? 2 : i - 1]));
				A_Triplet.push_back(T(2 * v_index + 1, 2 * nV + 2 * t,
					-cot[t][i] * diff_y[t][i] + cot[t][i == 0 ? 2 : i - 1] * diff_y[t][i == 0 ? 2 : i - 1]));
				A_Triplet.push_back(T(2 * v_index, 2 * nV + 2 * t + 1,
					-cot[t][i] * diff_y[t][i] + cot[t][i == 0 ? 2 : i - 1] * diff_y[t][i == 0 ? 2 : i - 1]));
				A_Triplet.push_back(T(2 * v_index + 1, 2 * nV + 2 * t + 1,
					cot[t][i] * diff_x[t][i] - cot[t][i == 0 ? 2 : i - 1] * diff_x[t][i == 0 ? 2 : i - 1]));
			}

			else if (v_index == index2)
			{
				b(2 * (size_t)nV + 2 * (size_t)t) -= (-cot[t][i] * diff_x[t][i] + cot[t][i == 0 ? 2 : i - 1] * diff_x[t][i == 0 ? 2 : i - 1]) * two_points[2];
				b(2 * (size_t)nV + 2 * (size_t)t) -= (-cot[t][i] * diff_y[t][i] + cot[t][i == 0 ? 2 : i - 1] * diff_y[t][i == 0 ? 2 : i - 1]) * two_points[3];
				b(2 * (size_t)nV + 2 * (size_t)t + 1) -= (-cot[t][i] * diff_y[t][i] + cot[t][i == 0 ? 2 : i - 1] * diff_y[t][i == 0 ? 2 : i - 1]) * two_points[2];
				b(2 * (size_t)nV + 2 * (size_t)t + 1) -= (cot[t][i] * diff_x[t][i] - cot[t][i == 0 ? 2 : i - 1] * diff_x[t][i == 0 ? 2 : i - 1]) * two_points[3];
			}
			i++;
		}
		A_Triplet.push_back(T(2 * nV + 2 * t, 2 * nV + 2 * t, tmp));
		A_Triplet.push_back(T(2 * nV + 2 * t + 1, 2 * nV + 2 * t + 1, tmp));
	}

	for (MyMesh::VertexIter v_it = ptr_object_->vertices_begin(); v_it != ptr_object_->vertices_end(); ++v_it) {
		int i = (*v_it).idx();
		double tmp = 0;

		if (i == index1 || i == index2)
		{
			A_Triplet.push_back(T(2 * i, 2 * i, 1));
			A_Triplet.push_back(T(2 * i + 1, 2 * i + 1, 1));
			b(2 * (size_t)i) = i == index1 ? two_points[0] : two_points[2];
			b(2 * (size_t)i + 1) = i == index1 ? two_points[1] : two_points[3];
			continue;
		}
		
		for (MyMesh::VertexOHalfedgeIter voh_it = ptr_object_->voh_begin(*v_it); voh_it != ptr_object_->voh_end(*v_it); ++voh_it) {
			double coe_ij = 0;
			// Left Triangle
			MyMesh::FaceHandle triangle = ptr_object_->face_handle(*voh_it);
			size_t t = 0;
			int index;
			if (!ptr_object_->is_boundary(*voh_it))
			{
				t = triangle.idx();
				index = 0;
				for (MyMesh::FaceVertexIter fv_it = ptr_object_->fv_iter(triangle); fv_it.is_valid(); ++fv_it) {
					if ((*fv_it).idx() == i) break;
					index++;
				}
				coe_ij += cot[t][index];
			}
			// Right Triangle
			triangle = ptr_object_->face_handle(ptr_object_->opposite_halfedge_handle(*voh_it));
			if (!ptr_object_->is_boundary(ptr_object_->opposite_halfedge_handle(*voh_it)))
			{
				t = triangle.idx();
				index = 0;
				for (MyMesh::FaceVertexIter fv_it = ptr_object_->fv_iter(triangle); fv_it.is_valid(); ++fv_it) {
					if ((*fv_it).idx() == i) break;
					index++;
				}
				coe_ij += cot[t][index == 0 ? 2 : index - 1];
			}
			// Set
			tmp += coe_ij;
			int v_index = ptr_object_->to_vertex_handle(*voh_it).idx();
			if (v_index != index1 && v_index != index2)
			{
				A_Triplet.push_back(T(i * 2, 2 * v_index, -coe_ij));
				A_Triplet.push_back(T(i * 2 + 1, 2 * v_index + 1, -coe_ij));
			}
			else if (v_index == index2)
			{
				b((size_t)i * 2) -= (-coe_ij) * two_points[2];
				b((size_t)i * 2 + 1) -= (-coe_ij) * two_points[3];
			}
		}
		A_Triplet.push_back(T(i * 2, i * 2, tmp));
		A_Triplet.push_back(T(i * 2 + 1, i * 2 + 1, tmp));
	}
	//cout << b << endl;

	// solve the equations
	SparseMatrixType A(2 * (size_t)nV + 2 * (size_t)nT, 2 * (size_t)nV + 2 * (size_t)nT);
	Solve_LLT solver;
	A.setFromTriplets(A_Triplet.begin(), A_Triplet.end());
	//cout << A << endl;
	solver.compute(A);

	Eigen::VectorXd V = solver.solve(b);
	//cout << V << endl;

	// Update the solutions
	//赋值
	float x_min = INT_MAX;
	float x_max = INT_MIN;
	float y_min = INT_MAX;
	float y_max = INT_MIN;
	for (auto v_it = ptr_object_->vertices_begin(); v_it != ptr_object_->vertices_end(); v_it++) {
		int i = (*v_it).idx();
		MyMesh::Point new_vert;
		if (i != index1 && i != index2) {
			new_vert = MyMesh::Point(V(2 * i), V(2 * i + 1), 0);
		}
		else if (i == index1) {
			new_vert = MyMesh::Point(two_points[0], two_points[1], 0);
		}
		else {
			new_vert = MyMesh::Point(two_points[2], two_points[3], 0);
		}

		if (new_vert[0] > x_max) x_max = new_vert[0];
		else if(new_vert[0] < x_min) x_min = new_vert[0];
		if (new_vert[1] > y_max) y_max = new_vert[1];
		else if (new_vert[1] < y_min) y_min = new_vert[1];

		ptr_object_->property(cogs, *v_it) = new_vert;
	}

	FitIn(x_min, x_max, y_min, y_max);
}

void Pt_ASAP::flatten()
{
	size_t nT = ptr_object_->n_faces();
	std::vector<std::vector<Eigen::Vector2f>> x_plane(nT);
	cot.resize(nT);

	for (MyMesh::FaceIter f_it = ptr_object_->faces_begin(); f_it != ptr_object_->faces_end(); ++f_it) {
		int i = (*f_it).idx();
		std::vector<MyMesh::VertexHandle> x_tmp;
		for (MyMesh::FaceVertexIter fv_it = ptr_object_->fv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
			x_tmp.push_back(*fv_it);
		}
		
		x_plane[i].push_back(Eigen::Vector2f(0, 0));
		x_plane[i].push_back(Eigen::Vector2f((ptr_object_->point(x_tmp[1]) - ptr_object_->point(x_tmp[0])).norm(), 0));
		float x = ((ptr_object_->point(x_tmp[2]) - ptr_object_->point(x_tmp[0])) | (ptr_object_->point(x_tmp[1]) - ptr_object_->point(x_tmp[0]))) / (ptr_object_->point(x_tmp[1]) - ptr_object_->point(x_tmp[0])).norm();
		x_plane[i].push_back(Eigen::Vector2f(x, pow(pow((ptr_object_->point(x_tmp[2]) - ptr_object_->point(x_tmp[0])).norm(), 2) - pow(x, 2), 0.5)));
		
		if (Parameterization::EqualVector(ptr_object_->point(x_tmp[0]), ptr_object_->point(x_tmp[1])) || Parameterization::EqualVector(ptr_object_->point(x_tmp[1]), ptr_object_->point(x_tmp[2])) || Parameterization::EqualVector(ptr_object_->point(x_tmp[2]), ptr_object_->point(x_tmp[0])))
		{
			cot[i] = { 1, 1, 1 };
			diff_x.push_back({ 1, 1, 1 });
			diff_y.push_back({ 1, 1, 1 });
			continue;
		}

		diff_x.push_back({ (double)x_plane[i][0][0] - x_plane[i][1][0], (double)x_plane[i][1][0] - x_plane[i][2][0], (double)x_plane[i][2][0] - x_plane[i][0][0] });
		diff_y.push_back({ (double)x_plane[i][0][1] - x_plane[i][1][1], (double)x_plane[i][1][1] - x_plane[i][2][1], (double)x_plane[i][2][1] - x_plane[i][0][1] });

		for (int j = 0; j < 3; j++) {
			Eigen::Vector2f a = (x_plane[i][j] - x_plane[i][j == 0 ? 2 : j - 1]);
			Eigen::Vector2f b = (x_plane[i][j == 2 ? 0 : j + 1] - x_plane[i][j == 0 ? 2 : j - 1]);
			
			float cos_theta = a.dot(b) / (a.norm() * b.norm());
			float abs_cot = pow(cos_theta * cos_theta / (1 - cos_theta * cos_theta), 0.5);

			cot[i].push_back(cos_theta > 0 ? abs_cot : -abs_cot);
		}
	}
}

Pt_ARAP::Pt_ARAP(MyMesh* ptr_mesh_, int dw, int dh, DifferentialGeometry* ddg)
{
	ptr_object_ = ptr_mesh_;
	ptr_object_->add_property(cogs);
	w = dw;
	h = dh;
	int nV = ptr_object_->n_vertices();
	int nF = ptr_object_->n_faces();
	setTimes(10);

	Pt_ASAP init_pt(ptr_mesh_, dw, dh, ddg);
	init_pt.parameterization();
	for (auto v_it = ptr_object_->vertices_begin(); v_it != ptr_object_->vertices_end(); v_it++) {
		int i = (*v_it).idx();
		MyMesh::Point new_vert = init_pt.ptr_object_->property(init_pt.cogs, *v_it);
		ptr_object_->property(cogs, *v_it) = new_vert;
	}

	flatten();
	LoadDDG(ddg);
	// set two fixed vertices
	if (nV < 2) {
		printf("ERROR::ARAP::Init:\n"
			"\t""need more vertices\n");
		return;
	}
	size_t a = DDG->Boundaries[0].size() / 2;
	auto v1 = ptr_object_->from_vertex_handle(DDG->Boundaries[0][0]);
	auto v2 = ptr_object_->to_vertex_handle(DDG->Boundaries[0][a]);
	fixed_vertices_.push_back(v1.idx());
	fixed_coords_.push_back(Eigen::Vector2f(0, 0));
	L.resize(nF);
}

Pt_ARAP::~Pt_ARAP()
{
}

void Pt_ARAP::parameterization()
{
	size_t nV = ptr_object_->n_vertices();
	SparseMatrixType A(nV, nV);
	GlobalMatrixA(A);
	//cout << A << endl;

	Solve_LLT solver;
	solver.compute(A);

	if (times < 0 || times > 100)
		times = 10;
	//std::cout << "time = " << times << std::endl;

	for (int count = 0; count < times; count++)
	{
		LocalSetL();
		GlobalSolveU(solver);
	}
}

void Pt_ARAP::flatten()
{
	size_t nT = ptr_object_->n_faces();
	x_plane.resize(nT);
	cot.resize(nT);
	for (MyMesh::FaceIter f_it = ptr_object_->faces_begin(); f_it != ptr_object_->faces_end(); ++f_it) {
		int i = (*f_it).idx();
		std::vector<MyMesh::VertexHandle> x_tmp;
		for (MyMesh::FaceVertexIter fv_it = ptr_object_->fv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
			x_tmp.push_back(*fv_it);
		}

		if (Parameterization::EqualVector(ptr_object_->point(x_tmp[0]), ptr_object_->point(x_tmp[1])) || Parameterization::EqualVector(ptr_object_->point(x_tmp[1]), ptr_object_->point(x_tmp[2])) || Parameterization::EqualVector(ptr_object_->point(x_tmp[2]), ptr_object_->point(x_tmp[0])))
		{
			cot[i] = { 1, 1, 1 };
			x_plane[i] = { Eigen::Vector2f(1, 1), Eigen::Vector2f(1, 1), Eigen::Vector2f(1, 1) };
			continue;
		}
		
		x_plane[i].push_back(Eigen::Vector2f(0, 0));
		x_plane[i].push_back(Eigen::Vector2f((ptr_object_->point(x_tmp[1]) - ptr_object_->point(x_tmp[0])).norm(), 0));
		float x = ((ptr_object_->point(x_tmp[2]) - ptr_object_->point(x_tmp[0])) | (ptr_object_->point(x_tmp[1]) - ptr_object_->point(x_tmp[0]))) / (ptr_object_->point(x_tmp[1]) - ptr_object_->point(x_tmp[0])).norm();
		x_plane[i].push_back(Eigen::Vector2f(x, pow(pow((ptr_object_->point(x_tmp[2]) - ptr_object_->point(x_tmp[0])).norm(), 2) - pow(x, 2), 0.5)));

		for (int j = 0; j < 3; j++) {
			Eigen::Vector2f a = (x_plane[i][j] - x_plane[i][j == 0 ? 2 : j - 1]);
			Eigen::Vector2f b = (x_plane[i][j == 2 ? 0 : j + 1] - x_plane[i][j == 0 ? 2 : j - 1]);

			float cos_theta = a.dot(b) / (a.norm() * b.norm());
			float abs_cot = pow(cos_theta * cos_theta / (1 - cos_theta * cos_theta), 0.5);

			cot[i].push_back(cos_theta > 0 ? abs_cot : -abs_cot);
		}
	}
}

void Pt_ARAP::setTimes(int time)
{
	times = time;
}

void Pt_ARAP::LocalSetL()
{
	for (MyMesh::FaceIter f_it = ptr_object_->faces_begin(); f_it != ptr_object_->faces_end(); ++f_it)
	{
		int i = (*f_it).idx();
		std::vector<MyMesh::VertexHandle> vertices;
		for (MyMesh::FaceVertexIter fv_it = ptr_object_->fv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
			vertices.push_back(*fv_it);
		}
		Eigen::Matrix2d X, U;
		X << (x_plane[i][0] - x_plane[i][1])[0], (x_plane[i][1] - x_plane[i][2])[0],
			(x_plane[i][0] - x_plane[i][1])[1], (x_plane[i][1] - x_plane[i][2])[1];
		U << (ptr_object_->property(cogs, vertices[0]) - ptr_object_->property(cogs, vertices[1]))[0], (ptr_object_->property(cogs, vertices[1]) - ptr_object_->property(cogs, vertices[2]))[0],
			(ptr_object_->property(cogs, vertices[0]) - ptr_object_->property(cogs, vertices[1]))[1], (ptr_object_->property(cogs, vertices[1]) - ptr_object_->property(cogs, vertices[2]))[1];
		Eigen::Matrix2d J = U * X.inverse();

		Eigen::JacobiSVD<Eigen::Matrix2d> svd(J, Eigen::ComputeFullU | Eigen::ComputeFullV);
		if (J.determinant() > 0)
			L[i] = svd.matrixU() * svd.matrixV().transpose();
		else
		{
			Eigen::Matrix2d D;
			D(0, 0) = 1; D(0, 1) = 0; D(1, 0) = 0; D(1, 1) = -1;
			L[i] = svd.matrixU() * D * svd.matrixV().transpose();
		}
		
		//cout << "L:" << triangle->L << endl;
	}
}

void Pt_ARAP::GlobalSolveU(Solve_LLT& solver)
{
	size_t nV = ptr_object_->n_vertices();
	Eigen::VectorXd b_x = Eigen::VectorXd::Zero(nV);
	Eigen::VectorXd b_y = Eigen::VectorXd::Zero(nV);

	for (auto v_it = ptr_object_->vertices_begin(); v_it != ptr_object_->vertices_end(); v_it++)
	{
		int i = (*v_it).idx();
		if (i == fixed_vertices_[0])
		{
			b_x(i) += fixed_coords_[0][0];
			b_y(i) += fixed_coords_[0][1];
			continue;
		}
		/*if (i == fixed_vertices_[1])
		{
			b_x(i) += fixed_coords_[1][0];
			b_y(i) += fixed_coords_[1][1];
			continue;
		}*/
		Eigen::Vector2d sum = Eigen::Vector2d::Zero();
		for (MyMesh::VertexOHalfedgeIter voh_it = ptr_object_->voh_begin(*v_it); voh_it != ptr_object_->voh_end(*v_it); ++voh_it)
		{
			// Left Triangle
			MyMesh::FaceHandle triangle = ptr_object_->face_handle(*voh_it);
			size_t t = 0;
			int index;
			double coe_ij = 0;
			if (!ptr_object_->is_boundary(*voh_it))
			{
				t = triangle.idx();
				index = 0;
				for (MyMesh::FaceVertexIter fv_it = ptr_object_->fv_iter(triangle); fv_it.is_valid(); ++fv_it) {
					if ((*fv_it).idx() == i) break;
					index++;
				}
				auto diffx = x_plane[t][index] - x_plane[t][index == 2 ? 0 : index + 1];
				coe_ij += cot[t][index];
				sum += cot[t][index] * L[t] * Eigen::Vector2d(diffx[0], diffx[1]);
			}
			// Right Triangle
			triangle = ptr_object_->face_handle(ptr_object_->opposite_halfedge_handle(*voh_it));
			if (!ptr_object_->is_boundary(ptr_object_->opposite_halfedge_handle(*voh_it)))
			{
				t = triangle.idx();
				index = 0;
				for (MyMesh::FaceVertexIter fv_it = ptr_object_->fv_iter(triangle); fv_it.is_valid(); ++fv_it) {
					if ((*fv_it).idx() == i) break;
					index++;
				}
				auto diffx = x_plane[t][index] - x_plane[t][index == 0 ? 2 : index - 1];
				coe_ij += cot[t][index == 0 ? 2 : index - 1];
				sum += cot[t][index == 0 ? 2 : index - 1] * L[t] * Eigen::Vector2d(diffx[0], diffx[1]);
			}
			// Set
			int j = ptr_object_->to_vertex_handle(*voh_it).idx();
			if (j == fixed_vertices_[0])
			{
				b_x(i) -= (-coe_ij) * fixed_coords_[0][0];
				b_y(i) -= (-coe_ij) * fixed_coords_[0][1];
			}
			/*else if (j == fixed_vertices_[1])
			{
				b_x(i) -= (-coe_ij) * fixed_coords_[1][0];
				b_y(i) -= (-coe_ij) * fixed_coords_[1][1];
			}*/

		}

		b_x(i) += sum(0);
		b_y(i) += sum(1);
	}
	//cout << b_x << endl << b_y << endl;
	Eigen::VectorXd u_x = solver.solve(b_x);
	Eigen::VectorXd u_y = solver.solve(b_y);
	//cout << u_x << endl << u_y << endl;
	float x_min = INT_MAX;
	float x_max = INT_MIN;
	float y_min = INT_MAX;
	float y_max = INT_MIN;
	for (auto v_it = ptr_object_->vertices_begin(); v_it != ptr_object_->vertices_end(); v_it++) {
		int i = (*v_it).idx();
		MyMesh::Point new_vert = MyMesh::Point(u_x(i), u_y(i), 0);

		if (new_vert[0] > x_max) x_max = new_vert[0];
		else if (new_vert[0] < x_min) x_min = new_vert[0];
		if (new_vert[1] > y_max) y_max = new_vert[1];
		else if (new_vert[1] < y_min) y_min = new_vert[1];

		ptr_object_->property(cogs, *v_it) = new_vert;
	}

	FitIn(x_min, x_max, y_min, y_max);
}

void Pt_ARAP::GlobalMatrixA(SparseMatrixType& A)
{
	std::vector<T> A_Triplet;
	int nV = ptr_object_->n_vertices();

	for (MyMesh::VertexIter v_it = ptr_object_->vertices_begin(); v_it != ptr_object_->vertices_end(); ++v_it) {
		int i = (*v_it).idx();
		if (i == fixed_vertices_[0])
		{
			A_Triplet.push_back(T(i, i, 1));
			continue;
		}
		double sum = 0;

		for (MyMesh::VertexOHalfedgeIter voh_it = ptr_object_->voh_begin(*v_it); voh_it != ptr_object_->voh_end(*v_it); ++voh_it) {
			double coe_ij = 0;
			// Left Triangle
			MyMesh::FaceHandle triangle = ptr_object_->face_handle(*voh_it);
			size_t t = 0;
			int index;
			if (!ptr_object_->is_boundary(*voh_it))
			{
				t = triangle.idx();
				index = 0;
				for (MyMesh::FaceVertexIter fv_it = ptr_object_->fv_iter(triangle); fv_it.is_valid(); ++fv_it) {
					if ((*fv_it).idx() == i) break;
					index++;
				}
				coe_ij += cot[t][index];
			}
			// Right Triangle
			triangle = ptr_object_->face_handle(ptr_object_->opposite_halfedge_handle(*voh_it));
			if (!ptr_object_->is_boundary(ptr_object_->opposite_halfedge_handle(*voh_it)))
			{
				t = triangle.idx();
				index = 0;
				for (MyMesh::FaceVertexIter fv_it = ptr_object_->fv_iter(triangle); fv_it.is_valid(); ++fv_it) {
					if ((*fv_it).idx() == i) break;
					index++;
				}
				coe_ij += cot[t][index == 0 ? 2 : index - 1];
			}
			// Set
			sum += coe_ij;
			int v_index = ptr_object_->to_vertex_handle(*voh_it).idx();
			if (v_index != fixed_vertices_[0])
				A_Triplet.push_back(T(i, v_index, -coe_ij));
		}
		A_Triplet.push_back(T(i, i, sum));
	}
	A.setFromTriplets(A_Triplet.begin(), A_Triplet.end());
}

Pt_SLIM::Pt_SLIM(MyMesh* ptr_mesh_, int dw, int dh, DifferentialGeometry* ddg)
{
	ptr_object_ = ptr_mesh_;
	ptr_object_->add_property(cogs);
	w = dw;
	h = dh;
	int nV = ptr_object_->n_vertices();
	int nF = ptr_object_->n_faces();
	setTimes(100);

	Pt_ARAP init_pt(ptr_mesh_, dw, dh, ddg);
	init_pt.parameterization();
	for (auto v_it = ptr_object_->vertices_begin(); v_it != ptr_object_->vertices_end(); v_it++) {
		int i = (*v_it).idx();
		MyMesh::Point new_vert = init_pt.ptr_object_->property(init_pt.cogs, *v_it);
		ptr_object_->property(cogs, *v_it) = new_vert;
	}

	flatten();
	LoadDDG(ddg);
	// set two fixed vertices
	if (nV < 2) {
		printf("ERROR::ARAP::Init:\n"
			"\t""need more vertices\n");
		return;
	}
	size_t a = DDG->Boundaries[0].size() / 2;
	auto v1 = ptr_object_->from_vertex_handle(DDG->Boundaries[0][0]);
	auto v2 = ptr_object_->to_vertex_handle(DDG->Boundaries[0][a]);
	fixed_vertices_.push_back(v1.idx());
	fixed_coords_.push_back(Eigen::Vector2f(0, 0));
	L.resize(nF);
	W.resize(nF);
	p.resize(nV);
}

Pt_SLIM::~Pt_SLIM()
{
}

void Pt_SLIM::parameterization()
{
	size_t nV = ptr_object_->n_vertices();
	SparseMatrixType A(nV * 2, nV * 2);
	Solve_LLT solver;

	if (times < 0 || times > 5000)
		times = 2500;
	//std::cout << "time = " << times << std::endl;

	for (int count = 0; count < times; count++)
	{
		//std::cout << "Iter: " << count << std::endl;
		LocalSetL();
		UpdateGlobalMatrixA(A);
		solver.compute(A);
		GlobalSolveU(solver);
		LineSearch(A);
	}
}

void Pt_SLIM::flatten()
{
	size_t nT = ptr_object_->n_faces();
	x_plane.resize(nT);
	cot.resize(nT);
	for (MyMesh::FaceIter f_it = ptr_object_->faces_begin(); f_it != ptr_object_->faces_end(); ++f_it) {
		int i = (*f_it).idx();
		std::vector<MyMesh::VertexHandle> x_tmp;
		for (MyMesh::FaceVertexIter fv_it = ptr_object_->fv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
			x_tmp.push_back(*fv_it);
		}

		if (Parameterization::EqualVector(ptr_object_->point(x_tmp[0]), ptr_object_->point(x_tmp[1])) || Parameterization::EqualVector(ptr_object_->point(x_tmp[1]), ptr_object_->point(x_tmp[2])) || Parameterization::EqualVector(ptr_object_->point(x_tmp[2]), ptr_object_->point(x_tmp[0])))
		{
			cot[i] = { 1, 1, 1 };
			x_plane[i] = { Eigen::Vector2f(1, 1), Eigen::Vector2f(1, 1), Eigen::Vector2f(1, 1) };
			continue;
		}

		x_plane[i].push_back(Eigen::Vector2f(0, 0));
		x_plane[i].push_back(Eigen::Vector2f((ptr_object_->point(x_tmp[1]) - ptr_object_->point(x_tmp[0])).norm(), 0));
		float x = ((ptr_object_->point(x_tmp[2]) - ptr_object_->point(x_tmp[0])) | (ptr_object_->point(x_tmp[1]) - ptr_object_->point(x_tmp[0]))) / (ptr_object_->point(x_tmp[1]) - ptr_object_->point(x_tmp[0])).norm();
		x_plane[i].push_back(Eigen::Vector2f(x, pow(pow((ptr_object_->point(x_tmp[2]) - ptr_object_->point(x_tmp[0])).norm(), 2) - pow(x, 2), 0.5)));

		for (int j = 0; j < 3; j++) {
			Eigen::Vector2f a = (x_plane[i][j] - x_plane[i][j == 0 ? 2 : j - 1]);
			Eigen::Vector2f b = (x_plane[i][j == 2 ? 0 : j + 1] - x_plane[i][j == 0 ? 2 : j - 1]);

			float cos_theta = a.dot(b) / (a.norm() * b.norm());
			float abs_cot = pow(cos_theta * cos_theta / (1 - cos_theta * cos_theta), 0.5);

			cot[i].push_back(cos_theta > 0 ? abs_cot : -abs_cot);
		}
	}
}

void Pt_SLIM::setTimes(int time)
{
	times = time;
}

void Pt_SLIM::LocalSetL()
{
	for (MyMesh::FaceIter f_it = ptr_object_->faces_begin(); f_it != ptr_object_->faces_end(); ++f_it)
	{
		int i = (*f_it).idx();
		std::vector<MyMesh::VertexHandle> vertices;
		for (MyMesh::FaceVertexIter fv_it = ptr_object_->fv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
			vertices.push_back(*fv_it);
		}
		Eigen::Matrix2d X, U;
		X << (x_plane[i][0] - x_plane[i][1])[0], (x_plane[i][1] - x_plane[i][2])[0],
			(x_plane[i][0] - x_plane[i][1])[1], (x_plane[i][1] - x_plane[i][2])[1];
		U << (ptr_object_->property(cogs, vertices[0]) - ptr_object_->property(cogs, vertices[1]))[0], (ptr_object_->property(cogs, vertices[1]) - ptr_object_->property(cogs, vertices[2]))[0],
			(ptr_object_->property(cogs, vertices[0]) - ptr_object_->property(cogs, vertices[1]))[1], (ptr_object_->property(cogs, vertices[1]) - ptr_object_->property(cogs, vertices[2]))[1];
		Eigen::Matrix2d J = U * X.inverse();

		Eigen::JacobiSVD<Eigen::Matrix2d> svd(J, Eigen::ComputeFullU | Eigen::ComputeFullV);
		Eigen::Matrix2d S = svd.matrixU().transpose() * J * svd.matrixV();
		Eigen::Matrix2d Sw;
		Sw <<	sqrt((S(0, 0) - pow(S(0, 0), -3)) / (S(0, 0) - 1)), 0, 
				0, sqrt((S(1, 1) - pow(S(1, 1), -3)) / (S(1, 1) - 1));

		L[i] = svd.matrixU() * svd.matrixV().transpose();
		
		W[i] = svd.matrixU().transpose() * Sw * svd.matrixU();
		if (S(0, 0) * S(1, 1) < 0) std::cout << "Flip!" << std::endl;
		//W[i] = Eigen::Matrix2d().Identity();
		//cout << "L:" << triangle->L << endl;
	}
}

void Pt_SLIM::GlobalSolveU(Solve_LLT& solver)
{
	size_t nV = ptr_object_->n_vertices();
	Eigen::VectorXd b = Eigen::VectorXd::Zero(nV * 2);
	
	for (auto v_it = ptr_object_->vertices_begin(); v_it != ptr_object_->vertices_end(); v_it++)
	{
		int i = (*v_it).idx();
		if (i == fixed_vertices_[0])
		{
			b(i) += fixed_coords_[0][0];
			b(i+nV) += fixed_coords_[0][1];
			continue;
		}
		/*if (i == fixed_vertices_[1])
		{
			b_x(i) += fixed_coords_[1][0];
			b_y(i) += fixed_coords_[1][1];
			continue;
		}*/
		Eigen::Vector2d sum = Eigen::Vector2d::Zero();
		for (MyMesh::VertexOHalfedgeIter voh_it = ptr_object_->voh_begin(*v_it); voh_it != ptr_object_->voh_end(*v_it); ++voh_it)
		{
			// Left Triangle
			MyMesh::FaceHandle triangle = ptr_object_->face_handle(*voh_it);
			size_t t = 0;
			int index;
			double coe_ijxx = 0;
			double coe_ijxy = 0;
			double coe_ijyx = 0;
			double coe_ijyy = 0;
			if (!ptr_object_->is_boundary(*voh_it))
			{
				t = triangle.idx();
				index = 0;
				for (MyMesh::FaceVertexIter fv_it = ptr_object_->fv_iter(triangle); fv_it.is_valid(); ++fv_it) {
					if ((*fv_it).idx() == i) break;
					index++;
				}
				auto diffx = x_plane[t][index] - x_plane[t][index == 2 ? 0 : index + 1];
				coe_ijxx += cot[t][index] * W[t](0, 0);
				coe_ijxy += cot[t][index] * W[t](0, 1);
				coe_ijyx += cot[t][index] * W[t](1, 0);
				coe_ijyy += cot[t][index] * W[t](1, 1);
				sum += cot[t][index] * W[t] * L[t] * Eigen::Vector2d(diffx[0], diffx[1]);
			}
			// Right Triangle
			triangle = ptr_object_->face_handle(ptr_object_->opposite_halfedge_handle(*voh_it));
			if (!ptr_object_->is_boundary(ptr_object_->opposite_halfedge_handle(*voh_it)))
			{
				t = triangle.idx();
				index = 0;
				for (MyMesh::FaceVertexIter fv_it = ptr_object_->fv_iter(triangle); fv_it.is_valid(); ++fv_it) {
					if ((*fv_it).idx() == i) break;
					index++;
				}
				auto diffx = x_plane[t][index] - x_plane[t][index == 0 ? 2 : index - 1];
				coe_ijxx += cot[t][index == 0 ? 2 : index - 1] * W[t](0, 0);
				coe_ijxy += cot[t][index == 0 ? 2 : index - 1] * W[t](0, 1);
				coe_ijyx += cot[t][index == 0 ? 2 : index - 1] * W[t](1, 0);
				coe_ijyy += cot[t][index == 0 ? 2 : index - 1] * W[t](1, 1);
				sum += cot[t][index == 0 ? 2 : index - 1] * W[t] * L[t] * Eigen::Vector2d(diffx[0], diffx[1]);
			}
			// Set
			int j = ptr_object_->to_vertex_handle(*voh_it).idx();
			if (j == fixed_vertices_[0])
			{
				b(i) -= (-coe_ijxx) * fixed_coords_[0][0] + (-coe_ijxy) * fixed_coords_[0][1];
				b(i+nV) -= (-coe_ijyx) * fixed_coords_[0][0] + (-coe_ijyy) * fixed_coords_[0][1];
			}
			/*else if (j == fixed_vertices_[1])
			{
				b_x(i) -= (-coe_ij) * fixed_coords_[1][0];
				b_y(i) -= (-coe_ij) * fixed_coords_[1][1];
			}*/

		}

		b(i) += sum(0);
		b(i+nV) += sum(1);
	}
	//cout << b_x << endl << b_y << endl;
	Eigen::VectorXd u = solver.solve(b);
	//cout << u_x << endl << u_y << endl;
	for (auto v_it = ptr_object_->vertices_begin(); v_it != ptr_object_->vertices_end(); v_it++) {
		int i = (*v_it).idx();
		MyMesh::Point new_vert = MyMesh::Point(u(i), u(i+nV), 0);
		p[(*v_it).idx()] = new_vert - ptr_object_->property(cogs, *v_it);
	}

	//FitIn(x_min, x_max, y_min, y_max);
}

void Pt_SLIM::UpdateGlobalMatrixA(SparseMatrixType& A)
{
	std::vector<T> A_Triplet;
	int nV = ptr_object_->n_vertices();

	for (MyMesh::VertexIter v_it = ptr_object_->vertices_begin(); v_it != ptr_object_->vertices_end(); ++v_it) {
		int i = (*v_it).idx();
		if (i == fixed_vertices_[0])
		{
			A_Triplet.push_back(T(i, i, 1));
			A_Triplet.push_back(T(i + nV, i + nV, 1));
			continue;
		}
		double sumxx = 0;
		double sumxy = 0;
		double sumyx = 0;
		double sumyy = 0;

		for (MyMesh::VertexOHalfedgeIter voh_it = ptr_object_->voh_begin(*v_it); voh_it != ptr_object_->voh_end(*v_it); ++voh_it) {
			double coe_ijxx = 0;
			double coe_ijxy = 0;
			double coe_ijyx = 0;
			double coe_ijyy = 0;
			// Left Triangle
			MyMesh::FaceHandle triangle = ptr_object_->face_handle(*voh_it);
			size_t t = 0;
			int index;
			if (!ptr_object_->is_boundary(*voh_it))
			{
				t = triangle.idx();
				index = 0;
				for (MyMesh::FaceVertexIter fv_it = ptr_object_->fv_iter(triangle); fv_it.is_valid(); ++fv_it) {
					if ((*fv_it).idx() == i) break;
					index++;
				}
				coe_ijxx += cot[t][index] * W[t](0, 0);
				coe_ijxy += cot[t][index] * W[t](0, 1);
				coe_ijyx += cot[t][index] * W[t](1, 0);
				coe_ijyy += cot[t][index] * W[t](1, 1);
			}
			// Right Triangle
			triangle = ptr_object_->face_handle(ptr_object_->opposite_halfedge_handle(*voh_it));
			if (!ptr_object_->is_boundary(ptr_object_->opposite_halfedge_handle(*voh_it)))
			{
				t = triangle.idx();
				index = 0;
				for (MyMesh::FaceVertexIter fv_it = ptr_object_->fv_iter(triangle); fv_it.is_valid(); ++fv_it) {
					if ((*fv_it).idx() == i) break;
					index++;
				}
				coe_ijxx += cot[t][index == 0 ? 2 : index - 1] * W[t](0, 0);
				coe_ijxy += cot[t][index == 0 ? 2 : index - 1] * W[t](0, 1);
				coe_ijyx += cot[t][index == 0 ? 2 : index - 1] * W[t](1, 0);
				coe_ijyy += cot[t][index == 0 ? 2 : index - 1] * W[t](1, 1);
			}
			// Set
			sumxx += coe_ijxx;
			sumxy += coe_ijxy;
			sumyx += coe_ijyx;
			sumyy += coe_ijyy;
			int v_index = ptr_object_->to_vertex_handle(*voh_it).idx();
			if (v_index != fixed_vertices_[0]) {
				A_Triplet.push_back(T(i, v_index, -coe_ijxx));
				A_Triplet.push_back(T(i, v_index+nV, -coe_ijxy));
				A_Triplet.push_back(T(i + nV, v_index, -coe_ijyx));
				A_Triplet.push_back(T(i + nV, v_index + nV, -coe_ijyy));

			}
		}
		A_Triplet.push_back(T(i, i, sumxx));
		A_Triplet.push_back(T(i, i + nV, sumxy));
		A_Triplet.push_back(T(i + nV, i, sumyx));
		A_Triplet.push_back(T(i + nV, i + nV, sumyy));
	}
	A.setFromTriplets(A_Triplet.begin(), A_Triplet.end());
}

void Pt_SLIM::LineSearch(SparseMatrixType& A)
{
	float alpha_max = INT_MAX;
	for (MyMesh::FaceIter f_it = ptr_object_->faces_begin(); f_it != ptr_object_->faces_end(); ++f_it)
	{
		int i = (*f_it).idx();
		std::vector<MyMesh::VertexHandle> vertices;
		for (MyMesh::FaceVertexIter fv_it = ptr_object_->fv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
			vertices.push_back(*fv_it);
		}
		auto e1 = (ptr_object_->property(cogs, vertices[1]) - ptr_object_->property(cogs, vertices[0]));
		auto e2 = (ptr_object_->property(cogs, vertices[2]) - ptr_object_->property(cogs, vertices[0]));
		auto p1 = p[vertices[1].idx()] - p[vertices[0].idx()];
		auto p2 = p[vertices[2].idx()] - p[vertices[0].idx()];

		float a = p1[0] * p2[1] - p1[1] * p2[0];
		float b = p1[0] * e2[1] + e1[0] * p2[1] - p2[0] * e1[1] - e2[0] * p1[1];
		float c = e1[0] * e2[1] - e1[1] * e2[0];
		float delta = b * b - 4 * a * c;
		if (delta < 0) continue;
		else {
			float r1 = (-b - sqrt(delta)) / (2 * a);
			float r2 = (-b + sqrt(delta)) / (2 * a);
			if (r1 > 0) {
				if (alpha_max > r1) {
					alpha_max = r1;
				}
				continue;
			}
			else if (r2 > 0) {
				if (alpha_max > r2) {
					alpha_max = r2;
				}
				continue;
			}
		}
	}
	//std::cout << alpha_max << std::endl;
	float alpha = std::min(alpha_max * 0.8, 1.0);

	size_t nV = ptr_object_->n_vertices();
	Eigen::VectorXd b = Eigen::VectorXd::Zero(nV * 2);
	Eigen::VectorXd xk = Eigen::VectorXd::Zero(nV * 2);
	Eigen::VectorXd pk = Eigen::VectorXd::Zero(nV * 2);

	for (auto v_it = ptr_object_->vertices_begin(); v_it != ptr_object_->vertices_end(); v_it++)
	{
		int i = (*v_it).idx();
		xk(i) = ptr_object_->property(cogs, *v_it)[0];
		xk(i + nV) = ptr_object_->property(cogs, *v_it)[1];
		pk(i) = p[i][0];
		pk(i + nV) = p[i][1];
		if (i == fixed_vertices_[0])
		{
			b(i) += fixed_coords_[0][0];
			b(i + nV) += fixed_coords_[0][1];
			continue;
		}
		/*if (i == fixed_vertices_[1])
		{
			b_x(i) += fixed_coords_[1][0];
			b_y(i) += fixed_coords_[1][1];
			continue;
		}*/
		Eigen::Vector2d sum = Eigen::Vector2d::Zero();
		for (MyMesh::VertexOHalfedgeIter voh_it = ptr_object_->voh_begin(*v_it); voh_it != ptr_object_->voh_end(*v_it); ++voh_it)
		{
			// Left Triangle
			MyMesh::FaceHandle triangle = ptr_object_->face_handle(*voh_it);
			size_t t = 0;
			int index;
			double coe_ijxx = 0;
			double coe_ijxy = 0;
			double coe_ijyx = 0;
			double coe_ijyy = 0;
			if (!ptr_object_->is_boundary(*voh_it))
			{
				t = triangle.idx();
				index = 0;
				for (MyMesh::FaceVertexIter fv_it = ptr_object_->fv_iter(triangle); fv_it.is_valid(); ++fv_it) {
					if ((*fv_it).idx() == i) break;
					index++;
				}
				auto diffx = x_plane[t][index] - x_plane[t][index == 2 ? 0 : index + 1];
				coe_ijxx += cot[t][index] * W[t](0, 0);
				coe_ijxy += cot[t][index] * W[t](0, 1);
				coe_ijyx += cot[t][index] * W[t](1, 0);
				coe_ijyy += cot[t][index] * W[t](1, 1);
				sum += cot[t][index] * W[t] * L[t] * Eigen::Vector2d(diffx[0], diffx[1]);
			}
			// Right Triangle
			triangle = ptr_object_->face_handle(ptr_object_->opposite_halfedge_handle(*voh_it));
			if (!ptr_object_->is_boundary(ptr_object_->opposite_halfedge_handle(*voh_it)))
			{
				t = triangle.idx();
				index = 0;
				for (MyMesh::FaceVertexIter fv_it = ptr_object_->fv_iter(triangle); fv_it.is_valid(); ++fv_it) {
					if ((*fv_it).idx() == i) break;
					index++;
				}
				auto diffx = x_plane[t][index] - x_plane[t][index == 0 ? 2 : index - 1];
				coe_ijxx += cot[t][index == 0 ? 2 : index - 1] * W[t](0, 0);
				coe_ijxy += cot[t][index == 0 ? 2 : index - 1] * W[t](0, 1);
				coe_ijyx += cot[t][index == 0 ? 2 : index - 1] * W[t](1, 0);
				coe_ijyy += cot[t][index == 0 ? 2 : index - 1] * W[t](1, 1);
				sum += cot[t][index == 0 ? 2 : index - 1] * W[t] * L[t] * Eigen::Vector2d(diffx[0], diffx[1]);
			}
			// Set
			int j = ptr_object_->to_vertex_handle(*voh_it).idx();
			if (j == fixed_vertices_[0])
			{
				b(i) -= (-coe_ijxx) * fixed_coords_[0][0] + (-coe_ijxy) * fixed_coords_[0][1];
				b(i + nV) -= (-coe_ijyx) * fixed_coords_[0][0] + (-coe_ijyy) * fixed_coords_[0][1];
			}
			/*else if (j == fixed_vertices_[1])
			{
				b_x(i) -= (-coe_ij) * fixed_coords_[1][0];
				b_y(i) -= (-coe_ij) * fixed_coords_[1][1];
			}*/

		}

		b(i) += sum(0);
		b(i + nV) += sum(1);
	}

	float al_a = pk.transpose() * A.transpose() * A * pk;
	float al_b_2 = (pk.transpose() * A.transpose() * (A * xk - b))(0,0);

	float alpha_t = -al_b_2 / al_a;
	if (alpha_t > alpha_max) alpha = alpha_max;
	else if (alpha_t > 0) alpha = alpha_t;

	float x_min = INT_MAX;
	float x_max = INT_MIN;
	float y_min = INT_MAX;
	float y_max = INT_MIN;
	for (auto v_it = ptr_object_->vertices_begin(); v_it != ptr_object_->vertices_end(); v_it++) {
		int i = (*v_it).idx();
		MyMesh::Point new_vert = MyMesh::Point(xk(i) + alpha * pk(i), xk(i + nV) + alpha * pk(i + nV), 0);

		if (new_vert[0] > x_max) x_max = new_vert[0];
		else if (new_vert[0] < x_min) x_min = new_vert[0];
		if (new_vert[1] > y_max) y_max = new_vert[1];
		else if (new_vert[1] < y_min) y_min = new_vert[1];

		ptr_object_->property(cogs, *v_it) = new_vert;
	}

	FitIn(x_min, x_max, y_min, y_max);
}
