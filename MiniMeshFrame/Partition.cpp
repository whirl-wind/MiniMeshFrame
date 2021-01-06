#include "Partition.h"
#include <gl/glut.h>
#include <gl/GLU.h>
#include "colorbar.h"

Partition::Partition()
{
	classnum = 0;
}

Partition::~Partition()
{
}

void Partition::LoadDDG(DifferentialGeometry* ddg)
{
	DDG = ddg;
}

int Partition::get_classnum()
{
	return classnum;
}

void Partition::partition()
{
}

void Partition::Draw()
{
}

Pa_VSA::Pa_VSA(MyMesh *ptr, DifferentialGeometry* ddg, int n)
{
	ptr_mesh_ = ptr;
	classnum = n;
	LoadDDG(ddg);
	ptr_mesh_->add_property(label);

	init();
	float eps = 1E-6;
	float loss = 0;
	float d = 1;
	int iter = 0;
	int iter_max = 100;
	while(d > eps) {
		partition();
		d = loss;
		loss = proxy();
		d = abs(loss - d);
		iter++;
		if (iter >= iter_max) break;
	}
	std::cout << "After " << iter << " iters with loss: " << loss << std::endl;
}

Pa_VSA::~Pa_VSA()
{
}

void Pa_VSA::init()
{
	//初步分类
	float max_a = 0;
	std::vector<Eigen::Vector3f> p;
	p.resize(classnum);
	conqueue.resize(ptr_mesh_->n_faces());
	
	int temp = 0;
	for (MyMesh::FaceIter f_it = ptr_mesh_->faces_begin(); f_it != ptr_mesh_->faces_end(); ++f_it) {
		ptr_mesh_->property(label, *f_it) = -1;
		if (DDG->FaceArea[(*f_it).idx()] > max_a) {
			temp = (*f_it).idx();
			max_a = DDG->FaceArea[(*f_it).idx()];
		}
		conqueue[(*f_it).idx()] = false;
	}
	ptr_mesh_->property(label, ptr_mesh_->face_handle(temp)) = 0;
	conqueue[temp] = true;
	p[0] = Eigen::Vector3f(ptr_mesh_->normal(ptr_mesh_->face_handle(temp)).data());

	std::vector<float> dist;
	dist.resize(ptr_mesh_->n_faces());
	float max_d = 0;
	for (MyMesh::FaceIter f_it = ptr_mesh_->faces_begin(); f_it != ptr_mesh_->faces_end(); ++f_it) {
		Eigen::Vector3f nf = Eigen::Vector3f(ptr_mesh_->normal(f_it).data());
		float d = (nf - p[0]).squaredNorm() * DDG->FaceArea[(*f_it).idx()];
		dist[(*f_it).idx()] = d;
		ptr_mesh_->property(label, *f_it) = 0;
		if (max_d <= d && !conqueue[(*f_it).idx()]) {
			max_d = d;
			temp = (*f_it).idx();
		}
	}
	ptr_mesh_->property(label, ptr_mesh_->face_handle(temp)) = 1;
	conqueue[temp] = true;
	p[1] = Eigen::Vector3f(ptr_mesh_->normal(ptr_mesh_->face_handle(temp)).data());
	for (int i = 1; i < classnum - 1; i++) {
		max_d = 0;
		for (MyMesh::FaceIter f_it = ptr_mesh_->faces_begin(); f_it != ptr_mesh_->faces_end(); ++f_it) {
			Eigen::Vector3f nf = Eigen::Vector3f(ptr_mesh_->normal(f_it).data());
			float d = (nf - p[i]).squaredNorm() * DDG->FaceArea[(*f_it).idx()];
			if (dist[(*f_it).idx()] > d) {
				dist[(*f_it).idx()] = d;
				ptr_mesh_->property(label, *f_it) = i;
			}

			if (max_d <= dist[(*f_it).idx()] && !conqueue[(*f_it).idx()]) {
				max_d = dist[(*f_it).idx()];
				temp = (*f_it).idx();
			}
		}
		ptr_mesh_->property(label, ptr_mesh_->face_handle(temp)) = i + 1;
		conqueue[temp] = true;
		p[i + 1] = Eigen::Vector3f(ptr_mesh_->normal(ptr_mesh_->face_handle(temp)).data());
	}
	for (MyMesh::FaceIter f_it = ptr_mesh_->faces_begin(); f_it != ptr_mesh_->faces_end(); ++f_it) {
		Eigen::Vector3f nf = Eigen::Vector3f(ptr_mesh_->normal(f_it).data());
		float d = (nf - p[classnum - 1]).squaredNorm() * DDG->FaceArea[(*f_it).idx()];
		if (dist[(*f_it).idx()] > d) {
			dist[(*f_it).idx()] = d;
			ptr_mesh_->property(label, *f_it) = classnum - 1;
		}
	}

	//同化
	bool s = true;
	while (s) {
		s = false;
		for (MyMesh::FaceIter f_it = ptr_mesh_->faces_begin(); f_it != ptr_mesh_->faces_end(); ++f_it) {
			int* a = new int[3];
			int idx = 0;
			for (MyMesh::FaceFaceIter ff_it = ptr_mesh_->ff_begin(*f_it); ff_it != ptr_mesh_->ff_end(*f_it); ++ff_it) {
				a[idx] = ptr_mesh_->property(label, *f_it) - ptr_mesh_->property(label, *ff_it);
				idx++;
			}
			if (a[0] + a[1] == 0 || a[1] + a[2] == 0 || a[2] + a[0] == 0) continue;
			else if (a[0] == a[1]) {
				ptr_mesh_->property(label, *f_it) -= a[0];
				s = true;
				break;
			}
			else if (a[1] == a[2]) {
				ptr_mesh_->property(label, *f_it) -= a[1];
				s = true;
				break;
			}
			else if (a[2] == a[0]) {
				ptr_mesh_->property(label, *f_it) -= a[2];
				s = true;
				break;
			}
		}
	}

	//真正的严格分类
	R.clear();
	for (int i = 0; i < ptr_mesh_->n_faces(); i++) {
		conqueue[i] = false;
	}
	s = true;
	int classidx = 0;
	while (s) {
		s = false;
		for (MyMesh::FaceIter f_it = ptr_mesh_->faces_begin(); f_it != ptr_mesh_->faces_end(); ++f_it) {
			if (!conqueue[(*f_it).idx()]) {
				R.resize(classidx + 1);
				std::queue<MyMesh::FaceHandle> oneclass;
				oneclass.push(*f_it);
				conqueue[(*f_it).idx()] = true;
				sort_into_oneclass(oneclass, classidx);
				s = true;
				classidx++;
				break;
			}
		}
	}
	classnum = classidx;
	std::cout << "Mesh was classified into " << classnum << " planes" << std::endl;
	P.resize(classnum);
	for (int i = 0; i < classnum; i++) {
		Eigen::Vector3f pi = Eigen::Vector3f(0, 0, 0);
		for (int j = 0; j < R[i].size(); j++) {
			ptr_mesh_->property(label, ptr_mesh_->face_handle(R[i][j])) = i;
			pi += DDG->FaceArea[R[i][j]] * Eigen::Vector3f(ptr_mesh_->normal(ptr_mesh_->face_handle(R[i][j])).data());
		}
		P[i] = pi.normalized();
	}
}

float Pa_VSA::proxy()
{
	float loss = 0.0;
	for (int i = 0; i < classnum; i++) {
		Eigen::Vector3f pi = Eigen::Vector3f(0, 0, 0);
		for (int j = 0; j < R[i].size(); j++) {
			pi += DDG->FaceArea[R[i][j]] * Eigen::Vector3f(ptr_mesh_->normal(ptr_mesh_->face_handle(R[i][j])).data());
		}
		loss += pi.squaredNorm();
		P[i] = pi.normalized();
	}
	return loss;
}

void Pa_VSA::partition()
{
	//找到区域初始种子点
	for (int i = 0; i < ptr_mesh_->n_faces(); i++) {
		conqueue[i] = false;
	}
	std::vector< std::tuple<MyMesh::FaceHandle, int> > sortclass;
	for (int i = 0; i < classnum; i++) {
		float min_d = INT_MAX;
		int temp;
		for (int j = 0; j < R[i].size(); j++) {
			Eigen::Vector3f nf = Eigen::Vector3f(ptr_mesh_->normal(ptr_mesh_->face_handle(R[i][j])).data());
			float d = (nf - P[i]).squaredNorm() * DDG->FaceArea[R[i][j]];
			if (min_d > d) {
				min_d = d;
				temp = R[i][j];
			}
		}
		R[i].clear();
		R[i].push_back(temp);
		ptr_mesh_->property(label, ptr_mesh_->face_handle(temp)) = i;
		conqueue[temp] = true;
		for (MyMesh::FaceFaceIter ff_it = ptr_mesh_->ff_begin(ptr_mesh_->face_handle(temp)); ff_it != ptr_mesh_->ff_end(ptr_mesh_->face_handle(temp)); ++ff_it) {
			sortclass.push_back(std::tuple<MyMesh::FaceHandle, int>(*ff_it, i));
		}
	}
	//扩张区域
	while (!sortclass.empty()) {
		float min_d = INT_MAX;
		int temp;
		MyMesh::FaceHandle f;
		int classidx;
		for (int i = 0; i < sortclass.size(); i++) {
			std::tie(f, classidx) = sortclass[i];
			Eigen::Vector3f nf = Eigen::Vector3f(ptr_mesh_->normal(f).data());
			float d = (nf - P[classidx]).squaredNorm() * DDG->FaceArea[f.idx()];
			if (min_d > d) {
				min_d = d;
				temp = i;
			}
		}
		std::tie(f, classidx) = sortclass[temp];
		if (temp != sortclass.size() - 1) {
			std::iter_swap(sortclass.begin() + temp, sortclass.end() - 1);
		}
		sortclass.pop_back();

		temp = f.idx();
		if (!conqueue[temp]) {
			conqueue[temp] = true;
			R[classidx].push_back(temp);
			ptr_mesh_->property(label, ptr_mesh_->face_handle(temp)) = classidx;
			for (MyMesh::FaceFaceIter ff_it = ptr_mesh_->ff_begin(f); ff_it != ptr_mesh_->ff_end(f); ++ff_it) {
				if (!conqueue[(*ff_it).idx()]) {
					sortclass.push_back(std::tuple<MyMesh::FaceHandle, int>(*ff_it, classidx));
				}
			}
		}
	}
	//同化
	bool s = true;
	while (s) {
		s = false;
		for (MyMesh::FaceIter f_it = ptr_mesh_->faces_begin(); f_it != ptr_mesh_->faces_end(); ++f_it) {
			int* a = new int[3];
			int idx = 0;
			for (MyMesh::FaceFaceIter ff_it = ptr_mesh_->ff_begin(*f_it); ff_it != ptr_mesh_->ff_end(*f_it); ++ff_it) {
				a[idx] = ptr_mesh_->property(label, *f_it) - ptr_mesh_->property(label, *ff_it);
				idx++;
			}
			if (a[0] + a[1] == 0 || a[1] + a[2] == 0 || a[2] + a[0] == 0) continue;
			else if (a[0] == a[1]) {
				ptr_mesh_->property(label, *f_it) -= a[0];
				s = true;
				break;
			}
			else if (a[1] == a[2]) {
				ptr_mesh_->property(label, *f_it) -= a[1];
				s = true;
				break;
			}
			else if (a[2] == a[0]) {
				ptr_mesh_->property(label, *f_it) -= a[2];
				s = true;
				break;
			}
		}
	}
}

void Pa_VSA::sort_into_oneclass(std::queue<MyMesh::FaceHandle> &oneclass, int classidx)
{
	if (oneclass.empty()) return;
	MyMesh::FaceHandle f = oneclass.front();
	oneclass.pop();
	R[classidx].push_back(f.idx());

	for (MyMesh::FaceFaceIter ff_it = ptr_mesh_->ff_begin(f); ff_it != ptr_mesh_->ff_end(f); ++ff_it) {
		if (ptr_mesh_->property(label, *ff_it) == ptr_mesh_->property(label, f) && !conqueue[(*ff_it).idx()]) {
			oneclass.push(*ff_it);
			conqueue[(*ff_it).idx()] = true;
		}
	}
	sort_into_oneclass(oneclass, classidx);
}

void Pa_VSA::Draw()
{
	ColorBar cb;
	cb.setRange(0, classnum);
	for (MyMesh::FaceIter f_it = ptr_mesh_->faces_begin(); f_it != ptr_mesh_->faces_end(); ++f_it) {
		QColor color = cb.getColor(classnum - ptr_mesh_->property(label, *f_it) - 0.5);
		glColor3f(color.redF(), color.greenF(), color.blueF());
		glBegin(GL_TRIANGLES);
		for (MyMesh::FaceVertexIter fv_it = ptr_mesh_->fv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
			if (ptr_mesh_->status(*fv_it).deleted()) continue;
			glNormal3fv(ptr_mesh_->normal(*fv_it).data());
			glVertex3fv(ptr_mesh_->point(*fv_it).data());
		}
		glEnd();
	}

	glColor3f(1.0, 1.0, 1.0);
}
