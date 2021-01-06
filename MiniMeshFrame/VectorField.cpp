#include "VectorField.h"
#include "Polynomial.h"
#include <gl/glut.h>
#include <gl/GLU.h>
#include "colorbar.h"

VectorField::VectorField()
{
}

VectorField::VectorField(MyMesh* ptr, int mm, int nn)
{
	if (mm <= 0) m = 1;
	else m = mm;
	if (nn <= 0) n = 1;
	else n = nn;

	UpdateMesh(ptr);
}

VectorField::~VectorField()
{
	ptr_mesh_->remove_property(vectorfield_theta);
	ptr_mesh_->remove_property(localcoord);
	ptr_mesh_->remove_property(vectorfield_idx);
}

void VectorField::BuildLocal()
{
	for (MyMesh::FaceIter f_it = ptr_mesh_->faces_begin(); f_it != ptr_mesh_->faces_end(); ++f_it) {

		MyMesh::FaceHalfedgeIter fh_it = ptr_mesh_->fh_iter(*f_it);
		auto v1 = ptr_mesh_->point(ptr_mesh_->from_vertex_handle(fh_it));
		auto v2 = ptr_mesh_->point(ptr_mesh_->to_vertex_handle(fh_it));
		std::complex<double> z{ 1.0, 0.0 };
		ptr_mesh_->property(localcoord, fh_it) = z;

		fh_it++;
		auto v3 = ptr_mesh_->point(ptr_mesh_->to_vertex_handle(fh_it));
		double x1 = (v2 - v1).norm();
		double s = ((v2 - v1) % (v3 - v1)).norm();
		double x2 = ((v3 - v1) | (v2 - v1).normalized());
		
		double a = x2 - x1;
		double b = s / x1;
		double c = sqrt(a * a + b * b);
		z.real(a / c);
		z.imag(b / c);
		ptr_mesh_->property(localcoord, fh_it) = z;

		fh_it++;
		a = -x2;
		b = -b;
		c = sqrt(a * a + b * b);
		z.real(a / c);
		z.imag(b / c);
		ptr_mesh_->property(localcoord, fh_it) = z;
	}
}

void VectorField::SetFix(vector<size_t> fidx, vector< vector<double> > fvectorfield)
{
	int k = fidx.size();
	if (k != fvectorfield.size()) {
		std::cout << "VectorField: wrong fixing size" << std::endl;
		return;
	}

	for (int i = 0; i < k; i++) {
		int kk = fvectorfield[i].size();
		if (kk < n) {
			std::cout << "VectorField: wrong fixing size" << std::endl;
			return;
		}
		ptr_mesh_->property(vectorfield_theta, ptr_mesh_->face_handle(fidx[i])).assign(fvectorfield[i].begin(), fvectorfield[i].begin() + n);
	}

	fix_idx.resize(k);
	fix_vectorfield.resize(k);
	fix_idx.assign(fidx.begin(), fidx.end());
	for (int i = 0; i < k; i++) {
		fix_vectorfield[i].assign(fvectorfield[i].begin(), fvectorfield[i].begin() + n);
	}
}

void VectorField::UpdateMesh(MyMesh* ptr)
{
	ptr_mesh_ = ptr;
	ptr_mesh_->add_property(vectorfield_theta);
	ptr_mesh_->add_property(localcoord);
	ptr_mesh_->add_property(vectorfield_idx);

	BuildLocal();
}

void VectorField::BuildVectorField()
{
}

void VectorField::SetColor()
{

	std::cout << "SetColor Begin" << std::endl;
	MyMesh::FaceHandle f_it;
	for (int k = 0; k < fix_idx.size(); k++) {
		f_it = ptr_mesh_->face_handle(fix_idx[k]);
		ptr_mesh_->property(vectorfield_idx, f_it).resize(n);
		for (int i = 0; i < n; i++) {
			ptr_mesh_->property(vectorfield_idx, f_it)[i] = i;
		}
	}
	
	
	for (MyMesh::FaceHalfedgeIter fh_it = ptr_mesh_->fh_begin(f_it); fh_it != ptr_mesh_->ff_end(f_it); fh_it++) {
		SetColorIter(f_it, fh_it);
	}

	std::cout << "Done\n" << std::endl;
}

void VectorField::SetColorIter(MyMesh::FaceHandle f_it, MyMesh::HalfedgeHandle fh_it)
{
	MyMesh::FaceHandle cf_it = ptr_mesh_->face_handle(ptr_mesh_->opposite_halfedge_handle(fh_it));
	if (!ptr_mesh_->property(vectorfield_idx, cf_it).empty()) return;

	auto ef_ = std::conj(-ptr_mesh_->property(localcoord, fh_it));
	auto eg_ = std::conj(ptr_mesh_->property(localcoord, ptr_mesh_->opposite_halfedge_handle(fh_it)));

	vector<size_t> array;
	array.resize(n);
	for (size_t i = 0; i < n; i++) {
		array[i] = i;
	}

	ptr_mesh_->property(vectorfield_idx, cf_it).resize(n);
	for (size_t i = 0; i < n; i++) {
		auto uf_theta = ptr_mesh_->property(vectorfield_theta, f_it)[i];
		std::complex<double> uf{ cos(uf_theta), sin(uf_theta) };
		auto ufef_ = uf * ef_;
		double min = INT_MAX;
		size_t idx = array[0];
		for (size_t ii = 0; ii < array.size(); ii++) {
			auto ug_theta = ptr_mesh_->property(vectorfield_theta, cf_it)[array[ii]];
			std::complex<double> ug{ cos(ug_theta), sin(ug_theta) };
			auto z = (ug * eg_ - ufef_);
			double l = pow(z.real(), 2) + pow(z.imag(), 2);
			if (l < min) {
				min = l;
				idx = ii;
			}
		}

		ptr_mesh_->property(vectorfield_idx, cf_it)[array[idx]] = ptr_mesh_->property(vectorfield_idx, f_it)[i];
		// 这里移除的算法复杂度是O(1)，将待删除元素与最后一个元素交换再pop_back
		if (idx == array.size() - 1)
		{
			array.pop_back();
		}
		else
		{
			iter_swap(array.begin() + idx, array.end() - 1);
			array.pop_back();
		}
	}
	array.clear();
	array.shrink_to_fit();

	for (MyMesh::FaceHalfedgeIter cfh_it = ptr_mesh_->fh_begin(cf_it); cfh_it != ptr_mesh_->ff_end(cf_it); cfh_it++) {
		SetColorIter(cf_it, cfh_it);
	}
}

void VectorField::Draw(float size)
{
	ColorBar cb;
	cb.setRange(0, n);

	for (MyMesh::FaceIter f_it = ptr_mesh_->faces_begin(); f_it != ptr_mesh_->faces_end(); ++f_it) {
		for (int i = 0; i < n; i++) {
			QColor color = cb.getColor(i + 0.5);
			glColor3f(color.redF(), color.greenF(), color.blueF());

			MyMesh::FaceHalfedgeIter fh_it = ptr_mesh_->fh_iter(*f_it);
			auto v1 = ptr_mesh_->point(ptr_mesh_->from_vertex_handle(fh_it));
			auto v2 = ptr_mesh_->point(ptr_mesh_->to_vertex_handle(fh_it));

			fh_it++;
			auto v3 = ptr_mesh_->point(ptr_mesh_->to_vertex_handle(fh_it));
			double x1 = (v2 - v1).norm();
			double s = ((v2 - v1) % (v3 - v1)).norm();
			double x2 = ((v3 - v1) | (v2 - v1).normalized());
			double h = s / x1;

			auto center = (v1 + v2 + v3) / 3.0;
			double c;
			if (ptr_mesh_->property(vectorfield_idx, f_it).empty()) {
				c = ptr_mesh_->property(vectorfield_theta, f_it)[i];
			}
			else {
				c = ptr_mesh_->property(vectorfield_theta, f_it)[ptr_mesh_->property(vectorfield_idx, f_it)[i]];
			}

			glBegin(GL_LINES);
			for (int j = 0; j < m; j++) {
				double ky = sin(c + 2 * M_PI * j / m) / h;
				double kx = (cos(c + 2 * M_PI * j / m) - ky * x2) / x1;

				auto vector = ky * (v3 - v1) + kx * (v2 - v1);
				vector.normalize();

				glNormal3fv(ptr_mesh_->normal(f_it).data());
				glVertex3f(center[0], center[1], center[2]);
				glVertex3f(center[0] + size * vector[0], center[1] + size * vector[1], center[2] + size * vector[2]);
			}
			glEnd();
		}
	}

	glColor3f(1.0, 1.0, 1.0);
}

NpolyVectorFields::NpolyVectorFields(MyMesh* ptr, int nn, int mm)
{
	if (mm <= 0) m = 1;
	else m = mm;
	if (nn <= 0) n = 1;
	else n = nn;

	UpdateMesh(ptr);
}

NpolyVectorFields::~NpolyVectorFields()
{
	ptr_mesh_->remove_property(vectorfield_theta);
	ptr_mesh_->remove_property(localcoord);
	ptr_mesh_->remove_property(vectorfield_idx);
}

void NpolyVectorFields::BuildLaplace()
{
	Lm.clear();
	Lm.resize(n + 1);
	//solver.clear();
	//solver.resize(n + 1);

	int nF = ptr_mesh_->n_faces();
	for (int i = 0; i < n + 1; i++) {
		Eigen::SparseMatrix< std::complex<double> > L(nF, nF);
		std::vector<  Eigen::Triplet< std::complex<double> > > L_Triplet;
	
		for (MyMesh::FaceIter f_it = ptr_mesh_->faces_begin(); f_it != ptr_mesh_->faces_end(); ++f_it) {
			int idx_f = (*f_it).idx();
			bool s = false;
			for (int j = 0; j < fix_idx.size(); j++) {
				if (idx_f == fix_idx[j]) {
					L_Triplet.push_back(Eigen::Triplet< std::complex<double> >(idx_f, idx_f, std::complex<double>{ 1, 0 }));
					s = true;
					break;
				}
			}
			if (s) continue;
			
			std::complex<double> ff{ 0,0 };
			for (MyMesh::FaceHalfedgeIter fh_it = ptr_mesh_->fh_iter(*f_it); fh_it.is_valid(); ++fh_it) {
				int idx_g = ptr_mesh_->face_handle(ptr_mesh_->opposite_halfedge_handle(fh_it)).idx();
				s = false;
				for (int j = 0; j < fix_idx.size(); j++) {
					if (idx_g == fix_idx[j]) {
						s = true;
						break;
					}
				}

				auto ef_ = std::conj(ptr_mesh_->property(localcoord, fh_it));
				auto eg_ = std::conj(-ptr_mesh_->property(localcoord, ptr_mesh_->opposite_halfedge_handle(fh_it)));

				//std::cout << ptr_mesh_->property(localcoord, fh_it) << " _ " << ef_ << std::endl;
				//std::cout << -ptr_mesh_->property(localcoord, ptr_mesh_->opposite_halfedge_handle(fh_it)) << " _ " << eg_ << std::endl;

				ff += pow(ef_, n - i);

				if (s) continue;
				L_Triplet.push_back(Eigen::Triplet< std::complex<double> >(idx_f, idx_g, -pow(eg_, n - i)));
				//L.coeffRef(idx_f, idx_g) = -pow(eg_, i);
			}
			L_Triplet.push_back(Eigen::Triplet< std::complex<double> >(idx_f, idx_f, ff));
			//L.coeffRef(idx_f, idx_f) = ff;
		}
		/*
		for (vector<Eigen::Triplet< std::complex<double> > >::iterator iter = L_Triplet.begin(); iter != L_Triplet.end(); ++iter)
		{
			std::cout << (*iter).col() << " " << (*iter).row() << " " << (*iter).value() <<" "; //使用 * 访问迭代器所指向的元素
		}
		std::cout << std::endl;
		*/
		L.setFromTriplets(L_Triplet.begin(), L_Triplet.end());
		L.makeCompressed();
		
		//std::cout << "\n" << i << ": \n"<< L << std::endl;

		Lm[i] = L;
		/*
		solver[i] = new Eigen::SparseLU< Eigen::SparseMatrix< std::complex<double> > >();
		solver[i]->compute(L);
		if (solver[i]->info() != Eigen::Success) {
			std::cout << "NpolyVectorFields: " << i << "th decomposition failed!" << std::endl;
		}*/
	}
}

void NpolyVectorFields::BuildVectorField()
{
	vector<  Eigen::VectorXcd* > b;
	vector<  Eigen::VectorXcd* > x;
	int nF = ptr_mesh_->n_faces();
	b.resize(n + 1);
	x.resize(n + 1);
	for (int i = 0; i < n + 1; i++) {
		b[i] = new Eigen::VectorXcd(nF);
		(*b[i]).setZero();
	}

	for (int j = 0; j < fix_idx.size(); j++) {
		vector<std::complex<double> > c;
		//c.push_back(std::complex<double>{1, 0});
		c.push_back(-std::complex<double>{cos(fix_vectorfield[j][0]), sin(fix_vectorfield[j][0])});
		c.push_back(std::complex<double>{1, 0});
		Polynomial<std::complex<double> > p(c.begin(), c.end());
		//std::cout << p << std::endl;
		//std::cout << p[1] << std::endl;
		for (int k = 1; k < n; k++) {
			vector<std::complex<double> > co;
			//co.push_back(std::complex<double>{1, 0});
			co.push_back(-std::complex<double>{cos(fix_vectorfield[j][k]), sin(fix_vectorfield[j][k])});
			co.push_back(std::complex<double>{1, 0});
			p *= Polynomial<std::complex<double> >(co.begin(), co.end());
		}
		
		int idx_g = fix_idx[j];
		for (int i = 0; i < n + 1; i++) {
			(*b[i])(idx_g) =  p[i];
		}
		for (MyMesh::FaceHalfedgeIter fh_it = ptr_mesh_->fh_iter(ptr_mesh_->face_handle(fix_idx[j])); fh_it.is_valid(); ++fh_it) {
			int idx_f = ptr_mesh_->face_handle(ptr_mesh_->opposite_halfedge_handle(fh_it)).idx();
			auto eg_ = std::conj(-ptr_mesh_->property(localcoord, fh_it));
			auto ef_ = std::conj(ptr_mesh_->property(localcoord, ptr_mesh_->opposite_halfedge_handle(fh_it)));

			for (int i = 0; i < n + 1; i++) {
				(*b[i])(idx_f) -= -pow(eg_, n - i) * p[i];
			}
		}
	}

	for (int i = 0; i < n + 1; i++) {
		//std::cout << i << std::endl;
		//std::cout << (*b[i]).transpose() << std::endl;
		Eigen::SparseLU< Eigen::SparseMatrix< std::complex<double> > > solver;
		//std::cout << Lm[i].toDense() << std::endl;
		solver.compute(Lm[i]);
		Eigen::VectorXcd xx = solver.solve(*b[i]);
		x[i] = new Eigen::VectorXcd(xx);
		//std::cout << *x[i] << std::endl;
		//std::cout << (Lm[i] * (*x[i]) - (*b[i])).transpose() << std::endl;
	}

	for (MyMesh::FaceIter f_it = ptr_mesh_->faces_begin(); f_it != ptr_mesh_->faces_end(); ++f_it) {
		int idx_f = (*f_it).idx();
		
		ptr_mesh_->property(vectorfield_idx, f_it).clear();

		bool s = false;
		for (int j = 0; j < fix_idx.size(); j++) {
			if (idx_f == fix_idx[j]) {
				s = true;
				break;
			}
		}
		if (s) continue;
		
		ptr_mesh_->property(vectorfield_theta, f_it).clear();
		vector<std::complex<double> > c;
		for (int i = 0; i < n + 1; i++) {
			c.push_back((*x[i])(idx_f));
		}
		Polynomial< std::complex<double> > p(c.begin(), c.end());
		int degree = p.Degree();
		//std::cout << p << std::endl;
		
		std::vector< std::complex<long double> > roots = p.FindRoots();

		for (int k = 0; k < degree; k++) {
			double a = roots[k].real();
			double b = roots[k].imag();
			double c = sqrt(a * a + b * b);

			double theta = acos(a / c);
			if (b < 0) theta *= -1.0;

			ptr_mesh_->property(vectorfield_theta, f_it).push_back(theta);
			//std::cout << roots[k] << ": " << theta << "    ";
		}
		//std::cout << std::endl;
	}

	std::cout << "Done\n" << std::endl;
}

void NpolyVectorFields::UpdateMesh(MyMesh* ptr)
{
	ptr_mesh_ = ptr;
	ptr_mesh_->add_property(vectorfield_theta);
	ptr_mesh_->add_property(vectorfield_idx);
	ptr_mesh_->add_property(localcoord);
	int nF = ptr_mesh_->n_faces();
	if (nF < 2) return;

	BuildLocal();

	fix_idx.clear();
	fix_vectorfield.clear();

	vector<size_t> fidx;
	vector<vector< double > > fvectorfield;
	fidx.push_back(0);
	fidx.push_back(nF / 2);

	vector< double > fix;
	for (int i = 0; i < n; i++) {
		fix.push_back(2 * M_PI / n * i);
	}
	
	fvectorfield.push_back(fix);
	fvectorfield.push_back(fix);

	SetFix(fidx, fvectorfield);
	BuildLaplace();
}
