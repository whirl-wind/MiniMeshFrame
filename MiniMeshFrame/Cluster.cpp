#include "Cluster.h"
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <stdlib.h>
#include <math.h>
#include <limits.h> //for INT_MIN INT_MAX
#include <time.h>
#include <iostream>
#include <gl/glut.h>
#include <gl/GLU.h>
#include "colorbar.h"

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



Cluster::Cluster()
{
	nonlinear_type_ = kDefault;

	ptr_object_ = new MyMesh();
	class_num = 1;
	triangles.clear();
	canonical_triangles.clear();
}

Cluster::Cluster(MyMesh * ptr_mesh_, int cla)
{
	nonlinear_type_ = kLM;

	ptr_object_ = ptr_mesh_;
	class_num = cla;

	if (class_num < 1) {
		std::cout << "Class number must be positive!" << std::endl;
		class_num = 1;
	}
	triangles.clear();
	canonical_triangles.clear();
	
	vector<Matrix3f> temp_triangles;
	temp_triangles.clear();
	for (MyMesh::FaceIter f_it = ptr_object_->faces_begin(); f_it != ptr_object_->faces_end(); f_it++) {
		Matrix3f temp;
		int i = 0;
		for (MyMesh::FaceVertexIter fv_it = ptr_object_->fv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
			if (i >= 3) {
				std::cout << "Find a face which is not a triangle!" << std::endl;
				break;
			}
			auto vert = ptr_object_->point(fv_it.handle());
			temp(0, i) = vert[0];
			temp(1, i) = vert[1];
			temp(2, i) = vert[2];
			i++;
		}
		triangles.push_back(temp);
		canonical_triangles.push_back(Into_Canonical_Triangle(temp));
	}
	std::cout << "Clustering, please wait..." << std::endl;
	clock_t start = clock();
	kmeans();
	clock_t end = clock();
	std::cout << "Finished,thank you! Lasting: " << end - start << " ms" << std::endl;
	std::cout << std::endl;
}


Cluster::~Cluster()
{
}

int Cluster::get_minIdx(int idx_f)
{
	return clusterAssment[idx_f].min_Index;
}

Vector3f Cluster::Into_Canonical_Triangle(Matrix3f t)
//(a, xb, h)
{
	Vector3f CanonicalTriangle;
	float a, b, c, temp, p;
	a = sqrt(((t.col(0) - t.col(1)).transpose()*(t.col(0) - t.col(1)))(0, 0));
	b = sqrt(((t.col(0) - t.col(2)).transpose()*(t.col(0) - t.col(2)))(0, 0));
	c = sqrt(((t.col(2) - t.col(1)).transpose()*(t.col(2) - t.col(1)))(0, 0));

	
	if (a < b) { temp = b; b = a; a = temp; }
	if (a < c) { temp = a; a = c; c = temp; }
	if (b < c) { temp = c; c = b; b = temp; }
	
	p = (a + b + c) / 2;
		
	CanonicalTriangle[0] = a;
	CanonicalTriangle[2] = 2 * sqrt(p*(p - a)*(p - b)*(p - c)) / a;
	CanonicalTriangle[1] = sqrt(b * b - CanonicalTriangle[2] * CanonicalTriangle[2]);

	return CanonicalTriangle;
}

Matrix3f Cluster::Into_Triangle(Vector3f t)
{
	Matrix3f Triangle;
	Triangle << 0,		0,		0,
				0,		0,		t(2),
				0,		t(0),	t(1);
	return Triangle;
}

void Cluster::Set_Canonical_Triangle()
{
	canonical_centroids_triangles.clear();
	vector< int > ::iterator it = centroids_index.begin();
	while (it != centroids_index.end()) {
		canonical_centroids_triangles.push_back(canonical_triangles[(*it)]);
		it++;
	}
}

void Cluster::Set_Rigid_Transformation()
{
	int it_idx = 0;
	vector< Matrix3f > ::iterator it = triangles.begin();
	rigid_transformations.clear();
	while (it != triangles.end()) {
		Matrix3f A = (*it);
		Matrix3f B = Into_Triangle(canonical_centroids_triangles[clusterAssment[it_idx].min_Index]);

		Matrix3f PA, PB, I, M, R0, RB, U, V, R;
		Vector3f T0, J, a_, b_, T;	// 3x1 float matrix.
		I = Matrix3f::Identity(3, 3);
		J = Vector3f::Ones(3);
		float min_Dist = INT_MAX;
			
		int cla = 0;
		int cla_temp = 0;
		for (int j = 0; j < 2; j++) {
			PA.col(0) = A.col((0 + j) % 2);
			PA.col(1) = A.col((1 + j) % 2);
			PA.col(2) = A.col(2);
			for (int i = 0; i < 3; i++) {
				PB.col(0) = B.col((0 + i) % 3);
				PB.col(1) = B.col((1 + i) % 3);
				PB.col(2) = B.col((2 + i) % 3);

				M = Matrix3f::Zero();

				a_ = PA * J / 3.0;
				b_ = PB * J / 3.0;

				for (int k = 0; k < 3; k++) {
					M += (PA.col(k) - a_) * (PB.col(k) - b_).transpose();
				}

				JacobiSVD<Eigen::MatrixXf> svd(M, ComputeThinU | ComputeThinV);
				U = svd.matrixU();
				V = svd.matrixV();

				R0 = U * V.transpose();
				if (R0.determinant() < 0) {
					Matrix3f  S = U.inverse() * M * V.transpose().inverse();
					for (int k = 0; k < 3; k++) {
						if (S(k, k) < 0) {
							V.row(k) *= -1;
							break;
						}
					}
				}
				R0 = U * V.transpose();
				T0 = a_ - R0 * b_;

				float temp_d = 0;
				RB = R0 * PB;
				for (int k = 0; k < 3; k++) {
					temp_d += ((RB.col(k) + T0 - PA.col(k)).transpose()*(RB.col(k) + T0 - PA.col(k)))(0, 0);
				}
				if (temp_d < min_Dist) {
					min_Dist = temp_d;
					R = R0;
					T = T0;
					cla = cla_temp;
				}
				cla_temp++;
			}
		}

		Matrix3f R_;
		switch (cla) {
		case 0: //cla 0 == 0 1 2
			R_ << 1, 0, 0,
				  0, 1, 0,
				  0, 0, 1;
			break;
		case 3: //cla 3 == 1 0 2
			R_ << 0, 1, 0,
				  1, 0, 0,
				  0, 0, 1;
			break;
		case 2: //cla 2 == 1 2 0
			R_ << 0, 1, 0,
				  0, 0, 1,
				  1, 0, 0;
			break;
		case 5: //cla 5 == 0 2 1
			R_ << 1, 0, 0,
				  0, 0, 1,
				  0, 1, 0;
			break;
		case 1: //cla 1 == 2 0 1
			R_ << 0, 0, 1,
				  1, 0, 0,
				  0, 1, 0;
			break;
		case 4: //cla 4 == 2 1 0
			R_ << 0, 0, 1,
				  0, 1, 0,
				  1, 0, 0;
			break;
		default:
			break;
		}
		rigid_transformations.push_back(transNode(T, R, R_));

		it++;
		it_idx++;
	}
}

bool Cluster::addCentroids()
{
	if (nonlinear_type_ == kDefault || nonlinear_type_ == kT) {
		int centroids_num = centroids_index.size();
		if (centroids_num == class_num) return false;

		float  max_Dist = INT_MIN;
		int it_idx = 0;
		int idx = 0;
		int min_index;
		bool s = false;
		vector< tNode > ::iterator it = clusterAssment.begin();
		while (it != clusterAssment.end())
		{	
			//find farthest triangle
			if ((*it).min_Dist > max_Dist && centroids_index[(*it).min_Index] != it_idx) {
				max_Dist = (*it).min_Dist;
				min_index = (*it).min_Index;
				idx = distance(clusterAssment.begin(), it);
				s = true;
			}
			it++;
			it_idx++;
		}
		if (!s) {
			std::cout << "class_num is too big!" << std::endl;
			return false;
		}
		clusterAssment[idx].min_Dist = 0;
		clusterAssment[idx].min_Index = centroids_num;
		centroids_index.push_back(idx);

	}
	else if (nonlinear_type_ == kLM) {
		int centroids_num = canonical_centroids_triangles.size();
		if (centroids_num == class_num) return false;

		float  max_Dist = INT_MIN;
		int it_idx = 0;
		int idx = 0;
		int min_index;
		bool s = false;
		vector< tNode > ::iterator it = clusterAssment.begin();
		while (it != clusterAssment.end())
		{
			if ((*it).min_Dist > max_Dist) {
				max_Dist = (*it).min_Dist;
				min_index = (*it).min_Index;
				idx = distance(clusterAssment.begin(), it);
				s = true;
			}
			it++;
			it_idx++;
		}
		if (!s) {
			std::cout << "class_num is too big!" << std::endl;
			return false;
		}

		int count = 0;
		it = clusterAssment.begin();
		while (it != clusterAssment.end())
		{
			if ((*it).min_Index == min_index) {
				count++;
			}
			it++;
		}
		if (count == 1) std::cout << "class_num  is a little big!" << std::endl;
		clusterAssment[idx].min_Dist = 0;
		clusterAssment[idx].min_Index = centroids_num;
		centroids_index.push_back(idx);
		canonical_centroids_triangles.push_back(canonical_triangles[idx]);
	}
	
	return true;
}

void Cluster::updateCluster()
{
	if (nonlinear_type_ == kDefault) {
		int centroids_num = centroids_index.size();
		vector<tNode>::iterator it_t = clusterAssment.begin();
		float  *min_Dist = new float[centroids_num];
		for (int i = 0; i < centroids_num; i++) {
			min_Dist[i] = 0;
		}
		int idx_t = 0;
		while (it_t != clusterAssment.end()) {
			if ((*it_t).min_Dist == 0) {
				it_t++;
				idx_t++;
				continue;
			}
			else {
				int idx_c = 0;
				vector<int>::iterator it_idx = centroids_index.begin();
				while (it_idx != centroids_index.end()) {
					min_Dist[idx_c] = distTriangles(canonical_triangles[(*it_idx)], canonical_triangles[idx_t]);
					it_idx++;
					idx_c++;
				}

				float min_d = INT_MAX;
				int min_idx = 0;
				for (int i = 0; i < centroids_num; i++) {
					if (min_Dist[i] < min_d) {
						min_d = min_Dist[i];
						min_idx = i;
					}
				}

				(*it_t).min_Dist = min_d;
				(*it_t).min_Index = min_idx;
			}
			it_t++;
			idx_t++;
		}
	}
	else if (nonlinear_type_ == kT) {
		int centroids_num = centroids_index.size();
		vector<tNode>::iterator it_t = clusterAssment.begin();
		float  *min_Dist = new float[centroids_num];
		for (int i = 0; i < centroids_num; i++) {
			min_Dist[i] = 0;
		}
		int idx_t = 0;
		while (it_t != clusterAssment.end()) {
			if ((*it_t).min_Dist == 0) {
				it_t++;
				idx_t++;
				continue;
			}
			else {
				int idx_c = 0;
				vector<int>::iterator it_idx = centroids_index.begin();
				while (it_idx != centroids_index.end()) {
					min_Dist[idx_c] = distTriangles(triangles[(*it_idx)], triangles[idx_t]);
					it_idx++;
					idx_c++;
				}

				float min_d = INT_MAX;
				int min_idx = 0;
				for (int i = 0; i < centroids_num; i++) {
					if (min_Dist[i] < min_d) {
						min_d = min_Dist[i];
						min_idx = i;
					}
				}

				(*it_t).min_Dist = min_d;
				(*it_t).min_Index = min_idx;
			}
			it_t++;
			idx_t++;
		}
	}
	else if (nonlinear_type_ == kLM) {
		int centroids_num = canonical_centroids_triangles.size();
		vector<tNode>::iterator it_t = clusterAssment.begin();
		float  *min_Dist = new float[centroids_num];
		for (int i = 0; i < centroids_num; i++) {
			min_Dist[i] = 0;
		}
		int idx_t = 0;
		while (it_t != clusterAssment.end()) {
			if ((*it_t).min_Dist == 0) {
				it_t++;
				idx_t++;
				continue;
			}
			else {
				int idx_c = 0;
				vector<Vector3f>::iterator it_c = canonical_centroids_triangles.begin();
				while (it_c != canonical_centroids_triangles.end()) {
					min_Dist[idx_c] = distTriangles((*it_c), canonical_triangles[idx_t]);
					it_c++;
					idx_c++;
				}

				float min_d = INT_MAX;
				int min_idx = 0;
				for (int i = 0; i < centroids_num; i++) {
					if (min_Dist[i] < min_d) {
						min_d = min_Dist[i];
						min_idx = i;
					}
				}

				(*it_t).min_Dist = min_d;
				(*it_t).min_Index = min_idx;
			}
			it_t++;
			idx_t++;
		}
	}
}

void Cluster::updateCentroids()
{
	if (nonlinear_type_ == kDefault) {
		int centroids_num = centroids_index.size();
		int  *new_centroids_index = new int[centroids_num];
		for (int i = 0; i < centroids_num; i++) {
			new_centroids_index[i] = centroids_index[i];
		}
		centroids_index.clear();

		vector<tNode>::iterator it_c = clusterAssment.begin();
		float  *min_Dist = new float[centroids_num];
		for (int i = 0; i < centroids_num; i++) {
			min_Dist[i] = INT_MAX;
		}

		int idx_c = 0;
		while (it_c != clusterAssment.end()) {
			int m_idx = (*it_c).min_Index;
			float d = 0;
			vector<float> d_12;
			vector<int> idx_v;
			d_12.clear();
			idx_v.clear();
			vector<tNode>::iterator it_t = clusterAssment.begin();
			int idx_t = 0;
			while (it_t != clusterAssment.end()) {
				if ((*it_t).min_Index == m_idx) {
					float d12 = distTriangles(canonical_triangles[idx_c], canonical_triangles[idx_t]);
					d += d12;
					d_12.push_back(d12);
					idx_v.push_back(idx_t);
				}
				it_t++;
				idx_t++;
			}

			if (d < min_Dist[m_idx]) {
				new_centroids_index[m_idx] = idx_c;
				min_Dist[m_idx] = d;
				vector<int>::iterator it_idv = idx_v.begin();
				int idx_d = 0;
				while (it_idv != idx_v.end()) {
					clusterAssment[(*it_idv)].min_Dist = d_12[idx_d];
					it_idv++;
					idx_d++;
				}
			}
			it_c++;
			idx_c++;
		}
		for (int i = 0; i < centroids_num; i++) {
			centroids_index.push_back(new_centroids_index[i]);//update centroids of cluster[it]
		}
	}
	else if (nonlinear_type_ == kT) {
		int centroids_num = centroids_index.size();
		int  *new_centroids_index = new int[centroids_num];
		for (int i = 0; i < centroids_num; i++) {
			new_centroids_index[i] = centroids_index[i];
		}
		centroids_index.clear();

		vector<tNode>::iterator it_c = clusterAssment.begin();
		float  *min_Dist = new float[centroids_num];
		for (int i = 0; i < centroids_num; i++) {
			min_Dist[i] = INT_MAX;
		}

		int idx_c = 0;
		while (it_c != clusterAssment.end()) {
			int m_idx = (*it_c).min_Index;
			float d = 0;
			vector<float> d_12;
			vector<int> idx_v;
			d_12.clear();
			idx_v.clear();
			vector<tNode>::iterator it_t = clusterAssment.begin();
			int idx_t = 0;
			while (it_t != clusterAssment.end()) {
				if ((*it_t).min_Index == m_idx) {
					float d12 = distTriangles(triangles[idx_c], triangles[idx_t]);
					d += d12;
					d_12.push_back(d12);
					idx_v.push_back(idx_t);
				}
				it_t++;
				idx_t++;
			}

			if (d < min_Dist[m_idx]) {
				new_centroids_index[m_idx] = idx_c;
				min_Dist[m_idx] = d;
				vector<int>::iterator it_idv = idx_v.begin();
				int idx_d = 0;
				while (it_idv != idx_v.end()) {
					clusterAssment[(*it_idv)].min_Dist = d_12[idx_d];
					it_idv++;
					idx_d++;
				}
			}
			it_c++;
			idx_c++;
		}
		for (int i = 0; i < centroids_num; i++) {
			centroids_index.push_back(new_centroids_index[i]);//update centroids of cluster[it]
		}
	}
	else if (nonlinear_type_ == kLM) {
		LM_NonLinearSearch();
	}
}

void Cluster::updateTriangles()
{
	triangles.clear();
	canonical_triangles.clear();
	vector<Matrix3f> temp_triangles;
	temp_triangles.clear();
	for (MyMesh::FaceIter f_it = ptr_object_->faces_begin(); f_it != ptr_object_->faces_end(); f_it++) {
		Matrix3f temp;
		int i = 0;
		for (MyMesh::FaceVertexIter fv_it = ptr_object_->fv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
			if (i >= 3) {
				std::cout << "Find a face which is not a triangle!" << std::endl;
				break;
			}
			auto vert = ptr_object_->point(fv_it.handle());
			temp(0, i) = vert[0];
			temp(1, i) = vert[1];
			temp(2, i) = vert[2];
			i++;
		}
		triangles.push_back(temp);
		canonical_triangles.push_back(Into_Canonical_Triangle(temp));
	}
}

struct xtradata {
	int min_idx;
	Cluster *c;
};

void cluster_distance(float *p, float *x, int m, int n, void *data)
{
	register int i;
	struct xtradata *dat;
	dat = (struct xtradata *)data;
	vector<Cluster::tNode>::iterator it_c = (dat->c)->clusterAssment.begin();
	int idx_c = 0;
	for (i = 0; i < n; ++i)
	{
		x[i] = 0;
	}
	while (it_c != (dat->c)->clusterAssment.end()) {
		if ((*it_c).min_Index == dat->min_idx) {
			float d = (dat->c)->distTriangles(Vector3f(p[0], p[1], p[2]), (dat->c)->triangles[idx_c]);
			for (i = 0; i < n; ++i)
			{
				x[i] += d;
			}
		}
		it_c++;
		idx_c++;
	}
	
}


void Cluster::LM_NonLinearSearch()
{
	int centroids_num = centroids_index.size();
	int m = 3, n = 3;     //m表示待求参数的维度，n表示测量值的维度
	float *p = new float[m];
	float *x = new float[n];
	float opts[LM_OPTS_SZ], info[LM_INFO_SZ];
	struct xtradata data;
	data.c = this;
	for (int k = 0; k < centroids_num; k++) {
		//待求参数的初值
		data.min_idx = k;
		p[0] = canonical_centroids_triangles[k](0);
		p[1] = canonical_centroids_triangles[k](1);
		p[2] = canonical_centroids_triangles[k](2);
		//模拟n次测量的结果，由于是最小值，设为零
		for (int i = 0; i < n; i++) x[i] = 0.0;
		opts[0] = LM_INIT_MU;
		// 调用迭代入口函数
		int ret = slevmar_dif(
			cluster_distance,      //描述测量值之间关系的函数指针
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

		canonical_centroids_triangles[k] = Vector3f(p);
	}
}

void Cluster::init()
{
	centroids_index.clear();
	canonical_centroids_triangles.clear();
	clusterAssment.clear();
	if(nonlinear_type_ == kDefault || nonlinear_type_ == kT) centroids_index.push_back(0);
	else if (nonlinear_type_ == kLM) canonical_centroids_triangles.push_back(canonical_triangles[0]);
	clusterAssment.push_back(tNode(0, 0));

	for (int i = 1; i < triangles.size(); i++) {
		float d12;
		if (nonlinear_type_ == kDefault) d12 = distTriangles(canonical_triangles[0], canonical_triangles[i]);
		else if (nonlinear_type_ == kT)  d12 = distTriangles(triangles[0], triangles[i]);
		else if (nonlinear_type_ == kLM) d12 = distTriangles(canonical_centroids_triangles[0], canonical_triangles[i]);;
		
		clusterAssment.push_back(tNode(0, d12));
	}
}

void Cluster::kmeans()
{
	if (triangles.size() == 0) return;
	init();
	updateCluster();
	updateCentroids();

	while (addCentroids()) {
		updateCluster();
		updateCentroids();
	}

	if (nonlinear_type_ == kDefault || nonlinear_type_ == kT) Set_Canonical_Triangle();
	Set_Rigid_Transformation();
	Set_Faces_Colors();
}


#define min3v(v1, v2, v3)   ((v1)>(v2)? ((v2)>(v3)?(v3):(v2)):((v1)>(v3)?(v3):(v2)))
#define max3v(v1, v2, v3)   ((v1)<(v2)? ((v2)<(v3)?(v3):(v2)):((v1)<(v3)?(v3):(v1)))
typedef struct
{
	int  red;              // [0,255]
	int  green;            // [0,255]
	int  blue;             // [0,255]
}COLOR_RGB;
typedef struct
{
	float hue;              // [0,360]
	float saturation;       // [0,100]
	float luminance;        // [0,100]
}COLOR_HSL;

// Converts RGB to HSL
static void RGBtoHSL(const COLOR_RGB *rgb, COLOR_HSL *hsl)
{
	float h = 0, s = 0, l = 0;
	// normalizes red-green-blue values
	float r = rgb->red / 255.0f;
	float g = rgb->green / 255.0f;
	float b = rgb->blue / 255.0f;
	float maxVal = max3v(r, g, b);
	float minVal = min3v(r, g, b);
	// hue
	if (maxVal == minVal)
	{
		h = 0; // undefined
	}
	else if (maxVal == r && g >= b)
	{
		h = 60.0f*(g - b) / (maxVal - minVal);
	}
	else if (maxVal == r && g < b)
	{
		h = 60.0f*(g - b) / (maxVal - minVal) + 360.0f;
	}
	else if (maxVal == g)
	{
		h = 60.0f*(b - r) / (maxVal - minVal) + 120.0f;
	}
	else if (maxVal == b)
	{
		h = 60.0f*(r - g) / (maxVal - minVal) + 240.0f;
	}
	// luminance
	l = (maxVal + minVal) / 2.0f;
	// saturation
	if (l == 0 || maxVal == minVal)
	{
		s = 0;
	}
	else if (0 < l && l <= 0.5f)
	{
		s = (maxVal - minVal) / (maxVal + minVal);
	}
	else if (l > 0.5f)
	{
		s = (maxVal - minVal) / (2 - (maxVal + minVal)); //(maxVal-minVal > 0)?
	}
	hsl->hue = (h > 360) ? 360 : ((h < 0) ? 0 : h);
	hsl->saturation = ((s > 1) ? 1 : ((s < 0) ? 0 : s)) * 100;
	hsl->luminance = ((l > 1) ? 1 : ((l < 0) ? 0 : l)) * 100;
}

// Converts HSL to RGB
static void HSLtoRGB(const COLOR_HSL *hsl, COLOR_RGB *rgb)
{
	float h = hsl->hue;                  // h must be [0, 360]
	float s = hsl->saturation / 100.f; // s must be [0, 1]
	float l = hsl->luminance / 100.f;      // l must be [0, 1]
	float R, G, B;
	if (hsl->saturation == 0)
	{
		// achromatic color (gray scale)
		R = G = B = l * 255.0f;
	}
	else
	{
		float q = (l < 0.5f) ? (l * (1.0f + s)) : (l + s - (l*s));
		float p = (2.0f * l) - q;
		float Hk = h / 360.0f;
		float T[3];
		T[0] = Hk + 0.3333333f; // Tr   0.3333333f=1.0/3.0
		T[1] = Hk;              // Tg
		T[2] = Hk - 0.3333333f; // Tb
		for (int i = 0; i < 3; i++)
		{
			if (T[i] < 0) T[i] += 1.0f;
			if (T[i] > 1) T[i] -= 1.0f;
			if ((T[i] * 6) < 1)
			{
				T[i] = p + ((q - p)*6.0f*T[i]);
			}
			else if ((T[i] * 2.0f) < 1) //(1.0/6.0)<=T[i] && T[i]<0.5
			{
				T[i] = q;
			}
			else if ((T[i] * 3.0f) < 2) // 0.5<=T[i] && T[i]<(2.0/3.0)
			{
				T[i] = p + (q - p) * ((2.0f / 3.0f) - T[i]) * 6.0f;
			}
			else T[i] = p;
		}
		R = T[0] * 255.0f;
		B = T[1] * 255.0f;
		G = T[2] * 255.0f;
	}
	rgb->red = (int)((R > 255) ? 255 : ((R < 0) ? 0 : R));
	rgb->green = (int)((G > 255) ? 255 : ((G < 0) ? 0 : G));
	rgb->blue = (int)((B > 255) ? 255 : ((B < 0) ? 0 : B));
}

void Cluster::Set_Faces_Colors()
{
	clusterNumber.clear();
	clusterNumber.reserve(class_num);
	vector< OpenMesh::Vec3uc > colors;
	ColorBar cb;
	cb.setRange(0, class_num);
	for (int i = 0; i < class_num; i++) {
		OpenMesh::Vec3uc color_i(cb.getColor(i + 0.5).red(), cb.getColor(i + 0.5).green(), cb.getColor(i + 0.5).blue());
		colors.push_back(color_i);
		clusterNumber.push_back(0);
	}
	
	int idx = 0;
	if (ptr_object_->has_face_colors()) {
		for (MyMesh::FaceIter f_it = ptr_object_->faces_begin(); f_it != ptr_object_->faces_end(); f_it++) {
			ptr_object_->set_color(f_it.handle(), colors[clusterAssment[idx].min_Index]);
			clusterNumber[clusterAssment[idx].min_Index]++;
			idx++;
		}
	}
	else {
		for (MyMesh::FaceIter f_it = ptr_object_->faces_begin(); f_it != ptr_object_->faces_end(); f_it++) {
			//ptr_object_->set_color(f_it.handle(), colors[clusterAssment[idx].min_Index]);
			clusterNumber[clusterAssment[idx].min_Index]++;
			idx++;
		}
	}

}

void Cluster::show_info()
{
	if ((nonlinear_type_ == kDefault || nonlinear_type_ == kT) && !centroids_index.empty()) {
		for (int i = 0; i < class_num; i++) {
			std::cout << centroids_index[i] << " :  " << canonical_triangles[centroids_index[i]].transpose() << std::endl;
			std::cout << "ClassNum: " << clusterNumber[i] << std::endl;
		}
	}
	else if (nonlinear_type_ == kLM) {
		for (int i = 0; i < class_num; i++) {
			std::cout << i + 1 << " :  " << canonical_centroids_triangles[i].transpose() << std::endl;
			std::cout << "ClassNum: " << clusterNumber[i] << std::endl;
		}
	}

	std::cout << std::endl;
}

void Cluster::Draw()
{
	glBegin(GL_LINES);
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(1, 1, 1);
	Vector3f J = Vector3f::Ones(3);
	//std::cout << this->rigid_transformations.size() << std::endl;
	for (int i = 0; i < this->rigid_transformations.size(); i++) {
		Matrix3f equivalant_triangle;
		equivalant_triangle = this->rigid_transformations[i].R * this->Into_Triangle(this->canonical_centroids_triangles[this->clusterAssment[i].min_Index]) * this->rigid_transformations[i].R_ + this->rigid_transformations[i].T * J.transpose();
		glVertex3f(equivalant_triangle(0, 0), equivalant_triangle(1, 0), equivalant_triangle(2, 0));
		glVertex3f(equivalant_triangle(0, 1), equivalant_triangle(1, 1), equivalant_triangle(2, 1));
		glVertex3f(equivalant_triangle(0, 1), equivalant_triangle(1, 1), equivalant_triangle(2, 1));
		glVertex3f(equivalant_triangle(0, 2), equivalant_triangle(1, 2), equivalant_triangle(2, 2));
		glVertex3f(equivalant_triangle(0, 2), equivalant_triangle(1, 2), equivalant_triangle(2, 2));
		glVertex3f(equivalant_triangle(0, 0), equivalant_triangle(1, 0), equivalant_triangle(2, 0));
	}
	glEnd();
	int idx = 0;
	glBegin(GL_TRIANGLES);
	for (MyMesh::FaceIter f_it = ptr_object_->faces_begin(); f_it != ptr_object_->faces_end(); ++f_it) {
		glClear(GL_COLOR_BUFFER_BIT);
		glColor3f((int(ptr_object_->color(f_it.handle()).data()[0]) - int('0')) / 255.0, (int(ptr_object_->color(f_it.handle()).data()[1]) - int('0')) / 255.0, (int(ptr_object_->color(f_it.handle()).data()[2]) - int('0')) / 255.0);
		Matrix3f equivalant_triangle;
		equivalant_triangle = this->rigid_transformations[idx].R * this->Into_Triangle(this->canonical_centroids_triangles[this->clusterAssment[idx].min_Index]) * this->rigid_transformations[idx].R_ + this->rigid_transformations[idx].T * J.transpose();
		glNormal3fv(ptr_object_->normal(*f_it).data());
		glVertex3f(equivalant_triangle(0, 0), equivalant_triangle(1, 0), equivalant_triangle(2, 0));
		glVertex3f(equivalant_triangle(0, 1), equivalant_triangle(1, 1), equivalant_triangle(2, 1));
		glVertex3f(equivalant_triangle(0, 2), equivalant_triangle(1, 2), equivalant_triangle(2, 2));
		idx++;
	}
	glEnd();
	glColor3f(1.0, 1.0, 1.0);
}

float Cluster::distTriangles(Matrix3f A, Matrix3f B)
{
	Matrix3f PA, PB, I, M, R, RB, U, V;
	Vector3f T, J, a_, b_;	// 3x1 float matrix.
	I = Matrix3f::Identity(3,3);
	J = Vector3f::Ones(3);
	float min_Dist = INT_MAX;

	for (int j = 0; j < 2; j++) {
		PA.col(0) = A.col((0 + j) % 2);
		PA.col(1) = A.col((1 + j) % 2);
		PA.col(2) = A.col(2);
		for (int i = 0; i < 3; i++) {
			PB.col(0) = B.col((0 + i) % 3);
			PB.col(1) = B.col((1 + i) % 3);
			PB.col(2) = B.col((2 + i) % 3);

			M = Matrix3f::Zero();

			a_ = PA * J / 3.0;
			b_ = PB * J / 3.0;

			for (int k = 0; k < 3; k++) {
				M += (PA.col(k) - a_) * (PB.col(k) - b_).transpose();
			}

			JacobiSVD<Eigen::MatrixXf> svd(M, ComputeThinU | ComputeThinV);
			U = svd.matrixU();
			V = svd.matrixV();

			R = U * V.transpose();
			if (R.determinant() < 0) {
				Matrix3f  S = U.inverse() * M * V.transpose().inverse();
				for (int k = 0; k < 3; k++) {
					if (S(k, k) < 0) {
						V.row(k) *= -1;
						break;
					}
				}
			}
			R = U * V.transpose();
			T = a_ - R * b_;

			float temp_d = 0;
			RB = R * PB;
			for (int k = 0; k < 3; k++) {
				temp_d += ((RB.col(k) + T - PA.col(k)).transpose()*(RB.col(k) + T - PA.col(k)))(0,0);
			}
			if (temp_d < min_Dist) {
				min_Dist = temp_d;
			}
		}
	}

	return min_Dist;
}

float Cluster::distTriangles(Vector3f A, Vector3f B)
{	
	float a1_a2 = A[0] - B[0];
	float b1_b2 = A[1] - B[1];
	float h1_h2 = A[2] - B[2];

	return a1_a2 * a1_a2 + b1_b2 * b1_b2 - a1_a2 * b1_b2 + h1_h2 * h1_h2;
}

float Cluster::distTriangles(Vector3f A, Matrix3f B)
{
	return distTriangles(A, Into_Canonical_Triangle(B));
}