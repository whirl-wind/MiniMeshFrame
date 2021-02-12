#include "renderingwidget.h"
#include <QKeyEvent>
#include <QColorDialog>
#include <QFileDialog>
#include <iostream>
#include <QtWidgets/QMenu>
#include <QtWidgets/QAction>
#include <QTextCodec>
#include <gl/GLU.h>
#include <gl/glut.h>
#include <algorithm>
#include "mainwindow.h"
#include "ArcBall.h"
#include "globalFunctions.h"
#include "MinimalSurface.h"

RenderingWidget::RenderingWidget(QWidget *parent, MainWindow* mainwindow)
	: QGLWidget(parent), ptr_mainwindow_(mainwindow), eye_distance_(5.0),
	has_lighting_(false), is_draw_point_(true), is_draw_edge_(false), is_draw_face_(false), is_draw_texture_(false)
{
	ptr_arcball_ = new CArcBall(width(), height());
	ptr_mesh_ = new MyMesh();
	ptr_backup_mesh_ = new MyMesh();
	ptr_origin_mesh_ = new MyMesh();
	ptr_target_mesh_ = new MyMesh();

	is_load_texture_ = false;
	is_draw_axes_ = false;
	is_draw_meancurvature_ = false;
	is_draw_absmeancurvature_ = false;
	is_draw_gaussiancurvature_ = false;
	is_draw_loss_ = false;
	is_draw_boundary_ = false;
	is_draw_interaction_ = false;
	is_set_fix_ = false;
	is_set_move_ = false;
	need_to_update_vectorfield_ = true;
	is_draw_vectorfield_ = false;
	is_rotate_ = false;
	mesh_is_clustered_ = false;
	is_draw_cluster_ = false;
	angle = 0;

	is_draw_plane_ = false;

	pa = new Partition();
	cb = new ColorBar(Qt::Horizontal, this);
	DDG = new DifferentialGeometry();
	DDG_backup_ = new DifferentialGeometry();
	DDG_target_ = new DifferentialGeometry();
	need_to_update_DDG = true;

	IXs.clear();
	select_ix = nullptr;

	eye_goal_[0] = eye_goal_[1] = eye_goal_[2] = 0.0;
	eye_direction_[0] = eye_direction_[1] = 0.0;
	eye_direction_[2] = 1.0;

	//m_timer = new QElapsedTimer;
	//m_timer->start();
}

RenderingWidget::~RenderingWidget()
{
	SafeDelete(ptr_arcball_);
	SafeDelete(ptr_mesh_);
	SafeDelete(ptr_backup_mesh_);
}

void RenderingWidget::initializeGL()
{
	glClearColor(0.2, 0.2, 0.2, 0.0);
	b_color = QColor(int(0.2 * 255), int(0.2 * 255), int(0.2 * 255));
	glShadeModel(GL_SMOOTH);

	glEnable(GL_DOUBLEBUFFER);
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_DEPTH_TEST);
	glClearDepth(1);

	SetLight();

}

void RenderingWidget::resizeGL(int w, int h)
{
	h = (h == 0) ? 1 : h;

	ptr_arcball_->reSetBound(w, h);

	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluPerspective(45.0, GLdouble(w) / GLdouble(h), 0.001, 1000);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void RenderingWidget::paintGL()
{
	glShadeModel(GL_SMOOTH);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	if (has_lighting_)
	{
		SetLight();
	}
	else
	{
		glDisable(GL_LIGHTING);
		glDisable(GL_LIGHT0);
	}

	DrawSelection(is_set_fix_ || is_set_move_);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	register vec eyepos = eye_distance_ * eye_direction_;
	gluLookAt(eyepos[0], eyepos[1], eyepos[2],
		eye_goal_[0], eye_goal_[1], eye_goal_[2],
		0.0, 1.0, 0.0);
	glPushMatrix();

	glRotatef(angle, 0, 1, 0);
	if (is_rotate_) {
		angle = (angle + 1) % 360;
	}
	glMultMatrixf(ptr_arcball_->GetBallMatrix());

	Render();
	glPopMatrix();
	//std::cout << 1 << std::endl;
}

void RenderingWidget::timerEvent(QTimerEvent * e)
{
	updateGL();
}

void RenderingWidget::unify(MyMesh *ptr_mesh_,float size)
{
	if (ptr_mesh_->n_vertices() < 2)return;
	float xmax_, xmin_, ymax_, ymin_, zmax_, zmin_;
	MyMesh::VertexIter v_it = ptr_mesh_->vertices_begin();
	MyMesh::Point vert = ptr_mesh_->point(v_it.handle());
	xmax_ = xmin_ = vert[0];
	ymax_ = ymin_ = vert[1];
	zmax_ = zmin_ = vert[2];
	for (v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); v_it++) {
		vert = ptr_mesh_->point(v_it.handle());
		if(vert[0] > xmax_) xmax_ = vert[0];
		else if(vert[0] < xmin_) xmin_ = vert[0];
		if (vert[1] > ymax_) ymax_ = vert[1];
		else if (vert[1] < ymin_) ymin_ = vert[1];
		if (vert[2] > zmax_) zmax_ = vert[2];
		else if (vert[2] < zmin_) zmin_ = vert[2];
	}
	
	float scaleX = xmax_ - xmin_;
	float scaleY = ymax_ - ymin_;
	float scaleZ = zmax_ - zmin_;
	float scaleMax;

	if (scaleX < scaleY)
	{
		scaleMax = scaleY;
	}
	else
	{
		scaleMax = scaleX;
	}
	if (scaleMax < scaleZ)
	{
		scaleMax = scaleZ;
	}
	float scaleV = size / scaleMax;
	MyMesh::Point centerPos((xmin_ + xmax_) / 2.f, (ymin_ + ymax_) / 2.f, (zmin_ + zmax_) / 2.f);
	for (v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); v_it++) 
	{
		vert = ptr_mesh_->point(v_it.handle());
		ptr_mesh_->set_point(*v_it, (vert - centerPos) * scaleV);
	}
}

void RenderingWidget::mousePressEvent(QMouseEvent *e)
{
	switch (e->button())
	{
	case Qt::LeftButton:
		ptr_arcball_->MouseDown(e->pos());
		break;
	case Qt::MidButton:
		current_position_ = e->pos();
		break;
	case Qt::RightButton:
		if (is_set_fix_ || is_set_move_) {
			start_point_ = end_point_ = e->pos();
		}
		else if (is_draw_interaction_) {
			MyMesh::Point near_point, far_point;
			std::tie(near_point, far_point) = GetSelectionRay(e->x(), e->y());
			select_ix = GetSelectedInteraction(near_point, far_point);
			current_position_ = e->pos();
			if (select_ix) {
				select_ix->update_selection();
				MyMesh::Point d = near_point + ((select_ix->center - near_point) | (far_point - near_point)) * (far_point - near_point) / (far_point - near_point).sqrnorm();
				select_ix->handle_point = d;
				df->update_IXs(IXs, select_ix);
				df->Init();
			}
		}
		break;
	default:
		break;
	}

	updateGL();
}
void RenderingWidget::mouseMoveEvent(QMouseEvent *e)
{
	switch (e->buttons())
	{
		setCursor(Qt::ClosedHandCursor);
	case Qt::LeftButton:
		ptr_arcball_->MouseMove(e->pos());
		break;
	case Qt::MidButton:
		eye_goal_[0] -= 4.0*GLfloat(e->x() - current_position_.x()) / GLfloat(width());
		eye_goal_[1] += 4.0*GLfloat(e->y() - current_position_.y()) / GLfloat(height());
		current_position_ = e->pos();
		break;
	case Qt::RightButton:
		if (is_set_fix_ || is_set_move_) {
			end_point_ = e->pos();
		}
		else if (is_draw_interaction_) {
			if (select_ix) {
				MyMesh::Point near_point, far_point;
				std::tie(near_point, far_point) = GetSelectionRay(e->x(), e->y());
				MyMesh::Point d = near_point + ((select_ix->center - near_point) | (far_point - near_point)) * (far_point - near_point) / (far_point - near_point).sqrnorm();
				//if (need_to_update_DDG) UpdateDDG();
				select_ix->Move(d - select_ix->handle_point);
				select_ix->handle_point = d;
				for (int i = 0; i < IXs.size(); i++) {
					if (IXs[i] != select_ix) IXs[i]->update_center();
				}
				df->deform();
				need_to_update_DDG = true;
				ptr_mesh_->update_normals();
			}
		}
		break;
	default:
		break;
	}

	updateGL();
}
void RenderingWidget::mouseDoubleClickEvent(QMouseEvent *e)
{
	switch (e->button())
	{
	case Qt::LeftButton:
		BackUp();
		break;
	default:
		break;
	}
	updateGL();
}
void RenderingWidget::mouseReleaseEvent(QMouseEvent *e)
{
	switch (e->button())
	{
	case Qt::LeftButton:
		ptr_arcball_->MouseUp(e->pos());
		setCursor(Qt::ArrowCursor);
		break;

	case Qt::RightButton:
		if (is_set_fix_) {
			select_ix = SetIX(start_point_, end_point_, true);
			IXs.push_back(select_ix);
			is_set_fix_ = false;
			is_draw_interaction_ = true;
		}
		else if (is_set_move_) {
			select_ix = SetIX(start_point_, end_point_, false);
			IXs.push_back(select_ix);
			is_set_move_ = false;
			is_draw_interaction_ = true;
		}
		else if (is_draw_interaction_) {
			if (select_ix) {
				MyMesh::Point near_point, far_point;
				std::tie(near_point, far_point) = GetSelectionRay(e->x(), e->y());
				MyMesh::Point d = near_point + ((select_ix->center - near_point) | (far_point - near_point)) * (far_point - near_point) / (far_point - near_point).sqrnorm();
				//if (need_to_update_DDG) UpdateDDG();
				select_ix->Move(d - select_ix->handle_point);
				select_ix->handle_point = d;
				for (int i = 0; i < IXs.size(); i++) {
					if (IXs[i] != select_ix) IXs[i]->update_center();
				}
				df->deform();
				need_to_update_DDG = true;
				ptr_mesh_->update_normals();
			}
		}
		break;
	default:
		break;
	}

	updateGL();
}

void RenderingWidget::wheelEvent(QWheelEvent *e)
{
	eye_distance_ += e->delta()*0.001;
	eye_distance_ = eye_distance_ < 0 ? 0 : eye_distance_;

	updateGL();
}

void RenderingWidget::keyPressEvent(QKeyEvent *e)
{
	switch (e->key())
	{
	case Qt::Key_A:
		break;
	default:
		break;
	}
}

void RenderingWidget::keyReleaseEvent(QKeyEvent *e)
{
	switch (e->key())
	{
	case Qt::Key_A:
		break;
	default:
		break;
	}
}

void RenderingWidget::Render()
{
	DrawAxes(is_draw_axes_);
	DrawInteraction(is_draw_interaction_);

	DrawCluster(is_draw_cluster_);

	DrawLoss(is_draw_loss_);
	DrawAbsMeanCurvature(is_draw_absmeancurvature_);
	DrawMeanCurvature(is_draw_meancurvature_);
	DrawGaussianCurvature(is_draw_gaussiancurvature_);

	DrawPlane(is_draw_plane_);

	DrawPoints(is_draw_point_);
	DrawBoundary(is_draw_boundary_);
	DrawVectorField(is_draw_vectorfield_);
	DrawEdge(is_draw_edge_);
	DrawFace(is_draw_face_);
	DrawTexture(is_draw_texture_);

	if (!(is_draw_absmeancurvature_ || is_draw_meancurvature_ || is_draw_gaussiancurvature_ || is_draw_loss_)) CloseColorBar();
}

void RenderingWidget::SetLight()
{
	static GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
	static GLfloat mat_shininess[] = { 50.0 };
	static GLfloat light_position[] = { 1.0, 1.0, 1.0, 0.0 };
	static GLfloat white_light[] = { 0.8, 0.8, 0.8, 1.0 };
	static GLfloat lmodel_ambient[] = { 0.3, 0.3, 0.3, 1.0 };

	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, white_light);
	glLightfv(GL_LIGHT0, GL_SPECULAR, white_light);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);

	glColorMaterial(GL_FRONT, GL_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);
}

void RenderingWidget::ClearFaceColor()
{
	for (MyMesh::FaceIter f_it = ptr_mesh_->faces_begin(); f_it != ptr_mesh_->faces_end(); ++f_it) {
		ptr_mesh_->set_color(f_it.handle(), OpenMesh::Vec3uc(225,225,225));
	}
	updateGL();
}


void RenderingWidget::Restore()
{
	ptr_mesh_->assign(*ptr_origin_mesh_);
	ptr_mesh_->update_normals();
	ptr_mesh_->request_vertex_status();
	ptr_mesh_->request_edge_status();
	ptr_mesh_->request_halfedge_status();
	ptr_mesh_->request_face_status();


	ClearFaceColor();
	if(need_to_update_DDG) UpdateDDG();
	select_ix = nullptr;
	IXs.clear();
	IXs.push_back(new Interaction(ptr_mesh_));
	BackUp();
	updateGL();
}

void RenderingWidget::SetBackground()
{
	QColor color = QColorDialog::getColor(Qt::white, this, tr("background color"));
	b_color = color;
	GLfloat r = (color.red()) / 255.0f;
	GLfloat g = (color.green()) / 255.0f;
	GLfloat b = (color.blue()) / 255.0f;
	GLfloat alpha = color.alpha() / 255.0f;
	glClearColor(r, g, b, alpha);
	updateGL();
}

void RenderingWidget::ReadMesh()
{
	QString filename = QFileDialog::
		getOpenFileName(this, tr("Read Mesh"),
			"D:\\Zz\\CG\\Materials", tr("Meshes (*.obj;*.off)"));

	if (filename.isEmpty())
	{
		emit(operatorInfo(QString("Read Mesh Failed!")));
		return;
	}
	//中文路径支持
	QTextCodec *code = QTextCodec::codecForName("gd18030");
	QTextCodec::setCodecForLocale(code);

	QByteArray byfilename = filename.toLocal8Bit();
	if (!OpenMesh::IO::read_mesh(*ptr_mesh_, byfilename.data()))
	{
		std::cerr << "Cannot Open mesh to file " << byfilename.data() << std::endl;
		return;
	}
	unify(ptr_mesh_, 2.0);
	ptr_mesh_->request_face_normals();
	ptr_mesh_->request_vertex_normals();
	ptr_mesh_->request_halfedge_normals();
	ptr_mesh_->update_normals();

	ptr_mesh_->request_vertex_status();
	ptr_mesh_->request_edge_status();
	ptr_mesh_->request_halfedge_status();
	ptr_mesh_->request_face_status();
	ptr_mesh_->request_face_colors();

	ptr_origin_mesh_->assign(*ptr_mesh_);
	UpdateDDG();
	BackUp();
	ClearFaceColor();

	select_ix = nullptr;
	IXs.clear();
	IXs.push_back(new Interaction(ptr_mesh_));

	df = new Df_ARAP(ptr_mesh_,IXs, select_ix, DDG_backup_);

	emit(operatorInfo(QString("Read Mesh from") + filename + QString(" Done")));
	emit(meshInfo(ptr_mesh_->n_vertices(), ptr_mesh_->n_edges(), ptr_mesh_->n_faces()));
	std::cout << "\nMesh: " << std::endl;
	std::cout << "vertices: " << ptr_mesh_->n_vertices() << std::endl;
	std::cout << "edges: " << ptr_mesh_->n_edges() << std::endl;
	std::cout << "faces: " << ptr_mesh_->n_faces() << std::endl;
	std::cout << std::endl;

	is_draw_vectorfield_ = false;
	need_to_update_vectorfield_ = true;
	is_draw_plane_ = false;
	mesh_is_clustered_ = false;
	is_draw_cluster_ = false;

	updateGL();
}

void RenderingWidget::MorphingTo()
{
	QString filename = QFileDialog::
		getOpenFileName(this, tr("Read Mesh"),
			"D:\\Zz\\CG\\Materials", tr("Meshes (*.obj;*.off)"));

	if (filename.isEmpty())
	{
		emit(operatorInfo(QString("Read Mesh Failed!")));
		return;
	}
	//中文路径支持
	QTextCodec* code = QTextCodec::codecForName("gd18030");
	QTextCodec::setCodecForLocale(code);

	QByteArray byfilename = filename.toLocal8Bit();
	if (!OpenMesh::IO::read_mesh(*ptr_target_mesh_, byfilename.data()))
	{
		std::cerr << "Cannot Open mesh to file " << byfilename.data() << std::endl;
		return;
	}
	unify(ptr_target_mesh_, 2.0);
	ptr_target_mesh_->request_face_normals();
	ptr_target_mesh_->request_vertex_normals();
	ptr_target_mesh_->request_halfedge_normals();
	ptr_target_mesh_->update_normals();

	DDG_target_->update_mesh(ptr_target_mesh_);

	emit(operatorInfo(QString("Read Target Mesh from") + filename + QString(" Done")));

	if (ptr_target_mesh_->n_vertices() != ptr_mesh_->n_vertices() || ptr_target_mesh_->n_edges() != ptr_mesh_->n_edges() || ptr_target_mesh_->n_faces() != ptr_mesh_->n_faces()) {
		std::cout << "Target OBJ dosen't Match!" << std::endl;
		std::cout << std::endl;
		return;
	}

	Mp_ARAP mp(ptr_mesh_, ptr_target_mesh_, DDG_target_);
	float dt = mp.t;
	int frame = int(1 / dt);
	std::cout << "Begin Morphing: " << std::endl;
	for (int i = 0; i < frame; i++) {
		float time = 1 / float(frame - i);
		//float time = i * dt;
		mp.morphing(time, i);
		UpdateDDG();
		updateGL();
	}
	UpdateDDG();
	std::cout << "Done" << std::endl;
}

void RenderingWidget::WriteMesh()
{
	if (ptr_mesh_->n_vertices() == 0)
	{
		emit(QString("The Mesh is Empty !"));
		return;
	}
	QString filename = QFileDialog::
		getSaveFileName(this, tr("Write Mesh"),
			"..", tr("Meshes (*.obj)"));

	if (filename.isEmpty())
		return;
	if (!OpenMesh::IO::write_mesh(*ptr_mesh_, filename.toLatin1().data()))
	{
		std::cerr << "Cannot Write mesh to file "<< filename.toLatin1().data() << std::endl;
		return;
	}

	emit(operatorInfo(QString("Write Mesh to ") + filename + QString(" Done")));
}

void RenderingWidget::LoadTexture()
{
	QString filename = QFileDialog::getOpenFileName(this, tr("Load Texture"),
		"D:\\Zz\\CG\\Materials", tr("Images(*.bmp *.jpg *.png *.jpeg)"));
	if (filename.isEmpty())
	{
		emit(operatorInfo(QString("Load Texture Failed!")));
		return;
	}

	texture_filename_ = filename;

	glGenTextures(1, &texture_[0]);
	QImage tex1, buf;
	if (!buf.load(filename))
	{
		//        QMessageBox::warning(this, tr("Load Fialed!"), tr("Cannot Load Image %1").arg(filenames.at(0)));
		emit(operatorInfo(QString("Load Texture Failed!")));
		return;
		/*
		QImage dummy(128, 128, QImage::Format_ARGB32);
		dummy.fill(Qt::green);
		buf = dummy;
		*/
	}
	texture_img_ = cv::imread(filename.toLatin1().data());
	cv::cvtColor(texture_img_, texture_img_, cv::COLOR_BGR2RGB);

	tex1 = QGLWidget::convertToGLFormat(buf);
	glBindTexture(GL_TEXTURE_2D, texture_[0]);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_NEAREST);
	gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGB, tex1.width(), tex1.height(),
		GL_RGBA, GL_UNSIGNED_BYTE, tex1.bits());

	is_load_texture_ = true;
	ptr_mainwindow_->action_pt_window->trigger();
	emit(operatorInfo(QString("Load Texture from ") + filename + QString(" Done")));
}

void RenderingWidget::UpdateDDG()
{
	ptr_mesh_->update_normals();
	DDG->update_mesh(ptr_mesh_);
	need_to_update_DDG = false;
}

void RenderingWidget::BackUp()
{
	//backup
	ptr_backup_mesh_->assign(*ptr_mesh_);
	DDG_backup_->update_mesh(ptr_backup_mesh_);
}

void RenderingWidget::CheckDrawPoint(bool bv)
{
	is_draw_point_ = bv;
	updateGL();
}
void RenderingWidget::CheckDrawEdge(bool bv)
{
	is_draw_edge_ = bv;
	updateGL();
}
void RenderingWidget::CheckDrawFace(bool bv)
{
	is_draw_face_ = bv;
	updateGL();
}
void RenderingWidget::CheckDrawBoundary(bool bv)
{
	is_draw_boundary_ = bv;
	updateGL();
}
void RenderingWidget::CheckLight(bool bv)
{
	has_lighting_ = bv;
	updateGL();
}
void RenderingWidget::CheckDrawTexture(bool bv)
{
	is_draw_texture_ = bv;
	if (is_draw_texture_)
		glEnable(GL_TEXTURE_2D);
	else
		glDisable(GL_TEXTURE_2D);

	updateGL();
}
void RenderingWidget::CheckDrawAxes(bool bv)
{
	is_draw_axes_ = bv;
	updateGL();
}

void RenderingWidget::CheckRotate(bool bv)
{
	is_rotate_ = bv;
	updateGL();
	if (is_rotate_) m_nTimerId = this->startTimer(50);
	else this->killTimer(m_nTimerId);
}

void RenderingWidget::CheckDrawInteraction(bool bv)
{
	is_draw_interaction_ = bv;
	updateGL();
}

void RenderingWidget::CheckDrawVectorField(bool bv)
{
	is_draw_vectorfield_ = bv;
	updateGL();
}

void RenderingWidget::CheckDrawLoss(bool bv)
{
	is_draw_loss_ = bv;
	updateGL();
}

void RenderingWidget::CheckDrawMeanCurvature(bool bv)
{
	is_draw_meancurvature_ = bv;
	updateGL();
}

void RenderingWidget::CheckDrawAbsMeanCurvature(bool bv)
{
	is_draw_absmeancurvature_ = bv;
	updateGL();
}

void RenderingWidget::CheckDrawGaussianCurvature(bool bv)
{
	is_draw_gaussiancurvature_ = bv;
	updateGL();
}

void RenderingWidget::CheckDrawPlane(bool bv)
{
	is_draw_plane_ = bv;
	if (is_draw_plane_) {
		if (ptr_mesh_ == NULL) return;
		if (ptr_mesh_->n_faces() == 0) return;
		if (need_to_update_DDG) UpdateDDG();
		int n;
		std::cout << "Please input the num of vector: (its just a reference, not the final classnum!)" << std::endl;
		std::cin >> n;
		if (n <= 0) return;
		if (n > ptr_mesh_->n_faces()) n = ptr_mesh_->n_faces();
		if(n!=pa->get_classnum()) pa = new Pa_VSA(ptr_mesh_, DDG, n);
	}
	updateGL();
}

void RenderingWidget::CheckDrawCluster(bool bv)
{
	is_draw_cluster_ = bv;
	if (is_draw_cluster_&& !mesh_is_clustered_) {
		ClusterMesh();
	}
	updateGL();
}

void RenderingWidget::MinimalSurface_Local()
{
	MinimalSurface ms(ptr_mesh_);
	ms.local_methods();
	need_to_update_DDG = true;
	updateGL();
}

void RenderingWidget::MinimalSurface_Global()
{
	MinimalSurface ms(ptr_mesh_);
	ms.global_methods();
	need_to_update_DDG = true;
	updateGL();
}

void RenderingWidget::AddNoise()
{
	if (need_to_update_DDG) UpdateDDG();

	double rate = 0.03;
	Smoothing* soomthing = new Smoothing();
	soomthing->LoadMesh(ptr_mesh_);
	double scale = ((DDG->bbox_x_max - DDG->bbox_x_min) + (DDG->bbox_y_max - DDG->bbox_y_min) + (DDG->bbox_z_max - DDG->bbox_z_min)) / 3.0 * rate;
	soomthing->AddNoise(scale);
	need_to_update_DDG = true;
	updateGL();
}

void RenderingWidget::GlobalSmooth()
{
	if (need_to_update_DDG) UpdateDDG();

	BilateralMeshDenoising* soomthing = new BilateralMeshDenoising();
	soomthing->LoadMesh(ptr_mesh_);
	soomthing->LoadDDG(DDG);
	soomthing->SetGlobal(true);
	soomthing->Set_GlobalParameter(0.3, 0.5, 0.4);
	soomthing->SmoothMesh();
	need_to_update_DDG = true;
	updateGL();
}

void RenderingWidget::LocalSmooth()
{
	if (need_to_update_DDG) UpdateDDG();

	BilateralMeshDenoising* soomthing = new BilateralMeshDenoising();
	soomthing->LoadMesh(ptr_mesh_);
	soomthing->LoadDDG(DDG);
	soomthing->SetGlobal(false);
	soomthing->Set_LocalParameter(20);
	soomthing->SmoothMesh();
	need_to_update_DDG = true;
	updateGL();
}

/*
void RenderingWidget::Simplification_QEM()
{
	if (simplification->Qualified())
	{
		simplification->simplify_one_step();
		updateGL();
	}

	std::cout << "n: " << ptr_mesh_->n_vertices()<< std::endl;
	std::cout << "cost: " << simplification->update_min_cost() << std::endl;
	std::cout << "qsize: " << simplification->q.size() << std::endl;
	//UpdateDDG();
}
*/
void RenderingWidget::Simplification_QEM()
{
	if (need_to_update_DDG) UpdateDDG();

	Sp_QEM* simplification = new Sp_QEM(ptr_mesh_, ptr_mesh_->n_vertices()/2, 10.0);
	int ccc = ptr_mesh_->n_faces();
	while (simplification->Qualified())
	{
		simplification->simplify_one_step();
		//ccc = 0;
		//for (MyMesh::FaceIter f_it = ptr_mesh_->faces_begin(); f_it != ptr_mesh_->faces_end(); ++f_it) {
		//	if (ptr_mesh_->status(*f_it).deleted()) continue;
		//	ccc++;
		//}
		updateGL();
	}

	ptr_mesh_->garbage_collection();
	UpdateDDG();
	updateGL();

	std::cout << "Sp_Mesh: " << std::endl;
	std::cout << "vertices: " << ptr_mesh_->n_vertices() << std::endl;
	std::cout << "edges: " << ptr_mesh_->n_edges() << std::endl;
	std::cout << "faces: " << ptr_mesh_->n_faces() << std::endl;
	std::cout << "sp_min_cost: " << simplification->update_min_cost() << std::endl;
	std::cout << std::endl;
}

void RenderingWidget::BuildNPolyVectorField()
{
	if (ptr_mesh_ == NULL) return;
	if (ptr_mesh_->n_faces() == 0) return;

	int n = 1, m = 4;
	std::cout << "Please input n: " << std::endl;
	std::cin >> n;
	std::cout << "Please input m: " << std::endl;
	std::cin >> m;
	vf = new NpolyVectorFields(ptr_mesh_, n, m);
	vf->BuildVectorField();
	//vf->SetColor();
	need_to_update_vectorfield_ = false;
	updateGL();
}

std::tuple<MyMesh::Point, MyMesh::Point> RenderingWidget::GetSelectionRay(int mouse_x, int mouse_y)
{
	// 获取 Model-View、Projection 矩阵 & 获取Viewport视区
	GLdouble modelview[16];
	GLdouble projection[16];
	GLint viewport[4];
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	register vec eyepos = eye_distance_ * eye_direction_;
	gluLookAt(eyepos[0], eyepos[1], eyepos[2],
		eye_goal_[0], eye_goal_[1], eye_goal_[2],
		0.0, 1.0, 0.0);
	glPushMatrix();
	glMultMatrixf(ptr_arcball_->GetBallMatrix());
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	glPopMatrix();
	glGetDoublev(GL_PROJECTION_MATRIX, projection);
	glGetIntegerv(GL_VIEWPORT, viewport); 
	GLfloat winX, winY, winZ;
	winX = (float)mouse_x;
	winY = (float)viewport[3] - (float)mouse_y;
	glReadPixels(int(winX), int(winY), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ);

	GLdouble world_x, world_y, world_z; 
	// 获取近裁剪面上的交点
	gluUnProject(winX, winY, 0.0,
		modelview, projection, viewport,
		&world_x, &world_y, &world_z);
	MyMesh::Point near_point(world_x, world_y, world_z); 
	// 获取远裁剪面上的交点
	gluUnProject(winX, winY, 1.0,
		modelview, projection, viewport,
		&world_x, &world_y, &world_z);
	MyMesh::Point far_point(world_x, world_y, world_z);

	return { near_point, far_point };
}

QPoint RenderingWidget::GetViewPoint(MyMesh::Point)
{
	return QPoint();
}

Interaction* RenderingWidget::GetSelectedInteraction(MyMesh::Point near_point, MyMesh::Point far_point)
{
	Interaction* selected;
	selected = nullptr;
	float dist = 1.0;
	for (int i = 0; i < IXs.size(); i++) {
		auto nor = ((IXs[i]->center - near_point) | (far_point - near_point))* (far_point - near_point)/ (far_point - near_point).sqrnorm();
		MyMesh::Point d = near_point + nor;
		float new_dist = (d - IXs[i]->center).norm();
		if (new_dist < dist) {
			dist = new_dist;
			selected = IXs[i];
			selected->handle_point = d;
		}
	}
	return selected;
}

Interaction* RenderingWidget::SetIX(QPoint start_point_, QPoint end_point_, bool fix)
{
	double x0, x1, y0, y1;
	if (start_point_.x() < end_point_.x()) {
		x0 = start_point_.x();
		x1 = end_point_.x();
	}
	else {
		x1 = start_point_.x();
		x0 = end_point_.x();
	}
	if (start_point_.y() < end_point_.y()) {
		y0 = start_point_.y();
		y1 = end_point_.y();
	}
	else {
		y1 = start_point_.y();
		y0 = end_point_.y();
	}

	// 获取 Model-View、Projection 矩阵 & 获取Viewport视区
	GLdouble modelview[16];
	GLdouble projection[16];
	GLint viewport[4];
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	register vec eyepos = eye_distance_ * eye_direction_;
	gluLookAt(eyepos[0], eyepos[1], eyepos[2],
		eye_goal_[0], eye_goal_[1], eye_goal_[2],
		0.0, 1.0, 0.0);
	glPushMatrix();
	glMultMatrixf(ptr_arcball_->GetBallMatrix());
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	glPopMatrix();
	glGetDoublev(GL_PROJECTION_MATRIX, projection);
	glGetIntegerv(GL_VIEWPORT, viewport);
	GLdouble winX, winY, winZ;
	GLdouble world_x, world_y, world_z;
	vector<MyMesh::VertexHandle> point_array;
	for (MyMesh::VertexIter v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); ++v_it) {
		world_x = ptr_mesh_->point(*v_it)[0];
		world_y = ptr_mesh_->point(*v_it)[1];
		world_z = ptr_mesh_->point(*v_it)[2];

		gluProject(world_x, world_y, world_z,
			modelview, projection, viewport,
			&winX, &winY, &winZ);

		if (winX > x0 && winX<x1 && (float)viewport[3] - winY > y0 && (float)viewport[3] - winY < y1)
			point_array.push_back(*v_it);
	}
	
	return new Interaction(ptr_mesh_,point_array,fix);
}

void RenderingWidget::Set_Fix()
{
	is_draw_interaction_ = true;
	is_set_fix_ = true;
	is_set_move_ = false;
	start_point_ = end_point_ = QPoint(0, 0);
}

void RenderingWidget::Set_Move()
{
	is_draw_interaction_ = true;
	is_set_move_ = true;
	is_set_fix_ = false;
	start_point_ = end_point_ = QPoint(0, 0);
}

void RenderingWidget::DeleteIX()
{
	if (IXs.size() > 0) {
		vector<Interaction*>::iterator it = std::find(IXs.begin(), IXs.end(), select_ix);
		if (it != IXs.end())
		{
			int idx = std::distance(IXs.begin(), it);
			Interaction* p = IXs[idx];
			IXs[idx] = nullptr;
			delete p;
			IXs.erase(IXs.begin() + idx);
			//IXs.pop_back();
		}
	}
	select_ix = nullptr;
	updateGL();
}

void RenderingWidget::UpdateRemeshing()
{
	UpdateDDG();
	IncrementalRemeshing rm(ptr_mesh_, DDG);
	rm.remesh(DDG->average_edgelength);
	UpdateDDG();
	updateGL();
}

void RenderingWidget::Triangulate()
{
	DelaunayTriangulation Tri(ptr_mesh_);
	Tri.triangulate();
	/*
	bool s = false;
	for (MyMesh::FaceIter f_it = ptr_mesh_->faces_begin(); f_it != ptr_mesh_->faces_end(); ++f_it) {
		for (MyMesh::FaceEdgeIter fe_it = ptr_mesh_->fe_begin(*f_it); fe_it != ptr_mesh_->fe_end(*f_it); fe_it++) {
			if (!ptr_mesh_->is_boundary(fe_it)) {
				auto he = ptr_mesh_->halfedge_handle(*fe_it, 0);
				MyMesh::Point v1 = ptr_mesh_->point(ptr_mesh_->from_vertex_handle(he));
				MyMesh::Point v2 = ptr_mesh_->point(ptr_mesh_->to_vertex_handle(he));
				MyMesh::Point v0 = ptr_mesh_->point(ptr_mesh_->to_vertex_handle(ptr_mesh_->next_halfedge_handle(he)));
				MyMesh::Point v3 = ptr_mesh_->point(ptr_mesh_->to_vertex_handle(ptr_mesh_->next_halfedge_handle(ptr_mesh_->opposite_halfedge_handle(he))));

				auto e11 = (v1 - v0).normalized();
				auto e12 = (v2 - v0).normalized();
				auto e21 = (v1 - v3).normalized();
				auto e22 = (v2 - v3).normalized();

				float cosa = e11[0] * e12[0] + e11[1] * e12[1];
				float sina = abs(e11[0] * e12[1] - e11[1] * e12[0]);
				float cosb = e21[0] * e22[0] + e21[1] * e22[1];
				float sinb = abs(e22[0] * e21[1] - e22[1] * e21[0]);
				
				if (cosa * sinb + sina * cosb < 0) {
					if (ptr_mesh_->is_flip_ok(fe_it)) {
						s = true;
						ptr_mesh_->flip(fe_it);
						ptr_mesh_->garbage_collection();
						break;
					}
				}
			}
		}
		if (s) break;
	}
	*/

	std::cout << "New Mesh: " << std::endl;
	std::cout << "vertices: " << ptr_mesh_->n_vertices() << std::endl;
	std::cout << "edges: " << ptr_mesh_->n_edges() << std::endl;
	std::cout << "faces: " << ptr_mesh_->n_faces() << std::endl;
	std::cout << std::endl;

	//UpdateDDG();
	updateGL();
}

void RenderingWidget::ClusterMesh()
{
	if (ptr_mesh_->n_faces() == 0) return;

	mesh_is_clustered_ = true;
	is_draw_cluster_ = true;
	int class_num = 1;
	std::cout << "Please input the class number: " << std::endl;
	std::cin >> class_num;
	cluster = new Cluster(ptr_mesh_, class_num);
	cluster->show_info();

	updateGL();
}

void RenderingWidget::ModifyMesh()
{
	if (!mesh_is_clustered_) {
		if (ptr_mesh_->n_faces() == 0) return;
		int class_num = 1;
		std::cout << "Please input the class number: " << std::endl;
		std::cin >> class_num;
		cluster = new Cluster(ptr_mesh_, class_num);
		mesh_is_clustered_ = true;
	}

	if (cluster->triangles.size() == 0) return;

	is_draw_cluster_ = true;

	float opt_err = 10.0;
	std::cout << "Please input the tolerance error: " << std::endl;
	std::cin >> opt_err;

	glo = new GlobalLinearOptimization(cluster,opt_err);
	
	///////////////
	std::cout << "Modifying... " << std::endl;
	int MAX_C = 2000;
	int MIN_C = 100;
	float error = 0.075;

	float lo = glo->get_loss();
	while (lo > glo->t_err) {
		std::cout << glo->get_cluster()->class_num << ") loss: " << lo << std::endl;
		glo->get_cluster()->addClass();
		glo->get_cluster()->updateInfo();
		lo = glo->get_loss();
		updateGL();
	}
	std::cout << glo->get_cluster()->class_num << ") loss: " << lo << std::endl;
	std::cout << "Cluster_OPT Finished!\n" << std::endl;

	vector<float> lo_mem;
	lo = INT_MAX;
	int count = 0;
	lo_mem.clear();
	while (lo > error && count < MAX_C) {

		lo = glo->One_Step();
		//glo->get_cluster()->show_info();
		updateGL();

		count++;
		std::cout << count << ") loss: " << lo << std::endl;

		lo_mem.push_back(lo);
		if (lo_mem.size() >= 10) {
			int count_lo = 0;
			for (int k = 0; k < lo_mem.size() - 1; k++) {
				if (lo_mem[k] <= lo) count_lo++;
				lo_mem[k] = lo_mem[k + 1];
			}
			if (count_lo > 8 && count > MIN_C) break;
			lo_mem.pop_back();
		}
	}

	glo->loss = lo;
	////////////////

	cluster = glo->get_cluster();
	cluster->updateTriangles();
	cluster->updateCentroids();
	cluster->updateInfo();

	std::cout << "------------------------------------------------- " << std::endl;
	glo->get_cluster()->show_info();
	std::cout << "Cluster Num: " << glo->get_cluster()->class_num << std::endl;
	std::cout << "Total Loss: " << glo->loss << std::endl;
	std::cout << "------------------------------------------------- \n" << std::endl;

	UpdateDDG();
	updateGL();
}

void RenderingWidget::DrawAxes(bool bV)
{
	if (!bV)
		return;
	//x axis
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	glVertex3f(0, 0, 0);
	glVertex3f(0.7, 0.0, 0.0);
	glEnd();
	glPushMatrix();
	glTranslatef(0.7, 0, 0);
	glRotatef(90, 0.0, 1.0, 0.0);
	//glutSolidCone(0.02, 0.06, 20, 10);
	glPopMatrix();

	//y axis
	glColor3f(0.0, 1.0, 0.0);
	glBegin(GL_LINES);
	glVertex3f(0, 0, 0);
	glVertex3f(0.0, 0.7, 0.0);
	glEnd();

	glPushMatrix();
	glTranslatef(0.0, 0.7, 0);
	glRotatef(90, -1.0, 0.0, 0.0);
	//glutSolidCone(0.02, 0.06, 20, 10);
	glPopMatrix();

	//z axis
	glColor3f(0.0, 0.0, 1.0);
	glBegin(GL_LINES);
	glVertex3f(0, 0, 0);
	glVertex3f(0.0, 0.0, 0.7);
	glEnd();
	glPushMatrix();
	glTranslatef(0.0, 0, 0.7);
	//glutSolidCone(0.02, 0.06, 20, 10);
	glPopMatrix();

	glColor3f(1.0, 1.0, 1.0);
}

void RenderingWidget::DrawPoints(bool bv)
{
	if (!bv || ptr_mesh_ == NULL) return;
	if (ptr_mesh_->n_vertices() == 0) return;
	if (!ptr_mesh_->has_vertex_normals()) {
		std::cout << "ERROR: Standard vertex property 'Normals' not available!" << std::endl;
	}
	
	glBegin(GL_POINTS);
	for (MyMesh::VertexIter v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); ++v_it) {
		if (ptr_mesh_->status(*v_it).deleted()) continue;
		glNormal3fv(ptr_mesh_->normal(*v_it).data());
		glVertex3fv(ptr_mesh_->point(*v_it).data());
	}
	glEnd();
}

void RenderingWidget::DrawEdge(bool bv)
{
	if (!bv || ptr_mesh_ == NULL) return;
	if (ptr_mesh_->n_edges() == 0) return;

	glBegin(GL_LINES);
	for (MyMesh::HalfedgeIter he_it = ptr_mesh_->halfedges_begin(); he_it != ptr_mesh_->halfedges_end(); ++he_it) {
		//glNormal3fv(ptr_mesh_->normal(ptr_mesh_->from_vertex_handle(*he_it)).data());
		if (ptr_mesh_->status(*he_it).deleted()) continue;
		glNormal3fv(ptr_mesh_->normal(*he_it).data());
		glVertex3fv(ptr_mesh_->point(ptr_mesh_->from_vertex_handle(*he_it)).data());
		glVertex3fv(ptr_mesh_->point(ptr_mesh_->to_vertex_handle(*he_it)).data());
	}
	glEnd();
}

void RenderingWidget::DrawFace(bool bv)
{
	if (!bv || ptr_mesh_ == NULL) return;
	if (ptr_mesh_->n_faces() == 0) return;
	if (!ptr_mesh_->has_vertex_normals()) {
		std::cout << "ERROR: Standard vertex property 'Normals' not available!" << std::endl;
	}
	glBegin(GL_TRIANGLES);
	for (MyMesh::FaceIter f_it = ptr_mesh_->faces_begin(); f_it != ptr_mesh_->faces_end(); ++f_it) {
			
		glClear(GL_COLOR_BUFFER_BIT);
		glColor3f((int(ptr_mesh_->color(f_it.handle()).data()[0]) - int('0')) / 255.0, (int(ptr_mesh_->color(f_it.handle()).data()[1]) - int('0')) / 255.0, (int(ptr_mesh_->color(f_it.handle()).data()[2]) - int('0')) / 255.0);
		for (MyMesh::FaceVertexIter fv_it = ptr_mesh_->fv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
			if (ptr_mesh_->status(*fv_it).deleted()) continue;
			glNormal3fv(ptr_mesh_->normal(*fv_it).data());
			glVertex3fv(ptr_mesh_->point(*fv_it).data());
		}			
	}
	glEnd();
	glColor3f(1.0, 1.0, 1.0);
}

void RenderingWidget::DrawBoundary(bool bv)
{
	if (!bv || ptr_mesh_ == NULL) return;
	if (ptr_mesh_->n_edges() == 0) return;
	if (need_to_update_DDG) UpdateDDG();

	glBegin(GL_LINES);
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(1.0, 0.0, 0.0);
	for (int i = 0; i < DDG->Boundaries.size(); i++) {
		for (int j = 0; j < DDG->Boundaries[i].size(); j++) {
			glNormal3fv(ptr_mesh_->normal(ptr_mesh_->opposite_halfedge_handle(DDG->Boundaries[i][j])).data());
			glVertex3fv(ptr_mesh_->point(ptr_mesh_->from_vertex_handle(DDG->Boundaries[i][j])).data());
			glVertex3fv(ptr_mesh_->point(ptr_mesh_->to_vertex_handle(DDG->Boundaries[i][j])).data());
		}
	}
	glEnd();

	glColor3f(1.0, 1.0, 1.0);
}

void RenderingWidget::DrawTexture(bool bv)
{
	if (!bv) return;
	if (ptr_mesh_->n_faces() == 0 || !is_load_texture_) return;
	
	//默认使用球面纹理映射，效果不好
	if (!ptr_mesh_->has_vertex_texcoords2D()) {
		std::cout << "ERROR: Standard vertex property 'Texcoords2D' not available!" << std::endl;
		ptr_mesh_->request_vertex_texcoords2D();
	}
	
	//glBindTexture(GL_TEXTURE_2D, texture_[0]);
	glBegin(GL_TRIANGLES);
	for (MyMesh::FaceIter f_it = ptr_mesh_->faces_begin(); f_it != ptr_mesh_->faces_end(); ++f_it)
	{
		for (auto fh_it = ptr_mesh_->fh_begin(f_it.handle()); fh_it != ptr_mesh_->fh_end(f_it.handle()); ++fh_it) {
			MyMesh::HalfedgeHandle pedge = fh_it.handle();
			glTexCoord2fv(ptr_mesh_->texcoord2D(ptr_mesh_->to_vertex_handle(pedge)).data());
			glNormal3fv(ptr_mesh_->normal(ptr_mesh_->to_vertex_handle(pedge)).data());
			glVertex3fv(ptr_mesh_->point(ptr_mesh_->to_vertex_handle(pedge)).data());
		}
	}
	glEnd();
}

void RenderingWidget::DrawInteraction(bool bv)
{
	if (!bv || ptr_mesh_ == NULL) return;
	if (ptr_mesh_->n_vertices() == 0) return;

	if (select_ix) select_ix->Draw_Select();

	for (int i = 0; i < IXs.size(); i++) {
		IXs[i]->Draw();
	}
}

void RenderingWidget::DrawSelection(bool bv)
{
	if (!bv) return;

	MyMesh::Point near_point_0, far_point_0;
	std::tie(near_point_0, far_point_0) = GetSelectionRay(start_point_.x(), start_point_.y());
	MyMesh::Point near_point_1, far_point_1;
	std::tie(near_point_1, far_point_1) = GetSelectionRay(end_point_.x(), end_point_.y());
	MyMesh::Point near_point_2, far_point_2;
	std::tie(near_point_2, far_point_2) = GetSelectionRay(start_point_.x(), end_point_.y());
	MyMesh::Point near_point_3, far_point_3;
	std::tie(near_point_3, far_point_3) = GetSelectionRay(end_point_.x(), start_point_.y());

	register vec eyepos = eye_distance_ * eye_direction_;
	gluLookAt(eyepos[0], eyepos[1], eyepos[2],
		eye_goal_[0], eye_goal_[1], eye_goal_[2],
		0.0, 1.0, 0.0);
	glPushMatrix();
	glMultMatrixf(ptr_arcball_->GetBallMatrix());

	MyMesh::Point a = ((near_point_0 + far_point_0) / 2.0);
	MyMesh::Point b = ((near_point_2 + far_point_2) / 2.0);
	MyMesh::Point c = ((near_point_1 + far_point_1) / 2.0);
	MyMesh::Point d = ((near_point_3 + far_point_3) / 2.0);

	if (is_set_fix_) glColor3f(1.0, 0.5, 0.0);
	else glColor3f(0.5, 1.0, 0.0);
	glPointSize(3);
	glBegin(GL_POINTS);
	glVertex3fv(a.data());
	glVertex3fv(b.data());
	glVertex3fv(c.data());
	glVertex3fv(d.data());
	glEnd();
	glPointSize(1);

	glBegin(GL_LINES);
	glVertex3fv(a.data());
	glVertex3fv(b.data());
	glVertex3fv(b.data());
	glVertex3fv(c.data());
	glVertex3fv(c.data());
	glVertex3fv(d.data());
	glVertex3fv(d.data());
	glVertex3fv(a.data());
	glEnd();

	glPopMatrix();

	glColor3f(1.0, 1.0, 1.0);
}

void RenderingWidget::DrawVectorField(bool bv)
{
	if (!bv || ptr_mesh_ == NULL) return;
	if (ptr_mesh_->n_faces() == 0) return;
	if (DDG->Boundaries.size() > 0) return;

	if (need_to_update_vectorfield_) BuildNPolyVectorField();
	
	vf->Draw(DDG->average_edgelength/3.0);
}

void RenderingWidget::DrawPlane(bool bv)
{
	if (!bv || pa->get_classnum() <= 0) return;
	pa->Draw();
}

void RenderingWidget::DrawCluster(bool bv)
{
	if (!bv || ptr_mesh_ == NULL) return;
	cluster->Draw();
}


void RenderingWidget::DrawLoss(bool bv)
{
	if (!bv || ptr_mesh_ == NULL) return;
	if (ptr_mesh_->n_faces() == 0) return;
	if (ptr_mesh_->n_vertices() != ptr_backup_mesh_->n_vertices())return;
	vector<double> loss;

	loss.resize(ptr_mesh_->n_vertices());
	for (MyMesh::VertexIter v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); ++v_it) {
		loss[(*v_it).idx()] = (ptr_mesh_->point(*v_it) - ptr_backup_mesh_->point(ptr_backup_mesh_->vertex_handle((*v_it).idx()))).length();
	}

	//vector<double>::iterator max = max_element(loss.begin(), loss.end());
	cb->setRange(0.0, 0.2);

	ShowColorBar();
	for (MyMesh::FaceIter f_it = ptr_mesh_->faces_begin(); f_it != ptr_mesh_->faces_end(); ++f_it) {
		glBegin(GL_TRIANGLES);
		glShadeModel(GL_SMOOTH);
		for (MyMesh::FaceVertexIter fv_it = ptr_mesh_->fv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
			QColor color = cb->getColor(loss[(*fv_it).idx()]);
			glColor3f(color.redF(), color.greenF(), color.blueF());
			glNormal3fv(ptr_mesh_->normal(*fv_it).data());
			glVertex3fv(ptr_mesh_->point(*fv_it).data());
		}
		glEnd();
	}
	glColor3f(1.0, 1.0, 1.0);
}


void RenderingWidget::DrawMeanCurvature(bool bv)
{
	if (!bv || ptr_mesh_ == NULL) return;
	if (ptr_mesh_->n_faces() == 0) return;
	if (need_to_update_DDG) UpdateDDG();

	DDG->Caculate_MeanCurvature();
	vector<double> mc;
	DDG->Get_MeanCurvature(mc);
	vector<double>::iterator max = max_element(mc.begin(), mc.end());
	vector<double>::iterator min = min_element(mc.begin(), mc.end());
	cb->setRange(*min, *max);

	ShowColorBar();
	for (MyMesh::FaceIter f_it = ptr_mesh_->faces_begin(); f_it != ptr_mesh_->faces_end(); ++f_it) {
		glBegin(GL_TRIANGLES);
		glShadeModel(GL_SMOOTH);
		for (MyMesh::FaceVertexIter fv_it = ptr_mesh_->fv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
			QColor color = cb->getColor(mc[(*fv_it).idx()]);
			glColor3f(color.redF(), color.greenF(), color.blueF());
			glNormal3fv(ptr_mesh_->normal(*fv_it).data());
			glVertex3fv(ptr_mesh_->point(*fv_it).data());
		}
		glEnd();
	}
	glColor3f(1.0, 1.0, 1.0);
}

void RenderingWidget::DrawAbsMeanCurvature(bool bv)
{
	if (!bv || ptr_mesh_ == NULL) return;
	if (ptr_mesh_->n_faces() == 0) return;
	if (need_to_update_DDG) UpdateDDG();

	DDG->Caculate_MeanCurvature();
	vector<double> mc;
	DDG->Get_MeanCurvature(mc);
	vector<double>::iterator max = max_element(mc.begin(), mc.end());
	vector<double>::iterator min = min_element(mc.begin(), mc.end());
	if (*max <= 0) {
		cb->setRange(-(*max), -(*min));
	}
	else if (*min <= 0) {
		cb->setRange(0, *max);
	}
	else cb->setRange(*min, *max);

	ShowColorBar();
	for (MyMesh::FaceIter f_it = ptr_mesh_->faces_begin(); f_it != ptr_mesh_->faces_end(); ++f_it) {
		glBegin(GL_TRIANGLES);
		glShadeModel(GL_SMOOTH);
		for (MyMesh::FaceVertexIter fv_it = ptr_mesh_->fv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
			QColor color = cb->getColor(abs(mc[(*fv_it).idx()]));
			glColor3f(color.redF(), color.greenF(), color.blueF());
			glNormal3fv(ptr_mesh_->normal(*fv_it).data());
			glVertex3fv(ptr_mesh_->point(*fv_it).data());
		}
		glEnd();
	}
	glColor3f(1.0, 1.0, 1.0);
}

void RenderingWidget::DrawGaussianCurvature(bool bv)
{
	if (!bv || ptr_mesh_ == NULL) return;
	if (ptr_mesh_->n_faces() == 0) return;
	if (need_to_update_DDG) UpdateDDG();

	DDG->Caculate_GaussianCurvature();
	vector<double> gc;
	DDG->Get_GaussianCurvature(gc);
	vector<double>::iterator max = max_element(gc.begin(), gc.end());
	vector<double>::iterator min = min_element(gc.begin(), gc.end());
	cb->setRange(*min, *max);

	ShowColorBar();
	for (MyMesh::FaceIter f_it = ptr_mesh_->faces_begin(); f_it != ptr_mesh_->faces_end(); ++f_it) {
		glBegin(GL_TRIANGLES);
		glShadeModel(GL_SMOOTH);
		for (MyMesh::FaceVertexIter fv_it = ptr_mesh_->fv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
			QColor color = cb->getColor(gc[(*fv_it).idx()]);
			glColor3f(color.redF(), color.greenF(), color.blueF());
			glNormal3fv(ptr_mesh_->normal(*fv_it).data());
			glVertex3fv(ptr_mesh_->point(*fv_it).data());
		}
		glEnd();
	}
	glColor3f(1.0, 1.0, 1.0);
}

void RenderingWidget::ShowColorBar()
{
	int w = this->width();
	cb->move(w - 100, 50);
	cb->setbackground(b_color);
	cb->resize(30 + 20 + 5, 343 + 50 + 5);
	cb->show();
}

void RenderingWidget::CloseColorBar()
{
	cb->close();
}

