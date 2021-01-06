#pragma once
#ifndef RENDERINGWIDGET_H
#define RENDERINGWIDGET_H

#include <QGLWidget>
#include <QEvent>
#include <QElapsedTimer>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Handles.hh>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/core/core.hpp>
#include "mainwindow.h"
#include "imagewindow.h"
#include "HE_mesh/Vec.h"
#include "Interaction.h"
#include "Deformation.h"
#include "DifferentialGeometry.h"
#include "colorbar.h"
#include "Smoothing.h"
#include "Morphing.h"
#include "Simplification.h"
#include "VectorField.h"
#include "Remeshing.h"
#include "Triangulation.h"
#include "Partition.h"
#include "Cluster.h"
#include "GlobalLinearOptimization.h"

using trimesh::vec;
using trimesh::point;

class MainWindow;
class CArcBall;

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;


class RenderingWidget : public QGLWidget
{
	Q_OBJECT

public:
	RenderingWidget(QWidget *parent, MainWindow* mainwindow = 0);
	~RenderingWidget();

protected:
	void initializeGL();
	void resizeGL(int w, int h);
	void paintGL();
	void timerEvent(QTimerEvent *e);
	void unify(MyMesh* ptr_mesh_, float size); //调节mesh适应窗口
	
	// mouse events
	void mousePressEvent(QMouseEvent *e);
	void mouseMoveEvent(QMouseEvent *e);
	void mouseReleaseEvent(QMouseEvent *e);
	void mouseDoubleClickEvent(QMouseEvent *e);
	void wheelEvent(QWheelEvent *e);

public:
	void keyPressEvent(QKeyEvent *e);
	void keyReleaseEvent(QKeyEvent *e);

signals:
	void meshInfo(int, int, int);
	void operatorInfo(QString);

private:
	void Render();
	void SetLight();

public slots:
	void Restore();
	void SetBackground();
	void ReadMesh();
	void WriteMesh();
	void LoadTexture();
	void UpdateDDG();
	void BackUp();
	void MorphingTo();

	void CheckDrawPoint(bool bv);
	void CheckDrawEdge(bool bv);
	void CheckDrawFace(bool bv);
	void CheckLight(bool bv);
	void CheckDrawTexture(bool bv);
	void CheckDrawAxes(bool bv);
	void CheckRotate(bool bv);

	void CheckDrawBoundary(bool bv);
	void CheckDrawInteraction(bool bv);
	void CheckDrawVectorField(bool bv);

	void CheckDrawLoss(bool bv);
	void CheckDrawMeanCurvature(bool bv);
	void CheckDrawAbsMeanCurvature(bool bv);
	void CheckDrawGaussianCurvature(bool bv);

	void CheckDrawPlane(bool bv);
	void CheckDrawCluster(bool bv);

	//GeometryProcess
	void MinimalSurface_Local();
	void MinimalSurface_Global();
	void AddNoise();
	void GlobalSmooth();
	void LocalSmooth();
	void Simplification_QEM();

	//VectorField
	void BuildNPolyVectorField();

	//Interaction
	std::tuple<MyMesh::Point, MyMesh::Point> GetSelectionRay(int mouse_x, int mouse_y);
	QPoint GetViewPoint(MyMesh::Point);
	Interaction*  GetSelectedInteraction(MyMesh::Point, MyMesh::Point);
	Interaction* SetIX(QPoint start_point_, QPoint end_point_, bool fix);
	void Set_Fix();
	void Set_Move();
	void DeleteIX();

	//Remeshing
	void UpdateRemeshing();

	//Triangulation
	void Triangulate();

	//ClusterMesh
	void ClusterMesh();
	void ModifyMesh();

private:
	void DrawAxes(bool);
	void DrawPoints(bool);
	void DrawEdge(bool);
	void DrawFace(bool);
	void DrawBoundary(bool);
	void DrawTexture(bool);
	void DrawInteraction(bool);
	void DrawSelection(bool);
	void DrawVectorField(bool);
	void DrawPlane(bool);
	void DrawCluster(bool);

	void DrawLoss(bool);

	void DrawMeanCurvature(bool);
	void DrawAbsMeanCurvature(bool);
	void DrawGaussianCurvature(bool);
	void ShowColorBar();
	void CloseColorBar();

public:
	MainWindow					*ptr_mainwindow_;
	CArcBall					*ptr_arcball_;
	MyMesh						*ptr_mesh_;

	// Texture
	cv::Mat						texture_img_;
	GLuint						texture_[1];
	bool						is_load_texture_;
	QString						texture_filename_;

	// eye
	GLfloat						eye_distance_;
	point						eye_goal_;
	vec							eye_direction_;
	QPoint						current_position_;

	// Render information
	bool						is_draw_point_;
	bool						is_draw_edge_;
	bool						is_draw_face_;
	bool						is_draw_texture_;
	bool						has_lighting_;
	bool						is_draw_axes_;
	bool						is_rotate_;

	bool						is_draw_loss_;
	bool						is_draw_meancurvature_;
	bool						is_draw_absmeancurvature_;
	bool						is_draw_gaussiancurvature_;

	bool						is_draw_boundary_;
	bool						is_draw_interaction_;
	bool						is_set_fix_;
	bool						is_set_move_;
	bool						is_draw_vectorfield_;
	bool						need_to_update_vectorfield_;
	bool						is_draw_controlmesh_;

	bool						is_draw_plane_;

	bool						mesh_is_clustered_;
	bool						is_draw_cluster_;

	DifferentialGeometry		*DDG;
private:
	int							m_nTimerId;
	int							angle;
	QPoint start_point_;
	QPoint end_point_;

	MyMesh						*ptr_target_mesh_;
	DifferentialGeometry		*DDG_target_;

	MyMesh						*ptr_origin_mesh_;
	MyMesh						*ptr_backup_mesh_;
	DifferentialGeometry		*DDG_backup_;
	ColorBar					*cb;
	QColor						b_color;
	bool						need_to_update_DDG;

	Deformation					*df;
	vector<Interaction*>		IXs;
	Interaction					*select_ix;
	VectorField					*vf;
	Partition					*pa;

	GlobalLinearOptimization*	glo;
	Cluster*					cluster;
};

#endif // RENDERINGWIDGET_H
