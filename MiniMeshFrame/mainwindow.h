#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtWidgets/QMainWindow>
#include "ui_mainwindow.h"

class QLabel;
class QPushButton;
class QCheckBox;
class QGroupBox;
class RenderingWidget;
class ImageWindow;

class MainWindow : public QMainWindow
{
	Q_OBJECT

public:
	MainWindow(QWidget *parent = 0);
	~MainWindow();

	RenderingWidget					*ptr_renderingwidget_;
	ImageWindow						*ptr_imagewindow_;

	//Parameterization
	QAction							*action_pt_window;

private:
	void CreateActions();
	void CreateMenus();
	void CreateToolBars();
	void CreateStatusBar();
	void CreateRenderGroup();
	void ImageWidget_Show();

protected:
	void keyPressEvent(QKeyEvent *e);
	void keyReleaseEvent(QKeyEvent *e);

public slots:
	void ShowMeshInfo(int npoint, int nedge, int nface);
	void OpenFile();
	void ShowAbout();

private:
	Ui::MainWindowClass ui;

	// Basic
	QMenu							*menu_file_;
	QMenu							*menu_edit_;
	QMenu							*menu_help_;
	QToolBar						*toolbar_file_;
	QToolBar						*toolbar_edit_;
	QToolBar						*toolbar_basic_;
	QAction							*action_new_;
	QAction							*action_open_;
	QAction							*action_save_;
	QAction							*action_saveas_;

	QAction							*action_aboutqt_;
	QAction							*action_about_;

	// Basic Operator Tool
	QAction							*action_loadmesh_;
	QAction							*action_loadtexture_;
	QAction							*action_background_;
	QAction							*action_restore_;

	// Render RadioButtons
	QCheckBox						*checkbox_point_;
	QCheckBox						*checkbox_edge_;
	QCheckBox						*checkbox_face_;
	QCheckBox						*checkbox_light_;
	QCheckBox						*checkbox_texture_;
	QCheckBox						*checkbox_axes_;
	QCheckBox						*checkbox_rotate_;

	QGroupBox						*groupbox_render_;

	// Differential Geometry
	QCheckBox						*checkbox_loss_;
	QCheckBox						*checkbox_meancurvature_;
	QCheckBox						*checkbox_absolute_meancurvature_;
	QCheckBox						*checkbox_gaussiancurvature_;

	QGroupBox						*groupbox_chart_;

	QCheckBox						*checkbox_boundary_;
	QCheckBox						*checkbox_interaction_;


	QToolBar						*toolbar_vectorfield_;
	QAction							*action_vectorfield_npoly_;
	QCheckBox						*checkbox_vectorfield_npoly_;

	QGroupBox						*groupbox_geometry_;

	// Information
	QLabel							*label_meshinfo_;
	QLabel							*label_operatorinfo_;

	//GeometryProcess
	QToolBar						*toolbar_geometryprocess_;
	QAction							*action_ms_globalmethods_;
	QAction							*action_smoothing_local_;
	QAction							*action_noising_;
	QAction							*action_sp_qem_;

	//Interaction
	QToolBar						*toolbar_interaction_;
	QAction							*action_fix_;
	QAction							*action_move_;
	QAction							*action_delete_;

	//Parameterization
	QToolBar						*toolbar_parameterization_;

	//Morphing
	QToolBar						*toolbar_morphing_;
	QAction							*action_morphingto_;

	//Remeshing
	QToolBar						*toolbar_remeshing_;
	QAction							*action_remeshing_;

	//Triangulation
	QToolBar						*toolbar_triangulation_;
	QAction							*action_triangulation_;

	QGroupBox*						groupbox_else_;

	//Partition
	QCheckBox						*checkbox_plane_;

	//ClusterMesh
	QCheckBox*						checkbox_cluster_;
	QToolBar*						toolbar_clustermesh_;
	QAction*						action_clustermesh;
	QAction*						action_modifymesh;
};

#endif // MAINWINDOW_H
