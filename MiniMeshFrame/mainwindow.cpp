#include "mainwindow.h"

#include <QtWidgets/QLabel>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QMessageBox>
#include <QKeyEvent>
#include "renderingwidget.h"
#include "imagewindow.h"
#include "imagewidget.h"

MainWindow::MainWindow(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);

	ptr_renderingwidget_ = new RenderingWidget(this,this);
	//	setCentralWidget(renderingwidget_);
	ptr_imagewindow_ = nullptr;
	setGeometry(300, 150, 800, 600);

	CreateActions();
	CreateMenus();
	CreateToolBars();
	CreateStatusBar();
	CreateRenderGroup();

//////////////////////////////////////////////////////////////////////////////////////////////////////////

	action_ms_globalmethods_ = new QAction(tr("MinimalSurface"), this);
	connect(action_ms_globalmethods_, SIGNAL(triggered()), ptr_renderingwidget_, SLOT(MinimalSurface_Global()));

	action_noising_ = new QAction(tr("AddNoise"), this);
	connect(action_noising_, SIGNAL(triggered()), ptr_renderingwidget_, SLOT(AddNoise()));

	action_smoothing_local_ = new QAction(tr("Smooth"), this);
	connect(action_smoothing_local_, SIGNAL(triggered()), ptr_renderingwidget_, SLOT(LocalSmooth()));

	action_sp_qem_ = new QAction(tr("Simplification"), this);
	connect(action_sp_qem_, SIGNAL(triggered()), ptr_renderingwidget_, SLOT(Simplification_QEM()));

	this->addToolBarBreak(Qt::TopToolBarArea);
	toolbar_geometryprocess_ = addToolBar(tr("GeometryProcess"));
	toolbar_geometryprocess_->addAction(action_ms_globalmethods_);
	toolbar_geometryprocess_->addAction(action_noising_);
	toolbar_geometryprocess_->addAction(action_smoothing_local_);
	toolbar_geometryprocess_->addAction(action_sp_qem_);

	action_fix_ = new QAction(tr("SetFix"), this);
	connect(action_fix_, SIGNAL(triggered()), ptr_renderingwidget_, SLOT(Set_Fix()));

	action_move_ = new QAction(tr("SetMove"), this);
	connect(action_move_, SIGNAL(triggered()), ptr_renderingwidget_, SLOT(Set_Move()));

	action_delete_ = new QAction(tr("DeleteIX"), this);
	connect(action_delete_, SIGNAL(triggered()), ptr_renderingwidget_, SLOT(DeleteIX()));

	toolbar_interaction_ = addToolBar(tr("Interaction"));
	toolbar_interaction_->addAction(action_fix_);
	toolbar_interaction_->addAction(action_move_);
	toolbar_interaction_->addAction(action_delete_);

//////////////////////////////////////////////////////////////////////////////////////////////////////////

	QVBoxLayout *layout_left = new QVBoxLayout;
	layout_left->addWidget(groupbox_render_);
	layout_left->addWidget(groupbox_geometry_);
	layout_left->addWidget(groupbox_chart_);
	layout_left->addWidget(groupbox_else_);
	layout_left->addStretch(4);

	QHBoxLayout *layout_main = new QHBoxLayout;

	layout_main->addLayout(layout_left);
	layout_main->addWidget(ptr_renderingwidget_);
	layout_main->setStretch(1, 1);
	this->centralWidget()->setLayout(layout_main);

	ui.mainToolBar->setVisible(false);
	toolbar_file_->setVisible(false);
}

MainWindow::~MainWindow()
{

}

void MainWindow::CreateActions()
{
	action_new_ = new QAction(QIcon(":/MiniMeshFrame/Resources/images/new.png"), tr("&New"), this);
	action_new_->setShortcut(QKeySequence::New);
	action_new_->setStatusTip(tr("Create a new file"));

	action_open_ = new QAction(QIcon(":/MiniMeshFrame/Resources/images/open.png"), tr("&Open..."), this);
	action_open_->setShortcuts(QKeySequence::Open);
	action_open_->setStatusTip(tr("Open an existing file"));
	connect(action_open_, SIGNAL(triggered()), ptr_renderingwidget_, SLOT(ReadMesh()));

	action_save_ = new QAction(QIcon(":/MiniMeshFrame/Resources/images/save.png"), tr("&Save"), this);
	action_save_->setShortcuts(QKeySequence::Save);
	action_save_->setStatusTip(tr("Save the document to disk"));
	connect(action_save_, SIGNAL(triggered()), ptr_renderingwidget_, SLOT(WriteMesh()));

	action_saveas_ = new QAction(tr("Save &As..."), this);
	action_saveas_->setShortcuts(QKeySequence::SaveAs);
	action_saveas_->setStatusTip(tr("Save the document under a new name"));
	//	connect(action_saveas_, SIGNAL(triggered()), imagewidget_, SLOT(SaveAs()));

	action_loadmesh_ = new QAction(tr("ReadOBJ"), this);
	action_loadtexture_ = new QAction(tr("LoadTexture"), this);
	action_background_ = new QAction(tr("ChangeBackground"), this);
	action_restore_ = new QAction(tr("Restore"), this);

	connect(action_loadmesh_, SIGNAL(triggered()), ptr_renderingwidget_, SLOT(ReadMesh()));
	connect(action_loadtexture_, SIGNAL(triggered()), ptr_renderingwidget_, SLOT(LoadTexture()));
	connect(action_background_, SIGNAL(triggered()), ptr_renderingwidget_, SLOT(SetBackground()));
	connect(action_restore_, SIGNAL(triggered()), ptr_renderingwidget_, SLOT(Restore()));

	action_pt_window = new QAction(tr("pt_window"), this);
	connect(action_pt_window, &QAction::triggered, this, &MainWindow::ImageWidget_Show);

	action_morphingto_ = new QAction(tr("morphing_to"), this);
	connect(action_morphingto_, SIGNAL(triggered()), ptr_renderingwidget_, SLOT(MorphingTo()));

	action_vectorfield_npoly_ = new QAction(tr("update_vectorfield"), this);
	connect(action_vectorfield_npoly_, SIGNAL(triggered()), ptr_renderingwidget_, SLOT(BuildNPolyVectorField()));

	action_remeshing_ = new QAction(tr("remeshing"), this);
	connect(action_remeshing_, SIGNAL(triggered()), ptr_renderingwidget_, SLOT(UpdateRemeshing()));

	action_triangulation_ = new QAction(tr("triangulation"), this);
	connect(action_triangulation_, SIGNAL(triggered()), ptr_renderingwidget_, SLOT(Triangulate()));

	action_clustermesh = new QAction(tr("cluster_mesh"), this);
	connect(action_clustermesh, SIGNAL(triggered()), ptr_renderingwidget_, SLOT(ClusterMesh()));
	action_modifymesh = new QAction(tr("modify_mesh"), this);
	connect(action_modifymesh, SIGNAL(triggered()), ptr_renderingwidget_, SLOT(ModifyMesh()));
}

void MainWindow::CreateMenus()
{
	menu_file_ = menuBar()->addMenu(tr("&File"));
	menu_file_->setStatusTip(tr("File menu"));
	menu_file_->addAction(action_new_);
	menu_file_->addAction(action_open_);
	menu_file_->addAction(action_save_);
	menu_file_->addAction(action_saveas_);
}

void MainWindow::CreateToolBars()
{
	toolbar_file_ = addToolBar(tr("File"));
	toolbar_file_->addAction(action_new_);
	toolbar_file_->addAction(action_open_);
	toolbar_file_->addAction(action_save_);

	toolbar_basic_ = addToolBar(tr("Basic"));
	toolbar_basic_->addAction(action_loadmesh_);
	toolbar_basic_->addAction(action_loadtexture_);
	toolbar_basic_->addAction(action_background_);
	toolbar_basic_->addAction(action_restore_);

	toolbar_parameterization_ = addToolBar(tr("toolbar_parameterization_"));
	toolbar_parameterization_->addAction(action_pt_window);

	toolbar_morphing_ = addToolBar(tr("toolbar_morphing_"));
	toolbar_morphing_->addAction(action_morphingto_);

	toolbar_vectorfield_ = addToolBar(tr("toolbar_vectorvield_"));
	toolbar_vectorfield_->addAction(action_vectorfield_npoly_);

	toolbar_remeshing_ = addToolBar(tr("toolbar_remeshing_"));
	toolbar_remeshing_->addAction(action_remeshing_);

	toolbar_triangulation_ = addToolBar(tr("toolbar_triangulation_"));
	toolbar_triangulation_->addAction(action_triangulation_);

	toolbar_clustermesh_ = addToolBar(tr("toolbar_clustermesh_"));
	toolbar_clustermesh_->addAction(action_clustermesh);
	toolbar_clustermesh_->addAction(action_modifymesh);
}

void MainWindow::CreateStatusBar()
{
	label_meshinfo_ = new QLabel(QString("MeshInfo: p: %1 e: %2 f: %3").arg(0).arg(0).arg(0));
	label_meshinfo_->setAlignment(Qt::AlignCenter);
	label_meshinfo_->setMinimumSize(label_meshinfo_->sizeHint());

	label_operatorinfo_ = new QLabel();
	label_operatorinfo_->setAlignment(Qt::AlignVCenter);


	statusBar()->addWidget(label_meshinfo_);
	connect(ptr_renderingwidget_, SIGNAL(meshInfo(int, int, int)), this, SLOT(ShowMeshInfo(int, int, int)));

	statusBar()->addWidget(label_operatorinfo_);
	connect(ptr_renderingwidget_, SIGNAL(operatorInfo(QString)), label_operatorinfo_, SLOT(setText(QString)));
}

void MainWindow::CreateRenderGroup()
{
	checkbox_point_ = new QCheckBox(tr("Point"), this);
	connect(checkbox_point_, SIGNAL(clicked(bool)), ptr_renderingwidget_, SLOT(CheckDrawPoint(bool)));
	checkbox_point_->setChecked(true);

	checkbox_edge_ = new QCheckBox(tr("Edge"), this);
	connect(checkbox_edge_, SIGNAL(clicked(bool)), ptr_renderingwidget_, SLOT(CheckDrawEdge(bool)));

	checkbox_face_ = new QCheckBox(tr("Face"), this);
	connect(checkbox_face_, SIGNAL(clicked(bool)), ptr_renderingwidget_, SLOT(CheckDrawFace(bool)));

	checkbox_light_ = new QCheckBox(tr("Light"), this);
	connect(checkbox_light_, SIGNAL(clicked(bool)), ptr_renderingwidget_, SLOT(CheckLight(bool)));

	checkbox_texture_ = new QCheckBox(tr("Texture"), this);
	connect(checkbox_texture_, SIGNAL(clicked(bool)), ptr_renderingwidget_, SLOT(CheckDrawTexture(bool)));

	checkbox_axes_ = new QCheckBox(tr("Axes"), this);
	connect(checkbox_axes_, SIGNAL(clicked(bool)), ptr_renderingwidget_, SLOT(CheckDrawAxes(bool)));

	checkbox_rotate_ = new QCheckBox(tr("Rotate"), this);
	connect(checkbox_rotate_, SIGNAL(clicked(bool)), ptr_renderingwidget_, SLOT(CheckRotate(bool)));

	checkbox_boundary_ = new QCheckBox(tr("Boundary"), this);
	connect(checkbox_boundary_, SIGNAL(clicked(bool)), ptr_renderingwidget_, SLOT(CheckDrawBoundary(bool)));

	checkbox_interaction_ = new QCheckBox(tr("Interaction"), this);
	connect(checkbox_interaction_, SIGNAL(clicked(bool)), ptr_renderingwidget_, SLOT(CheckDrawInteraction(bool)));

	checkbox_meancurvature_ = new QCheckBox(tr("MeanCurvature"), this);
	connect(checkbox_meancurvature_, SIGNAL(clicked(bool)), ptr_renderingwidget_, SLOT(CheckDrawMeanCurvature(bool)));

	checkbox_absolute_meancurvature_ = new QCheckBox(tr("AbsMeanCurvature"), this);
	connect(checkbox_absolute_meancurvature_, SIGNAL(clicked(bool)), ptr_renderingwidget_, SLOT(CheckDrawAbsMeanCurvature(bool)));

	checkbox_gaussiancurvature_ = new QCheckBox(tr("GaussianCurvature"), this);
	connect(checkbox_gaussiancurvature_, SIGNAL(clicked(bool)), ptr_renderingwidget_, SLOT(CheckDrawGaussianCurvature(bool)));

	checkbox_loss_ = new QCheckBox(tr("Loss"), this);
	connect(checkbox_loss_, SIGNAL(clicked(bool)), ptr_renderingwidget_, SLOT(CheckDrawLoss(bool)));

	checkbox_vectorfield_npoly_ = new QCheckBox(tr("VectorField"), this);
	connect(checkbox_vectorfield_npoly_, SIGNAL(clicked(bool)), ptr_renderingwidget_, SLOT(CheckDrawVectorField(bool)));

	checkbox_plane_ = new QCheckBox(tr("Plane"), this);
	connect(checkbox_plane_, SIGNAL(clicked(bool)), ptr_renderingwidget_, SLOT(CheckDrawPlane(bool)));

	checkbox_cluster_ = new QCheckBox(tr("Cluster"), this);
	connect(checkbox_cluster_, SIGNAL(clicked(bool)), ptr_renderingwidget_, SLOT(CheckDrawCluster(bool)));

	groupbox_render_ = new QGroupBox(tr("Render"), this);
	QVBoxLayout* render_layout = new QVBoxLayout(groupbox_render_);
	render_layout->addWidget(checkbox_point_);
	render_layout->addWidget(checkbox_edge_);
	render_layout->addWidget(checkbox_face_);
	render_layout->addWidget(checkbox_texture_);
	render_layout->addWidget(checkbox_light_);
	render_layout->addWidget(checkbox_axes_);
	render_layout->addWidget(checkbox_rotate_);

	groupbox_chart_ = new QGroupBox(tr("Chart"), this);
	QVBoxLayout* chart_layout = new QVBoxLayout(groupbox_chart_);
	chart_layout->addWidget(checkbox_loss_);
	chart_layout->addWidget(checkbox_meancurvature_);
	chart_layout->addWidget(checkbox_absolute_meancurvature_);
	chart_layout->addWidget(checkbox_gaussiancurvature_);

	groupbox_geometry_ = new QGroupBox(tr("Geometry"), this);
	QVBoxLayout* geometry_layout = new QVBoxLayout(groupbox_geometry_);
	geometry_layout->addWidget(checkbox_boundary_);
	geometry_layout->addWidget(checkbox_interaction_);
	geometry_layout->addWidget(checkbox_vectorfield_npoly_);

	groupbox_else_ = new QGroupBox(tr("Else"), this);
	QVBoxLayout* else_layout = new QVBoxLayout(groupbox_else_);
	else_layout->addWidget(checkbox_plane_);
	else_layout->addWidget(checkbox_cluster_);
}

void MainWindow::ImageWidget_Show()
{
	if (ptr_renderingwidget_->ptr_mesh_->n_vertices() < 3) return;
	if (!ptr_renderingwidget_->is_load_texture_) return;
	if (!ptr_imagewindow_) ptr_imagewindow_ = new ImageWindow(this);
	else ptr_imagewindow_->imagewidget_->Open(ptr_renderingwidget_->texture_img_);
	ptr_imagewindow_->show();
}

void MainWindow::keyPressEvent(QKeyEvent *e)
{

}

void MainWindow::keyReleaseEvent(QKeyEvent *e)
{

}

void MainWindow::ShowMeshInfo(int npoint, int nedge, int nface)
{
	label_meshinfo_->setText(QString("MeshInfo: p: %1 e: %2 f: %3").arg(npoint).arg(nedge).arg(nface));
}

void MainWindow::OpenFile()
{

}

void MainWindow::ShowAbout()
{
	QMessageBox::information(this, "About QtMeshFrame-1.0.1",

		QString("<h3>This MeshFrame provides some operations about *.obj files sunch as") +
		" IO, render with points , edges, triangles or textures and some interactions with mouse."
		" A fix light source is provided for you."
		"This is a basic and raw frame for handling meshes. The mesh is of half_edge struct.\n"
		"Please contact" "<font color=blue> wkcagd@mail.ustc.edu.cn<\font><font color=black>, Kang Wang if you has any questions.<\font><\h3>"
		,
		QMessageBox::Ok);
}