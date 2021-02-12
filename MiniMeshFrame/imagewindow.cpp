#include "imagewindow.h"
#include <QtWidgets>
#include <QImage>
#include <QPainter>
#include "imagewidget.h"


ImageWindow::ImageWindow(MainWindow *mainwindow)
	: ptr_mainwindow_(mainwindow)
{
	setGeometry(mainwindow->pos().x()+ mainwindow->size().width()/4*3-8.5, mainwindow->pos().y()+128, max(mainwindow->size().width()/4,300), max(mainwindow->size().width()/4,300));

	ptr_imagewidget_ = new ImageWidget(this);
	setCentralWidget(ptr_imagewidget_);
	setWindowIcon(QIcon("./Resources/images/Icon.jpg"));

	CreateActions();
	CreateMenus();
	CreateToolBars();
	CreateStatusBar();

//////////////////////////////////////////////////////////////////////////////////////////////////////////
	action_pt_uniform = new QAction(tr("pt_uniform"), this);
	action_pt_weightedleastsquares = new QAction(tr("pt_weightedleastsquares"), this);
	action_pt_shapepreserving = new QAction(tr("pt_shapepreserving"), this);
	action_pt_asap = new QAction(tr("pt_asap"), this);
	action_pt_arap = new QAction(tr("pt_arap"), this);
	action_pt_slim = new QAction(tr("pt_slim"), this);
	connect(action_pt_uniform, SIGNAL(triggered()), ptr_imagewidget_, SLOT(Pt_Uniform()));
	connect(action_pt_weightedleastsquares, SIGNAL(triggered()), ptr_imagewidget_, SLOT(Pt_Weightedleastsquares()));
	connect(action_pt_shapepreserving, SIGNAL(triggered()), ptr_imagewidget_, SLOT(Pt_Shapepreserving()));
	connect(action_pt_asap, SIGNAL(triggered()), ptr_imagewidget_, SLOT(Pt_Asap()));
	connect(action_pt_arap, SIGNAL(triggered()), ptr_imagewidget_, SLOT(Pt_Arap()));
	connect(action_pt_slim, SIGNAL(triggered()), ptr_imagewidget_, SLOT(Pt_Slim()));

	this->addToolBarBreak(Qt::TopToolBarArea);
	toolbar_parameterization_ = addToolBar(tr("toolbar_parameterization_"));
	toolbar_parameterization_->addAction(action_pt_asap);
	toolbar_parameterization_->addAction(action_pt_arap);
	toolbar_parameterization_->addAction(action_pt_slim);
	toolbar_parameterization_->addAction(action_pt_uniform);
	toolbar_parameterization_->addAction(action_pt_weightedleastsquares);
	toolbar_parameterization_->addAction(action_pt_shapepreserving);

//////////////////////////////////////////////////////////////////////////////////////////////////////////
}

ImageWindow::~ImageWindow()
{

}

void ImageWindow::closeEvent(QCloseEvent *e)
{

}

void ImageWindow::paintEvent(QPaintEvent* paintevent)
{

}

void ImageWindow::AboutImageWarping()
{
	QMessageBox::about(this, tr("ImageWarping"), tr("This is the CG homework made by PB17000123"));
}

void ImageWindow::AboutGCL()
{
	QMessageBox::about(this, tr("GCL"), tr("Graphics&Geometric Computing Laboratory (GCL) "
		"at the University of Science and Technology of China (USTC) "
		"primarily focuses on the mathematical foundations of <b>computer graphics, geometric modeling,</b> "
		"<b>image processing, scientific computing, and information visualization.</b> "
		"The lab was built in 2012 with faculties from School of Mathematical Sciences, "
		"School of Information Science and Technology, and School of Computer Science and Technology. "
		"The USTC graphics and geometric computing research efforts are led by <b>Prof. Ligang Liu.</b>"));
}

void ImageWindow::CreateActions()
{
	minidraw_about_action_ = new QAction(tr("&About MiniDraw..."), this);
	minidraw_about_action_->setStatusTip(tr("About MiniDraw..."));
	connect(minidraw_about_action_, &QAction::triggered, this, &ImageWindow::AboutImageWarping);

	gcl_about_action_ = new QAction(tr("&About GCL..."), this);
	gcl_about_action_->setStatusTip(tr("About GCL..."));
	connect(gcl_about_action_, &QAction::triggered, this, &ImageWindow::AboutGCL);

	QDir::setCurrent(QCoreApplication::applicationDirPath());

	action_new_ = new QAction(QIcon(":/MiniMeshFrame/Resources/images/new.png"), tr("&New"), this);
	action_new_->setShortcut(QKeySequence::New);
	action_new_->setStatusTip(tr("Create a new file"));

	action_open_ = new QAction(QIcon(":/MiniMeshFrame/Resources/images/open.png"), tr("&Open..."), this);
	action_open_->setShortcuts(QKeySequence::Open);
	action_open_->setStatusTip(tr("Open an existing file"));
	//connect(action_open_, SIGNAL(triggered()), ptr_imagewidget_, SLOT(Open()));
	connect(action_open_, SIGNAL(triggered()), ptr_mainwindow_->ptr_renderingwidget_, SLOT(LoadTexture()));
	
	action_save_ = new QAction(QIcon(":/MiniMeshFrame/Resources/images/save.png"), tr("&Save"), this);
	action_save_->setShortcuts(QKeySequence::Save);
	action_save_->setStatusTip(tr("Save the document to disk"));
	connect(action_save_, SIGNAL(triggered()), ptr_imagewidget_, SLOT(Save()));

	action_saveas_ = new QAction(tr("Save &As..."), this);
	//action_saveas_->setShortcuts(QKeySequence::SaveAs);
	action_saveas_->setStatusTip(tr("Save the document under a new name"));
	connect(action_saveas_, SIGNAL(triggered()), ptr_imagewidget_, SLOT(SaveAs()));

	action_restore_ = new QAction(tr("Restore"), this);
	action_restore_->setStatusTip(tr("Show origin image"));
	connect(action_restore_, SIGNAL(triggered()), ptr_imagewidget_, SLOT(Restore()));

	action_exit_ = new QAction(tr("Exit"), this);
	action_exit_->setStatusTip(tr("Exit"));
	connect(action_exit_, &QAction::triggered, this, &QApplication::exit);
}

void ImageWindow::CreateMenus()
{
	menu_file_ = menuBar()->addMenu(tr("&File"));
	menu_file_->setStatusTip(tr("File menu"));
	menu_file_->addAction(action_new_);
	menu_file_->addAction(action_open_);
	menu_file_->addAction(action_save_);
	menu_file_->addAction(action_saveas_);
	menu_file_->addAction(action_exit_);

	menu_edit_ = menuBar()->addMenu(tr("&Edit"));
	menu_edit_->setStatusTip(tr("Edit menu"));
	menu_edit_->addAction(action_restore_);

	about_menu_ = menuBar()->addMenu(tr("&About"));
	about_menu_->addAction(minidraw_about_action_);
	about_menu_->addAction(gcl_about_action_);
}

void ImageWindow::CreateToolBars()
{
	toolbar_file_ = addToolBar(tr("File"));
	toolbar_file_->addAction(action_new_);
	toolbar_file_->addAction(action_open_);
	toolbar_file_->addAction(action_save_);
	toolbar_file_->addAction(action_saveas_);

	// Add separator in toolbar 
	toolbar_file_->addSeparator();
	toolbar_file_->addAction(action_restore_);
}

void ImageWindow::CreateStatusBar()
{
	statusBar()->showMessage(tr("Ready"));
}
