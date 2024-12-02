#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QKeyEvent>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

int MainWindow::parseArgs()
{
    return ui->widget->parse_command_line();
}


void MainWindow::keyPressEvent(QKeyEvent *event)
{
    if (event->key()==Qt::Key_X &&
            (event->modifiers() & Qt::ControlModifier)){
        close();
        event->accept();
    }
}

