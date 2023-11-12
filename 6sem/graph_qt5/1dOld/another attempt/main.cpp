#include "widget.h"
#include <QApplication>
#include <QtWidgets/QMainWindow>
#include <QScreen>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QAction>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QMessageBox>
#include <string>
#include <functional>

int main(int argc, char *argv[]){
    QApplication a(argc, argv);

    QMainWindow *window = new QMainWindow;
    QMenuBar *tool_bar = new QMenuBar (window);
    Widget *graph_area = new Widget (window);
    QAction *action;
    // Widget w;
    if (graph_area->parse_command_line (argc, argv)) {
        QMessageBox::warning (0, "Wrong input arguments!", "Wrong input arguments!");
        return -1;
    }
    graph_area->plot(graph_area);
    action = tool_bar->addAction ("&Increase n", graph_area, SLOT (incr_n ()));
    action->setShortcut (QString ("Ctrl+C"));
    action = tool_bar->addAction ("&Decrease n", graph_area, SLOT (decr_n ()));
    action->setShortcut (QString ("Ctrl+T"));
    action = tool_bar->addAction ("&Next function", graph_area, SLOT (incrFuncNum ()));
    action->setShortcut (QString ("Ctrl+Y"));
    action = tool_bar->addAction ("&Previos function", graph_area, SLOT (decrFuncNum ()));
    action->setShortcut (QString ("Ctrl+U"));

    action = tool_bar->addAction ("E&xit", window, SLOT (close ()));
    action->setShortcut (QString ("Ctrl+X"));

    window->setMenuBar (tool_bar);
    window->setCentralWidget (graph_area);
    window->setWindowTitle ("Graph");

    QScreen *screen = QGuiApplication::primaryScreen();
    QRect screenGeometry = screen->geometry();
    int height = screenGeometry.height();
    int width = screenGeometry.width();

    window->resize(width, height);
    window->show ();
    a.exec ();
    delete window;
    return 0;

    // graph_area->show();
    //
    // return a.exec();
}
