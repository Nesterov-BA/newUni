#include "mainwindow.h"
#include <QApplication>

int main(int argc, char *argv[]){
    QApplication a(argc, argv);
    MainWindow w;
    if (w.parseArgs()){
        qWarning("Wrong input arguments!, error code: %d", w.parseArgs());
        return -1;
    }
    w.show();
    return a.exec();
}
