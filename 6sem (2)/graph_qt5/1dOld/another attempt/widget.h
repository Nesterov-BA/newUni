#include <functional>

#ifndef WIDGET_H
#define WIDGET_H

#include <QWidget>

namespace Ui {
    class Widget;
}

class Widget : public QWidget {
    Q_OBJECT

public:
    double a, b;
    int func_id = 0;
    int n;

    int parse_command_line (int argc, char *argv[]);
    std::string f_name;
    std::function <double(double)> f;

    void plot(QWidget *parent);

    explicit Widget(QWidget *parent = 0);
    // explicit Widget(QWidget *parent = 0, int a = 0);
    ~Widget();

private:
    Ui::Widget *ui;
public slots:
    void incr_n();
    void decr_n();
    void incrFuncNum();
    void decrFuncNum();
};


#endif // WIDGET_H
