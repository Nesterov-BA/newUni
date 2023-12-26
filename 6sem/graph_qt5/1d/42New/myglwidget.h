#ifndef MYGLWIDGET_H
#define MYGLWIDGET_H
#include <QGLWidget>

class MyGLWidget : public QGLWidget{
    Q_OBJECT
    int n; //number of points by Y
    int k; //id of the approxiamted function
    double step;
    int nInt = 20;
    double a; //left end by X
    double b; //right end by X
    double min_z; //function's minimum
    double max_z; //function's maximum
    double absmax; //function's absolute maximum
    const char *f_name;
    double(*f)(double); //function
    double(*d2f)(double); //2nd derivative
    int view_id=0; //what to be viewed
    int s=0; //squeeze-stretch
    int p=0; //disturbance
    double d2Start, d2End;
    double *F=nullptr; //matrix of values
    double *apprVal=nullptr;
    double *intPoints=nullptr; //matrix of interpolation coefficients
    double *coeff=nullptr; //intermediate matrix
    double *c=nullptr; //array of points by X
    double *ksi=nullptr;
    double extr[2]; //factic extrema
public:

    MyGLWidget(QWidget* parent=nullptr);
    virtual ~MyGLWidget();
    int parse_command_line();
    void allocate();
    void print_console();
    void change_func();
    void extrema_hunt();
    void draw_area();
    void func_graph();
    void appr_graph();
    void err_graph();
    void press0();
    void press23();
    void press45();
    void press67();
    void printwindow();
protected:
    void initializeGL();
    void paintGL();
    void resizeGL(int width, int height);
    QSize minimumSizeHint() const;
    QSize sizeHint() const;
    void keyPressEvent(QKeyEvent *e);
    void setProjection();

public slots:
    // slots for xyz-rotation slider
    void setXRotation(int angle);
    void setYRotation(int angle);
    void setZRotation(int angle);
signals:
    // signaling rotation from mouse movement
    void xRotationChanged(int angle);
    void yRotationChanged(int angle);
    void zRotationChanged(int angle);
private:
    //void drawPyramid();
    int xRot;
    int yRot;
    int zRot;
    QPoint lastPos;
};

#endif // MYGLWIDGET_H
