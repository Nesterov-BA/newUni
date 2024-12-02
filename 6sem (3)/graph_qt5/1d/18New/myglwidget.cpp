#include <QtWidgets>
#include <QtOpenGL>
#include <GL/glu.h>
#include "myglwidget.h"
#include "help.hpp"
#include "chebyshev.hpp"

MyGLWidget::MyGLWidget(QWidget *parent)
    : QGLWidget(QGLFormat(QGL::SampleBuffers), parent){

}

MyGLWidget::~MyGLWidget(){
    free(F);
    free(c);
    free(coeff);
    free(apprVal);
    free(intPoints);
}

int MyGLWidget::parse_command_line(){
    //FILE *ff;
    QStringList args=QApplication::arguments();
    if(args.size()!=5)
        return -3;
    /*QByteArray byteFilename = args.at(1).toLocal8Bit();
    const char* fileName = byteFilename.data();
    qWarning("%c\n", *fileName);
    ff=fopen("obl.txt", "r");
    if(!ff)
        return -1;

    fscanf(ff, "%lf", &a);
    fscanf(ff, "%lf", &c);
    fscanf(ff, "%lf", &b);
    fscanf(ff, "%lf", &d);
    */
    n = args.at(1).toInt();
    a = args.at(2).toInt();
    b = args.at(3).toInt();
    k = args.at(4).toInt();

    if(k<0 || k>7)
        return -9;
    step =(b-a)/n;
    allocate();
    points(c, n, a, b);
    points(intPoints, (n-1)*nInt+1, a, b);
    change_func();
    coefficents_eval(c, F, n, coeff, d2Start, d2End);
    evalArr(apprVal, nInt, n, coeff, c, intPoints);
    extrema_hunt();
    print_console();
    return 0;
}

void MyGLWidget::allocate()
{
    if(F)
        free(F);
    if(c)
        free(c);
   if(coeff)
        free(coeff);
   if(apprVal)
        free(apprVal);
   if(intPoints)
       free(intPoints);

    F=(double*)malloc((n)*sizeof(double));
    c=(double*)malloc((n)*sizeof(double));
    coeff = (double*)malloc(((n-1)*4)*sizeof(double));
    apprVal = (double*)malloc(((n-1)*nInt+1)*sizeof(double));
    intPoints = (double*)malloc(((n-1)*nInt+1)*sizeof(double));
}

void MyGLWidget::print_console(){
    printf("format: %d\n", view_id);
    printf("area: [%lf;%lf]\n", a,b);
    printf("squeeze-stretch: %d\n", s);
    printf("points: %d \n", n);
    printf("disturbance: %d\n", p);
    printf("function absmax: %10.3e\n", absmax);
    printf("factic absmax: %10.3e\n", max(fabs(extr[0]),fabs(extr[1])));
    printf("extrema: %10.3e;%10.3e\n", extr[0], extr[1]);
}

void MyGLWidget::change_func(){
    switch(k){
        case 0:
            f_name="k=0 f(x,y)=1";
            f=f0;
            d2f = d2f0;
            break;
        case 1:
            f_name="k=1 f(x,y)=x";
            f=f1;
            d2f = d2f1;
            break;
        case 2:
            f_name="k=2 f(x,y)=x^2";
            f=f2;
            d2f = d2f2;
            break;
        case 3:
            f_name="k=3 f(x,y)=x^3";
            f=f3;
            d2f = d2f3;
            break;
        case 4:
            f_name="k=4 f(x,y)=x^4";
            f=f4;
            d2f = d2f4;
            break;
        case 5:
            f_name="k=5 f(x,y)=exp(x)";
            f=f5;
            d2f = d2f5;
            break;
        case 6:
            f_name="k=6 f(x,y)=1/(25*x^2+1)";
            f=f6;
            d2f = d2f6;
            break;
    }
    Fill_F(F, n, c, f);
    d2Start = d2f(a);
    d2End = d2f(b);
    max_z=max_matr(F,n);
    min_z=min_matr(F,n);
    absmax=max(fabs(max_z),fabs(min_z));
//    if(p)
//        F[ny*(nx/2)+ny/2]+=(0.1*absmax*p);
    /*if(p){
        F[(ny+2)*(nx/2+1)+ny/2+1]+=(0.1*absmax*p);
        if(F[(ny+2)*(nx/2+1)+ny/2+1]>max_z)
            max_z=F[(ny+2)*(nx/2+1)+ny/2+1];
        if(F[(ny+2)*(nx/2+1)+ny/2+1]<min_z)
            min_z=F[(ny+2)*(nx/2+1)+ny/2+1];
    }*/
}

void MyGLWidget::extrema_hunt(){
    if(view_id==0){
        extr[1]=max_z;
        extr[0]=min_z;
        if(p>0 && F[n/2]>max_z)
            extr[1]=F[n/2];
        if(p<0 && F[n/2]<min_z)
            extr[0]=F[n/2];
    }
    if(view_id==1){
        extr[1]=max_matr(apprVal,(n-1)*nInt + 1);
        extr[0]=min_matr(apprVal,(n-1)*nInt + 1);
    }
    if(view_id==2){
        extr[0]=F[0]-apprVal[0];
        extr[1]=F[0]-apprVal[0];
        for(int i=0; i<n; ++i){
            if(F[i] - apprVal[nInt*i]>extr[1])
                extr[1]=F[i] - apprVal[nInt*i];
            if(F[i] - apprVal[nInt*i]<extr[0])
                extr[0]=F[i] - apprVal[nInt*i];
        }
    }
}

void MyGLWidget::draw_area(){
    glColor3f(0.0,0.0,0.0);
    glBegin(GL_LINE_STRIP);
    for(int j=0; j<n; ++j)
        glVertex2f(c[j],0);
    glEnd();
}

void MyGLWidget::func_graph(){
    glColor3f(1.0,1.0,0.0);
    glBegin(GL_LINE_STRIP);
    for(int j=0; j<n; ++j)
        glVertex2f(c[j],F[j]);
    glEnd();
}

void MyGLWidget::appr_graph(){
    glColor3f(1.0,0.0,1.0);
    glBegin(GL_LINE_STRIP);
    for(int j=0; j<(n-1)*nInt+1; ++j)
        glVertex2f(intPoints[j],apprVal[j]);
    glEnd();
}

void MyGLWidget::err_graph(){
    glColor3f(0.0,1.0,1.0);
    glBegin(GL_LINE_STRIP);
    for(int j=0; j<n; ++j)
        glVertex2f(c[j],0);
    glEnd();
}

void MyGLWidget::press0(){ // change function
    allocate();
    points(c, n, a, b);
    points(intPoints, (n-1)*nInt+1, a, b);
    change_func();
    coefficents_eval(c, F, n, coeff, d2Start, d2End);
    evalArr(apprVal, nInt, n, coeff, c, intPoints);
}

void MyGLWidget::press23(){ // change size of area
    allocate();
    points(c, n, a, b);
    points(intPoints, (n-1)*nInt+1, a, b);
    Fill_F(F, n, c, f);
    d2Start = d2f(a);
    d2End = d2f(b);
    max_z=max_matr(F,n);
    min_z=min_matr(F,n);
    absmax=max(fabs(max_z),fabs(min_z));
    coefficents_eval(c, F, n, coeff, d2Start, d2End);
    evalArr(apprVal, nInt, n, coeff, c, intPoints);
    extrema_hunt();
    print_console();
    //transponse(   apprVal, nx*ny);
}

void MyGLWidget::press45(){ // change number of interpolation points
    allocate();
    points(c, n, a, b);
    points(intPoints, (n-1)*nInt+1, a, b);
    Fill_F(F, n, c, f);
    d2Start = d2f(a);
    d2End = d2f(b);
    max_z=max_matr(F,n);
    min_z=min_matr(F,n);
    absmax=max(fabs(max_z),fabs(min_z));
    coefficents_eval(c, F, n, coeff, d2Start, d2End);
    evalArr(apprVal, nInt, n, coeff, c, intPoints);
    extrema_hunt();
    print_console();
    //transponse(   apprVal, nx*ny);
}

void MyGLWidget::press67(){ // change disturbance
    Fill_F(F, n, c, f);
    d2Start = d2f(a);
    d2End = d2f(b);
    max_z=max_matr(F,n);
    min_z=min_matr(F,n);
    absmax=max(fabs(max_z),fabs(min_z));
    coefficents_eval(c, F, n, coeff, d2Start, d2End);
    evalArr(apprVal, nInt, n, coeff, c, intPoints);
    extrema_hunt();
    print_console();
}

void MyGLWidget::printwindow(){
    qglColor(Qt::black);
//    if(nx<5){
//        renderText(10, 30, QString("TOO FEW POINTS by X"));
//    }
//    if(nx>50){
//        renderText(10, 30, QString("TOO MUCH POINTS by X"));
//    }
//    if(ny<5){
//        renderText(10, 60, QString("TOO FEW POINTS by Y"));
//    }
//    if(ny>50){
//        renderText(10, 60, QString("TOO MUCH POINTS by Y"));
//    }
    //if(nx<=50 && nx>=5 && ny<=50 && ny>=5){
        renderText(0, 15, f_name);
        renderText(10, 30, QString("format: %1").arg(view_id));
        renderText(10, 60, QString("points: %1").arg(n));
        renderText(10, 75, QString("p: %1").arg(p));
        renderText(10, 90, QString("absmax(fact): %1( %2)").arg(absmax).arg(max(fabs(extr[0]), fabs(extr[1]))));
    //}
}

QSize MyGLWidget::minimumSizeHint() const{
    return QSize(50, 50);
}

QSize MyGLWidget::sizeHint() const{
    return QSize(400, 400);
}

void MyGLWidget::setXRotation(int angle){
    if (angle != xRot) {
        xRot = angle;
        angle = angle>180 ? angle-360 : angle;
        emit xRotationChanged(angle);
        updateGL();
    }
}

void MyGLWidget::setYRotation(int angle){
    if (angle != yRot) {
        yRot = angle;
        angle = angle>180 ? angle-360 : angle;
        emit yRotationChanged(angle);
        updateGL();
    }
}

void MyGLWidget::setZRotation(int angle){
    if (angle != zRot) {
        zRot = angle;
        angle = angle>180 ? angle-360 : angle;
        emit zRotationChanged(angle);
        updateGL();
    }
}

void MyGLWidget::initializeGL(){
    qglClearColor(Qt::white);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);
    static GLfloat lightPosition[4] = { 0, 0, 10, 1.0 };
    glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
    glEnable(GL_LIGHT0);
}

void MyGLWidget::setProjection() {
    double modelM[16],projM[16];
    int view[4]={-1,-1,2,2};
    double verts[8][3];
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glGetDoublev(GL_MODELVIEW_MATRIX,modelM);
    glGetDoublev(GL_PROJECTION_MATRIX,projM);
    // ищем проекции вершин куба [a,b]x[c,d]x[minz,maxz] на плоскость экрана
    gluProject(a,min(extr[0],0),min(extr[0],0),modelM,projM,view,&verts[0][0],&verts[0][1],&verts[0][2]);
    gluProject(b,min(extr[0],0),min(extr[0],0),modelM,projM,view,&verts[1][0],&verts[1][1],&verts[1][2]);
    gluProject(a,max(extr[1],0),min(extr[0],0),modelM,projM,view,&verts[2][0],&verts[2][1],&verts[2][2]);
    gluProject(b,max(extr[1],0),min(extr[0],0),modelM,projM,view,&verts[3][0],&verts[3][1],&verts[3][2]);
    gluProject(a,min(extr[0],0),max(extr[1],0),modelM,projM,view,&verts[4][0],&verts[4][1],&verts[4][2]);
    gluProject(b,min(extr[0],0),max(extr[1],0),modelM,projM,view,&verts[5][0],&verts[5][1],&verts[5][2]);
    gluProject(a,max(extr[1],0),max(extr[1],0),modelM,projM,view,&verts[6][0],&verts[6][1],&verts[6][2]);
    gluProject(b,max(extr[1],0),max(extr[1],0),modelM,projM,view,&verts[7][0],&verts[7][1],&verts[7][2]);
    // ищем квадрат в плоскости экрана куда они вписываются
    double minX=verts[0][0], maxX=verts[0][0],
            minY=verts[0][1], maxY=verts[0][1],
            minZ=verts[0][2], maxZ=verts[0][2];
    for (int i=1; i<8; i++) {
        if (verts[i][0]<minX) minX=verts[i][0];
        if (verts[i][1]<minY) minY=verts[i][1];
        if (verts[i][2]<minZ) minZ=verts[i][2];
        if (verts[i][0]>maxX) maxX=verts[i][0];
        if (verts[i][1]>maxY) maxY=verts[i][1];
        if (verts[i][2]>maxZ) maxZ=verts[i][2];
    }
    double sz=maxX-minX;
    if (maxY-minY>sz) sz=maxY-minY;
    sz*=1.1;
    // "центр масс"
    double sx=(maxX+minX)/2,sy=(maxY+minY)/2;
    glOrtho(sx-sz/2, sx+sz/2, sy-sz/2, sy+sz/2, -20, +20);
    glMatrixMode(GL_MODELVIEW);
}

void MyGLWidget::paintGL(){
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    setProjection();
   // if(nx>=5 && nx<=50 && ny<=50 && ny>=5){
        glLineWidth(2.0);
        glBegin(GL_LINES);
        glColor3d(1.0,0.0,0);
        glVertex2d(a,0);
        glVertex2d(b,0);
        glColor3d(0,0.0,1.0);
        glVertex2d(0,min(0,extr[0]));
        glVertex2d(0,max(0,extr[1]));
        glEnd();
        if(view_id!=2)
            draw_area();
        if(view_id==0)
            func_graph();
        if(view_id==1)
            appr_graph();
        if(view_id==2)
            err_graph();
   // }
    printwindow();
}

void MyGLWidget::resizeGL(int width, int height){
    glViewport(0, 0, width, height);
}


void MyGLWidget::keyPressEvent(QKeyEvent* e){
    switch (e->key()){
        case Qt::Key_0:
            k=(k+1)%7;
            press0();
            break;
        case Qt::Key_1:
            view_id=(view_id+1)%3;
            break;
        case Qt::Key_2:
            { double w=(b-a)/2;
            a=a-w;
            b=b+w;
            ++s;
            press23();}
            break;
        case Qt::Key_3:
            { double w=(b-a)/4;
            a=a+w;
            b=b-w;
            --s;
            press23();}
            break;
        case Qt::Key_4:
            n*=2;
            press45();
            break;
        case Qt::Key_5:
            n/=2;
            press45();
            break;
        case Qt::Key_6:
            ++p;
            F[n/2]+=(0.1*absmax);
            press67();
            break;
        case Qt::Key_7:
            --p;
            F[n/2]-=(0.1*absmax);
            press67();
            break;
    }
    if (e->key()==Qt::Key_X &&
            (e->modifiers() & Qt::ControlModifier)){
        close();
    }
    e->accept();
    extrema_hunt();
    print_console();
    updateGL();
}
