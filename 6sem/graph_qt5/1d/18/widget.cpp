#include "widget.h"
#include "ui_widget.h"
#include <QLogValueAxis>
#include <QLineSeries>
#include <QValueAxis>
#include <QChart>
#include <QPen>
#include <QChartView>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <QLineF>
#include <vector>
#include <coefficents.hpp>
// #include "qcustomplot.h"

using namespace QtCharts;
using namespace std;

bool is_number(std::string s){
    //Eliminate obvious irritants that could spoil the party
    //Handle special cases here, e.g. return true for "+", "-", "" if they are acceptable as numbers to you
    if (s == "" || s == "." || s == "+" || s == "-" || s == "+." || s == "-.") return false;

    //Remove leading / trailing spaces **IF** they are acceptable to you
    while (s.size() > 0 && s[0] == ' ') s = s.substr(1, s.size() - 1);
    while (s.size() > 0 && s[s.size() - 1] == ' ') s = s.substr(0, s.size() - 1);


    //Remove any leading + or - sign
    if (s[0] == '+' || s[0] == '-')
        s = s.substr(1, s.size() - 1);

    //Remove decimal points
    long prevLength = s.size();

    size_t start_pos = 0;
    while((start_pos = s.find(".", start_pos)) != std::string::npos)
        s.replace(start_pos, 1, "");

    //If the string had more than 2 decimal points, return false.
    if (prevLength > s.size() + 1) return false;

    //Check that you are left with numbers only!!
    //Courtesy selected answer by Charles Salvia above
    std::string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}


int Widget::parse_command_line (int argc, char *argv[])
{
  if (argc != 5 || !is_number(std::string(argv[1])) || !is_number(std::string(argv[2])) || !is_number(std::string(argv[3])) || !is_number(std::string(argv[4])))
    return -1;

  a = std::stod(argv[1]); //начало отрезка
  b = std::stod(argv[2]); //конец отрезка
  n = std::stoi(argv[3]); //кол-во точек интерполяции
  func_id = std::stoi(argv[4]);
  if (b - a < 1.e-6 || n <= 0)
    return -1;

  return 0;
}

static int randomBetween(int low, int high, int seed)
{
    qsrand(seed); // Установка базового числа для отсчёта рандома в qrand
    return (qrand() % ((high + 1) - low) + low);
}

static double f_0 (double x){
  return sin(x);
}

void Widget::incr_n (){
  n = n * 2;
  plot(this);
}

void Widget::decr_n (){
    if (n > 3){
        n = n / 2;
    }
  plot(this);

}

void Widget::incrFuncNum(){
    if (func_id >-1 && func_id<6){
        func_id++;
    }
  plot(this);

}


void Widget::decrFuncNum(){
    if (func_id >0 && func_id<7){
        func_id--;
    }
    plot(this);

}


void Widget::plot(QWidget *parent) {


    // Создаём представление графика
    QChartView *chartView = new QChartView();
    // Добавляем его в горизонтальный Layout
    ui->horizontalLayout->addWidget(chartView);
    int k = 3;

    double h = (b - a) / (n - 1);
    double* x_ = new double[n]; // вектор точек интерполирования
    for (int j = 0; j < n ; j++){
        x_[j] = a + h * j;
    }
    double* fValues = new double[n];
    for (int j = 0; j < n; j++){
        fValues[j] = fun(func_id, x_[j]);
    }
    double* coeff = new double[4*n];
    coefficents_eval(x_, fValues, n, func_id, coeff);

    double* evalPoints = new double[k*(n-1)+1];
    double* evalValues = new double[k*(n-1)+1];
    h = h/k;
    for (int j = 0; j < k * (n-1)+1; j++){
        evalPoints[j] = a + h*j;
    }
    for (int j = 0; j < n-1; j++){
        for (int i = 0; i < k; i++)
        {
            evalValues[k*j+i] = approx_eval(evalPoints[k*j+i], coeff, j, x_);
        }
    }
    evalValues[k*(n-1)] = approx_eval(evalPoints[k*(n-1)], coeff, n-2, x_);


    //=====================================================================//



    QLineSeries *series = new QLineSeries();
    for (int j = 0; j < k*(n-1)+1; j++){
        *series << QPointF(evalPoints[j], evalValues[j]);

    }

    QLineSeries *series2 = new QLineSeries();
    for (int j = 0; j < k*(n-1)+1; j++){
        int j1 = j/k;
        *series2 << QPointF(evalPoints[j], fValues[j1]);

    }

    QPen pen;
    pen.setWidth(3);
    pen.setBrush(Qt::black);

    series->setPen(pen);

    QPen pen1;
    pen1.setWidth(3);
    pen1.setBrush(Qt::blue);

    series2->setPen(pen1);

    // Создаём график и добаляем в него синусоиду
    QChart *chart = new QChart();
    chartView->setChart(chart);
    chart->addSeries(series);

    //chart->addSeries(series2);
    chart->legend()->hide();

    QFont font;
    font.setPixelSize(18);
    chart->setTitleFont(font);
    chart->setTitleBrush(QBrush(Qt::black));
    chart->setTitle("f(x) = sin(x)");

    // Добавим всплывающую подсказку для графика
    chart->setToolTip(QString("f(x) = sin(x), a = %1, b = %2, n = %3").arg(a)
                      .arg(b)
                      .arg(n));

    // Настройка осей графика
    int k1;
    if (n < 12){
        k1 = n + 1;
    }
    else {
        k1 = 13;
    }
    QValueAxis *axisX = new QValueAxis();
    axisX->setTitleText("x");
    axisX->setTitleFont(font);
    axisX->setTitleBrush(QBrush(Qt::black));
    axisX->setLabelFormat("%.2f");
    axisX->setTickCount(k1);
    axisX->setRange(a, b);
    chart->addAxis(axisX, Qt::AlignTop);
    series->attachAxis(axisX);

    QValueAxis *axisY = new QValueAxis();
    axisY->setTitleText("f(x)");
    axisY->setTitleFont(font);
    axisY->setTitleBrush(QBrush(Qt::black));
    axisY->setLabelFormat("%.2f");
    axisY->setTickCount(13);
    axisY->setRange(fmin(fValues[0], 0), fmax(fValues[n-1],2));
    chart->addAxis(axisY, Qt::AlignLeft);
    series->attachAxis(axisY);


    QPen axisPen(Qt::black);
    axisPen.setWidth(2);
    axisX->setLinePen(axisPen);
    axisY->setLinePen(axisPen);

    QBrush axisBrush(Qt::black);
    axisX->setLabelsBrush(axisBrush);
    axisY->setLabelsBrush(axisBrush);



    // Устанавливаем график в представление
    system("pause");

    ui->horizontalLayout->removeWidget(chartView);
}


Widget::Widget(QWidget *parent) :
    // QWidget(parent),
    ui(new Ui::Widget)
{

    ui->setupUi(this);

    //plot(this);
}

Widget::~Widget()
{
    delete ui;
}
