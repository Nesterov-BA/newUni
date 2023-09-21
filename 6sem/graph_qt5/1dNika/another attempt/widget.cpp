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

  a = std::stod(argv[1]);
  b = std::stod(argv[2]);
  n = std::stoi(argv[3]);
  func_id = 0;

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

void Widget::plot(QWidget *parent) {
    a = 0;
    b = M_PI;

    // Создаём представление графика
    QChartView *chartView = new QChartView();
    // Добавляем его в горизонтальный Layout
    ui->horizontalLayout->addWidget(chartView);


    double h = (b - a) / (n - 1);
    vector<double> x_(n, 0);
    for (int j = 1; j < n + 1; j++){
        x_[j - 1] = a + h * (j - 1);
    }

    vector<double> d(n, 0);
    for (int j = 1; j < n - 1; j++){
        d[j] = (sin(x_[j + 1]) - sin(x_[j - 1])) / (2 * h);
    }
    d[0] = (3 * (sin(x_[1]) - sin(x_[0])) - d[1]) / 2;
    d[n - 1] = (3 * (sin(x_[n - 1]) - sin(x_[n - 2])) - d[n - 2]) / 2;

    vector<double> x(4 * n, 0);
    for (int j = 0; j < 4 * n; j++){
        if (j % 4 == 0){
            x[j] = x_[j / 4];
        }
        else {
            x[j] = x[j / 4] + (b - a) / (4 * n) * (j % 4);
        }
    }



    //=====================================================================//




    QLineSeries *series = new QLineSeries();
    for (int j = 0; j < n; j++){
        for (int i = 0; i < 4; i+=4){
            *series << QPointF(x[j * 4 + i],
                sin(x_[j]) +
                 d[j] * (x[j * 4 + i] - x_[j]) +
                 (3 * (sin(x_[j + 1]) - sin(x_[j])) / (x_[j + 1] - x_[j]) - 2 * d[j] - 2 * d[j + 1]) / (x_[j + 1] - x_[j]) * (x[j * 4 + i] - x_[j]) * (x[j * 4 + i] - x_[j]) +
                 (d[j] + d[j + 1] - 2 * ((sin(x_[j + 1]) - sin(x_[j])) / (x_[j + 1] - x_[j]))) / ((x_[j + 1] - x_[j])*(x_[j + 1] - x_[j])) * (x[j * 4 + i] - x_[j]) * (x[j * 4 + i] - x_[j]) * (x[j * 4 + i] - x_[j]));
        }
    }

    QPen pen;
    pen.setWidth(3);
    pen.setBrush(Qt::black);

    series->setPen(pen);

    // Создаём график и добаляем в него синусоиду
    QChart *chart = new QChart();
    chartView->setChart(chart);
    chart->addSeries(series);
    chart->legend()->hide();

    QFont font;
    font.setPixelSize(18);
    chart->setTitleFont(font);
    chart->setTitleBrush(QBrush(Qt::black));
    chart->setTitle("f(x) = sin(x)");

    // Добавим всплывающую подсказку для графика
    chart->setToolTip(QString("f(x) = sin(x), n = %1").arg(n));

    // Настройка осей графика
    int k;
    if (n < 12){
        k = n + 1;
    }
    else {
        k = 13;
    }
    QValueAxis *axisX = new QValueAxis();
    axisX->setTitleText("x");
    axisX->setTitleFont(font);
    axisX->setTitleBrush(QBrush(Qt::black));
    axisX->setLabelFormat("%.2f");
    axisX->setTickCount(k);
    axisX->setRange(0, 4);
    chart->addAxis(axisX, Qt::AlignBottom);
    series->attachAxis(axisX);

    QValueAxis *axisY = new QValueAxis();
    axisY->setTitleText("f(x)");
    axisY->setTitleFont(font);
    axisY->setTitleBrush(QBrush(Qt::black));
    axisY->setLabelFormat("%.2f");
    axisY->setTickCount(5);
    axisY->setRange(0, 1.5);
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

    plot(this);
}

Widget::~Widget()
{
    delete ui;
}
