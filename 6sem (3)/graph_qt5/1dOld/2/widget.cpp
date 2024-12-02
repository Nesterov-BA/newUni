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

// ================== Soliton of SLAU ================== //

vector<double> solve_tridiagonal_system(vector<double> a, vector<double> b, vector<double> c, vector<double> d) {
    int n = a.size();
    vector<double> x(n);

    // Forward elimination
    for (int i = 1; i < n; i++) {
        double m = a[i] / b[i-1];
        b[i] = b[i] - m * c[i-1];
        d[i] = d[i] - m * d[i-1];
    }

    // Backward substitution
    x[n-1] = d[n-1] / b[n-1];
    for (int i = n - 2; i >= 0; i--) {
        x[i] = (d[i] - c[i] * x[i+1]) / b[i];
    }

    return x;
}

//  int main(void){
//      vector<double> a = {0, 1, 1};
// vector<double> b = {2, 3, 4};
// vector<double> c = {1, 1, 0};
// vector<double> d = {1, 2, 3};
// vector<double> x = solve_tridiagonal_system(a, b, c, d);

// ======================================================= //

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


    double h = M_PI / (n - 1);

    vector<double> x(n, 0);
    for (int j = 0; j < n; j++) {
        x[j] = h * j;
    }

    vector<double> ksi(n + 1, 0);
    ksi[0] = -h / 2;
    ksi[n] = M_PI + h / 2;
    for (int j = 1; j < n; j++) {
        ksi[j] = (x[j - 1] + x[j]) / 2;
    }

    vector<double> a(n + 1, 0);
    a[0] = 0; a[1] = 2 / h;
    vector<double> b(n + 1, 0);
    b[0] = 2 / h; b[n] = 2 / h;
    vector<double> c(n + 1, 0);
    c[n] = 0; c[n - 1] = 2 / h;
    vector<double> d(n + 1, 0);
    d[0] = 0; // sin(0) = 0
    d[n] = 0; // sin(M_PI) = 0
    cout << a[0] << " " << a[1] << " " << b[0] << " " << b[n] << " " << c[n] << " " << c[n - 1] << endl;

    for (int i = 2; i < n + 1; i++){
        a[i] = 1 / h;
        b[i - 1] = 6 / h;
        c[i - 2] = 1 / h;
        d[i - 1] = 4 / h * (sin(x[i - 2]) + sin(x[i - 1]));
        // cout << a[i] << " " << b[i - 1] << " " << c[i - 2] << endl;
    }

    vector<double> v = solve_tridiagonal_system(a, b, c, d);


    //=====================================================================//
    vector<double> pnts(3 * n, 0);
    for (int i = 0; i < 3 * n; i++){
        pnts[i] = i * M_PI / (3 * n - 1);
    }

    QLineSeries *series = new QLineSeries();
    for (int i = 0; i < 3 * n; i++) {
        int i_ = i / 3;
        double a2 = 2 * (sin(x[i_]) - v[i_]) / h;
        double a3 = 2 * (v[i_ + 1] + v[i_] - 2 * sin(x[i_])) / (h * h);
        // cout << pnts[i] << " " << ksi[i_] << " " << x[i_] << " " << v[i_] + a2 * (pnts[i] - ksi[i_]) + a3 * (pnts[i] - ksi[i]) * (pnts[i] - x[i_]) << endl;
        *series << QPointF(pnts[i], v[i_] + a2 * (pnts[i] - ksi[i_]) + a3 * (pnts[i] - ksi[i]) * (pnts[i] - x[i_]));
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
