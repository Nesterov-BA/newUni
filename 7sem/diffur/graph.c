#include "stdio.h"

int numberOfPoints;
double* errorLog;
double* t1;
double* x1;
double* y2;  


void getArray()
{
  FILE *number = fopen("numberOfPoints.txt", "r");
  FILE *errors =fopen("errorLog.txt", "r");
  FILE *data = fopen("data.txt", "r");
  int num1 = fscanf(number, "%d", &numberOfPoints);
  t1 = new double[numberOfPoints];
  x1 = new double[numberOfPoints];
  y2 = new double[numberOfPoints];
  int num;
  for (int i = 0; i<numberOfPoints; i++)
  {
    num =fscanf(data, "%lf %lf %lf", &x1[i], &y2[i], &t1[i]);
    printf("%d, %lf, %lf, %lf \n", num, t1[i], x1[i], y2[i]);
  }
}


void graph() {
   TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,500);

   c1->SetGrid();
   getArray();
   const int n = 20;
   double x[n], y[n];
   for (int i=0;i<n;i++) {
     x[i] = i*0.1;
     y[i] = 10*sin(x[i]+0.2);
//     printf(" i %i %f %f \n",i,x[i],y[i]);
   }
   TGraph *gr = new TGraph(numberOfPoints,x1,y2);
   gr->SetLineColor(2);
   gr->SetLineWidth(4);
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);
   gr->SetTitle("a simple graph");
   gr->GetXaxis()->SetTitle("X title");
   gr->GetYaxis()->SetTitle("Y title");
   gr->Draw("ACP");

   // TCanvas::Update() draws the frame, after which one can change it
   c1->Update();
   c1->GetFrame()->SetBorderSize(12);
   c1->Modified();
}
