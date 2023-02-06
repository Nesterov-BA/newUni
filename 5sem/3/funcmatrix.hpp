#ifndef MYFUNC
#define MYFUNC


double func(int n, int i, int j)
{
    i++;
    j++;
     double dn = double(n);
    int max = i > j ? i : j;
    int min = i < j ? i : j;
    if(i == j){
        if (i == 1){
            return -1;
        } else if (i == n){
            return -(dn-1)/dn;
        } else {
            return -2;
        }
    } else if (i==j+1 || j == i+1){
        return 1;
    } else {
        return 0;
    }
}

#endif