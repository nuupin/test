#include <stdio.h>
#include <math.h>

int main(void) {

    double f, f1, f2, dx, x, k1, k2, k3, k4, k;

    f1 = 1;
    dx = 0.01;

    for (int i=0; i<=100; i++) {
        x = i*dx;
        k1 = dx*x*f1;
        k2 = dx*(x + dx/2.)*(f1 + k1/2.);
        k3 = dx*(x + dx/2.)*(f1 + k2/2.);
        k4 = dx*(x + dx)*(f1 + k3);
        k = (k1 + 2*k2 + 2*k3 + k4)/6.;

        f2 = f1 + k;
        f = exp(x*x/2.);

        printf("x = %f, 計算値-> f1 = %f, 正確な値-> f = %f \n", x, f1, f);
        printf("こんにちは\n", );
        f1 = f2;
    }
    return 0;
}
