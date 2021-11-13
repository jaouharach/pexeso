#include<stdio.h>
#include<math.h>

int main(int argc, char const *argv[])
{
    float x = 2, y = 3, z = 2;
    float nc = (x*y)/(z*z);
    if(fmod(x*y, (z*z)) == 0)
        printf("Num cells = %.2f, OK!\n", nc);
    else   
        printf("Num cells = %.2f, Bad choice!\n", nc);
    
    return 0;
}
