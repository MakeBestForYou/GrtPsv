#include<stdio.h>
#include<stdlib.h>
#include"dispersion_io.h"

void nextline(FILE* fp,int n_line)
{
    char buf[1024] = { 0 };
    for(int i = 0;i<n_line;i++)
    {
        fgets(buf, sizeof(buf), fp);
    }
}
