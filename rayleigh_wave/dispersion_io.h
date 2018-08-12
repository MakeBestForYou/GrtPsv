#ifndef DISPERSION_IO_H
#define DISPERSION_IO_H

float** model_read(char* filename,int *model_layer);

void nextline(FILE* fp,int n_line);

#endif // DISPERSION_IO_H
