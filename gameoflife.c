#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <endian.h>
#include <mpi.h>

#define calcIndex(width, x,y)  ((y)*(width) + (x))

void show(unsigned* currentfield, int w, int h) {
  printf("\033[H");
  for (int y = 0; y < h; y++) {
    for (int x = 0; x < w; x++) printf(currentfield[calcIndex(w, x,y)] ? "\033[07m  \033[m" : "  ");
    printf("\033[E");
  }
  fflush(stdout);
}


float convert2BigEndian( const float inFloat )
{
  float retVal;
  char *floatToConvert = ( char* ) & inFloat;
  char *returnFloat    = ( char* ) & retVal;

  // swap the bytes into a temporary buffer
  returnFloat[0] = floatToConvert[3];
  returnFloat[1] = floatToConvert[2];
  returnFloat[2] = floatToConvert[1];
  returnFloat[3] = floatToConvert[0];

  return retVal;
}

void writeVTK(unsigned* currentfield, int w, int h, int t, char* prefix) {
  char name[1024] = "\0";
  sprintf(name, "%s_%d.vtk", prefix, t);
  FILE* outfile = fopen(name, "w");

  /*Write vtk header */
  fprintf(outfile,"# vtk DataFile Version 3.0\n");
  fprintf(outfile,"frame %d\n", t);
  fprintf(outfile,"BINARY\n");
  fprintf(outfile,"DATASET STRUCTURED_POINTS\n");
  fprintf(outfile,"DIMENSIONS %d %d %d \n", w, h, 1);
  fprintf(outfile,"SPACING 1.0 1.0 1.0\n");//or ASPECT_RATIO
  fprintf(outfile,"ORIGIN 0 0 0\n");
  fprintf(outfile,"POINT_DATA %d\n", h*w);
  fprintf(outfile,"SCALARS data float 1\n");
  fprintf(outfile,"LOOKUP_TABLE default\n");

  for (int y = 0; y < h; y++) {
    for (int x = 0; x < w; x++) {
      float value = currentfield[calcIndex(w, x,y)]; // != 0.0 ? 1.0:0.0;
      value = convert2BigEndian(value);
      fwrite(&value, 1, sizeof(float), outfile);
    }
  }
  fclose(outfile);
}

int getNumberOfAliveCells(unsigned* currentfield, int w, int h, int x, int y) {
  int aliveFields = 0;
  for (int j = y - 1; j <= y + 1; j = j + 1) {
    for (int i = x - 1; i <= x + 1; i = i + 1) {
      if (currentfield[calcIndex(w, i,j)]) {
        aliveFields += 1;
      }
    }
  }
  return aliveFields;
}



int evolve(unsigned* currentfield, unsigned* newfield, int w, int h, MPI_Status status, MPI_Comm comm, int rank, int size) {
  int changes = 0;
  /*
   for (int y = 0; y < h; y++) {
   for (int x = 0; x < w; x++) {
   */
  //TODO add game of life rules

  int sendcount = 1;
  int *sendbuf = calloc(sendcount, sizeof(unsigned));;
  MPI_Datatype sendtype = MPI_UNSIGNED;
  int dest = (rank + 1) % size;
  int sendtag = 1;

  int recvcount = 1;
  int *recvbuf = calloc(recvcount, sizeof(unsigned));
  MPI_Datatype recvtype = MPI_UNSIGNED;
  int source = (rank + size - 1) % size;
  int recvtag = 2;
  MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount, recvtype, source, recvtag, comm, &status);

  //TODO if changes == 0, the time loop will not run!
  return changes;
}

void filling(unsigned* currentfield, int w, int h) {
  for (int i = 0; i < h*w; i++) {
    currentfield[i] = (rand() < RAND_MAX / 10) ? 1 : 0; ///< init domain randomly
  }
}

void game(int w, int h, int timesteps, MPI_Status status, MPI_Comm comm, int rank, int size) {
  unsigned *currentfield = calloc(w*h, sizeof(unsigned));
  unsigned *newfield     = calloc(w*h, sizeof(unsigned));

  MPI_Comm newComm;
  int ndims = 1;
  int *dims = calloc(ndims, sizeof(unsigned));
  dims[0] = size;
  int *periods = calloc(ndims, sizeof(int));
  periods[0] = 1;
  int reorder = 1;

  MPI_Cart_create(comm, ndims, dims, periods, reorder, &newComm);

  int *myCoords = calloc(ndims, sizeof(int));
  MPI_Cart_coords(newComm, rank, ndims, myCoords);

  /*
   int *neighbourCoords = calloc(ndims, sizeof(int));
   neighbourCoords[0] = myCoords[0] - 1;
   int *negNeighbour = calloc (ndims, sizeof(int));
   MPI_Cart_rank(newComm, neighbourCoords, negNeighbour);
   int *posNeighbour = calloc (ndims, sizeof(int));
   neighbourCoords[0] = myCoords[0] + 1;
   MPI_Cart_rank(newComm, neighbourCoords, posNeighbour);
   printf("DEBUG: Thread No: [%d] Thread Size: [%d] My Neighbours are: [%d] and: [%d]\n",rank, size, negNeighbour[0], posNeighbour[0]);
   */

  int *sourceProcess = calloc (ndims, sizeof(int));
  int *destProcess = calloc (ndims, sizeof(int));
  MPI_Cart_shift(newComm, ndims, 1, sourceProcess, destProcess);
  printf("DEBUG: Thread No: [%d] Thread Size: [%d] My Neighbours are: [%d] and: [%d]\n",rank, size, sourceProcess[0], destProcess[0]);

  /*
   filling(currentfield, w, h);
   for (int t = 0; t < timesteps; t++) {
   // TODO consol output
   //     show(currentfield, w, h);
   //writeVTK(currentfield, w, h, t, "output");
   int changes = evolve(currentfield, newfield, w, h, status, comm, rank, size);
   if (changes == 0) {
   sleep(3);
   break;
   }

   // usleep(200000);

   //SWAP
   unsigned *temp = currentfield;
   currentfield = newfield;
   newfield = temp;
   }
   */

  free(currentfield);
  free(newfield);
}

int main(int c, char **v) {
  int rank, size;
  MPI_Status status;
  MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Init(&c, &v);
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  int w = 0, h = 0, timesteps = 10;
  if (c > 1) w = atoi(v[1]); ///< read width
  if (c > 2) h = atoi(v[2]); ///< read height
  if (c > 3) timesteps = atoi(v[3]);
  if (w <= 0) w = 30; ///< default width
  if (h <= 0) h = 30; ///< default height

  int length = w * h;
  w = (length / size) * (rank + 1);
  h = 1;

  /*
   int posNeighbour = (rank + 1) % size;
   int negNeighbour = (rank - 1 + size) % size;
   printf("DEBUG: Thread No: [%d] Thread Size: [%d] My Neighbours are: [%d] and: [%d]\n",rank, size, negNeighbour, posNeighbour);#
   */
  game(w, h, timesteps, status, comm, rank, size);
  
  MPI_Finalize();
  return 0;
}
