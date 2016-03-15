#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <endian.h>
#include <mpi.h>
#include <time.h>

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
    sprintf(name, "out/%s_%d.vtk", prefix, t);
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


int evolve(unsigned* currentfield, unsigned* newfield, int w, int h) {
    int changes = 0;


    /*
     for (int y = 0; y < h; y++) {
     for (int x = 0; x < w; x++) {
     int numberOfAliveCells = getNumberOfAliveCells(currentfield, w, h, x, y);
     int currentfieldValue = currentfield[calcIndex(w, x,y)];
     //Dead Field and 3 Living Fields -> resurrect Field       //Living Field with 2 or 3 living neighbours stays alive
     newfield[calcIndex(w, x,y)] = ((!currentfieldValue && numberOfAliveCells == 3) || (currentfieldValue && (numberOfAliveCells == 3 || numberOfAliveCells == 4)));
     if (newfield[calcIndex(w, x,y)] != currentfield[calcIndex(w, x,y)])
     {
     changed = 1;
     }

     }
     }
     */
    //TODO if changes == 0, the time loop will not run!
    return changes;
}

void filling(unsigned* currentfield, int w, int h, int rank) {
    clock_t t = time(NULL);
    srand((int) t * rank);
    for (int i = 0; i < h*w; i++) {
        currentfield[i] = (rand() < RAND_MAX / 10) ? 1 : 0; ///< init domain randomly
    }
}

void calculateNeighbourProcesses(int *myCoords, MPI_Comm comm, int *rightProcess, int *leftProcess)
{
    int *queryCoords = calloc(2, sizeof(unsigned));

    queryCoords[0] = myCoords[0];
    queryCoords[1] = myCoords[1] + 1;
    MPI_Cart_rank(comm, queryCoords, rightProcess);

    queryCoords[0] = myCoords[0];
    queryCoords[1] = myCoords[1] - 1;
    MPI_Cart_rank(comm, queryCoords, leftProcess);

    free(queryCoords);
}

void game(int w, int h, int timesteps, int rank, int size, MPI_Comm comm, MPI_Status status) {
    w = w / size;
    unsigned *currentfield = calloc(w*h, sizeof(unsigned));
    unsigned *newfield     = calloc(w*h, sizeof(unsigned));

    MPI_Comm newComm;
    int ndims = 2;
    int *dims = calloc(ndims, sizeof(unsigned));
    dims[0] = 1;
    dims[1] = size;
    int *periods = calloc(ndims, sizeof(int));
    periods[0] = 1;
    periods[1] = 1;
    int reorder = 0;

    MPI_Cart_create(comm, ndims, dims, periods, reorder, &newComm);
    int *myCoords = calloc(ndims, sizeof(int));
    MPI_Cart_coords(newComm, rank, ndims, myCoords);

    int *rightProcess = calloc (1, sizeof(int));
    int *leftProcess = calloc (1, sizeof(int));

    calculateNeighbourProcesses(myCoords, newComm, rightProcess, leftProcess);

    printf("DEBUG: Hello, my ID is [%d],  my right neighbours ID is [%d] my left neighbours ID is [%d]\n", rank, rightProcess[0], leftProcess[0]);

    filling(currentfield, w, h, rank);
    /*
     for (int t = 0; t < timesteps; t++) {

     writeVTK(currentfield, w, h, t, "output");
     int changes = evolve(currentfield, newfield, w, h);
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

    free(rightProcess);
    free(leftProcess);

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
    if (w <= 0) w = 32; ///< default width
    if (h <= 0) h = 32; ///< default height
    game(w, h, timesteps, rank, size, comm, status);
    MPI_Finalize();
    return 0;
}
